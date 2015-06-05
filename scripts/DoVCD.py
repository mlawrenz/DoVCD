import pylab, numpy
import operator
import heapq
import glob
import optparse


def check_size(array1, array2):
    if len(array1)!=len(array2):
        print "WARNING: size mismatch"
        print "trimming last vals of bigger array"
        val=min([len(array1), len(array2)])
        array1=array1[:val]
        array2=array2[:val]
    else:
        pass
    return array1, array2

def lorentzian(freqs, ir):
    vX=freqs
    vY=ir
    dRes=0.1
    dLorHalfWidth=100
    dGamma=9
    dSpectRoundUp=100
    bImagePersist=1
    dPlotOffset=20
    dResRecip=1.0/dRes # get freqs at resolution of 0.1
    vX=[(round(x*dResRecip))/dResRecip for x in vX]
    yMax=((numpy.ceil(numpy.max(vX)/dSpectRoundUp))*dSpectRoundUp) + dSpectRoundUp
    # create lorentzian to convolve
    vXTemp=numpy.arange(-dLorHalfWidth, dLorHalfWidth, dRes)
    x0=0
    vL=1.0/((numpy.pi*dGamma)*(1+((vXTemp-x0)/dGamma)**2))  
    vL=vL/max(vL)
    vSin=numpy.sin(numpy.arange(0,numpy.pi, numpy.pi/(2*dLorHalfWidth/dRes)))
    vL1=vL*vSin
    # do convolve
    vInConvA=numpy.zeros((yMax/dRes))
    for i in range(0, len(vX)):
        iIdx=round(vX[i]*10)
        vInConvA[iIdx]=vY[i]
    vXOut=numpy.arange(dRes, yMax+dRes, dRes)
    nYOut=numpy.convolve(vInConvA, vL1, 'same')
    size=len(nYOut)
    return 0.96*vXOut[:size], nYOut[:size]


def parse(spectra, file):
    data = spectra.setdefault(file,{'ir':[],'vcd':[], 'freqs':[]})
    fhandle=open(file)
    n=0
    complete=False
    for line in fhandle.readlines():
        if 'Sum of electronic and thermal Free Energies' in line:
            data['free_energy']=float(line.split()[7])
        if 'SCF Done' in line:
            data['energy']=float(line.split()[4]) # should end up final val
        if 'Frequencies' in line:
            data['freqs'].append(float(line.split()[2]))
            data['freqs'].append(float(line.split()[3]))
            data['freqs'].append(float(line.split()[4]))
        if 'IR Inten' in line:
            data['ir'].append(float(line.split()[3]))
            data['ir'].append(float(line.split()[4]))
            data['ir'].append(float(line.split()[5]))
        if 'Rot. str.' in line:
            data['vcd'].append(float(line.split()[3]))
            data['vcd'].append(float(line.split()[4]))
            data['vcd'].append(float(line.split()[5]))
        if '[Alpha]D (static)' in line:
            data['or']=float(line.split()[8])
        if 'Normal termination' in line:
            complete=True
    if complete==False:
        spectra.pop(file)
        print "---%s INCOMPLETE---" % file
    return spectra

def sort_kcal(data, type):
    num=len(data.keys())
    energies=dict()
    energies[type]=dict()
    for_sort=[(i,data[i][type]) for i in data.keys()]
    sorted_files=heapq.nsmallest(num,for_sort, operator.itemgetter(1))
    energies[type]['files']=[i[0] for i in sorted_files]
    min=data[sorted_files[0][0]][type]
    for file in energies[type]['files']:
        new=627.5*(data[file][type]-min)
        energies[type][file]=new
    return energies 

def get_weights(free_energies):
    sum=0
    boltz=dict()
    for file in free_energies['files']: 
        sum+=numpy.exp(-free_energies[file]/0.6)
    for file in free_energies['files']: 
        boltz[file]=numpy.exp(-free_energies[file]/0.6)/sum
    return boltz

def filter(spectra):
    print "---FILTERING DUPLICATE CONFS---"
    unique=[]
    d_degree=1.0
    d_ene=0.001
    for file in spectra.keys():
        pair=(round(spectra[file]['energy'], 5), round(spectra[file]['or'],2))
        duplicate=False
        if len(unique)!=0:
            for entry in unique:
                if abs(entry[0]-pair[0]) < d_ene:
                    if abs(entry[1]-pair[1]) < d_degree:
                        print pair, "matched", entry
                        duplicate=True
                        spectra.pop(file)
                        break
        if duplicate==False:
            unique.append(pair)
    return spectra

def main(prefix, type):
    files=glob.glob('%s*.log' % prefix)
    spectra=dict()
    specs=['ir', 'vcd']
    for file in files:
        spectra=parse(spectra, file)
    filter(spectra)
    for file in spectra.keys():
        numpy.savetxt('%s.freq.scaled' % (file.split('.log')[0]), [i*0.96 for i in spectra[file]['freqs']])
        for spec in specs:
            ohandle=open('%s.%s' % (file.split('.log')[0], spec), 'w')
            for (f,s) in zip([i*0.96 for i in spectra[file]['freqs']], spectra[file][spec]):
                ohandle.write('%0.4f\t%0.4f\n' % (f,s))
            ohandle.close()
    energies=sort_kcal(spectra, type)
    numpy.savetxt('sort_free_energies.txt', [(i, energies[type][i]) for i in energies[type]['files']], fmt='%s')
    print "---OUTPUT FILES---"
    print "sorted QM kcal/mol vals in sort_free_energies.txt"
    for file in energies[type]['files']:
        energy_val=round(energies[type][file], 3)
        or_val=round(spectra[file]['or'], 3)
        for i in energies[type]['files']:
            if i!=file:
                if round(energies[type][i], 3)==energy_val:
                    if round(spectra[i]['or'], 3)==or_val:
                        energies[type].pop(i)
                        spectra.pop(i)
                        index=energies[type]['files'].index(i)
                        energies[type]['files'].pop(index)
            else:
                pass
    boltz=get_weights(energies[type])
    specs=['vcd', 'ir', 'or']
    or_sum=0
    or_file=open('boltz_or_values.txt', 'w')
    print "unique confs with boltz wt and OR value in boltz_or_values.txt"
    for spec in specs:
        final_spectra=numpy.zeros((20)) # initial zero array
        if spec!='or':
            pylab.figure()
        for file in energies[type]['files']:
            if spec=='or':
                or_sum+=spectra[file]['or']*boltz[file]
                or_file.write('%s\t%0.4f\t%0.4f\n' % (file, round(boltz[file],5), spectra[file]['or']))
                #print "%s Boltz %s:  OR %s deg." % (file, round(boltz[file],5), spectra[file]['or'])
            else:
                spectra[file]['new_freqs'], spectra[file]['new_%s' % spec]=lorentzian(spectra[file]['freqs'], spectra[file][spec])
                if sum(final_spectra)==0:
                    final_spectra=boltz[file]*numpy.array(spectra[file]['new_%s' % spec])
                else:
                    final_spectra, spectra[file]['new_%s' % spec]=check_size(final_spectra, spectra[file]['new_%s' % spec])
                    final_spectra+=boltz[file]*numpy.array(spectra[file]['new_%s' % spec])
                #pylab.vlines(x=spectra[file]['freqs'], ymin=[0]*len(spectra[file][spec]),ymax=spectra[file][spec])
                if boltz[file] > 0.06:
                    spectra[file]['new_freqs'], spectra[file]['new_%s' % spec]=check_size(spectra[file]['new_freqs'], spectra[file]['new_%s' % spec])
                    pylab.plot(spectra[file]['new_freqs'], spectra[file]['new_%s' % spec], linewidth=boltz[file]*10, label=round(boltz[file],3))
                    ofile=open('weight_%s_%s_%s.txt' % (spec, round(boltz[file],3), spec), 'w')
                    for (f, s) in zip(spectra[file]['new_freqs'], spectra[file]['new_%s' % spec]):
                        ofile.write('%0.2f\t%0.8f\n' % (f,s))
                    ofile.close()
        if spec=='or':
            print "OR Rotation: %s deg." % round(or_sum,2)
        else:
            spectra[file]['new_freqs'], final_spectra=check_size(spectra[file]['new_freqs'], final_spectra)
            pylab.plot(spectra[file]['new_freqs'], final_spectra, color='k', linewidth=4)
            ofile=open('boltz_%s.txt' % spec, 'w')
            for (f, s) in zip(spectra[file]['new_freqs'],  final_spectra):
                ofile.write('%0.2f\t%0.8f\n' % (f,s))
            ofile.close()
            print "boltz averaged %s spectra and freqs in boltz_%s.txt" % (spec, spec)
            #pylab.xlim(0, int(max(spectra[file]['new_freqs'])))
            pylab.xlim(1100, 1725)
            pylab.xlabel('wavenumbers (cm$^{-1}$')
            pylab.ylabel('%s Intensity' % spec)
            pylab.legend()
            pylab.title(prefix)
            pylab.savefig('%s.png' % spec, dpi=300)
    or_file.close()
    pylab.show()

def parse_cmdln():
    parser=optparse.OptionParser()
    parser.add_option('-p','--prefix',dest='prefix',type='string', help='prefix for gaussian log files')
    parser.add_option('-e', action="store_true", dest="energy", help="perform calcs using energies, default uses free energies")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":
    (options,args)=parse_cmdln()
    if options.energy==True:
        type='energy'
    else:
        type='free_energy'
    main(prefix=options.prefix, type=type)

