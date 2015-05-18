import pylab, numpy
import operator
import heapq
import glob
import optparse


def check_size(array1, array2):
    if len(array1)!=len(array2):
        print "size mismatch"
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
    print "---sorted %s---" % type
    energies[type]=dict()
    for_sort=[(i,data[i][type]) for i in data.keys()]
    sorted_files=heapq.nsmallest(num,for_sort, operator.itemgetter(1))
    energies[type]['files']=[i[0] for i in sorted_files]
    min=data[sorted_files[0][0]][type]
    for file in energies[type]['files']:
        new=627.5*(data[file][type]-min)
        energies[type][file]=new
    #    print file, new
    return energies 

def get_weights(free_energies):
    sum=0
    boltz=dict()
    for file in free_energies['files']: 
        sum+=numpy.exp(-free_energies[file]/0.6)
    for file in free_energies['files']: 
        boltz[file]=numpy.exp(-free_energies[file]/0.6)/sum
    return boltz

def main(prefix, type):
    files=glob.glob('%s*.log' % prefix)
    spectra=dict()
    for file in files:
        spectra=parse(spectra, file)
    energies=sort_kcal(spectra, type)
    numpy.savetxt('sort_%s_%s.dat' % (prefix, type), [(i, energies[type][i]) for i in energies[type]['files']], fmt='%s')
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
    for file in energies[type]['files']:
        print file, boltz[file], energies[type][file]
    specs=['ir', 'vcd', 'or']
    or_sum=0
    for spec in specs:
        final_spectra=numpy.zeros((20)) # initial zero array
        if spec!='vcd':
            pylab.figure()
        for file in energies[type]['files']:
            if spec=='or':
                or_sum+=spectra[file]['or']*boltz[file]
                if boltz[file] > 0.05:
                    print "OR Boltz %s:  %s deg." % (boltz[file], spectra[file]['or'])
            else:
                spectra[file]['new_freqs'], spectra[file]['new_%s' % spec]=lorentzian(spectra[file]['freqs'], spectra[file][spec])
                if sum(final_spectra)==0:
                    final_spectra=boltz[file]*numpy.array(spectra[file]['new_%s' % spec])
                else:
                    final_spectra+=boltz[file]*numpy.array(spectra[file]['new_%s' %
spec])
                pylab.vlines(x=spectra[file]['freqs'], ymin=[0]*len(spectra[file][spec]),ymax=spectra[file][spec])
                if boltz[file] > 0.05:
                    pylab.plot(spectra[file]['new_freqs'], spectra[file]['new_%s' % spec], linewidth=boltz[file]*4)
        if spec=='or':
            print "OR Rotation: %s deg." % round(or_sum,2)
        else:
            pylab.plot(spectra[file]['new_freqs'], final_spectra, linewidth=4)
            pylab.xlim(0, int(max(spectra[file]['new_freqs'])))
            pylab.xlabel('wavenumbers (cm$^{-1}$')
            pylab.ylabel('%s Intensity' % spec)
            pylab.savefig('%s.png' % spec, dpi=300)
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

