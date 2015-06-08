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

def lorentzian(freqs, ir, gamma, res):
    vX=freqs
    vY=ir
    dLorHalfWidth=100
    dSpectRoundUp=100
    bImagePersist=1
    dPlotOffset=20
    resRecip=1.0/res # get freqs at resolution of 0.1
    vX=[(round(x*resRecip))/resRecip for x in vX]
    yMax=((numpy.ceil(numpy.max(vX)/dSpectRoundUp))*dSpectRoundUp) + dSpectRoundUp
    # create lorentzian to convolve
    vXTemp=numpy.arange(-dLorHalfWidth, dLorHalfWidth, res)
    x0=0
    vL=1.0/((numpy.pi*gamma)*(1+((vXTemp-x0)/gamma)**2))  
    vL=vL/max(vL)
    vSin=numpy.sin(numpy.arange(0,numpy.pi, numpy.pi/(2*dLorHalfWidth/res)))
    vL1=vL*vSin
    # do convolve
    vInConvA=numpy.zeros((yMax/res))
    for i in range(0, len(vX)):
        iIdx=round(vX[i]*10)
        vInConvA[iIdx]=vY[i]
    vXOut=numpy.arange(res, yMax+res, res)
    nYOut=numpy.convolve(vInConvA, vL1, 'same')
    size=len(nYOut)
    return vXOut[:size], nYOut[:size]


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

def apply_energy_window(window, spectra, energies, type):
    specs=['or', 'ir', 'vcd']
    import copy
    files=copy.copy(energies[type]['files'])
    for file in files:
        if energies[type][file] >=window:
            print file, energies[type][file], 'REMOVE'
            energies[type].pop(file)
            index=energies[type]['files'].index(file)
            energies[type]['files'].pop(index)
            spectra.pop(file)
        else:
            pass
    return spectra, energies
        

def get_weights(free_energies):
    sum=0
    boltz=dict()
    for file in free_energies['files']: 
        sum+=numpy.exp(-free_energies[file]/0.6)
    for file in free_energies['files']: 
        boltz[file]=numpy.exp(-free_energies[file]/0.6)/sum
    return boltz

def check_for_duplicates(type, energies, spectra, remove=False):
    unique=[]
    d_degree=1.0
    d_ene=0.5
    import copy
    files=copy.copy(energies[type]['files'])
    for file in files:
        pair=(round(energies[type][file], 5), round(spectra[file]['or'],2))
        duplicate=False
        if len(unique)!=0:
            for entry in unique:
                if abs(entry[0]-pair[0]) < d_ene:
                    if abs(entry[1]-pair[1]) < d_degree:
                        print pair, "(kcal/mol, OR) matched", entry
                        duplicate=True
                        if remove==True:
                            print "REMOVED DUPLICATE", file
                            spectra.pop(file)
                            energies[type].pop(file)
                            index=energies[type]['files'].index(file)
                            energies[type]['files'].pop(index)
                            break
        else:
            unique.append(pair)
            continue
        if duplicate==False:
            unique.append(pair)
    return spectra

def write_parsed_output(spectra, scale_factor):
    specs=['ir', 'vcd']
    for file in spectra.keys():
        numpy.savetxt('%s.freq' % (file.split('.log')[0]),spectra[file]['freqs'])
        numpy.savetxt('%s.freq.scaled' % (file.split('.log')[0]),[i*scale_factor for i in spectra[file]['freqs']])
        for spec in specs:
            ohandle=open('%s.%s' % (file.split('.log')[0], spec), 'w')
            for (f,s) in zip([i*scale_factor for i in spectra[file]['freqs']], spectra[file][spec]):
                ohandle.write('%0.4f\t%0.4f\n' % (f,s))
            ohandle.close()
    return

def get_optical_rotation(type, energies, spectra, boltz):
    or_sum=0
    or_file=open('sort_or_values.txt', 'w')
    print "sorted unique QM energy (kcal/mol) + OR values in sort_or_values.txt" 
    for file in energies[type]['files']:
        or_sum+=spectra[file]['or']*boltz[file]
        or_file.write('%s\t%0.4f\t%0.4f\n' % (file, round(energies[type][file],5), spectra[file]['or']))
    return or_sum

def plot_spectra(weighted_spectra, boltz):
    specs=['vcd', 'ir']
    for spec in specs:
        weights=[boltz[file] for file in boltz.keys()]
        weights=numpy.array(weights)
        order=numpy.argsort(weights)
        hi_to_low_weights=weights[order][::-1]
        print "SUM OF WEIGHTS", numpy.sum(weights)
        for (n, i) in enumerate(hi_to_low_weights):
            if numpy.sum(hi_to_low_weights[:n+1]) >= 0.90:
                max_weight_for_plot=i
                break
        pylab.figure()
        for weight in hi_to_low_weights:
            for file in boltz.keys():
                if boltz[file]==weight:
                    if weight >= max_weight_for_plot:
                        print "plotting", weight
                        pylab.plot(weighted_spectra[spec]['scaled_all_freqs'], weighted_spectra[spec][file], linewidth=weight*10, label=round(weight,3))
        pylab.plot(weighted_spectra[spec]['scaled_all_freqs'], weighted_spectra[spec]['final'], color='k', linewidth=4)
        pylab.xlim(1100, 1725)
        pylab.xlabel('wavenumbers (cm$^{-1}$)')
        pylab.ylabel('%s Intensity' % spec)
        pylab.legend()
        pylab.savefig('%s.png' % spec, dpi=300)
    pylab.show()

#############################################################

def main(prefix,pop,scale_factor, gamma, res, window, plot=False, removedup=False):
    if pop=='dg':
        type='free_energy'
    else:
        type='energy'
    float(res)
    gamma=float(gamma)
    scale_factor=float(scale_factor)
    files=glob.glob('%s*.log' % prefix)
    spectra=dict()
    #parse Gaussian log files for spectra info
    for file in files:
        spectra=parse(spectra, file)
    energies=sort_kcal(spectra, type)  #sort by energy
    if removedup==True:    # check for removedup
        print "CHECKING FOR DUPLICATES, WITH REMOVAL"
        check_for_duplicates(type, energies, spectra,remove=True)
    else:
        print "CHECKING FOR DUPLICATES, NO REMOVAL"
        check_for_duplicates(type, energies, spectra)
    if window!=None:
        window=float(window)
        print "REMOVING CONFS ABOVE %s %s" % (window, type)
        spectra, energies=apply_energy_window(window, spectra, energies, type)
    write_parsed_output(spectra, scale_factor) #write parsed output
    numpy.savetxt('sort_%s.txt' % type, [(i, energies[type][i]) for i in energies[type]['files']], fmt='%s')
    print "sorted unique QM (kcal/mol) values in sort_%s.txt" % type
    # get boltzmann weights from energy
    boltz=get_weights(energies[type])
    ohandle=open('weight.txt', 'w')
    for pair in sorted(boltz.iteritems(), key=operator.itemgetter(1)):
        ohandle.write('%s\t%s\n' % (pair[0], pair[1])) # write conf, boltz wt
    # get OR
    or_sum=get_optical_rotation(type, energies, spectra, boltz)
    print "OR Rotation: %s deg." % round(or_sum,2)
    weighted_spectra=dict()
    specs=['vcd', 'ir']
    for spec in specs:
        weighted_spectra[spec]=dict()
        weighted_spectra[spec]['final']=numpy.zeros((20)) # initial zero array
        for file in energies[type]['files']:
            unscaled_new_freq, convolved_spec=lorentzian(spectra[file]['freqs'], spectra[file][spec], gamma, res)
            unscaled_freq, convolved_spec=check_size(unscaled_new_freq, convolved_spec)
            weighted_spectra[spec][file]=convolved_spec
            ofile=open('%s.%s_convolved.txt' % (file.split('.log'),spec), 'w')
            scaled_new_freq=[i*scale_factor for i in unscaled_new_freq]
            for (f, s) in zip(scaled_new_freq, convolved_spec):
                ofile.write('%0.2f\t%0.8f\n' % (f,s))
            ofile.close()
            if sum(weighted_spectra[spec]['final'])==0:
                weighted_spectra[spec]['final']=boltz[file]*numpy.array(convolved_spec)
                weighted_spectra[spec]['scaled_all_freqs']=[i*scale_factor for i in unscaled_new_freq]
            else:
                weighted_spectra[spec]['final'], convolved_spec=check_size(weighted_spectra[spec]['final'], convolved_spec)
                weighted_spectra[spec]['final']+=boltz[file]*numpy.array(convolved_spec)
        ofile=open('final_weighted_%s.txt' % spec, 'w')
        for (f, s) in zip(weighted_spectra[spec]['scaled_all_freqs'], weighted_spectra[spec]['final']):
            ofile.write('%0.2f\t%0.8f\n' % (f,s))
        ofile.close()
        print "Boltzmann averaged %s spectra and scaled freqs in final_weighted_%s.txt" % (spec, spec)
    if plot==True:
        plot_spectra(weighted_spectra, boltz)

def parse_cmdln():
    parser=optparse.OptionParser()
    parser.add_option('--prefix',dest='prefix',type='string', help='Prefix for Gaussian log files.')
    parser.add_option('--scale',dest='scale_factor',type='string', help='Frequency scalar depending on theory level. Default=0.96 (for B3LYP-631G*).', default=0.96)
    parser.add_option('--pop',dest='pop',type='string',help='Selects total energies (de) or entropy-adjusted free energies (dg) to determine population. Default is dg.', default='dg')
    parser.add_option('--gamma',dest='gamma',type='string', help='Scales Lorenztian for peak reproduction. Default=4.', default=4.0)
    parser.add_option('--res',dest='resolution',type='string', help='Frequency resolution of Lorenztian. Default=0.1', default=0.1)
    parser.add_option('--window',dest='window',type='string',help='de/dg (depending on pop flag) cutoff in kcal/mol for conformations in prefix*.log to use. Default is no cutoff.', default=None)
    parser.add_option('--removedup', action="store_true", dest="removedup", help="If this flag is used, will remove duplicate conformations based on pop and OR. Without the flag, a check is run for duplicates but no removal.")
    parser.add_option('--plot', action="store_true", dest="plot", help="If this flag is used, spectra will be plotted.")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":
    (options,args)=parse_cmdln()
    # defaults
    plot=False
    removedup=False
    if options.removedup==True:
        removedup=True
    if options.plot==True:
        plot=True
    main(prefix=options.prefix, pop=options.pop, removedup=removedup,
scale_factor=options.scale_factor, gamma=options.gamma, res=options.resolution,
window=options.window, plot=plot)

