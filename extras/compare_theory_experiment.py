import numpy, pylab
import optparse


def main(spec_type, exp1, theory1, exp2, theory2=None):
    freqs=dict()
    spec=dict()
    spec['theory1']=numpy.loadtxt(theory1, usecols=(1,)) 
    freqs['theory1']=numpy.loadtxt(theory1, usecols=(0,)) 
    spec['exp1']=numpy.loadtxt(exp1, usecols=(1,)) 
    freqs['exp1']=numpy.loadtxt(exp1, usecols=(0,)) 
    spec['exp2']=numpy.loadtxt(exp2, usecols=(1,)) 
    freqs['exp2']=numpy.loadtxt(exp2, usecols=(0,)) 
    if theory2!=None: 
        spec['theory2']=numpy.loadtxt(theory2, usecols=(1,)) 
        freqs['theory2']=numpy.loadtxt(theory2, usecols=(0,)) 
    if spec_type=='ir':
        pylab.figure(figsize=(12,10))
        pylab.plot(freqs['exp1'], spec['exp1'], label=exp1.split('fix-')[1].split('AV')[0])
        pylab.ylim(0,0.8)
        pylab.xlim(1100,1725)
        pylab.ylabel('A')
        pylab.xlabel('Wavenumbers cm$^{-1}$')
        pylab.title(spec_type)
        pylab.legend()
        pylab.savefig('exp_ir.png', dpi=300)
        pylab.figure(figsize=(12,10))
        pylab.plot(freqs['theory1'], spec['theory1'], label='theory1')
        if theory2!=None: 
            pylab.plot(freqs['theory2'], spec['theory2'], label='theory2')
        pylab.xlim(1100,1725)
        pylab.ylabel('A')
        pylab.xlabel('Wavenumbers cm$^{-1}$')
        pylab.title(spec_type)
        pylab.legend()
        pylab.savefig('theory_ir.png', dpi=300)
    if spec_type=='vcd':
        pylab.figure(figsize=(12,10))
        pylab.plot(freqs['exp1'], spec['exp1'], label=exp1.split('fix-')[1].split('AV')[0])
        pylab.plot(freqs['exp2'], spec['exp2'], label=exp2.split('fix-')[1].split('AV')[0])
        pylab.xlim(1100,1725)
        pylab.ylabel('dA')
        pylab.xlabel('Wavenumbers cm$^{-1}$')
        pylab.title(spec_type)
        pylab.legend()
        pylab.savefig('exp_vcd.png', dpi=300)
        pylab.figure(figsize=(12,10))
        pylab.plot(freqs['theory1'], [i*1 for i in spec['theory1']], label='theory1=AM1944')
        if theory2!=None: 
            pylab.plot(freqs['theory2'], [i*1 for i in spec['theory2']],label='theory2=AMG1943')
        pylab.xlim(1100,1725)
        pylab.ylabel('dA')
        pylab.xlabel('Wavenumbers cm$^{-1}$')
        pylab.title(spec_type)
        pylab.legend()
        pylab.savefig('theory_vcd.png', dpi=300)
    pylab.show()

def parse_cmdln():
    parser=optparse.OptionParser()
    parser.add_option('--spec',dest='spec',type='string')
    parser.add_option('--exp1',dest='exp1',type='string')
    parser.add_option('--exp2',dest='exp2',type='string')
    parser.add_option('--theory1',dest='theory1',type='string')
    parser.add_option('--theory2',dest='theory2',type='string')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":	
    (options,args)=parse_cmdln()
    main(spec_type=options.spec, exp1=options.exp1, theory1=options.theory1, exp2=options.exp2, theory2=options.theory2)





