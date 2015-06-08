import numpy, pylab
import optparse


def main(input):
    fhandle=open(input)
    ncolumns=2
    columns=dict()
    for x in range(0, ncolumns):
        columns[x]=[]
    for line in fhandle.readlines():
        if len(line.split()) < ncolumns:
            print "skipping", line
            continue
        for x in range(0, ncolumns):
            tmp=line.split()[x]
            if tmp[0]=='.':
                columns[x].append('0'+tmp)
                continue
            if tmp[0:1]=='-.':
                columns[x].append('-0.'+tmp[2:])
                continue
            columns[x].append(tmp)
    fhandle.close()
    ohandle=open('fix-%s' % input, 'w')
    for (x,y) in zip(columns[0], columns[1]):
        ohandle.write('%0.2f\t%0.8f\n' % (float(x), float(y)))
    ohandle.close()

def parse_cmdln():
    parser=optparse.OptionParser()
    parser.add_option('-i','--input',dest='input',type='string')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":	
    (options,args)=parse_cmdln()
    main(input=options.input)

