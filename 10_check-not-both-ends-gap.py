import genutils
from optparse import  OptionParser


#############################################################################

USAGE = """10_check-not-both-ends-gap.py --in <selected longest fragment> --outdirbase <base directory for assembly output>

"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inFileName', help = 'sel longest fragment')
parser.add_option('--outdirbase',dest='outDirBase', help = 'base directory for assembly output')


(options,args)=parser.parse_args()

if options.inFileName is None:
    parser.error('file name for parsed assembly longest')
if options.outDirBase is None:
    parser.error('assembly base dir not given')

    
#############################################################################

delta = 30

assemDirBase = options.outDirBase

#print 'num is',options.num
inFile = open(options.inFileName)
outFileName = options.inFileName + '.notbothgap'
print 'checking both gap, output going to',outFileName

outFile = open(outFileName,'w')

bad = 0
good = 0
for line in inFile:
    ol = line
    line = line.rstrip()
    line = line.split()
    
#    print line

    k = line[0]
    c = k.split('_')
    c = c[0:-1]
    c = '_'.join(c)
    assemDir = assemDirBase + '/' + c + '/' + k + '/'
    scaffFA = assemDir + 'combined.scaffolds.fa'
    
    seqs = genutils.read_fasta_file_to_list(scaffFA)
    
    scaffName = line[2].split(':')[0]
    repeatPos = line[2].split(':')[-1]
    repeatPos = repeatPos.split('-')
    b = int(repeatPos[0])
    e = int(repeatPos[1])
    
    seq = seqs[scaffName]['seq']
#    print seq
#    print scaffName
#    print b,e
    l_b = b-1- delta + 1
    l_e = b-1
    r_b = e + 1
    r_e = r_b + delta - 1
    lSeq = seq[l_b-1:l_e]
    rSeq = seq[r_b-1:r_e]
#    print lSeq
#    print rSeq

    if lSeq == 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' and rSeq == 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN':
        #print 'skipping, ends in NN'
        bad += 1
    else:
        outFile.write(ol)
        good += 1

    
    
    
    
#    break


inFile.close()
outFile.close()

print 'good',good,'bad',bad