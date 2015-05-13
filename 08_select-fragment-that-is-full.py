import sys
from optparse import  OptionParser

##########################################################################################
def check_hasMin(hits):
    for i in hits:
        if i[3] >= minRMsize:
            return True
    return False
##########################################################################################
def check_hasMin_OneSide(hits):
    for i in hits:
        if i[3] >= minRMsize and ((i[1] >= minSide) or (i[2] >= minSide)):
            return True
    return False
##########################################################################################
def check_hasMin_TwoSides(hits):
    for i in hits:
        if i[3] >= minRMsize and (i[1] >= minSide) and (i[2] >= minSide):
            return True
    return False
##########################################################################################

#############################################################################

USAGE = """08_select-fragment-that-is-full.py --in <30 both file name>

"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inFileName', help = '30 both file name')


(options,args)=parser.parse_args()

if options.inFileName is None:
    parser.error('file name for parsed assembly table not given')
    
#############################################################################



minRMsize = 30
minSide = 30

outFileName = options.inFileName + '.sel'
print 'writing output to',outFileName

inFile = open(options.inFileName,'r')
outFile = open(outFileName,'w')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    id = line[0]
    contigs = line[3:]
    goodContigs = []
    for contig in contigs:
        i = contig.split(':')
        np = []
        np = [i[0],int(i[1]),int(i[2])]
        eB = int(i[3].split('-')[0])
        eE = int(i[3].split('-')[1])
        eS = eE - eB + 1
        np.append(eS)        
        hits=[]
        hits.append(np)
        if check_hasMin_TwoSides(hits) is True:
            goodContigs.append(contig)
#    if len(goodContigs) != 1:
#        print line
#        print goodContigs
#        sys.exit()
    nl = [id,str(len(goodContigs))]
    nl.extend(goodContigs)
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)

outFile.close()
inFile.close()