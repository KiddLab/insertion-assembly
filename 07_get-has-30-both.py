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

USAGE = """07_get-has-30-both.py --in <assembly summary file name>

      requires 30bp on both ends of call

"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inFileName', help = 'assembly output summary file name')


(options,args)=parser.parse_args()

if options.inFileName is None:
    parser.error('file name for parsed assembly table not given')
    
#############################################################################

outFileName = options.inFileName + '.30both'
print 'Writing output to',outFileName


minRMsize = 30
minSide = 30


inFile = open(options.inFileName,'r')
outFile = open(outFileName,'w')
for line in inFile:
    ol = line
    line = line.rstrip()
    line = line.split()
    if line[1] == 'HAS_SEQ':
        if int(line[2]) > 0:
            parts = line[3:]
            hits = []
            for i in parts:
                i = i.split(':')
                np = []
                np = [i[0],int(i[1]),int(i[2])]
                eB = int(i[3].split('-')[0])
                eE = int(i[3].split('-')[1])
                eS = eE - eB + 1                
                np.append(eS)
                hits.append(np)
            if check_hasMin_TwoSides(hits) is True:
                outFile.write(ol)


outFile.close()
inFile.close()