# check levels in VCF, output 6,7,8 calls only

from optparse import  OptionParser

###############################################################################
USAGE = """
select-retroseq-level678.py --in <input VCF> --out <output VCF>

Counts calls in each filter levle in RetroSeq VCF and output new VCF
consiting only of calls in level 6,7,8.


"""

parser = OptionParser(USAGE)
parser.add_option('--in',dest='vcfIn', help = 'vcf infput file')
parser.add_option('--out',dest='vcfOut', help = 'vcf out file')

(options, args) = parser.parse_args()

if options.vcfIn is None:
    parser.error('vcf input not given')
if options.vcfOut is None:
    parser.error('vcf output not given')

###############################################################################

#setup count bins for each level
counts = []
for i in range(0,9):
    counts.append(0)


inFile = open(options.vcfIn,'r')
outFile = open(options.vcfOut,'w')
for line in inFile:
    if line[0] == '#':
        continue
    ol = line

    line = line.rstrip()
    line = line.split()
    g = line[-1]
    fl = g.split(':')[2]
    fl = int(fl)
    counts[fl] += 1
    if fl >= 6:
        outFile.write(ol)
        
inFile.close()

sum = 0
for i in range(1,9):
    sum += counts[i]
    print 'flag %i\t%i' % (i,counts[i])
print 'Total\t%i' % sum

outFile.close()