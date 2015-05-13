# Jeff Kidd
# genutils
# This version does not include fastq functions that require Cython
# this makes it easier to distribute to others

import subprocess
import sys
import os
import signal
import tempfile
import numpy as np


#########################################################################################
# merges bed file intervals, uses call to bed tools
def merge_beds(bedInts):
    #change this to required temp dir you have access to
    tmpDir = '/home/jmkidd/kidd-lab-scratch/jmkidd-projects/tmp/'
    tmpDir = '.'
    
    tmpFile = tempfile.NamedTemporaryFile(dir=tmpDir,suffix='.bed',delete=False)
    tmpFileName = tmpFile.name
    for i in bedInts:
        i = [str(j) for j in i]
        nl = '\t'.join(i) + '\n'
        tmpFile.write(nl)
    tmpFile.close()
    sortedName = tmpFileName + '.sorted'
    cmd = 'sortBed -i %s > %s ' % (tmpFileName,sortedName)
    runCMD(cmd)
    mergedName= sortedName + '.merged'
    cmd = 'mergeBed -i %s > %s' % (sortedName,mergedName)
    runCMD(cmd)
    mergedInts = []
    inFile = open(mergedName,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        mergedInts.append([line[0],int(line[1]),int(line[2])] )
    inFile.close()
    # cleanup
    os.unlink(tmpFileName)
    os.unlink(sortedName)
    os.unlink(mergedName)
    return mergedInts
#########################################################################################    


##############################################################################
# Reads param file, and returns dictioanry of values
def read_params(inFileName):
    myParams = {}
    inFile = open(inFileName,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split('\t')
        myParams[line[0]] = line[1]
    inFile.close()
    return myParams
##############################################################################
# Returns complement of a bp.  If not ACGT then return same char
def complement(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'a':
        return 't'
    if c == 't':
        return 'a'
    if c == 'c':
        return 'g'
    if c == 'g':
        return 'c'
    # If not ACTGactg simply return same character
    return c   
##############################################################################
# Returns the reverse compliment of sequence 
def revcomp(seq):
    c = ''    
    seq = seq[::-1] #reverse
    # Note, this good be greatly sped up using list operations
    seq = [complement(i) for i in seq]
    c = ''.join(seq)
    return c
##############################################################################
# Takes ambiguity code, and returns list of two alleles 
def ambig2Nucs(c):
    if c == 'K':
        return ['G','T']
    if c == 'M':
        return ['A','C']
    if c == 'R':
        return  ['A','G']
    if c == 'Y':
        return ['C','T']
    if c == 'S':
        return ['C','G']
    if c == 'W':
        return ['A','T']
    if c == 'A':
        return ['A','A']
    if c == 'C':
        return ['C','C']
    if c == 'G':
        return ['G','G']
    if c == 'T':
        return ['T','T']
    # N is unknown
    if c == 'N':
        return ['N','N']
    # If get to here, then there is an error and the gentoype is not known
    print 'error in ambig2Nucs, what is: ',c
    exit(1)

##############################################################################
def is_het(c):
    if c  in ['A','C','T','G','a','c','t','g']:
        return False

    if c  in ['n','N']:
        return False
        
    if c  in ['K','M','R','Y','S','W']:
        return True
    
    print 'Error in is_het, what is=',c
    sys.exit(1)
##############################################################################

# Takes two nucleotides, and returns ambiguity code
# Returns 'N'  for unknowns
def nucs2Ambig(n1,n2):
    n1 = n1.upper()
    n2 = n2.upper()
    
    # fix for fq phased files that have het sites in them already
    if n1 in ['K','M','R','Y','S','W']:
        return 'N'
    if n2 in ['K','M','R','Y','S','W']:
        return 'N'

    
    if n1 in ['0','.','N','-','X'] :
        return 'N'
    if n2 in ['0','.','N','-','X'] :
        return 'N'
    if n1 == n2 :
        return n1
        
    if n1 in ['G','T'] and n2 in ['G','T'] :
        return 'K'
    if n1 in ['A','C'] and n2 in ['A','C'] :
        return 'M'        
    if n1 in ['A','G'] and n2 in ['A','G'] :
        return 'R'        
    if n1 in ['C','T'] and n2 in ['C','T'] :
        return 'Y'
    if n1 in ['C','G'] and n2 in ['C','G'] :
        return 'S'
    if n1 in ['A','T'] and n2 in ['A','T'] :
        return 'W'

    # If get to here, then something is wrong
    print 'error in nucs2Ambig, what is: ',n1,n2
    exit(1)

##############################################################################
# Takes 'A,B' and 'C' returns B if C == A, returns A if C == B, else exits
# Usefull for getting the 'other' SNP allele not found in a read
def otherAllele(gen,c):
    a,b = gen.split(',')
    if c == a:
        return b
    if c == b:
        return a
    print 'Eror in otherAllele: gen = ',gen, 'c = ',c
    exit(1)

##############################################################################
# Takes ambiguity code and allele, returns true if allele is represented in the
# ambiguity code, returns false otherwise
def matchAmbig(gen,c):
    gen = ambig2Nucs(gen)
    (a,b) = gen[0],gen[1]
    if c == a or c == b :
        return True
    return False
##############################################################################
# A and B are the SNP alleles, classifies the change
def classify_transition_transversion(a,b):
    a = a.upper()
    b = b.upper()
    if a == b:
        print 'not a change!',a,b
        sys.exit(1)
    if a not in ['A','C','G','T'] :
        print a,'is not in ACGT'
        sys.exit(1)
    if b not in ['A','C','G','T'] :
        print b,'is not in ACGT'
        sys.exit(1)

    nucType = {}
    nucType['A'] = 'purine'
    nucType['G'] = 'purine'
    nucType['T'] = 'pyrimidine'
    nucType['C'] = 'pyrimidine'

    if nucType[a] == nucType[b]:
        return 'transition'
    else:
        return 'transversion'

##############################################################################
# get non-ref allele count.  Assume that the site is biallelic, with onyl 1 alt allele
# returns tuple of (ref,nonref) int counts
def get_nonref_allele_count(line):
    refCount = 0
    nonRefCount = 0
    
    # Assume that first field of sample is the GT field
    for i in range(9,len(line)):
        sampleGen = line[i].split(':')[0]
        if sampleGen[0] == '.':
            continue
        if sampleGen[0] == '0':
            refCount += 1
        elif sampleGen[0] == '1':
            nonRefCount += 1
        if len(sampleGen) == 1:
            continue
        if sampleGen[2] == '0':
            refCount += 1
        elif sampleGen[2] == '1':
            nonRefCount += 1
    return (refCount,nonRefCount)
##############################################################################

# Returns true if chromosome is an autosome and not randomon or alternative
# False if chrm contains '_', 'chrM', 'chrY', 'chrX'
def isAuto(c):
    if c == 'chrM' :
        return False
    if c == 'chrX' :
        return False
    if c == 'chrY' :
        return False
    if '_' in c :
        return False
    return True

# Returns true if chromosome is an autosome or chrX
# False is chrm contains '_', 'chrM', 'chrY', 
def isAutoOrX(c):
    if c == 'chrM' :
        return False
    if c == 'chrX' :
        return True
    if c == 'chrY' :
        return False
    if '_' in c :
        return False
    return True

###############################################################################
# Function for filling out concordancy matrix
# matrix is 4 X 3 (last row is for doesn't match alleles)
# Returns the row,col index to be incremented 

#         Array
#       0HomRef 1Het 2HomAlt 3nocall
# 0 HomRef
# 1 Het
# 2 HomAlt
# 3 NoCall
# 4 AllelesMisMatch

def classify_concordance(refA, array1, array2, seq1, seq2):            
    row = -1
    col = -1
    # Array is 'N'
    if array1 =='N' or array2 == 'N':
        col = 3
        if seq1 == seq2 and seq1 == refA:
            row = 0
        elif seq1 == seq2 and seq1 != refA and seq1 != 'N' :
            row = 2
        elif seq1 == seq2 and seq1 == 'N' :
            row = 3
        elif seq1 != seq2:
            row = 1        

    # Array is homoz ref
    elif array1 == array2 and array1 == refA:
        col = 0
        if seq1 == seq2 and seq1 == refA:
            row = 0
        elif seq1 == seq2 and seq1 != refA and seq1 != 'N' :
            row = 2
        elif seq1 == seq2 and seq1 == 'N' :
            row = 3
        elif seq1 != seq2:
            row = 1

    # Array is homoz alt (but not 'N'!)
    elif array1 == array2 and array1 != refA:
        col = 2
        if seq1 == seq2 and seq1 == refA:
            row = 0
        elif seq1 == seq2 and seq1 != refA and seq1 != 'N' :
            row = 2
        elif seq1 == seq2 and seq1 == 'N' :
            row = 3
        elif seq1 != seq2:
            row = 1

    # Array is het 
    elif array1 != array2:
        col = 1
        if seq1 == seq2 and seq1 == refA:
            row = 0
        elif seq1 == seq2 and seq1 != refA and seq1 != 'N' :
            row = 2
        elif seq1 == seq2 and seq1 == 'N' :
            row = 3
        elif seq1 != seq2:
            if seq1 in [array1,array2] and seq2 in [array1,array2]:
                row = 1
            else :
                row = 4

    if row == -1 or col == -1 :
        print 'Do not know how to classify: ',refA,array1,array2,seq1,seq2
        print row, col
        exit(1)
    return (row,col)                    
###############################################################################
# Make string for printing out concordance matrix
def build_concordance_matrix_string(stats):
    nl = '\t\tArray Genotype\n'
    nl = nl + '\t' + '\t'.join(['Hom Ref','Het','Hom Alt','No call']) + '\n'
    nl = nl + 'Seq Hom Ref\t' + '\t'.join([str(i) for i in stats[0] ]) + '\n'
    nl = nl + 'Seq Het\t' + '\t'.join([str(i) for i in stats[1] ]) + '\n'
    nl = nl + 'Seq Hom Alt\t' + '\t'.join([str(i) for i in stats[2] ]) + '\n'
    nl = nl + 'Seq Nocall\t' + '\t'.join([str(i) for i in stats[3] ]) + '\n'
    nl = nl + 'Allele Mismatch\t' + '\t'.join([str(i) for i in stats[4] ]) + '\n'
    return nl    
###############################################################################
def calculate_total_concordance(stats):
    tot = 0
    agree = 0
    for r in range(3):
        for c in range(3):
            tot += stats[r][c]
            if r == c :
                agree += stats[r][c]
    if tot == 0 :
        return 0
    perAgree = 100.0*float(agree)/tot
    return perAgree
###############################################################################
def calculate_hom_alt_concordance(stats):
    tot = 0
    agree = 0
    c = 2
    for r in range(3):
        tot += stats[r][c]
        if r == c :
            agree += stats[r][c]
    if tot == 0 :
        return 0
    perAgree = 100.0*float(agree)/tot
    return perAgree
###############################################################################    
def calculate_het_concordance(stats):
    tot = 0
    agree = 0
    c = 1
    for r in range(3):
        tot += stats[r][c]
        if r == c :
            agree += stats[r][c]
    if tot == 0 :
        return 0
    perAgree = 100.0*float(agree)/tot
    return perAgree
###############################################################################
# Function for determining whether or not two bed intervals overlap
# returns True if intervals overlap, otherwise returns False
# options are chrom1, s1,e1, chrom2, s2,e2
# function assumes intervals are UCSC bed-style (0-based half open...)
# assumes that starts and ends are integers and that chroms are strings
def overlap_bed(chrom1, start1, end1, chrom2, start2, end2):
    if chrom1 != chrom2:
        return False
    # now, make sure things are in order, so put in left_s and left_e
    # and switch over to 1 based coordinates
    if start1 <= start2:
        left_start = start1 + 1
        left_end = end1
        right_start = start2 + 1
        right_end = end2
    else:
        left_start = start2 + 1
        left_end = end2
        right_start = start1 + 1
        right_end = end1
    
    # in order now, so just need to see if right_s is 'inside' of the left interval
    if right_start >= left_start and right_start <= left_end:
        return True
    else :
        return False

###############################################################################
def get_intersect_interval(targetC,targetB,targetE,intervalList):
    res = []
    for i in intervalList:
        overlap = overlap_bed(targetC,targetB,targetE,i[0],i[1],i[2])
        if overlap is True:
            res.append(i)
    return res
###############################################################################
# Function for chaning from genome names to plink names (chr1 --> 1, etc)
# does not deal with pseduo autosomal regions.  And, only works for human
def chrom_to_plink(chrom):
    numPart = chrom[3:]
    if numPart == 'X' :
        return '23'
    if numPart == 'Y' :
        return '24'
    if numPart == 'M' :
        return '26'
    # check to see if numPart is an integer, should be!
    intPart = int(numPart) #fails else
    return numPart
###############################################################################
# Function for chaning from genome names to plink names (chr1 --> 1, etc)
# does not deal with pseduo autosomal regions.  And, only works for human
def plink_to_chrom(numPart):
    if numPart == '23' :
        return 'chrX'
    if numPart == '25' :
        return 'chrX'
    if numPart == '24' :
        return 'chrY'
    if numPart == '26' :
        return 'chrM'
    return 'chr' + numPart

###############################################################################
# Helper function to run commands,
# doesn't check return or print to log.  Use for 'grep' so that doesn't
# fail if no items found
def runCMDNoFail(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
###############################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print 'command failed'
        print cmd
        sys.exit(1)
###############################################################################
# Helper function to run commands, handle return values and print to log file
def runCMD_output(cmd):
    val = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE)
    resLines = []
    for i in val.stdout:
       i = i.rstrip()
       resLines.append(i)
    return resLines
#############################################################################        
# Helper function to read in information from genome .fai file and return
# a dictionary containing chrom names and lengths
def read_chrom_len(faiFileName):
    chromLens = {}
    inFile = open(faiFileName,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        chromLens[line[0]] = int(line[1])
    inFile.close()
    return chromLens    
#############################################################################        
def read_chrom_len_list(faiFileName):
    chromLens = []
    inFile = open(faiFileName,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        chromLens.append([line[0],int(line[1])])
    inFile.close()
    return chromLens        
#############################################################################        
# Helper function to read in information from genome cm file and return
# a dictionary containing chrom names and CM lengths
def read_chrom_len_cm(geneticMapFile):
    chromLens = {}
    inFile = open(geneticMapFile,'r')
    for line in inFile:
        line = line.rstrip()
        if line[0] == '#' or line[:5] == 'chrom':
            continue
        line = line.split()
        chrom = line[0]
        cm = float(line[3])
        
        if chrom not in chromLens:
            chromLens[chrom] = cm
        else:
            if cm > chromLens[chrom]:
                chromLens[chrom] = cm
    inFile.close()
    return chromLens    
#############################################################################    
# Makes a dictionary of the info field in a vcf file
# returns the dictionary
#example: DP=5;AF1=1;CI95=0.5,1;DP4=0,0,4,1;MQ=51
def parse_vcf_info(infoField):
    info = {}
    infoList = infoField.split(';')
    for field in infoList:
       if field.count('=') == 1:
           (name,vals) = field.split('=')[0:2]
           vals = vals.split(',')
           if vals == 'true' or vals == 'True':
               vals = True
           if vals == 'false' or vals == 'False':
               vals = False
           info[name] = vals
       else:
           info[field] = 'PRESENT'
    return info
#############################################################################                 
# Return index of largest number in list
def find_largest_index(tmpList,listLen):
    largest = 0
    for i in range(1,listLen):
        if (tmpList[i] > tmpList[largest]):
            largest = i
    return largest        
#############################################################################         
# Return index of smallest number in list
def find_smallest_index(tmpList,listLen):
    smallest = 0
    for i in range(1,listLen):
        if (tmpList[i] < tmpList[smallest]):
            smallest = i
    return smallest        
#############################################################################         
# Return list of genotype possibilities
def get_gen_order_PL(alleleList):
    numAlleles = len(alleleList)
    PLlist = []
    for i in range(numAlleles):
        for j in range(i,numAlleles):
            PLlist.append(nucs2Ambig(alleleList[i],alleleList[j]))
    return PLlist    

#    if numAlleles == 1:
#        PLlist.append(alleleList[0])
#        return PLlist
#    if numAlleles == 2:
#        PLlist.append(alleleList[0])
#        PLlist.append(nucs2Ambig(alleleList[0],alleleList[1]))
#        PLlist.append(alleleList[1])
#        return PLlist
#    if numAlleles == 3:
# Note: I think the description on the mpileup webpage is wrong
# this order makes more sense with some examples, and is easier to figure out
# algorithmically.
#        PLlist.append(alleleList[0])
#        PLlist.append(nucs2Ambig(alleleList[0],alleleList[1]))
#        PLlist.append(nucs2Ambig(alleleList[0],alleleList[2]))
#        PLlist.append(alleleList[1])
#        PLlist.append(nucs2Ambig(alleleList[1],alleleList[2]))
#        PLlist.append(alleleList[2])
#        PLlist.append(alleleList[0])
#        PLlist.append(nucs2Ambig(alleleList[0],alleleList[1]))
#        PLlist.append(alleleList[1])
#        PLlist.append(nucs2Ambig(alleleList[0],alleleList[2]))
#        PLlist.append(nucs2Ambig(alleleList[1],alleleList[2]))
#        PLlist.append(alleleList[2])
#        return PLlist
#    print 'ERROR: there are',numAlleles
#    print 'How to put them together?  I am not sure of the order!'
#    sys.exit(1)
#############################################################################         

#############################################################################         
# Return index of param
def get_param_index(target,myList):
    for i in range(len(myList)):
        if myList[i] == target:
            return i
    return -1
#############################################################################         
#############################################################################         
def get_max_index(myList):
    max = myList[0]
    max_i = 0
    for i in range(len(myList)):
        if myList[i] > max:
            max = myList[i]
            max_i = i
    return max_i
#############################################################################         
def count_nonrefalleles_in_genotype(gen):
    nonRef = 0
    myGen = gen.split('/')
    for i in myGen:
        if i != '0' and i != '.':
            nonRef += 1
    return nonRef    
#####################################################################
def open_gzip_read(fileName):
    gc = 'gunzip -c ' + fileName
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # To deal with fact that might close file before reading all
    try:
        inFile = os.popen(gc, 'r')
    except:
        print "ERROR!! Couldn't open the file " + fileName + " using gunzip -c\n"
        sys.exit(1)
    return inFile
#####################################################################
def open_gzip_write(fileName):
    try:
        gc = 'gzip > ' + fileName
        outFile = os.popen(gc, 'w')
    except:
        print "ERROR!! Couldn't open the output file " + fileName+ " (with gzip)\n"
        sys.exit(1)
    return outFile
#####################################################################
def open_bam_read(fileName,reg=''):
    if reg == '':
        cmd = 'samtools view  ' + fileName
    else:
        cmd = 'samtools view  ' + fileName + ' ' + reg    
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # To deal with fact that might close file before reading all
    try:
        inFile = os.popen(cmd, 'r')
    except:
        print "ERROR!! Couldn't open the file " + fileName + " using samtools view -c\n"
        sys.exit(1)
    return inFile
#####################################################################
def get_readgroups_from_bam(bamFileName):
    cmd = 'samtools view -H ' + bamFileName
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # To deal with fact that might close file before reading all
    try:
        inFile = os.popen(cmd, 'r')
    except:
        print "ERROR!! Couldn't open the file " + fileName + " using samtools view -H\n"
        sys.exit(1)
    
    readGroups = []
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        if line[0] != '@RG':
            continue
        id = line[1]
        id = id.split(':')[1]
        readGroups.append(id)
    return readGroups
#####################################################################

def get_allele_counts_vcf_line(line):
    hasMissing = False
    counts = [0,0] # assume is only two alleles...
    for i in range(9,len(line)):
        gen = line[i].split(':')[0]
        if gen[0] == '.':
            hasMissing = True
            continue
        a1 = int(gen[0])
        counts[a1] += 1
        if len(gen) == 3:
            a2 = int(gen[2])
            counts[a2] += 1
    tot = counts[0] + counts[1]
    if counts[0] > counts[1]:
        maf = float(counts[1])/tot
    else:
        maf = float(counts[0])/tot
    counts.append(maf)
    counts.append(hasMissing)
    return counts
#####################################################################
def get_allele_counts_vcf_line_4options(line):
    hasMissing = False
    counts = [0,0,0,0] # assume is only 4 alleles...
    for i in range(9,len(line)):
        gen = line[i].split(':')[0]
        if gen[0] == '.':
            hasMissing = True
            continue
        a1 = int(gen[0])
        counts[a1] += 1
        if len(gen) == 3:
            a2 = int(gen[2])
            counts[a2] += 1
    tot = counts[0] + counts[1] + counts[2] + counts[3]
    counts.append(hasMissing)
    return counts
#####################################################################
def get_dps_vcf_line(line):
    DPpersample = []
    desc = line[8]
    desc = desc.split(':')
    DPindex = get_param_index('DP',desc)
    if DPindex == -1:
#        print 'NO DP INDEX'
#        print line
#        sys.exit()
         return [0]
    for i in range(9,len(line)):
        if line[i][0] == '.':
            dp = 0
        else:        
            sample = line[i].split(':')
            dp = sample[DPindex]
            if dp == '.': # because of BEAGLE imputation
                dp = 0
            else:
                dp = int(dp)

        DPpersample.append(dp)
    return DPpersample
#####################################################################
# returns N50, N90, Nk etc, assumes that data is reverse sorted, k is integer
# so will do k/100
def get_nx(data,totBp,k):
    targetBp = (k/100.0) * totBp
    sum = 0
    for i in range(0,len(data)):
        sum += data[i]
        if sum >= targetBp:
            break
    return data[i]
#####################################################################
# take list of data, will return mean, median, total, sum, and N50
def get_len_stats(data):
    myRes = {}
    #first, sort the data, in reverse
    data.sort(reverse = True)
    numRecords = 0
    totLen = 0
    for i in data:
        numRecords += 1
        totLen += i
    myRes['numRecords'] = numRecords
    myRes['totLen'] = totLen
    
    myRes['median'] = np.median(data)
    myRes['mean'] = np.mean(data)
    myRes['min'] = data[-1]
    myRes['max'] = data[0]
    
    myRes['n50'] = get_nx(data,myRes['totLen'],50)
    myRes['n10'] = get_nx(data,myRes['totLen'],10)
    myRes['n90'] = get_nx(data,myRes['totLen'],90)

    return myRes
#####################################################################
def get_axt_record(myFile):
    # is an array, data[0] is header (split already), then 
    # rest is sequence as string
    # if last item, then data will be empty
    data = []    
    
    
    header = myFile.readline()
    if header == '':
        return data
    while header[0] == '#':
        header = myFile.readline()
    header = header.rstrip()
    header = header.split()
    data.append(header)
    while True:
        line = myFile.readline()
        if line == '\n':
           break
        line = line.rstrip()
        data.append(line)
    return data
#####################################################################
def get_4l_record_len(myFile):
    #fastq style file...
    # just return sequence len
    # -1 if last record
    myLine1 = myFile.readline()
    if myLine1 == '':
        return -1
    myLine2 = myFile.readline()
    myLine3 = myFile.readline()
    myLine4 = myFile.readline()
    myLine2 = myLine2.rstrip()
    return len(myLine2)
#####################################################################
def get_4l_record_seq(myFile):
    #fastq style file...
    myLine1 = myFile.readline()
    if myLine1 == '':
        return ''
    myLine2 = myFile.readline()
    myLine3 = myFile.readline()
    myLine4 = myFile.readline()
    myLine1 = myLine1.rstrip()
    myLine1 = myLine1.split()[0]
    myLine1 = myLine1.replace('@','')
    myLine2 = myLine2.rstrip()
    return [myLine1,myLine2]
#####################################################################
def get_4l_record(myFile):
    #fastq style file...
    # just return sequence len
    # -1 if last record
    myLine1 = myFile.readline()
    if myLine1 == '':
        return ''
    myLine2 = myFile.readline()
    myLine3 = myFile.readline()
    myLine4 = myFile.readline()
    return [myLine1,myLine2,myLine3,myLine4]
#####################################################################
def get_2l_record(myFile):
    #2 line fasta style file...
    myLine1 = myFile.readline()
    if myLine1 == '':
        return ''
    myLine2 = myFile.readline()
    return [myLine1,myLine2]
#####################################################################
def get_2l_record_len(myFile):
    #fasta one line style file...
    # just return sequence len
    # -1 if last record
    myLine1 = myFile.readline()
    if myLine1 == '':
        return -1
    myLine2 = myFile.readline()
    myLine2 = myLine2.rstrip()
    return len(myLine2)
#####################################################################
# Calculate Vst as described by Redon et al 2006
# includes weights by sample size
# input is a list of lists, output is single value of vst
# use sample variances (n-1) for the calculation
def calc_vst(myData):
    ddVal = 1 # to do n-1 for the calculation
    numPops = len(myData)
#    print 'Have %i populations' % numPops
    totData = []
    for i in range(len(myData)):
        totData.extend(myData[i])
    numTot = len(totData)
    varTot = np.var(totData,ddof = ddVal)
    VsList = []
    aveVar = 0.0
    for i in range(len(myData)):       
        v = np.var(myData[i],ddof = ddVal)
        VsList.append(v)
        aveVar += (v*len(myData[i]))
    aveVar = aveVar / float(numTot)
    Vst = (varTot - aveVar) / varTot
    return Vst
#####################################################################
# calculates Vst, sets value < 0 to 0
def calc_vst_clip0(myData):
    vst = calc_vst(myData)
    if vst < 0.0:
        vst = 0.0
    return vst
#####################################################################
def print_sam_dict(myRec):
    keys = myRec.keys()
    keys.sort()
    for k in keys:
        print k,myRec[k]
#####################################################################
def parse_sam_line(myLine):
    res = {}
    res['seqName'] = myLine[0]
    res['flag'] = int(myLine[1])
    res['chrom'] = myLine[2]
    res['chromPos'] = int(myLine[3])
    res['mapQ'] = int(myLine[4])
    res['cigar'] = myLine[5]
    res['seq'] = myLine[9]
    res['seqLen'] = len(myLine[9])
    
    
    res['cigarExpand'] = expand_cigar(res['cigar'])
    res['qual'] = myLine[10]
    res['mateChrom'] = myLine[6]
    res['fragLen'] = int(myLine[8])
    
    res['cigarCounts']={}
    res['cigarCounts']['M'] = 0
    res['cigarCounts']['D'] = 0
    res['cigarCounts']['I'] = 0
    res['cigarCounts']['S'] = 0
    res['cigarCounts']['H'] = 0

    
    if res['flag'] & 0x10 != 0:
        res['reverseStrand'] = True
    else:
        res['reverseStrand'] = False

    if res['flag'] & 0x4 != 0:
        res['unMapped'] = True
    else:
        res['unMapped'] = False

    if res['flag'] & 0x400 != 0:
        res['isDuplicate'] = True
    else:
        res['isDuplicate'] = False

    if res['flag'] & 0x100 != 0:
        res['notPrimaryAlignment'] = True
    else:
        res['notPrimaryAlignment'] = False



    if res['flag'] & 0x1 != 0:
        res['isPaired'] = True
    else:
        res['isPaired'] = False


    if res['flag'] & 0x8 != 0:
        res['mateUnmapped'] = True
    else:
        res['mateUnmapped'] = False


    if res['flag'] & 0x40 != 0:
        res['isFirst'] = True
    else:
        res['isFirst'] = False



    
    
    for i in res['cigarExpand']:
        res['cigarCounts'][i[1]] += i[0]
        
    return res
#####################################################################
#returns lists of [int,flag]
def expand_cigar(cigar):
    res = []
    if cigar == '*':
        return res
    digits = ['0','1','2','3','4','5','6','7','8','9']
    accumulate = ''
    i = 0
    while True:
        if i == len(cigar):
            break
        if cigar[i] in digits:
            accumulate += cigar[i]
            i += 1
        else:
            d = int(accumulate)
            res.append([d,cigar[i]])
            i += 1
            accumulate = ''
    return res
#####################################################################
def read_list_from_expanded_cigar(eC):
    readList = []
    for i in eC:
        if i[1] != 'D':
            for j in range(i[0]):
                readList.append(i[1])
    return readList
#####################################################################
def get_rg_from_sam_line(myLine):
    rg = ''
    for i in range(10,len(myLine)):
        if myLine[i][0:2] == 'RG':
            rg = myLine[i].split(':')[-1]
            return rg
    return rg
#####################################################################
    
# assumes that both are the same length
def get_seq_mismatches(seq1,seq2):
    n = 0
    for (a,b) in zip(seq1,seq2):
        if a != b:
            n += 1
    return n       
#####################################################################
def init_blank_list(listLen,element):
    myList = []
    for i in range(listLen):
        myList.append(element)
    return myList
#####################################################################
def read_fasta_to_string(fastaFile):
    myString = ''
    inFile = open(fastaFile,'r')
    for line in inFile:
        line = line.rstrip()
        if line[0] == '>':
            continue
        myString += line
    inFile.close()
    return myString
###############################################################################
def read_fasta_file_to_list(fastaFile):
    myDict = {}
    inFile = open(fastaFile,'r')
    line = inFile.readline()
    line = line.rstrip()
    if line[0] != '>':
        print 'ERROR, FILE DOESNNOT START WITH >'
        sys.exit()
    myName = line[1:]
    myDict[myName] = {}
    myDict[myName]['seq'] = ''
    myDict[myName]['seqLen'] = 0    
    mySeq = ''
    while True:
        line = inFile.readline()
        if line == '':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            break
        line = line.rstrip()
        if line[0] == '>':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            myName = line[1:]
            myDict[myName] = {}
            myDict[myName]['seq'] = ''
            myDict[myName]['seqLen'] = 0    
            mySeq = ''
            continue
        mySeq += line
    inFile.close()
    return myDict
###############################################################################
def read_rm_file(rmFileName):
    rmLines = []
    inFile = open(rmFileName,'r')
    for line in inFile:
        if line == '\n':
            continue
        line = line.rstrip()
        line = line.split()
        if line[0] == 'There':
            return []
        if line[0] == 'SW':
            continue
        if line[0] == 'score':
            continue
        rmLines.append(line)
    inFile.close()
    return rmLines
###############################################################################
def add_breaks_to_line(seq,n=50):
    myList = []
    myList = [i for i in seq]
    newList = []
    c = 0
    for i in myList:
        newList.append(i)
        c += 1
        if c % n == 0 and c != (len(myList)):
            newList.append('\n')
    myStr = ''.join(newList)
    return myStr    
###############################################################################
def make_genotype_dictionary(vcf):
    vcfDict = {}
    inFile = open(vcf,'r')
    for line in inFile:
        line = line.rstrip()
        if line[0:2] == '##':
            continue
        line = line.split('\t')
        if line[0] == '#CHROM':
            header = line
            continue
        c = line[0]
        p = line[1]
        ref = line[3]
        alt = line[4]
        alt = alt.split(',')
        alleles = []
        alleles.append(ref)
        alleles.extend(alt)
        # assume is all PASS...
        for i in range(9,len(header)):
            sn = header[i]
            gen = line[i]
            a1 = (gen[0])
            a2 = (gen[2])
            if a1 != '.' and a2 != '.':
                a1 = alleles[int(a1)]
                a2 = alleles[int(a2)]
            genCode = nucs2Ambig(a1,a2)
            vcfDict[(sn,c,p)] = genCode
    inFile.close()
    return vcfDict
###############################################################################
def read_fq_like(myFile):
    myRec = {}
    line = myFile.readline()
    line = line.rstrip()
    myRec['name'] = line
    mySeq = ''
    while True:
        line = myFile.readline()
        if line == '+\n':
            break
        line = line.rstrip()
        mySeq += line
    myRec['seq'] = mySeq
    myRec['seqLen'] = len(myRec['seq'])
    myQual = ''
    while True:
        line = myFile.readline()
        if line == '':
            break
        line = line.rstrip()
        myQual += line
    myRec['qual'] = myQual
    myRec['qualLen'] = len(myRec['qual'])
    return myRec
###############################################################################
def print_dictionary(d):
    k = d.keys()
    k.sort()
    for i in k:
        print i,d[i]
###############################################################################
def print_dictionary_keys(d):
    k = d.keys()
    k.sort()
    for i in k:
        print i
###############################################################################
###############################################################################
def get_fq_like_multiple(myFile):
    myRec = {}
    line = myFile.readline()
    if line == '':
        return ''
    line = line.rstrip()
    myRec['name'] = line[1:]
    mySeq = ''
    numLines = 0
    while True:
        line = myFile.readline()
        if line == '+\n':
            break
        numLines += 1
        line = line.rstrip()
        mySeq += line
    myRec['seq'] = mySeq
    myRec['seqLen'] = len(myRec['seq'])
    myQual = ''
    for n in range(numLines):
        line = myFile.readline()
        line = line.rstrip()
        myQual += line
    myRec['qual'] = myQual
    myRec['qualLen'] = len(myRec['qual'])
    return myRec
###############################################################################
def determine_encoding_from_list_of_quals(qualList):
    qualCounts = []
    for i in range(127):
        qualCounts.append(0)
    
    for qualStr in qualList:
        for c in qualStr:
            qualCounts[ord(c)] += 1
    sangerCount = 0
    illuminaCount = 0
    for i in range(127):
        if i >=33 and i<= 58:
            sangerCount += qualCounts[i]
        if i >=75:
            illuminaCount += qualCounts[i]
#    print sangerCount, illuminaCount
    if sangerCount > 0 and illuminaCount == 0:
        return 'Sanger-33'
    if sangerCount == 0 and illuminaCount > 0:
        return 'Illumina-64'
#    print sangerCount,illuminaCount
#    for i in range(127):
#        print i,qualCounts[i]

    return 'Unknown'            
###############################################################################
def read_params_file_to_dict(myFile):
    myParams = {}
    inFile = open(myFile)
    for line in inFile:
        line = line.rstrip()
        line = line.split('\t')
        name = line[0]
        val = line[2]
        ty = line[1]
        if ty == 'str':
            myParams[name] = val
        
        else:
            print 'unknown param type'
            print ty
            print line
            sys.exit()
    inFile.close()
    return myParams
###############################################################################

