import subprocess
import sys
import signal
import os
import math


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
def get_site_intervals_from_table(bpIntervalTable):
    calls = []
    inFile = open(bpIntervalTable,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        id = line[0]
        if id == 'siteID':
             continue        
        c = line[1]
        
#        if id not in ['chr5_78331579']:
#            continue
        
#        c = c.replace('chr','') #for 1kg   # this should be an option, for now keep 'chr'
        b = int(line[4])
        e = int(line[5])
        p = int(line[2])
        calls.append([id,[c,b,e,p]])
    inFile.close()    
    # now I want to sort the list
    calls.sort(key = lambda x: x[1][1])
    calls.sort(key = lambda x: chromname_to_number(x[1][0])) # sort by chrom name
    return calls
###############################################################################
def chromname_to_number(c):
    if 'chr' == c[0:3]:
        cnum = c[3:]
    else:
        cnum = c
    if cnum.isdigit() is True:
        return int(cnum)
    if cnum == 'X':
        return 50
    if 'Y' in cnum:
        return 51
    if 'M' in cnum:
        return 52
    if 'Un' in cnum:
        return 53
    return cnum   # just return it...
###############################################################################
def open_bam_read(fileName,reg=''):
    # note -- skip over non-primary alignments...
    if reg == '':
        cmd = 'samtools view -F 256 ' + fileName
    else:
        cmd = 'samtools view -F 256 ' + fileName + ' ' + reg    
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # To deal with fact that might close file before reading all
    try:
        inFile = os.popen(cmd, 'r')
    except:
        print "ERROR!! Couldn't open the file " + fileName + " using samtools view -c\n"
        sys.exit(1)
    return inFile
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
    res['matePos'] = myLine[7]    
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
    # check for proper seqlen to update, 2015-05-05
    if myLine[9] == '*':  #not actually sequence present in SAM line
        res['seqLen'] = res['cigarCounts']['M'] + res['cigarCounts']['I']  + res['cigarCounts']['S'] + res['cigarCounts']['H']             
    return res
#####################################################################
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
###############################################################################
def make_seq_reads_file(samFileName,outFileName):
    inFile = open(samFileName,'r')
    outSeq = open(outFileName,'w')
    for line in inFile:
        line = line.rstrip()
        line = line.split('\t')
        samRec = parse_sam_line(line)
        seq = samRec['seq']
        qual = samRec['qual']
        if samRec['reverseStrand'] is True:
            qual = qual[::-1]
            seq = revcomp(seq)
        if samRec['isFirst'] is True:
            i = 0
        else:
            i = 1
        nl = [samRec['seqName'],str(i),seq,qual]
        nl = '\t'.join(nl) + '\n'
        outSeq.write(nl)
    inFile.close()
    outSeq.close()
###############################################################################
###############################################################################
def match_intervals_to_reads(outputReadFile):
    intervalsToReads = {}
    inFile = open(outputReadFile,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        id = line[0]
        readName = line[1]
        if id not in intervalsToReads:
            intervalsToReads[id] = {}
        intervalsToReads[id][readName] = 1
    inFile.close()
    return intervalsToReads
###############################################################################
def make_name_to_seq_dictionary(seqFileName):
    nameToSeq = {}
    inFile = open(seqFileName,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split('\t')
        rN = line[0]
        rNum = int(line[1])
        seq = line[2]
        qual = line[3]
        nameToSeq[(rN,rNum)] = [seq,qual]
    inFile.close()
    return nameToSeq
###############################################################################
###############################################################################
def write_fastq_for_site(myData,siteData):
    out1 = open(siteData['fq1'],'w')
    out2 = open(siteData['fq2'],'w')
    if siteData['siteID'] not in myData['intervalsToReads']:
        siteData['hasReads'] = False
        return
    else:
        siteData['hasReads'] = True
    
    readNames = myData['intervalsToReads'][siteData['siteID']]
    readNames = readNames.keys()
    readNames.sort()
    print 'Have %i read names' % len(readNames)
    for rN in readNames:
        out1.write('@%s\n%s\n+\n%s\n' % (rN,myData['nameToSeq'][(rN,0)][0],myData['nameToSeq'][(rN,0)][1]))
        out2.write('@%s\n%s\n+\n%s\n' % (rN,myData['nameToSeq'][(rN,1)][0],myData['nameToSeq'][(rN,1)][1]))
    out1.close()
    out2.close()
###############################################################################
###############################################################################
# do BWA alignment, in order to get the counts for downstream work.
# final result will be a BAM/SAM file
def align_to_alts_bwa(myData,siteData):
    if siteData['hasReads'] is False:
        #make fake empty files so subsequent steps will run
        siteData['outSAM'] = siteData['mappingOutDir'] + 'mapped.sam'
        siteData['outSamFilter'] = siteData['outSAM'] + '.filter'
        siteData['outSamSel'] = siteData['outSamFilter'] + '.sel'        

        outFile = open(siteData['outSamSel'],'w')
        outFile.close()
        siteData['outSamFilter'] = siteData['outSAM'] + '.filter'
        outFile = open(siteData['outSamSel'],'w')
        outFile.close()
        
        
                
        return
    
    siteData['targetFA'] = myData['alleleBase'] + 'locusAlleles/' + siteData['siteID'] + '/alleles.fa'
    print siteData['targetFA']
    siteData['outSAM'] = siteData['mappingOutDir'] + 'mapped.sam'


    # hard coded in -- this is for BWA-mem
    if False:
        cmd = myData['bwa'] + ' mem  -I 400,40,900,100 -M  ' +  siteData['targetFA'] + ' ' + siteData['fq1'] + ' ' + siteData['fq2'] + ' > ' + siteData['outSAM']
        print cmd
        runCMD(cmd)
    else:
        sai1 = siteData['fq1'] + '.sai'
        cmd = myData['bwa'] + ' aln -q 15 ' + siteData['targetFA'] + ' ' + siteData['fq1'] + ' > ' + sai1
        print cmd
        runCMD(cmd)
        sai2 = siteData['fq2'] + '.sai'
        cmd = myData['bwa'] + ' aln -q 15 ' + siteData['targetFA'] + ' ' + siteData['fq2'] + ' > ' + sai2
        print cmd
        runCMD(cmd)
        
        cmd = myData['bwa'] + ' sampe -A ' + siteData['targetFA'] + ' ' + sai1 + ' ' + sai2 + ' ' + siteData['fq1'] + ' ' + siteData['fq2'] + '> ' + siteData['outSAM']
        print cmd
        runCMD(cmd)
                        
        
    select_hits_from_sam(siteData)
###############################################################################
def select_hits_from_sam(siteData):
    siteData['outSamFilter'] = siteData['outSAM'] + '.filter'
    cmd = 'samtools view -S -F 256 %s > %s' % (siteData['outSAM'],siteData['outSamFilter'])
    print cmd
    runCMD(cmd)
    
    siteData['outSamSel'] = siteData['outSamFilter'] + '.sel'
    inFile = open(siteData['outSamFilter'],'r')
    outFile = open(siteData['outSamSel'],'w')
    while True:
        line = inFile.readline()
        if line == '':
            break
        if line[0] == '@':
            outFile.write(line)
            continue
        # should be order of pairs
        ol1 = line
        line = line.rstrip()
        line = line.split()
        samParse1 = parse_sam_line(line)  
        line = inFile.readline()
        ol2 = line
        line = line.rstrip()
        line = line.split()
        samParse2 = parse_sam_line(line)  
        
        if samParse1['seqName'] != samParse2['seqName'] :
            print 'Names not equal'
            print samParse1['seqName'],samParse2['seqName']
            sys.exit()

        if (samParse1['mapQ'] == 0) or (samParse2['mapQ'] == 0):
            continue
        
        if (samParse1['unMapped'] is True) or (samParse2['unMapped'] is True):
            continue
        outFile.write(ol1)
        outFile.write(ol2)        
    inFile.close()
    outFile.close()
###############################################################################
# functions for dealing with genotype likelihoods...
def read_samsel_hits(siteData):
   inFile = open(siteData['outSamSel'],'r')
   siteData['refQuals'] = []
   siteData['altQuals'] = []

   while True:
        line = inFile.readline()
        if line == '':
            break
        if line[0] == '@':
            continue
        # should be order of pairs
        ol1 = line
        line = line.rstrip()
        line = line.split()
        samParse1 = parse_sam_line(line)  
        line = inFile.readline()
        ol2 = line
        line = line.rstrip()
        line = line.split()
        samParse2 = parse_sam_line(line)
        # do some checks
        if samParse1['seqName'] != samParse2['seqName']:
            print 'seqnames do not match',data['siteID']
            sys.exit()
        if samParse1['chrom'] != samParse2['chrom']:
            print 'seq maped chroms not match',siteData['siteID']
            continue
        
        minMapQ = samParse1['mapQ']
        if samParse2['mapQ'] < minMapQ:
            minMapQ = samParse2['mapQ']
        
        if 'genome' in samParse1['chrom']:
            siteData['refQuals'].append(minMapQ)
        else:
            siteData['altQuals'].append(minMapQ)
   inFile.close()
   siteData['numRefFrag'] = len(siteData['refQuals'])
   siteData['numAltFrag'] = len(siteData['altQuals'])
   siteData['totFrags'] = siteData['numRefFrag']  + siteData['numAltFrag']
   siteData['refErrorProbs'] = phred_to_error_prop(siteData['refQuals'])
   siteData['altErrorProbs'] = phred_to_error_prop(siteData['altQuals'])
   
   # we will impose a max number of reads to prevent underflow issues
   # this tends to happen at cases where we have duplication/centromere stuff
   # going on anyway
   maxReads = 100
   if len(siteData['refErrorProbs']) > maxReads:
       siteData['refErrorProbs'] = siteData['refErrorProbs'][0:maxReads]
       print 'trimmed down ref counts to',len(siteData['refErrorProbs'])

   if len(siteData['altErrorProbs']) > maxReads:
       siteData['altErrorProbs'] = siteData['altErrorProbs'][0:maxReads]
       print 'trimmed down alt counts to',len(siteData['altErrorProbs'])


   siteData['numRefFrag'] = len(siteData['refErrorProbs'])
   siteData['numAltFrag'] = len(siteData['altErrorProbs'])
   siteData['totFrags'] = siteData['numRefFrag']  + siteData['numAltFrag']
   
###############################################################################
def phred_to_error_prop(myList):
    pList = []
    for i in myList:
        i = (-1)*i
        i = i/(10.0)
        i = 10**i
        pList.append(i)
    return pList
###############################################################################
def calc_gen_likelihood(siteData):
# based on equation from Heng Li paper
# we will do it VCF style, 0 --> hom ref, 1 --> het, 2 --> hom alt
    # het
    siteData['gl_1'] = 1.0/(2**siteData['totFrags'])

    # for now, we will do it in regular space, can convert to log space later if we need to
    # due to precision issues at higher coverage
    # hom ref
    p = 1.0
    for i in siteData['refErrorProbs']:
        p = p * (1.0-i)
    for i in siteData['altErrorProbs']:
        p = p * i    
    siteData['gl_0'] = p
    # hom alt
    p = 1.0
    for i in siteData['altErrorProbs']:
        p = p * (1.0-i)

    for i in siteData['refErrorProbs']:
        p = p * i    
    siteData['gl_2'] = p

    #put in a check for underflow to 0...
    minVal = 1e-200
    if siteData['gl_0'] <= minVal:
        siteData['gl_0'] = minVal

    if siteData['gl_1'] <= minVal:
        siteData['gl_1'] = minVal

    if siteData['gl_2'] <= minVal:
        siteData['gl_2'] = minVal
    make_scaled_likelihoods(siteData)
###############################################################################
def make_scaled_likelihoods(data):    
    maxLikelihood = max(data['gl_0'],data['gl_1'],data['gl_2'])
    scaledLikelihoods = [data['gl_0']/maxLikelihood,data['gl_1']/maxLikelihood,data['gl_2']/maxLikelihood]    
    print data['gl_0'],data['gl_1'],data['gl_2']
    for i in range(3):
        scaledLikelihoods[i] = math.log10(scaledLikelihoods[i])
        scaledLikelihoods[i] = -10.0 * scaledLikelihoods[i]
        scaledLikelihoods[i] = int(round(scaledLikelihoods[i]))
    data['scaledLikelihoods'] = scaledLikelihoods
###############################################################################
###############################################################################
def calc_gq(gLikeList,i):
    totP = sum(gLikeList)
    v = gLikeList[i]/totP
    v = 1.0 - v
    if v == 0.0 or v < 0.000001:
        v = 0.000001

    v = -10.0*math.log10(v)
    return v
###############################################################################



