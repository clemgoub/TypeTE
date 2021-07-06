import subprocess
import sys
import signal
import os
import math
import pysam
import age_parse



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
# Helper function to run commands, handle return values and print to log file
def runCMD_noFail(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print 'return output'
        print val
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
def open_bam_read(fileName,reference,reg=''):
    # reference for CRAM
    # is FALSE if not given
    
    if reference is False:
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
    else:
        # note -- skip over non-primary alignments...
        if reg == '':
            cmd = 'samtools view -F 256 -T ' + reference + ' ' + fileName
        else:
            cmd = 'samtools view -F 256 -T ' + reference + ' ' + fileName + ' ' + reg    
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

    if res['flag'] & 0x2048 != 0:
        res['SupplementaryAlignment'] = True
    else:
        res['SupplementaryAlignment'] = False


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
        siteData['outSamSel'] = siteData['outSAM'] + '.filter.sel.sam'

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
        
#        cmd = myData['bwa'] + ' sampe -A ' + siteData['targetFA'] + ' ' + sai1 + ' ' + sai2 + ' ' + siteData['fq1'] + ' ' + siteData['fq2'] + '> ' + siteData['outSAM']
# test of a = 1000 for larger insert sizes.....
        cmd = myData['bwa'] + ' sampe -a 1000 -A ' + siteData['targetFA'] + ' ' + sai1 + ' ' + sai2 + ' ' + siteData['fq1'] + ' ' + siteData['fq2'] + '> ' + siteData['outSAM']

        print cmd
        runCMD(cmd)
                        
        
    select_hits_from_sam(siteData,myData)
###############################################################################
def select_hits_from_sam(siteData,myData):
    siteData['outBAMFilter'] = siteData['outSAM'] + '.filter.bam'
    cmd = 'samtools view -b -F 256 %s > %s' % (siteData['outSAM'],siteData['outBAMFilter'])
    print cmd
    runCMD(cmd)
    
    
    siteData['outSamSel'] = siteData['outSAM'] + '.filter.sel.sam'

    inBamFile = pysam.AlignmentFile(siteData['outBAMFilter'],'rb')
    outSamFile = pysam.AlignmentFile(siteData['outSamSel'],'w',template=inBamFile)

    readCache = []
    for read in inBamFile:
        if read.is_secondary is True:
            continue
        if read.is_supplementary is True:
            continue
        if len(readCache) == 0:
            readCache.append(read)        
        elif len(readCache) == 1:
            if read.query_name == readCache[0].query_name:
                r1 = readCache[0]
                r2 = read
                readCache.pop()
                res = process_read_hits(r1,r2,myData)  
                if res is True:
                    outSamFile.write(r1)
                    outSamFile.write(r2)
            else:
                print 'CACHE not match!!'
                print readCache[0]
                print read
                print 'CACHE not match!!'
                sys.exit()
    inBamFile.close()
    outSamFile.close()        
###############################################################################
# do the filtering of reads for output
def process_read_hits(r1,r2,myData):
    #check if either is unmapped    
    if r1.is_unmapped is True:
        return False        
    if r2.is_unmapped is True:
        return False

    if r1.mapping_quality == 0:
        return False
    if r2.mapping_quality == 0:
        return False
        
    if r1.reference_name != r2.reference_name:
        return False

    # check to see if mapped entirely within exclusion region
        
            
    chromName = r1.reference_name
    totBpInExclude = 0
    numR1NotInExclude = 0
    for i in r1.get_reference_positions():
        p = i +1
        if (chromName,p) not in myData['excludeDict']:
            numR1NotInExclude += 1
        else:
            totBpInExclude += 1

    chromName = r2.reference_name
    numR2NotInExclude = 0
    for i in r2.get_reference_positions():
        p = i +1
        if (chromName,p) not in myData['excludeDict']:
            numR2NotInExclude += 1
        else:
            totBpInExclude += 1
            
#    print 'Number of bp not in exclude regions',numR1NotInExclude,numR2NotInExclude
    t = numR1NotInExclude + numR2NotInExclude
    if t < 10: # min number bp outside of exclusion region...
        return False

    # require at least 5bp mapped in excluded region if there are any
    # exclusions on that allele.  This presents false calls of hets due to ~1bp mismatch
    
    if chromName in myData['excludeDict'] and totBpInExclude < 5:
        return False

    return True    
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
# setup regions that will be excluded if both read pairs map entirely within them
def setup_exclusion(myData):
     myData['excludeDict'] = {}
     
     if myData['excludeFileName'] is None:
         print 'No exclusion regions to check'
         return
     inFile = open(myData['excludeFileName'],'r')
     for line in inFile:
         line = line.rstrip()
         line = line.split()
         c = line[0]
         b = int(line[1])
         e = int(line[2])
         myData['excludeDict'][c] = 'yes'
         for i in range(b,e+1):
              myData['excludeDict'][(c,i)] = 1     
     inFile.close()
     print 'Set up %i bp to exclude' % len(myData['excludeDict'])
   
###############################################################################   

#### functions for running the split read module

#####################################################################
# check to see if program is in PATH
# copied from https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
# is easier in Python 3 apparently...
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None
#####################################################################
def setup_locus_for_split(myData,siteID,siteChrom,sitePos):
    print 'Processign siteID',siteID
    origLocusAllelesFasta = myData['origLocusDir'] + siteID + '/alleles.fa'
    print origLocusAllelesFasta
    newLocusOutDir = myData['newLocusDir'] + siteID
    if os.path.isdir(newLocusOutDir) is False:
        print 'making dir',newLocusOutDir
        cmd = 'mkdir ' + newLocusOutDir
        print cmd
        runCMD(cmd)
    newLocusOutDir += '/'
    newLocusAllelesFasta = newLocusOutDir + 'alleles.fa'
    cmd = 'cp %s %s' % (origLocusAllelesFasta,newLocusAllelesFasta)
    print 'cp to new dir'
    print cmd
    runCMD(cmd)
    
    # indes for bwa mem
    cmd = myData['bwa'] + ' index ' + newLocusAllelesFasta
    print cmd
    runCMD(cmd)

    # index for geting len    
    cmd = 'samtools faidx ' + newLocusAllelesFasta
    print cmd
    runCMD(cmd)
    
    # get seq lens and figure out what the are
    tmp = []
    inFile = open(newLocusAllelesFasta + '.fai','r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        n = line[0]
        l = int(line[1])
        tmp.append([n,l])
    inFile.close()
    if len(tmp) != 2:
        print 'len not 2 -- something wrong!'
        print tmp
        sys.exit()
    if tmp[0][1] < tmp[1][1]:
        emptySeqName = tmp[0][0]
        insSeqName = tmp[1][0]
    else:
        emptySeqName = tmp[1][0]
        insSeqName = tmp[0][0]
    print 'Empty seq: %s  insertion seq %s' % (emptySeqName,insSeqName)
    emptSeqFA = newLocusAllelesFasta + '.empty.fa'
    insSeqFA = newLocusAllelesFasta + '.ins.fa'
    cmd = 'samtools faidx %s %s > %s' % (newLocusAllelesFasta,emptySeqName,emptSeqFA)
    print cmd
    runCMD(cmd)
    cmd = 'samtools faidx %s %s > %s' % (newLocusAllelesFasta,insSeqName,insSeqFA)
    print cmd
    runCMD(cmd)
    
    # run age_align to get the brkpnt and TSD info
    ageOutName = newLocusAllelesFasta + '.ageout'
    cmd = 'age_align -indel %s %s > %s' % (emptSeqFA,insSeqFA,ageOutName)
    print cmd
    # because returns 1 not 0 on completion
    runCMD_noFail(cmd)
    
    res = age_parse.process_align_file(ageOutName)
    
    # get the info we want.
    s1Name = res['seq1Name']
    s1BrkpntLeftEnd = res['align']['s1Left'][1]
    s1BrkpntRightStart = res['align']['s1Right'][0]
    
    print s1Name,s1BrkpntLeftEnd,s1BrkpntRightStart
    
    s2Name = res['seq2Name']
    s2BrkpntLeftEnd = res['align']['s2Left'][1]
    s2BrkpntRightStart = res['align']['s2Right'][0]

    print s2Name,s2BrkpntLeftEnd,s2BrkpntRightStart
    
    # id at brkpnts -- is te TSD
    tsdLen = res['IDat']['s2Len']
    print 'tsdlen',tsdLen
    emptyTSDStart = s1BrkpntRightStart
    emptyTSDEnd = s1BrkpntRightStart + tsdLen - 1
    
    insLeftTSDStart = s2BrkpntLeftEnd + 1
    insLeftTSDEnd = insLeftTSDStart + tsdLen - 1
    
    insRightTSDStart = s2BrkpntRightStart
    insRightTSDEnd =  insRightTSDStart + tsdLen - 1
    
    if tsdLen == 0:
        emptyTSDEnd= emptyTSDStart
        insLeftTSDEnd = insLeftTSDStart
        insRightTSDEnd = insRightTSDStart
        
    
    # info to print out
    nl = [siteID,tsdLen,s1Name,emptyTSDStart,emptyTSDEnd,s2Name,insLeftTSDStart,insLeftTSDEnd,insRightTSDStart,insRightTSDEnd,siteChrom,sitePos]
    print nl
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    myData['brkpntOutFile'].write(nl)
    
#####################################################################
def read_brkpntFile(myData):
    myData['brkpntIntervals'] = []
    inFile = open(myData['brkpntFile'],'r')
    for line in inFile:
        if line[0] == '#':
            continue
        line = line.rstrip()
        line = line.split()
        for i in [1,3,4,6,7,8,9]:
            line[i] = int(line[i])
        myData['brkpntIntervals'].append(line)
    inFile.close()
#####################################################################
def align_to_alts_split_mem(myData,siteData):
    siteData['outSAM'] = siteData['mappingOutDir'] + 'mapped.sam'
    siteData['targetFA'] = myData['newLocusDir'] + siteData['siteID'] + '/alleles.fa'    
    # run bwa mem, with defaults...
    cmd = myData['bwa'] + ' mem  -M  ' +  siteData['targetFA'] + ' ' + siteData['fq1'] + ' ' + siteData['fq2'] + ' > ' + siteData['outSAM']
    print cmd
    runCMD(cmd)        
    select_hits_from_sam_split_mem(siteData,myData)
#####################################################################
def get_cigar_counts(cig_string):    
    cigCounts = {}
    cigCounts['M'] = 0
    cigCounts['D'] = 0
    cigCounts['I'] = 0
    cigCounts['S'] = 0
    cigCounts['H'] = 0                    
    cigExpand = expand_cigar(cig_string)
    for i in cigExpand:
        cigCounts[i[1]] += i[0]        
    return cigCounts
#####################################################################



def select_hits_from_sam_split_mem(siteData,myData):

    siteData['outSAMmatch'] = siteData['outSAM'] + '.matchinfo'
    outMatch = open(siteData['outSAMmatch'],'w')
    # setup split read dict
    print siteData['siteInterval']
    print myData['tsdExtendWinLen']
    
    
    
    splitDict = {}
    splitDict['empty'] = {}
    for i in range( siteData['siteInterval'][3]-myData['tsdExtendWinLen'],siteData['siteInterval'][4] + myData['tsdExtendWinLen'] + 1):
        splitDict['empty'][(siteData['siteInterval'][2],i)] = 1

    splitDict['insLeft'] = {}
    for i in range( siteData['siteInterval'][6]-myData['tsdExtendWinLen'],siteData['siteInterval'][7] + myData['tsdExtendWinLen'] + 1):
        splitDict['insLeft'][(siteData['siteInterval'][5],i)] = 1

    splitDict['insRight'] = {}
    for i in range( siteData['siteInterval'][8]-myData['tsdExtendWinLen'],siteData['siteInterval'][9] + myData['tsdExtendWinLen'] + 1):
        splitDict['insRight'][(siteData['siteInterval'][5],i)] = 1

    print 'Lens for each split dict'
    print len(splitDict['empty'])
    print len(splitDict['insLeft'])
    print len(splitDict['insRight'])
    
    
    readMatchEmpty = 0
    readMatchLeft = 0
    readMatchRight = 0
    
    
    
    inBamFile = pysam.AlignmentFile(siteData['outSAM'],'r')
    # go through read by read..
    for read in  inBamFile:
        if read.is_secondary is True:
            continue
        if read.is_supplementary is True:
            continue
        if read.is_unmapped is True:
            continue
        if read.mapping_quality == 0:
            continue

        if read.mapping_quality < 20:  # be sure is really a better match to that allele...
            continue


       # remove cig hitsN
        cigCounts = get_cigar_counts(read.cigarstring)
        if cigCounts['S'] >= 10:
            print 'Skipping -- has at least 10 soft clips, so I question the alignment'
            continue
        chromName = read.reference_name
        countEmpty = 0
        countInsLeft = 0
        countInsRight = 0
        
        for i in read.get_reference_positions():
            p = i + 1
            if (chromName,p) in splitDict['empty']:
                countEmpty += 1        
            if (chromName,p) in splitDict['insLeft']:
                countInsLeft += 1        
            if (chromName,p) in splitDict['insRight']:
                countInsRight += 1        
                
        print 'Counts are:', countEmpty,countInsLeft,countInsRight
        # output ones non zero
        if max(countEmpty,countInsLeft,countInsRight) > 0:
            nl = [read.query_name,read.reference_name,read.cigarstring,countEmpty,countInsLeft,countInsRight]
            nl = [str(j) for j in nl]
            nl = '\t'.join(nl) + '\n'
            outMatch.write(nl)
        
        
        if countEmpty == len(splitDict['empty']):
            readMatchEmpty += 1
        if countInsLeft == len(splitDict['insLeft']):
            readMatchLeft += 1
        if countInsRight == len(splitDict['insRight']):
            readMatchRight += 1
    inBamFile.close()
    outMatch.close()
    # done with site
    print 'Done with site -- here are the counts for reads across each junction'
    print readMatchEmpty,readMatchLeft,readMatchRight
    nl = [siteData['siteID'],myData['sampleName'],readMatchEmpty,readMatchLeft,readMatchRight]
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    myData['splitSummaryFile'].write(nl)
#####################################################################
