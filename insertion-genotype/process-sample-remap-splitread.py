import sys
import os.path
import brkptgen
from optparse import OptionParser


###############################################################################
USAGE = """
python process-sample-remap-splitread.py
                           --brkpntfile <file of new breakpoint info>
                           --new_locus_dir <new directory of locus alleles files>
                           --samplename <name of sample>
                           --old_sample_gen_dir <directory of old sample genotype, contains *intervalreads, sam, seq>
                           --new_gen_base_dir <bas directory for new genotypes>
                           --bwamem <path/cmd for bwa mem [i.e, v 0.7.15 or newer] >





"""
parser = OptionParser(USAGE)
parser.add_option('--brkpntfile',dest='brkpntFile', help = 'file name of allele seqs to process')
parser.add_option('--new_locus_dir',dest='newLocusDir', help = 'base dir name for new locus dir')
parser.add_option('--samplename',dest='sampleName', help = 'sample name')
parser.add_option('--old_sample_gen_dir',dest='oldSampleGenDir', help = 'base dir name for old original sample gen dir')
parser.add_option('--new_gen_base_dir',dest='newGenBaseDir', help = 'base dir name for new genotypes')
parser.add_option('--bwamem',dest='bwaMem', help = 'path/cmd  bwa mem [i.e, v 0.7.15 or newer]')


(options, args) = parser.parse_args()

if options.brkpntFile is None:
    parser.error('brkpntFile name not given')
if options.newLocusDir is None:
    parser.error('newLocusDir name not given')
if options.sampleName is None:
    parser.error('sampleName  not given')
if options.oldSampleGenDir is None:
    parser.error('oldSampleGenDir not given')
if options.newGenBaseDir is None:
    parser.error('newGenBaseDir not given')
if options.bwaMem is None:
    parser.error('bwaMem path not given')
###############################################################################

#setup some info in dictionary
myData = {}
myData['bwa'] = options.bwaMem
myData['sampleName'] = options.sampleName
myData['brkpntFile'] = options.brkpntFile 
 
myData['tsdExtendWinLen'] = 5  # go +/- 5 bp on each side of each TSD, require there to be align to each
 
myData['newLocusDir'] = options.newLocusDir
if myData['newLocusDir'][-1] != '/' :
    myData['newLocusDir'] += '/'

myData['newGenBaseDir'] = options.newGenBaseDir
if myData['newGenBaseDir'][-1] != '/' :
    myData['newGenBaseDir'] += '/'

myData['oldSampleGenDir'] = options.oldSampleGenDir
if myData['oldSampleGenDir'][-1] != '/' :
    myData['oldSampleGenDir'] += '/'



##### check to see if dirs exist #####
if os.path.isdir(myData['newLocusDir']) is False:
    print 'Erro! Dir does not exist!'
    print myData['newLocusDir']
    sys.exit()

if os.path.isdir(myData['newGenBaseDir']) is False:
    print 'Erro! Dir does not exist!'
    print myData['newGenBaseDir']
    sys.exit()

if os.path.isdir(myData['oldSampleGenDir']) is False:
    print 'Erro! Dir does not exist!'
    print myData['oldSampleGenDir']
    sys.exit()
###############################################################################

samplesBase = myData['newGenBaseDir'] + '/' + myData['sampleName']
if os.path.isdir(samplesBase) is False:
    print 'Making samples base',samplesBase
    cmd = 'mkdir ' + samplesBase
    brkptgen.runCMD(cmd)
else:
    print samplesBase,'exists!'
myData['samplesBase'] = samplesBase + '/'

myData['mappingDirBase'] = myData['samplesBase'] + 'mapping'
if os.path.isdir(myData['mappingDirBase']) is False:
    print 'Making samples base',myData['mappingDirBase']
    cmd = 'mkdir ' + myData['mappingDirBase']
    brkptgen.runCMD(cmd)
else:
    print myData['mappingDirBase'],'exists!'

myData['mappingDirBase'] += '/'





# read data from initial process-sample
myData['outPutReadFile'] = myData['oldSampleGenDir'] + '%s.intervalreads' %  myData['sampleName']
myData['outputSAMFile'] = myData['oldSampleGenDir'] + '%s.intervalreads.sam' %  myData['sampleName']
myData['seqFileName'] = myData['outPutReadFile'] + '.seq'

if os.path.isfile(myData['outPutReadFile']) is False:
    print 'ERROR: file %s does not exist' % myData['outPutReadFile']
    sys.exit()

if os.path.isfile(myData['outputSAMFile']) is False:
    print 'ERROR: file %s does not exist' % myData['outputSAMFile']
    sys.exit()

if os.path.isfile(myData['seqFileName']) is False:
    brkptgen.make_seq_reads_file(myData['outputSAMFile'],myData['seqFileName'])
    print 'made seq file',myData['seqFileName']
else:
    print myData['seqFileName'],'seems to already exist'



###############################################################################
# Ready to start processing

brkptgen.read_brkpntFile(myData)
print 'Have %i intervals to process' % len(myData['brkpntIntervals'])

# assume reads already extracted, so get them setup
# make dictionary matching sequence reads to intervals
myData['intervalsToReads'] = brkptgen.match_intervals_to_reads(myData['outPutReadFile'])
print 'Read in matching info for %i intervals' % len(myData['intervalsToReads'])    

myData['nameToSeq'] = brkptgen.make_name_to_seq_dictionary(myData['seqFileName'])
print 'read in %i sequences' % len(myData['nameToSeq'])

myData['splitSummaryFileName'] =  myData['samplesBase'] + myData['sampleName'] + '.brkpnt-match.table.txt'
myData['splitSummaryFile'] = open(myData['splitSummaryFileName'],'w')
print 'Split summary output table name is',myData['splitSummaryFileName']
myData['splitSummaryFile'].write('#siteID\tsampleName\tnumMatchEmpty\tnumMatchInsLeft\tnumMatchInsRight\n')


siteInfoDict = {}

for interval_i in range(len(myData['brkpntIntervals'])):
    siteInterval = myData['brkpntIntervals'][interval_i]
    print interval_i, siteInterval
        
    siteID = siteInterval[0]
    print siteInterval
    
    siteInfoDict[siteID] = {}
    siteInfoDict[siteID]['chrom'] = siteInterval[10]
    siteInfoDict[siteID]['pos'] = siteInterval[11]
    if '_genome' in siteInterval[2]:
        siteInfoDict[siteID]['refType'] = 'refEmpty'
    else:
        siteInfoDict[siteID]['refType'] = 'refIns'
    
    
    
    siteData = {}
    siteData['siteID'] =  siteID
    siteData['siteInterval'] = siteInterval
    siteData['mappingOutDir'] = myData['mappingDirBase'] + siteID
    if os.path.isdir(siteData['mappingOutDir']) is False:
        cmd = 'mkdir ' + siteData['mappingOutDir']
        brkptgen.runCMD(cmd)
    siteData['mappingOutDir'] += '/'
            
    # make the two fastq files of reads
    siteData['fq1'] = siteData['mappingOutDir'] + 'read1.fq'
    siteData['fq2'] = siteData['mappingOutDir'] + 'read2.fq'
    brkptgen.write_fastq_for_site(myData,siteData)
    brkptgen.align_to_alts_split_mem(myData,siteData)    

myData['splitSummaryFile'].close()    
    

# make VCF
myData['outVCFName'] = myData['samplesBase'] + myData['sampleName'] + '.vcf'
print 'Writting VCF to',myData['outVCFName']

outVCF = open(myData['outVCFName'],'w')
outVCF.write('##fileformat=VCFv4.1\n')
outVCF.write('##fileDate=20190519\n')
outVCF.write('##reference=canFam3.1\n')
outVCF.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the empty and ins">\n')
outVCF.write('##FORMAT=<ID=BD,Number=1,Type=Integer,Description="Read Depth for breakpoint regions empty, ins left, ins right">\n')
outVCF.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

outVCF.write('##ALT=<ID=INS:ME:SINEC,Description="Insertion of SINEC element">\n')
outVCF.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Filtered Depth">\n')
header = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
header.append(myData['sampleName'])
nl = '\t'.join(header) + '\n'
outVCF.write(nl)


inFile = open(myData['splitSummaryFileName'],'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    if line[0] == '#siteID':
        continue
    
    siteID = line[0]
    # put in 'N' for ref... not right, but ok for now?
    # also <MEI> isn't exaclt right
    nl = [siteInfoDict[siteID]['chrom'],siteInfoDict[siteID]['pos'],siteID,'N','<MEI>','.','.']
    totDepth = int(line[2]) + int(line[3]) + int(line[4])
    i = 'DP=%i' % totDepth
    i += ';' + siteInfoDict[siteID]['refType']
    nl.append(i)   
    f = 'GT:AD:BD'
    nl.append(f)
    
    me = int(line[2])
    ml = int(line[3])
    mr = int(line[4])
    mAlt = ml + mr
    
    if me == 0 and mAlt == 0:
       naiveGen = '.'
    elif me == mAlt:
       naiveGen = '0/1'
    elif me > mAlt:
       naiveGen = '0/0'
    else:
       naiveGen = '1/1'
    


    gen = naiveGen + ':' + '%s' % me + ',' + '%s' % mAlt + ':' + '%s,%s,%s' % (me,ml,mr)
    nl.append(gen)
    nl = '\t'.join(nl) + '\n'
    outVCF.write(nl)
outVCF.close()
print 'Wrote out VCF file !  All done'











    
    


