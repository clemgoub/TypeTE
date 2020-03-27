import sys
import os.path
import brkptgen
from optparse import OptionParser

###############################################################################
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
        c = c.replace('chr','') #for 1kg
        b = int(line[4])
        e = int(line[5])
        calls.append([id,[c,b,e]])
    inFile.close()
    
    
    # now I want to sort the list
    calls.sort(key = lambda x: x[1][1])
#    calls.sort(key = lambda x: genutils.chrom_to_plink(x[1][0]))
    calls.sort(key = lambda x: x[1][0])
    return calls
###############################################################################


###############################################################################
USAGE = """
python process-sample-remap.py   --allelefile <allele file in proper format>
                           --allelebase <base directory for allele information>
                           --excludefile <file of coorindates to exclude contained mapping (1 based)>
                           --bwa <path/cmd for bwa 0.5.9
                           --samplename <name of sample>

Version assumes reads extracted already from the BAM, only does that mapping steps


"""
parser = OptionParser(USAGE)
parser.add_option('--allelefile',dest='alleleFile', help = 'file name of allele seqs to process')
parser.add_option('--allelebase',dest='alleleBase', help = 'base dir name for outputs')
parser.add_option('--bwa',dest='bwa', help = 'path/cmd for bwa 0.7.15')
parser.add_option('--samplename',dest='sampleName', help = 'sample name')
parser.add_option('--excludefile',dest='excludeFile', help = 'file of coordinates to exclude contained mappings')

(options, args) = parser.parse_args()

if options.alleleFile is None:
    parser.error('alleleFile name not given')
if options.alleleBase is None:
    parser.error('alleleBase name not given')
if options.bwa is None:
    parser.error('bwa cmd given')
if options.sampleName is None:
    parser.error('sample name not given')


###############################################################################
# this is messy and redundant because it is a merge of separate programs into
# one step...

#setup some info in dictionary
myData = {}
myData['bwa'] = options.bwa
myData['sampleName'] = options.sampleName
myData['alleleBase'] = options.alleleBase
if myData['alleleBase'][-1] != '/' :
    myData['alleleBase'] += '/'


myData['siteIntervals'] = brkptgen.get_site_intervals_from_table(options.alleleFile)

# for testing, only do first 10....
#myData['siteIntervals'] = myData['siteIntervals'][0:10]  

print 'Found %i siteIntervals' % len(myData['siteIntervals'])

myData['excludeFileName'] = options.excludeFile
brkptgen.setup_exclusion(myData)

samplesBase = myData['alleleBase'] + 'samples'
if os.path.isdir(samplesBase) is False:
    print 'Making samples base',samplesBase
    cmd = 'mkdir ' + samplesBase
    brkptgen.runCMD(cmd)
else:
    print samplesBase,'exists!'


myData['sampleBase'] = samplesBase + '/' + myData['sampleName']

if os.path.isdir(myData['sampleBase']) is False:
    print 'Making samples base',myData['sampleBase']
    cmd = 'mkdir ' + myData['sampleBase']
    brkptgen.runCMD(cmd)
else:
    print myData['sampleBase'],'exists!'

myData['sampleBase'] += '/'
myData['window_size'] = 0
myData['min_map_q'] = 20

myData['outPutReadFile'] = myData['sampleBase'] + '%s.intervalreads' %  myData['sampleName']
myData['outputSAMFile'] = myData['sampleBase'] + '%s.intervalreads.sam' %  myData['sampleName']


# will just add readlen to the starts to get intervals, this will be off by 1, but will match
# the output of RetroSeq

print 'Have %i intervals to process' % len(myData['siteIntervals'])

# at this point, we have extracted all of the reads for the sample!
print 'Assuming have extracted all reads for the sample!'

print 'Begin map to alt alleles steps!'

myData['outPutReadFile'] = myData['sampleBase'] + '%s.intervalreads' %  myData['sampleName']
myData['seqFileName'] = myData['outPutReadFile'] + '.seq'

print 'Have %i intervals to process' % len(myData['siteIntervals'])

if os.path.isfile(myData['seqFileName']) is False:
    brkptgen.make_seq_reads_file(myData['outputSAMFile'],myData['seqFileName'])
    print 'made seq file',myData['seqFileName']
else:
    print myData['seqFileName'],'seems to already exist'


# make dictionary matching sequence reads to intervals
myData['intervalsToReads'] = brkptgen.match_intervals_to_reads(myData['outPutReadFile'])
print 'Read in matching info for %i intervals' % len(myData['intervalsToReads'])    

myData['nameToSeq'] = brkptgen.make_name_to_seq_dictionary(myData['seqFileName'])
print 'read in %i sequences' % len(myData['nameToSeq'])

myData['mappingDirBase'] = myData['sampleBase'] + 'mapping'
if os.path.isdir(myData['mappingDirBase']) is False:
    print 'Making samples base',myData['mappingDirBase']
    cmd = 'mkdir ' + myData['mappingDirBase']
    brkptgen.runCMD(cmd)
else:
    print myData['mappingDirBase'],'exists!'

myData['mappingDirBase'] += '/'


# now we have to do the mapping for each region
for siteInterval in myData['siteIntervals']:
    siteID = siteInterval[0]
    print siteInterval
    siteData = {}
    siteData['siteID'] =  siteID
    siteData['mappingOutDir'] = myData['mappingDirBase'] + siteID
    if os.path.isdir(siteData['mappingOutDir']) is False:
        cmd = 'mkdir ' + siteData['mappingOutDir']
        brkptgen.runCMD(cmd)
    siteData['mappingOutDir'] += '/'
    
        
    # make the two fastq files of reads
    siteData['fq1'] = siteData['mappingOutDir'] + 'read1.fq'
    siteData['fq2'] = siteData['mappingOutDir'] + 'read2.fq'

    brkptgen.write_fastq_for_site(myData,siteData)
    brkptgen.align_to_alts_bwa(myData,siteData)    

print '*****\nAll Mapping Done\n******\n'
print 'Calculating genotype likelihoods\n'


myData['lhSummaryFileName'] =  myData['sampleBase'] + myData['sampleName'] + '.lksummary'


lhSummaryOut = open(myData['lhSummaryFileName'],'w')
print 'Summary output file name is',myData['lhSummaryFileName']
lhSummaryOut.write('siteID\tsampleName\ttotFrags\tnumRef\tnumAlt\tgl_00\tgl_01\tgl_11\tscaledGLs\n')



for siteInterval in myData['siteIntervals']:
    siteID = siteInterval[0]
    siteData = {}
    siteData['siteID'] =  siteID
    siteData['mappingOutDir'] = myData['mappingDirBase'] + siteID + '/'    
    siteData['outSAM'] = siteData['mappingOutDir'] + 'mapped.sam'
    siteData['outSamFilter'] = siteData['outSAM'] + '.filter'
    siteData['outSamSel'] = siteData['outSAM'] + '.filter.sel.sam'



    brkptgen.read_samsel_hits(siteData)
    brkptgen.calc_gen_likelihood(siteData)
        

    nl = [siteData['siteID']]
    nl.append(myData['sampleName'])
    nl.append(siteData['totFrags'])
    nl.append(siteData['numRefFrag'])
    nl.append(siteData['numAltFrag'])
    nl.append(siteData['gl_0'])
    nl.append(siteData['gl_1'])
    nl.append(siteData['gl_2'])
    sl = '%i,%i,%i' % (siteData['scaledLikelihoods'][0],siteData['scaledLikelihoods'][1],siteData['scaledLikelihoods'][2])
    nl.append(sl)
    
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    lhSummaryOut.write(nl)
    
    print nl

    
lhSummaryOut.close()




# make VCF
myData['lhSummaryFileName'] =  myData['sampleBase'] + myData['sampleName'] + '.lksummary'


myData['outVCFName'] = myData['sampleBase'] + myData['sampleName'] + '.vcf'
outVCF = open(myData['outVCFName'],'w')
outVCF.write('##fileformat=VCFv4.1\n')
outVCF.write('##fileDate=20170807\n')
outVCF.write('##reference=canFam3.1\n')
outVCF.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
outVCF.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">\n')
outVCF.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
outVCF.write('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">\n')
outVCF.write('##FORMAT=<ID=PL,Number=3,Type=Float,Description="Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt">\n')

outVCF.write('##ALT=<ID=INS:ME:CfERV,Description="Insertion of HERVK element">\n')
outVCF.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Filtered Depth">\n')
header = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']

header.append(myData['sampleName'])
nl = '\t'.join(header) + '\n'
outVCF.write(nl)

# get if ref or alter is smaller....

alleleSeqs = {}
inFile = open(options.alleleFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    if line[0] == 'siteID':
        continue
    siteID = line[0]
    ref = line[6]
    alt = line[7]
    if len(ref) < len(alt):
         alleleSeqs[siteID] = 'refEmpty'     
    else:
         alleleSeqs[siteID] = 'refIns'     
inFile.close()

# set siteID to pos....
siteIDToPos = {}
for siteInterval in myData['siteIntervals']:
    siteID = siteInterval[0]
    c = siteInterval[1][0]
    pos = str(siteInterval[1][3])
    siteIDToPos[siteID] = [c,pos]

inFile = open(myData['lhSummaryFileName'],'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    if line[0] == 'siteID':
        continue
    
    siteID = line[0]
    # put in 'N' for ref... not right, but ok for now?
    # also <MEI> isn't exaclt right
    nl = [siteIDToPos[siteID][0],siteIDToPos[siteID][1],siteID,'N','<MEI>','.','.']
    totDepth = int(line[2])
    i = 'DP=%i' % totDepth
    i += ';' + alleleSeqs[siteID]
    nl.append(i)   
    f = 'GT:AD:GQ:PL'
    nl.append(f)
    pl = line[8]
    naiveGen = '?'
    GQ = '?'
    gLikeList = [float(line[5]),float(line[6]),float(line[7])]
    plParts = pl.split(',')
    if pl == '0,0,0' :
        naiveGen  = '.'
        GQ = 0.0
    elif plParts[0] == '0':
        naiveGen  = '0/0'
        GQ = brkptgen.calc_gq(gLikeList,0)		
    elif plParts[1] == '0':
        naiveGen  = '0/1'
        GQ = brkptgen.calc_gq(gLikeList,1)
    elif plParts[2] == '0':
        naiveGen  = '1/1'
        GQ = brkptgen.calc_gq(gLikeList,2) 
    gen = naiveGen + ':' + line[3] + ',' + line[4]  + ':' + '%.2f' % GQ + ':' + pl
    nl.append(gen)
    nl = '\t'.join(nl) + '\n'
    outVCF.write(nl)
inFile.close()
outVCF.close()
print 'DONE!'
print 'VCF written to',myData['outVCFName']








