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
python process-sample.py   --allelefile <allele file in proper format>
                           --allelebase <base directory for allele information>
                           --bwa <path/cmd for bwa 0.7.15
                           --samplename <name of sample>
                           --bam <indexed bam of sample>


"""
parser = OptionParser(USAGE)
parser.add_option('--allelefile',dest='alleleFile', help = 'file name of allele seqs to process')
parser.add_option('--allelebase',dest='alleleBase', help = 'base dir name for outputs')
parser.add_option('--bwa',dest='bwa', help = 'path/cmd for bwa 0.7.15')
parser.add_option('--samplename',dest='sampleName', help = 'sample name')
parser.add_option('--bam',dest='bam', help = 'bam of sample')


(options, args) = parser.parse_args()

if options.alleleFile is None:
    parser.error('alleleFile name not given')
if options.alleleBase is None:
    parser.error('alleleBase name not given')
if options.bwa is None:
    parser.error('bwa cmd given')
if options.sampleName is None:
    parser.error('sample name not given')
if options.bam is None:
    parser.error('bam  not  given')



###############################################################################
# this is messy and redundant because it is a merge of separate programs into
# one step...

#setup some info in dictionary
myData = {}
myData['bwa'] = options.bwa
myData['sampleName'] = options.sampleName
myData['bam'] = options.bam
myData['alleleBase'] = options.alleleBase
if myData['alleleBase'][-1] != '/' :
    myData['alleleBase'] += '/'


myData['siteIntervals'] = brkptgen.get_site_intervals_from_table(options.alleleFile)
print 'Found %i siteIntervals' % len(myData['siteIntervals'])


print 'Will do sample name %s' % myData['sampleName']
print 'Bam file: %s' % myData['bam']
if os.path.isfile(myData['bam']) is False:
    print 'ERROR! BAM FILE Does not exist'
    print myData['bam']
    sys.exit()



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


myData['siteIntervals'] = myData['siteIntervals'][0:10]

print 'Have %i intervals to process' % len(myData['siteIntervals'])

# here, go through each site interval and get all reads that overlap....
outFile = open(myData['outPutReadFile'],'w')
outSAM = open(myData['outputSAMFile'],'w')
rN = 0
for call in myData['siteIntervals']:
    samLines = {}
    rpToDo = {}
    rpToFind = {}
    
    rN += 1
    if rN % 20 == 0:
        print '\n*****\nDoing %i of %i ...\n*****\n' % (rN,len(myData['siteIntervals']))
    id = call[0]
    c = call[1][0]
    b = call[1][1]
    e = call[1][2]
    reg = '%s:%i-%i' % (c,b,e)
    print reg
    numReads = 0
    numSC = 0
    bamIn = brkptgen.open_bam_read(myData['bam'],reg)
    for line in bamIn:
        numReads += 1
        ol = line
        line = line.rstrip()
        line = line.split()
        samParse = brkptgen.parse_sam_line(line)
        if samParse['isFirst'] is True:
            rNum = '1'
        else:
            rNum = '2'
        samLines[(samParse['seqName'],rNum)] = ol

        #get rid of things that are not good
        if samParse['unMapped'] is True:
            continue
        if samParse['isDuplicate'] is True:
            continue
        if samParse['notPrimaryAlignment'] is True:
            continue
        if samParse['mapQ'] < myData['min_map_q']:
            continue
        if samParse['isPaired'] is False:
            continue

        numSC += 1
        outFile.write('%s\t%s\t%i\t%i\t%s\n' % (id,samParse['seqName'],samParse['flag'],samParse['mapQ'],samParse['cigar']))
        rpToDo[samParse['seqName']] = 1
    bamIn.close()
    print 'Reads considered: %i, found: %i' % (numReads,numSC)
    print 'Found %i rp to do' % len(rpToDo)
    print 'initial sam dict contains %i reads' % len(samLines)

    #now want to go through list and print out sam lines for ones that we need
    numWritten = 0
    for rn in rpToDo:
        for rNum in ['1','2']:
            if (rn,rNum) in samLines:
                samLine = samLines[(rn,rNum)]
                outSAM.write(samLine)
                numWritten += 1
            else:
                rpToFind[(rn,rNum)] = 1
    print 'After first pass, we have %i reads to find' % len(rpToFind)
    for (rn,rNum) in rpToFind:
        # check since we had updated
        if (rn,rNum) in samLines:
            samLine = samLines[(rn,rNum)]
            outSAM.write(samLine)
            numWritten += 1
            continue
        
        if rNum == '1':
            otherNum = '2'
        elif rNum == '2':
            otherNum = '1'
            
        samLine = samLines[(rn,otherNum)]
        samLine = samLine.rstrip()
        samLine = samLine.split()
        samParse = brkptgen.parse_sam_line(samLine)
        #genutils.print_dictionary_keys(samParse)
        if samParse['mateChrom'] == '=':
            newReg = samParse['chrom'] + ':' + samParse['matePos'] + '-' + samParse['matePos']
        else:
            newReg = samParse['mateChrom'] + ':' + samParse['matePos'] + '-' + samParse['matePos']
         
        if samParse['mateUnmapped'] is True:
            print 'Mate is unmapped, not sure what to do, trying mapped region'
#            print samLine
            newReg = samParse['chrom'] + ':' + str(samParse['chromPos']) + '-' + str(samParse['chromPos'])


        bamIn = brkptgen.open_bam_read(myData['bam'],newReg)
        for line in bamIn:
            numReads += 1
            ol = line
            line = line.rstrip()
            line = line.split()
            samParse = brkptgen.parse_sam_line(line)
            if samParse['isFirst'] is True:
                r = '1'
            else:
                r = '2'
            samLines[(samParse['seqName'],r)] = ol
        bamIn.close()
        
        #now do the check
        if (rn,rNum) in samLines:
            samLine = samLines[(rn,rNum)]
            outSAM.write(samLine)
            numWritten += 1
        else:
            print 'cannot find!',rn,rNum
            sys.exit()

    print 'Wrote out',numWritten
outFile.close()
outSAM.close()


# at this point, we have extracted all of the reads for the sample!
print 'Have extracted all reads for the sample!'

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
    siteData['outSamSel'] = siteData['outSamFilter'] + '.sel'

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












