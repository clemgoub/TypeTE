import sys
import os
import brkptgen

from optparse import OptionParser
###############################################################################
USAGE = """
python create-alternative-alleles.py  --allelefile <allele file in proper format>
                                      --allelebase <base directory for allele information>
                                      --bwa <path/cmd for bwa 0.7.15


"""
parser = OptionParser(USAGE)
parser.add_option('--allelefile',dest='alleleFile', help = 'file name of allele seqs to process')
parser.add_option('--allelebase',dest='alleleBase', help = 'base dir name for outputs')
parser.add_option('--bwa',dest='bwa', help = 'path/cmd for bwa 0.7.15')


(options, args) = parser.parse_args()

if options.alleleFile is None:
    parser.error('alleleFile name not given')
if options.alleleBase is None:
    parser.error('alleleBase name not given')
if options.bwa is None:
    parser.error('bwa cmd given')
###############################################################################



alleleSeqs = {}
inFile = open(options.alleleFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()

    siteID = line[0]
    if siteID == 'siteID':
        continue
    print siteID,len(line)
    alleleSeqs[siteID] = {}
    alleleSeqs[siteID]['ref'] = line[6]
    alleleSeqs[siteID]['alt'] = line[7]
inFile.close()

print 'Read in sequences for %i sites from %s' % (len(alleleSeqs),options.alleleFile)


#check to see if options.alleleBase exists
if os.path.isdir(options.alleleBase) is True:
    print '%s exists, continuing!' % options.alleleBase
    options.alleleBase += '/'
else:
    print '%s does not exist, exiting.  Please fix!' % options.alleleBase
    sys.exit()
    

locusSeqDir = options.alleleBase  + 'locusAlleles'
cmd = 'mkdir -p ' + locusSeqDir
print cmd
brkptgen.runCMD(cmd)


for siteID in alleleSeqs:
    alleleDir = locusSeqDir + '/' + siteID
    cmd = 'mkdir -p ' + alleleDir
    brkptgen.runCMD(cmd)
    alleleDir += '/'
    alleleFa = alleleDir + 'alleles.fa'
    outFile = open(alleleFa,'w')
    outFile.write('>%s\n' % (siteID+'_genome'))
    seq = alleleSeqs[siteID]['ref']
    seq = brkptgen.add_breaks_to_line(seq)
    seq += '\n'
    outFile.write(seq)

    outFile.write('>%s\n' % (siteID+'_alternative'))
    seq = alleleSeqs[siteID]['alt']
    seq = brkptgen.add_breaks_to_line(seq)
    seq += '\n'    
    outFile.write(seq)
    outFile.close()
    cmd = '%s index %s' % (options.bwa, alleleFa)
    print cmd
    brkptgen.runCMD(cmd)
    
