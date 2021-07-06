import sys
import os.path
import brkptgen
from optparse import OptionParser


###############################################################################
USAGE = """
python setup-for-splitread.py   --allelefile <allele file to process in proper format>
                           --brkpntfile <file for new breakpoint info>
                           --orig_locusdir <original directory of locus alleles files>
                           --new_locusdir <new directory of lcous alleles files>
                           --bwamem <path/cmd for bwa mem [i.e, v 0.7.15 or newer] >



"""
parser = OptionParser(USAGE)
parser.add_option('--allelefile',dest='alleleFile', help = 'file name of allele seqs to process')
parser.add_option('--brkpntfile',dest='brkpntFile', help = 'file name of allele seqs to process')

parser.add_option('--orig_locusdir',dest='origLocusDir', help = 'base dir name for original locus dir')
parser.add_option('--new_locusdir',dest='newLocusDir', help = 'base dir name for new locus dir')
parser.add_option('--bwamem',dest='bwaMem', help = 'path/cmd  bwa mem [i.e, v 0.7.15 or newer]')


(options, args) = parser.parse_args()

if options.alleleFile is None:
    parser.error('alleleFile name not given')
if options.brkpntFile is None:
    parser.error('brkpntFile name not given')
if options.origLocusDir is None:
    parser.error('origLocusDir name not given')
if options.newLocusDir is None:
    parser.error('newLocusDir name not given')
if options.bwaMem is None:
    parser.error('bwaMem cmd given')
###############################################################################


print 'Hello!'


#setup some info in dictionary
myData = {}
myData['bwa'] = options.bwaMem

myData['origLocusDir'] = options.origLocusDir
if myData['origLocusDir'][-1] != '/' :
    myData['origLocusDir'] += '/'

myData['newLocusDir'] = options.newLocusDir
if myData['newLocusDir'][-1] != '/' :
    myData['newLocusDir'] += '/'



##### check to see if dirs exist #####
if os.path.isdir(myData['origLocusDir']) is False:
    print 'Erro! Dir does not exist!'
    print myData['origLocusDir']
    sys.exit()

if os.path.isdir(myData['newLocusDir']) is False:
    print 'Erro! Dir does not exist!'
    print myData['newLocusDir']
    sys.exit()

###############################################################################
#### check paths for programs that are needed

print 'Checking bwa...'
if brkptgen.which(myData['bwa']) is None:
	print 'ERROR bwa not found in path! please fix (module add?)'
	sys.exit()

print 'Checking samtools...'
if brkptgen.which('samtools') is None:
	print 'ERROR samtools not found in path! please fix (module add?)'
	sys.exit()

print 'Checking age_align...'
if brkptgen.which('age_align') is None:
	print 'ERROR age_align not found in path! please fix (module add?)'
	print 'ERROR need program age for alignment! check path!'
	sys.exit()

###############################################################################




# read in site intervals
myData['siteIntervals'] = brkptgen.get_site_intervals_from_table(options.alleleFile)
print 'Found %i siteIntervals' % len(myData['siteIntervals'])

myData['brkpntOutFile'] = open(options.brkpntFile,'w')
nl = ['#siteID','tsdLen','emptyName','emptyTSDStart','emptyTSDEnd','insName','insLeftTSDStart','insLeftTSDEnd','insRightTSDStart','insRightTSDEnd','chrom','pos']
nl = '\t'.join(nl) + '\n'
myData['brkpntOutFile'].write(nl)

for siteInterval in myData['siteIntervals']:
    siteID = siteInterval[0]
    print siteInterval
    brkptgen.setup_locus_for_split(myData,siteID,siteInterval[1][0],siteInterval[1][3])


myData['brkpntOutFile'].close()


