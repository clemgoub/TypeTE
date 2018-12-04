#! /bin/bash

###########################################
# reGenotypeTE - parameterfile.init       #
#                                         #
# configuration file for the pipeline     #
#                                         #
# Author: Clement Goubert                 #
# Date: 04/2018                           #
# Version: 1.0                            #
###########################################

### MAIN PARAMETERS
# input
VCF="" #Path to MELT vcf (.vcf or .vcf.gz)
BAMPATH="" # Path to the bams' folder
BAMFILE="" # indiv_name bam_name (table 2 fields)

# output
PROJECT="" # Name of the project
OUTDIR="" # Path to place the output directory (will be named after PROJECT); OUTDIR must exist

# options
individual_nb="" # number of individual per job
CPU="" # number of CPU
GENOME="" # Paht the the reference genome sequence
SCRIPT_PATH=$(pwd)
MAP="YES" #OR NO

# DEPENDENCIES PATH