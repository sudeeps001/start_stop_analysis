#!/bin/bash
set -e

###################################################################
# Shell script calling samtools, bedtools and python script for start and stop codon coverage analysis.
# The pipeline used here follows part of the method adopted by RiboTaper workflow for start and stop codon analysis.
# RiboTaper workflow: https://ohlerlab.mdc-berlin.de/software/RiboTaper_126/ 
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#######################################################################

# this shell script looks for the python file in the same root folder 
# edit the lines below to change this behaviour
SRC=$(dirname "$(readlink -f "$0")")
pythonscript="start_stop_analysis.py"
if ! [[ -f $SRC/$pythonscript ]]; then
	echo -e "\nERROR! cannot find python helper file ${pythonscript} in script source folder ${SRC}\nPlease check you files...\n"
	exit 1
fi
# check for samtools
SAM=`which samtools`
if [ -z "$SAM" ]; then
	echo -e "\nERROR! samtools not found in '$PATH'.\nPlease install samtools...\n"
	exit 1
fi
#check for bedtools
BEDTOOLS=`which bedtools`
if [ -z "$BEDTOOLS" ]; then
	echo -e "\nERROR! bedtools not found in '$PATH'.\nPlease install bedtools...\n"
	exit 1
fi
# shell script name
bn=$(basename $0)
# usage and stuff
if [ $# -lt 3 ]; then
	echo -e "\n Usage: ${bn} <RPF.bam> <bedfile> <name> <downsample_size>"
	echo -e " RPF.bam: Ribosome profiling bam file\n bedfile: Bed file with start and stop co-ordinates"
	echo -e " name: Output files name prefix\n downsample_size (optional): If number of alignments is > downsample_size, bam file will be downsampled (default: 10,000,000).\n"
	echo -e " For a full description of the included python helper script, call 'python ${SRC}/${pythonscript} -h'\n"
	exit 1
fi
RPF=$1
bed=$2
basen=$3
downsamp=$4
if ! [[ -f "$RPF" ]]; then
    echo -e "\nERROR! <RPF.bam> file not found\n"
    exit 1
fi
if ! [[ -f "$bed" ]]; then
     echo -e "\nERROR! <start_stop_bed> file not found\n"
     exit 1
fi
if [ -z "$basen" ]; then
	echo -e "\nERROR! output file name prefix MUST be given...\n"
	exit 1
fi
if [ -z "$downsamp" ]; then
	echo -e "\nDownsample size not given, using default value: 10,000,000...\n"
	downsamp=10000000
fi
echo "Checking number of mapped reads..."
mapped=$(samtools flagstat $RPF|perl -ne 'print $1 if /^(\d+)\s{1,}\+\s{1,}\d+\s{1,}mapped/')
if (( $downsamp < $mapped )); then
	echo "Number of aligned reads ${mapped} > ${downsamp}, downsampling to ${downsamp}..."
	cat <(samtools view -H $RPF ) <(samtools view -S $RPF | shuf -n ${downsamp}) | samtools view - -bS > ${basen}_samp.bam
	echo "Intersecting alignments with start/stop sites..."
	bamToBed -i ${basen}_samp.bam -bed12 -split | windowBed -w 100 -sm -b stdin -a $bed | cut -f 7- | sort -k1,1 -k2,2g | closestBed -s -t "last" -a stdin -b $bed > ${basen}_analysis.csv
else
	echo -e "Number of aligned reads: ${mapped} <= ${downsamp}, using all the aligned reads...\nIntersecting alignments with start/stop sites..."
	bamToBed -i $RPF -bed12 -split | windowBed -w 100 -sm -b stdin -a $bed | cut -f 7- | sort -k1,1 -k2,2g | closestBed -s -t "last" -a stdin -b $bed > ${basen}_analysis.csv
fi
if !  [[ -s ${basen}_analysis.csv ]]; then
    echo -e "\nERROR! no intersections found, check your input files\n"
    exit 1
fi
python ${SRC}/${pythonscript} --f ${basen}_analysis.csv --name $basen --c 20 --plot html
# some house keeping
gzip ${basen}_analysis.csv
echo -e "\n***Done***\n"
