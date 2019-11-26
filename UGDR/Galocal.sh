#!/bin/bash
#####################################################################################
#####################################################################################
# -*- coding: utf8 -*-
# This file is part of UGDR-BedratAmina
# This script extract informations from VCF files to give a human readable
# files resuming alleles variations, and difining region of recombinaition
# Copyright (C) 2016 Bedrat Amina - Meiogenix - Institut Curie.
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#####################################################################################
#####################################################################################
HELPME()
{
echo ""
echo "Usage: $0 -a Fasta -b FQfolder -c RESFolder"
echo -e "\t-a S288 or other fasta sequence"
echo -e "\t-b Folder where the reads (R_1.fq and R_2.fq) are stocked"
echo -e "\t-c Where to put the results"
echo -e ""
exit 1 # Exit script after printing help
}

while getopts "a:b:c:" opt
do
case "$opt" in
a ) Fasta="$OPTARG" ;;
b ) FQfolder="$OPTARG" ;;
c ) RESFolder="$OPTARG" ;;
? ) HELPME ;; # Print helpFunction in case parameter is non-existent
esac
done

# Print helpFunction in case parameters are empty
if [ -z "$Fasta" ] || [ -z "$FQfolder" ] || [ -z "$RESFolder" ]
then
echo ""
echo "Galogal is pipeline written in bash. this set of tools is used to map reads and call variants";
echo "The tools BWA, samtools, picard and GATK  are used and need to be installed before running Galocal. ";
echo "Importantly at the end of this step you will have a VCF file and a base level depth of coverage  file.";

HELPME
fi

# Begin script.
#-----------------------------------------------------------------------------
# Index the fast files  undo this part if you need to index your fasta
#-----------------------------------------------------------------------------

#samtools faidx /PATH/TO/YOUR/file.fasta
#java -jar ~/picard-tools-2.1.0/picard.jar CreateSequenceDictionary R=/PATH/TO/YOUR/file.fasta O=/OUTPUT/PATH/TO/YOUR/file.dict
#bwa index /PATH/TO/YOUR/file.fasta

#-------------------------------------------
#Path to your User repository from NGS group
#-------------------------------------------
#folder name is similar to read names
PathToF=$Fasta
PathToFq=$FQfolder/
#-------------------------------------------
#Path to your User repository from NGS group
#-------------------------------------------
PathToRBam=$RESFolder/
PathToRVCF=$RESFolder/
PathToRDofC=$RESFolder/

#######################
cd $PathToFq
i=0
for ReadsFq in $(ls -R)
do
if
[ -d $ReadsFq ]
then
cd $ReadsFq
echo $ReadsFq

fastafile=$PathToF
echo ${file%%_*}
Rreads=$PathToFq$ReadsFq/${ReadsFq%%.*}"_1.fastq"
Freads=$PathToFq$ReadsFq/${ReadsFq%%.*}"_2.fastq"
DoCovfile=$PathToRBam${ReadsFq%%.*}".txt"
vcffile=$PathToRVCF${ReadsFq%%.*}".vcf"
BamFile=$PathToRDofC${ReadsFq%%.*}".bam"

echo "##################"
echo "##################"
echo "Your PATH to Ref Fasta: " $fastafile
echo "Your PATH to R1 Reads : " $Rreads
echo "Your PATH to R2 Reads : " $Freads
echo $BamFile
echo "##################"
echo "##################"
########################

#bwa index $fastafile
#java -jar ~/picard-tools-2.1.0/picard.jar CreateSequenceDictionary R=$fastafile O=$fastafile.dict
#java -jar ~/picard-tools-2.1.0/picard.jar CreateSequenceDictionary R=/Users/MoiMeamina/Downloads/G4-Hunter-master/Mitochondria_NC_012920_1.fasta O=/Users/MoiMeamina/Downloads/G4-Hunter-master/Mitochondria_NC_012920_1.dict

#Aligne reads on the fastafile genome
bwa mem $fastafile $Rreads $Freads > bamfile.bam
########################
#Filter the reads to eliminates the unmapped reads
########################
samtools view -bS bamfile.bam | samtools view -b -F 4 | samtools sort -o tmpbamfile-S.bam

########################
#Workflow Picard
########################
# add replace groups

java -jar ~/picard-tools-2.1.0/picard.jar AddOrReplaceReadGroups I=tmpbamfile-S.bam O=tmpbamfile-ARgroups.bam VALIDATION_STRINGENCY=LENIENT RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
########################
#reorder

java -jar ~/picard-tools-2.1.0/picard.jar ReorderSam I=tmpbamfile-ARgroups.bam O=tmpbamfile-Reor.bam R=$fastafile CREATE_INDEX=TRUE
########################
#Marque and Eliminate duplicats

java -jar ~/picard-tools-2.1.0/picard.jar MarkDuplicates REMOVE_DUPLICATES= True I=tmpbamfile-Reor.bam O=$BamFile VALIDATION_STRINGENCY=LENIENT METRICS_FILE=marked_dup_metrics.txt
######################
#Index bam File

samtools index $BamFile

######################
#Variant calling
freebayes -f $fastafile $BamFile > $vcffile

########################
#Depth of covrage
java -jar ~/GATK-3.5/GenomeAnalysisTK.jar -T DepthOfCoverage -R $fastafile -o $DoCovfile -I $BamFile
#[-geneList refSeq.sorted.txt] #[-pt readgroup] #[-ct 4 -ct 6 -ct 10] #[-L my_capture_genes.interval_list]

########################
# All this tools generate secondary files that would be stocked at each Result repository
# we dele
########################
rm bamfile.bam
rm tmpbamfile-ARgroups.bam
rm tmpbamfile-Reor.bam
rm tmpbamfile-S.bam
rm tmpbamfile-Reor.bai
rm marked_dup_metrics.txt

((i++))
cd ..
fi
done

########################
# Some more secondary files
########################
cd $PathToRBam

rm *.sample_cumulative_coverage_counts
rm *.sample_cumulative_coverage_proportions
rm *.sample_interval_statistics
rm *.sample_interval_summary
rm *.sample_statistics
rm *.sample_summary

