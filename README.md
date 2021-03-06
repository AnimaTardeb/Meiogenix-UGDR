# UGDR #

This repository provides the UGDR pipeline that allows users to detect recombination between a reference/parent and a recombined yeast.
This script extracts information from VCF files, defines regions of recombination, Calculate normalized depth of coverage and plot recombination and depth of coverage.

## Basic Usage ##

<p align="center">
<img src="https://github.com/AnimaTardeb/Meiogenix-UGDR/blob/master/DATA/figure1.jpg" alt="" width="500" height="550">
</p>

<p align="center">
  <sub>The pipline is build by Amina Bedrat (2016) at Curie Institut and Meiogenix.
  </a>
</p>

![R](https://img.shields.io/badge/R-v3.6.0+-orange.svg)
![Python](https://img.shields.io/badge/python-v2.7+-blue.svg)
[![License](https://img.shields.io/badge/license-GNU_v3-green.svg)](https://www.gnu.org/licenses/gpl-3.0.fr.html)


## Requirements & Obtaining ##

* Python 2.7 and higher. 
* Samtools, BWA, Picard and freebayes
* Varient Calling Files (VCF) ([More informations about VCF file format](http://www.1000genomes.org/wiki/Analysis/vcf4.0)).
* To download UGDR, please use git to clone or download via:
```
git clone https://github.com/AnimaTardeb/Meiogenix-UGDR.git

```

## Table of Contents
- [Galocal](#Galocal)
- [UGDR](#UGDR)
- [Normalized Depth of Coverage](#NDoC)

------------

## Galocal ##
Galogal is a set of tools is used to map reads to the reference and call variants. The tools BWA, samtools, picard and GATK are used and need to be installed before running Galocal. 
Importantly at the end of this step you will have a VCF file and a base level depth of coverage file.

### Run Galocal ###

```
$./Galocal.sh

Usage: 
 ./Galocal.sh -a Fasta -b FQ_folder -c Result_folder
	-a S288 or other fasta sequence
	-b Folder where the reads (R_1.fq and R_2.fq) are stocked
	-c Folder to out put the results
 
```

## UGDR ##

UGDR a tool to  analyze alleles variation and identify regions ofrecombination in yeasts. This script compares two vcf files and extracts and plot recombination profile.

### Run UGDR ###

```
$ cd UGDR
$ python UGDR.py
Usage: 
python UGDR.py [options] REC_vcf_file [options] Reference_vcf_file [options] output_Repository 

Options:
  --version       show program's version number and exit
  -h, --help      show this help message and exit
  -i, --REC_rep   Folder of more than one recombinant to test
  -I, --REC_file  One recombinant (VCF file) to test (REC)
  -j, --par_file  Reference (VCF) file
  -o, --out_dir   Results folder 
```

> To use more then one recombinant_VCFfile you should gather them in one repository and choose the option -I 

### Description of output file  ###

* -ParentalAlleles.txt : File summarize all the alleles present in the VCF file of the reference strain (reference_VCFfile). 
* -RTGAlleles.txt : File summarize all the alleles presents in the VCF file of the RTG strain (RTG_VCFfile).
* -InvarAllelels.txt : File with alleles that the ratio does not vary between RTGs and reference strain.  
* -VarAlleles.txt : File with alleles that the ratio had varied between RTGs and reference strain.
* -RR.txt : Recombination region file, it extracts all the possible recombined region.

## Normalized Depth of Coverage ##

Normalized Depth of Coverage (NDoC), a tool to analyze alleles coverage and identify regions of coverage variation. The inputs are two depth of coverage files generated for UGDR. The per chromosome depth profile (strain__normalizedXKb.txt) is generated at the end.

> default XKB = 1KB: Normalized by 1kb.

### Run NDoC ###

```
$ python NDoC.py
Usage: 
python NDoC.py [options] Referenc_Depth_of_coverage_file [options] RTG_Depth_of_coverage_file [options] output_Repository 

Options:
  --version      show program's version number and exit
  -h, --help     show this help message and exit
  -i, --ref_f    Depth of coverage file of the reference strain
  -j, --RTG_f    Depth of coverage file of Recombined strain
  -o, --out_dir  Results repository

```

### Description of output file  ###

* 
*
*
 

### Bug Reports & Requests

Contact : bedrat.amina@gmail.com.

Please use the [issue tracker](https://github.com/AnimaTardeb/Meiogenix-UGDR/issues) to report any bugs or file feature requests.

