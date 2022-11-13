# UGDR #

UGDR is a pipeline to detect recombination in complex hybrid/recombined yeast genomes (<=4n). It is based on alleles ratio comparative analyses  and read depth coverage obtained from NGS reads.

## Basic Usage ##

<p align="center">
<img src="https://github.com/AnimaTardeb/Meiogenix-UGDR/blob/master/DATA/figure1.jpg" alt="" width="500" height="550">
</p>

<p align="center">
  <sub> UGDR workflow: a) Galocal, b) UGDR, c) NDoC.
  </a>
</p>

![R](https://img.shields.io/badge/R-v3.6.0+-orange.svg)
![Python](https://img.shields.io/badge/python-v2.7+-blue.svg)
[![License](https://img.shields.io/badge/license-GNU_v3-green.svg)](https://www.gnu.org/licenses/gpl-3.0.fr.html)


## Requirements ##

* Python 2.7 and higher. 
* Samtools, BWA, GATK and freebayes
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
Galocal is a bash script used to map reads to a reference genome, call variants and calculate the depth of covrage. A prior installation of those tools in needed before running Galocal [See here ](#Requirements). 

Importantly at the end of this step, Galocal you will output in the result folder : a VCF file and base level depth of coverage file.

If VCF files and/or depth of coverage files are already present the user directly run [UGDR ](#UGDR). 

### Run Galocal ###

```
$./Galocal.sh

Usage: 
 ./Galocal.sh -a Fasta -b Fq_folder -c Result_folder
	-a S288 or other fasta sequence
	-b Folder where the reads (R_1.fq and R_2.fq) are stocked
	-c Folder to out put the results
 
```

## UGDR ##

UGDR analyzes the alleles variation and identify regions of recombination in yeasts. This script compares two vcf files and plot the recombination profile.

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
  -I, --REC_file  **One** recombinant (VCF file) to test (recombined strain)
  -j, --par_file  Reference (VCF) file (reference strain)
  -o, --out_dir   Results folder 
```

> To run an example: 

```
$./UGDR.py -I ../DATA/Tetradploids/SRR3265444/SRR3265444.vcf -j ../DATA/Tetradploids/SRR3265445/SRR3265445.vcf -o /DATA 
```

> For analysing more than one recombinant_VCFfile choose the option -I .

### Description of the output  ###

* -ParentalAlleles.txt : A file to summarize all the alleles present in the VCF file of the reference strain (reference_VCFfile). 
* -RTGAlleles.txt : A file to summarize all the alleles presents in the VCF file of the RTG strain (RTG_VCFfile).
* -InvarAllelels.txt : A file of alleles which ratio does not vary between reference and RTG strains.  
* -VarAlleles.txt : A file of alleles which ratio varied between reference and RTG strains.
* -RR.txt : Recombination region file, it extracts all possible recombined regions.

## Normalized Depth of Coverage ##

Normalized Depth of Coverage (NDoC) analyzes alleles coverage and identifies regions of coverage variation. The inputs are two depth of coverage files generated by **gatk DepthOfCoverage** from  UGDR. The output is a per chromosome depth profile file (strain__normalizedxKb.txt).

> default xKB = 1KB: Normalized by 1kb.

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

* ref.txt and rec.txt : dept of coverage per 1Kb for each reference and recombined strains 
* strain__normalizedxKb.txt : Normalized depth of covrage used to plot the depth of coverage
* results of the depth of coverage in a pdf file 
 

### Bug Reports & Requests

Contact : bedrat.amina@gmail.com

Please use the [issue tracker](https://github.com/AnimaTardeb/Meiogenix-UGDR/issues) to report any bugs or file feature requests.

