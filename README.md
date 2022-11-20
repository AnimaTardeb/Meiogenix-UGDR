# UGDR #

UGDR is a pipeline that analyzes variations in SNP frequency to detect recombination and foster visualization of recombination tracks in natural and/or constructed hybrid yeast genomes. The detection of recombination in haploid (tetrads), diploid, aneuploidy and polyploid yeasts (limited to 4n) is achieved regardless of the ploidy of the yeast and the parental genome. UGDR also addresses the challenge of the continuous variation of the allele frequencies along with chromosomal loss and gain.

## Basic Usage ##

<p align="center">
<img src="https://github.com/AnimaTardeb/Meiogenix-UGDR/blob/master/DATA/figure1-2.jpg" alt="" width="800" height="650">
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

If VCF files and/or depth of coverage files are already present the user might skip this part and directly run [UGDR ](#UGDR). 

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
> Note: the depth of the base across the samples is now a parameter that the user can change to adapt it to their sequencing depth.
In this example we are running UGDR with a DP = 200 (-c 200)

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
  -c NUM, --DPQUAL=NUM  Filter applied on DP for 80X use 200

```

To run an example: 

```

$./UGDR.py -I ../../DATA/Tetradploids/SRR3265444/SRR3265444.vcf -j ../../DATA/Tetradploids/SRR3265445/SRR3265445.vcf -o ../../DATA -c 200

########################################################################
#        Results-SRR3265444 directory Created                              #
########################################################################

Reference file is proceeded. Alleles are summarized in -ReferencelAllele.txt file 
snp  :  75118
del  :  2120
complex  :  2598
ins  :  2041
mnp  :  916
_______________________________________________________________

Recombined file (REC) is proceeded. Alleles are summarized in -RECAllele.txt file

snp  :  79345
del  :  2217
complex  :  2673
ins  :  2116
mnp  :  967
_______________________________________________________________

 Computing Ratio with Qual >= 200 and AO >= 20 (you can modify the Qual and AO in SNPdistribution.py lines 234, 235 and 259).

_______________________________________________________________

 Extracting Recombined regions  
/DATA/Results-SRR3265444/
```

> See result ([here](https://github.com/AnimaTardeb/Meiogenix-UGDR/blob/master/DATA/Tetradploids/Results-SRR3265445-ref44/SRR3265445.pdf)
> For analysing more than one recombinant_VCFfile choose the option -I .

### Description of the output  ###

* -ParentalAlleles.txt : A file to summarize all the alleles present in the VCF file of the reference strain (reference_VCFfile). 
* -RTGAlleles.txt : A file to summarize all the alleles presents in the VCF file of the RTG strain (RTG_VCFfile).
* -InvarAllelels.txt : A file of alleles which ratio does not vary between reference and RTG strains.  
* -VarAlleles.txt : A file of alleles which ratio varied between reference and RTG strains.
* -RR.txt : Recombination region file, it extracts all possible recombined regions.
* the PDF output representing the recombination

## Normalized Depth of Coverage ##

Normalized Depth of Coverage (NDoC) analyzes alleles coverage and identifies regions of coverage variation. The inputs are two depth of coverage files generated by **gatk DepthOfCoverage** from  UGDR. The output is a per chromosome depth profile file (strain__normalizedxKb.txt).

> the -c argument will help you choose the window of the mean Best use = 1000.

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
  -c NUM         Mean depth of covrage per k kb


$ ./NDoC.py -i ../../DATA/NDof-Files/DoCALpha.txt -j ../../DATA/NDof-Files/DoCBeta.txt -o /Users/amyqueen/Downloads/Meiogenix-UGDR-master-2/DATA/NDof-Files -c 1000

########################################################################
#        DepthOfCov-Results-DoCBeta directory Created                              #
# Calculate the mean Depth of coverage per 1000 KB         
########################################################################
Reference DOC file:  DoCALpha
Number of base : 12015280 
 Total Mean of Coverage:  77.8140244755

Recombined DOC file:  DoCBeta
Number of base : 12015280 
 Total Mean of Coverage:  124.255924207

```

### Description of output file  ###

* ref.txt and rec.txt : dept of coverage per 1Kb for each reference and recombined strains 
* strain__normalizedxKb.txt : Normalized depth of covrage used to plot the depth of coverage
* results of the depth of coverage in a pdf file 
*  PDF : plot of depth of coverage  

> See result ([here](https://github.com/AnimaTardeb/Meiogenix-UGDR/blob/master/DATA/NDof-Files/DepthOfCov-Results-DoCBeta/DoCBeta.pdf)

### Plot Recombination ###

In addition to the plots created by UGDR and NDoC you can use to plot Mother and Daughter or the four tetrads in one figure. 

```
$ cd UGDR 
$ cat Alleles-rep-tetrad.R | /usr/local/bin/R --slave --args ../../DATA/VCF-tetrads/Tetrads-S288c-Results/ ../../DATA/VCF-tetrads/Tetrads-S288c-Results/PlotTetradsS288c.pdf 

[1] "Reading : "
[1] "../../DATA/VCF-tetrads/Tetrads-S288c-Results//Results-SRR2984915"
[1] "../../DATA/VCF-tetrads/Tetrads-S288c-Results//Results-SRR2984916"
[1] "../../DATA/VCF-tetrads/Tetrads-S288c-Results//Results-SRR2984917"
[1] "../../DATA/VCF-tetrads/Tetrads-S288c-Results//Results-SRR2984918"
[1] " "
[1] "PDF plot generated in "
[1] "../../DATA/VCF-tetrads/Tetrads-S288c-Results/PlotTetradsS288c.pdf"
```


The result is ([here](https://github.com/AnimaTardeb/Meiogenix-UGDR/blob/master/DATA/VCF-tetrads/Tetrads-S288c-Results/PlotTetradsS288c.pdf)) 


To plot Two yeasts (a, b), the Alleles-Rep-MD-DoC.R needs : 
* UGDR result folder for the yeast a
* UGDR result folder for the yeast b
* Normalized depth of coverage file for the yeast a
* Normalized depth of coverage file for the yeast b

>a, b could be mother and dauther, Parent1 and parent 2 or Two different cells e.g. triploid and tetraploid.

```
$ cat Alleles-Rep-MD-DoC.R | /usr/local/bin/R --slave --args ../../DATA/Triploids/Triploids-S288C/Results-SRR3265371 ../../DATA/Tetradploids/Results-SRR3265445-ref44   ../../DATA/Triploids/Triploids-S288C/DepthOfCov-Results-SRR3265371/SRR3265371__normalized10kb.txt ../../DATA/Tetradploids/DepthOfCov-Results-SRR3265445/SRR3265445__normalized1kb.txt ../../DATA/Tripltetraploids.pdf


[1] "Reading : "
[1] "../../DATA/Triploids/Triploids-S288C/Results-SRR3265371"
[1] "../../DATA/Tetradploids/Results-SRR3265445-ref44"
[1] " PDF plot generated in "
[1] "../../DATA/Tripltetraploids.pdf"

```
The result is ([here](https://github.com/AnimaTardeb/Meiogenix-UGDR/blob/master/DATA/Tripltetraploids.pdf)) 


### Bug Reports & Requests

Contact : bedrat.amina@gmail.com

Please use the [issue tracker](https://github.com/AnimaTardeb/Meiogenix-UGDR/issues) to report any bugs or file feature requests.

