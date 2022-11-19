#!/usr/local/bin/R
#####################################################################################
#####################################################################################
# -*- coding: utf8 -*-
# This file is part of UGDR-BA
# This script extract informations from VCF files to give a human readable
# files resuming alleles variations, and defining region of recombination.
# Copyright (C) 2016 Bedrat Amina - Institut Curie - Meiogenix
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
#cat Alleles-rep-allrep.R | /usr/local/bin/R --slave --args PATH2INPUTREP PATH2OUTPUT.pdf &
#NOTE parametre chr_info[,3] and chr_info[,4 ]pour S288C  et chr_info[,5] and chr_info[,6 ] pour SK1
#cat Alleles-rep-tetrad.R | /usr/local/bin/R --slave --args /Tetrads-S288c-Results /plotTetradP1S288c.pdf
#################################################################################################################################################################
args = commandArgs()

if (length(args)!=5){
    stop("Arguments missing. \
You need the folder containing the results from UGDR for the 4 tetrads and a pather to the output.pdf")
} else if (length(args)==5) {

    PATH2INPUT=args[4] #repository
    PATH2OUTPUT=args[5]
}


#################################################################################################################################################################
#####Importer les fichiers#####
chr_info <- read.table("Resources/chr_nb.txt",header=F) #Chr#0#length#centromer

pdf(file = PATH2OUTPUT, onefile = TRUE, width=8, height=11)
pathtoseg <- PATH2INPUT #
#################################################################################################################################################################;

listdirs<-list.dirs(pathtoseg)
#print (length(listdirs))
table <- data.frame("RTG" = character(), "Pvalue" = numeric(),"PV" = numeric(),stringsAsFactors=FALSE)
k<-0
####################################################################################################################################################
#plot invisible points to label Y axe
####################################################################################################################################################
plot(x=chr_info[,2],y=chr_info[,1],pch=2, col = "white", xlab="Coordinate (kb)", ylab="Chromosomes", ylim=rev(c(1,17)),xlim=c(0,1550000),cex.lab=1, cex.main=1,frame=FALSE, yaxt="n", xaxt="n", cex.axis=1)
redefinition=colors()[1:657]
palette(redefinition)
#####label axes
axis(side=2, at=c(1:16), labels=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"), tick=FALSE, line=-1, pos=NA, cex.axis=1.3,las=1)
axis(side=1, at=c(0,250000,500000,750000,1000000,1250000,1500000), labels=c(0,250,500,750,1000,1250,1500), tick=TRUE, line=-2, pos=NA, cex.axis=1.6,tcl=0.5)

#####legend
legend(900000, 8, c("Centromere."), cex= 0.7, col=c(153),pch=21,bty="n",pt.cex=1,pt.lwd=1);
legend(900000,8.25, c("Parent 1"), cex= 0.7, col=c(554),pch= 124,bty="n",pt.cex=1,pt.lwd=1);#red
legend(900000,8.5, c("Parent 1"), cex= 0.7, col=c(618),pch= 124,bty="n",pt.cex=1,pt.lwd=1);#blue


#####plot black line to define chromosomes
#####Chromosome length should vary from one strain to another

for (i in c(0.1,-0.1, -0.3, -0.5)){
    z= chr_info[,3]
    #print(z[1])
    points(x=c(0,z[1]), y=c(1+i,1+i), type="l"); points(x=c(0,z[2]), y=c(2+i,2+i), type="l"); points(x=c(0,z[3]), y=c(3+i,3+i), type="l"); points(x=c(0,z[4]), y=c(4+i,4+i), type="l");
    points(x=c(0,z[5]), y=c(5+i,5+i), type="l"); points(x=c(0,z[6]), y=c(6+i,6+i), type="l"); points(x=c(0,z[7]), y=c(7+i,7+i), type="l"); points(x=c(0,z[8]), y=c(8+i,8+i), type="l");
    points(x=c(0,z[9]), y=c(9+i,9+i), type="l"); points(x=c(0,z[10]), y=c(10+i,10+i), type="l"); points(x=c(0,z[11]), y=c(11+i,11+i), type="l"); points(x=c(0,z[12]), y=c(12+i,12+i),
    type="l"); points(x=c(0,z[13]), y=c(13+i,13+i), type="l"); points(x=c(0,z[14]), y=c(14+i,14+i), type="l"); points(x=c(0,z[15]), y=c(15+i,15+i), type="l"); points(x=c(0,z[16]), y=c(16+i,16+i), type="l")}
print("Reading : ")
for (t in 2:length(listdirs)){
    
    filespath<-listdirs[t]
    print (filespath)
    filesnameparental<-list.files(filespath, pattern = c("-InvarAllelels.txt"))
    filesnameVH<-list.files(filespath, pattern = c("-VarAlleles.txt" ))
    filesall<-list.files(filespath, pattern = c("-ReferenceAllele.txt"))
    
    RTG <-basename(filespath) #name of the RTG according to the file name
    RTGP<- read.table(paste(filespath,"/",filesnameparental,sep=""),fill=TRUE,sep="\t")
    RTGVH <- read.table(paste(filespath,"/",filesnameVH,sep=""),fill=TRUE,sep="\t")
    RTGparALL<-read.table(paste(filespath,"/", filesall,sep=""),fill=TRUE,sep="\t")
    
    #mtext(RTG)
    #####RTG alleles parental
    #print (t)
    #print(k)
    #print(1.5+k-t)
    
    points(x=subset(RTGP$V2, RTGP$V10<=0.8),y=(subset(RTGP$V1, RTGP$V10<=0.8)+(1.5+k-t)),pch=124, col = 314, lwd = 1.5, cex=0.75)
    #####RTG variable alleles
    points(x=subset(RTGVH$V2, RTGVH$V10==0.0),y=(subset(RTGVH$V1, RTGVH$V10==0.0)+(1.5+k-t)),pch=124, col = 554, lwd = 1.5, cex=0.75)#rouge
    points(x=subset(RTGVH$V2, RTGVH$V10>0.8),y=(subset(RTGVH$V1, RTGVH$V10>0.8)+(1.5+k-t)),pch=124, col = 618, lwd = 1.5, cex=0.75)#bleu
    points(x=subset(RTGVH$V2, 0.4<=RTGVH$V10 & RTGVH$V10<=0.55),y=(subset(RTGVH$V1, 0.4<=RTGVH$V10 & RTGVH$V10<=0.55)+(1.5+k-t)),pch=124, col = 155, lwd = 1.5, cex=0.75) #Noir
    k<-k+1.2
    
    #####plot ctm
    points(x=chr_info[,4],y=(chr_info[,1]-0.1),pch=21, col = 153, cex=1.25,lwd=1)
    points(x=chr_info[,4],y=(chr_info[,1]-0.3),pch=21, col = 153, cex=1.25,lwd=1)
    points(x=chr_info[,4],y=(chr_info[,1]-0.5),pch=21, col = 153, cex=1.25,lwd=1)
    points(x=chr_info[,4],y=(chr_info[,1]+0.1),pch=21, col = 153, cex=1.25,lwd=1)
    
    #####close
}

print("  \
PDF plot generated in ")
print(PATH2OUTPUT)

garbage <- dev.off()
