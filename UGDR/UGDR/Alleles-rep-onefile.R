#!/usr/local/bin/R
#####################################################################################
#####################################################################################
# -*- coding: utf8 -*-
# This file is part of R2D2-BA
# This script plot variable and invariable Alleles of an RTG
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
#cat Alleles-rep-onerep.R | /usr/local/bin/R --slave --args PATH2INPUTfile PATH2OUTPUT.pdf
#####################################################################################
args = commandArgs()
PATH2INPUT=args[4] #Repository of results obtained by R2D2BA
PATH2OUTPUT=args[5] #Path to a pdf output
#####################################################################################
#####Import fils#####
chr_info <- read.table("Resources/chr_nb.txt",header=F) #Chr#0#length#centromer

pdf(file = PATH2OUTPUT, onefile = TRUE, width=10, height=13)
pathtoseg <- PATH2INPUT
#####################################################################################
listdirs<-list.dirs(pathtoseg)

table <- data.frame("RTG" = character(), "Pvalue" = numeric(),"PV" = numeric(),stringsAsFactors=FALSE)
for (t in 1:length(listdirs)){
    filespath<-listdirs[t]
    #print (filespath)
    
    #filesnameInvAll<-list.files(filespath, pattern = c("-Parental.txt"))
    #filesnameVarAll<-list.files(filespath, pattern = c("-VH.txt" ))
    #filesRefAll<-list.files(filespath, pattern = c("-ParentalAllele.txt"))
    
    filesnameInvAll<-list.files(filespath, pattern = c("-InvarAllelels.txt"))
    filesnameVarAll<-list.files(filespath, pattern = c("-VarAlleles.txt" ))
    filesRefAll<-list.files(filespath, pattern = c("-ReferenceAllele.txt"))
    
    RTG <-basename(filespath) #name of the RTG according to the file name
    
    RTGP<- read.table(paste(filespath,"/",filesnameInvAll,sep=""),fill=TRUE,sep="\t")
    RTGVH <- read.table(paste(filespath,"/",filesnameVarAll,sep=""),fill=TRUE,sep="\t")
    RTGparALL<-read.table(paste(filespath,"/", filesRefAll,sep=""),fill=TRUE,sep="\t")
    
    ####################################################################################################################################################
    #plot invisible points to label Y axe
    ####################################################################################################################################################
    plot(x=chr_info[,2],y=chr_info[,1],pch=10, col = "white", xlab="Coordinate (kb)", ylab="Chromosomes", ylim=rev(c(1,17)),xlim=c(0,1550000),cex.lab=1, cex.main=1,frame=FALSE, yaxt="n", xaxt="n", cex.axis=2)
    redefinition=colors()[1:657]
    palette(redefinition)
    
    #####label axes
    axis(side=2, at=c(1:16), labels=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"), tick=FALSE, line=-2, pos=NA, cex.axis=1.5,las=1)
    axis(side=1, at=c(0,250000,500000,750000,1000000,1250000,1500000), labels=c(0,250,500,750,1000,1250,1500), tick=TRUE, line=-2, pos=NA, cex.axis=1.6,tcl=0.5)
    
    #####legend
    legend(900000,7.5, c("Centromere."), cex= 0.7, col=c(153),pch=21,bty="n",pt.cex=1,pt.lwd=1);
    legend(900000,8, c("Heterozygote parental markers. "), cex= 0.7, col=c(314),pch= 124,bty="n",pt.cex=1,pt.lwd=1);
    legend(900000,8.5, c("RTG invariable markers."), cex= 0.7, col=c(620),pch= 124,bty="n",pt.cex=1,pt.lwd=1);
    legend(900000,9.5, c("RTG Other varying markers ( allelic ratio != 0 and !=1)"), cex=0.7, bty="n",col=c(155),pch= 124,pt.cex=1,pt.lwd=1);
    legend(900000,9, c("RTG LOH (allelic ratio =0)"), cex= 0.7, col=c(552),pch= 124,bty="n",pt.cex=1,pt.lwd=1);
    legend(900000,10, c("RTG LOH (allelic ratio =1)"), cex= 0.7, col=c(490),pch= 124,bty="n",pt.cex=1,pt.lwd=1);


    mtext(RTG)
    
    #####plot black line to define chromosomes
    #####Chromosome length should vary from one strain to another
    
    for (i in c(-0.1, -0.3)){
        points(x=c(0,230218), y=c(1+i,1+i), type="l"); points(x=c(0,813184), y=c(2+i,2+i), type="l"); points(x=c(0,316620), y=c(3+i,3+i), type="l"); points(x=c(0,1531933), y=c(4+i,4+i), type="l");
        points(x=c(0,576874), y=c(5+i,5+i), type="l"); points(x=c(0,270161), y=c(6+i,6+i), type="l"); points(x=c(0,1090940), y=c(7+i,7+i), type="l"); points(x=c(0,562643), y=c(8+i,8+i), type="l");
        points(x=c(0,439888), y=c(9+i,9+i), type="l"); points(x=c(0,745751), y=c(10+i,10+i), type="l"); points(x=c(0,666816), y=c(11+i,11+i), type="l"); points(x=c(0,1078177), y=c(12+i,12+i),
        type="l"); points(x=c(0,924431), y=c(13+i,13+i), type="l"); points(x=c(0,784333), y=c(14+i,14+i), type="l"); points(x=c(0,1091291), y=c(15+i,15+i), type="l"); points(x=c(0,948066), y=c(16+i,16+i), type="l")}
    
    #####Reference
    points(x=subset(RTGparALL$V3,RTGparALL$V6<=0.8 & RTGparALL$V2=="snp"),y=(subset(RTGparALL$V1,RTGparALL$V6<=0.8 & RTGparALL$V2=="snp")-0.3),pch=124, col = 314, lwd = 7, cex=0.8)
    
    ####Invariable RTG alleles
    points(x=subset(RTGP$V2, RTGP$V10<=1.2),y=(subset(RTGP$V1, RTGP$V10<=1.2)-0.1),pch=124, col = 620, lwd = 1, cex=0.8)
    
    #####Variable RTG alleles
    points(x=subset(RTGVH$V2, RTGVH$V10==0.0),y=(subset(RTGVH$V1, RTGVH$V10==0.0)-0.1),pch=124, col = 552, lwd = 1, cex=0.8)#rouge
    points(x=subset(RTGVH$V2, RTGVH$V10>0.9),y=(subset(RTGVH$V1, RTGVH$V10>0.9)-0.1),pch=124, col = 490, lwd = 1, cex=0.8) #bleu
    points(x=subset(RTGVH$V2, 0.29<=RTGVH$V10 & RTGVH$V10<=0.67),y=(subset(RTGVH$V1, 0.29<=RTGVH$V10 & RTGVH$V10<=0.67)-0.1),pch=124, col = 155, lwd = 1, cex=0.8)#Noir
    
    #####plot centromere
    points(x=chr_info[,4],y=(chr_info[,1]-0.1),pch=21, col = 153, cex=1.5,lwd=1)
    points(x=chr_info[,4],y=(chr_info[,1]-0.3),pch=21, col = 153, cex=1.5,lwd=1)
    #####close
}
dev.off()