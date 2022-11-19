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
#cat Alleles-rep-onerep.R | /usr/local/bin/R --slave --args PATH2INPUTMrep PATHTODOFCOVFILE PATH2OUTPUT.pdf
#cat Alleles-Rep-MD-DoC.R | /usr/local/bin/R --slave --args S288c-MD/Results-A143R14 S288c-MD/Results-SRR2984813 /S288c-MD/DepthOfCov-Results-A143R14/A143R14__normalized1kb.txt /DepthOfCov-Results-SRR2984813/SRR2984813__normalized1kb.txt /DATA/DIploids_MD/MDREP.pdf
#####################################################################################
args = commandArgs()

if (length(args)!=8){
    stop("At least one argument must be supplied (input file)")
} else if (length(args)==8) {

    PATH2INPUTM=args[4] #Repository of results obtained by R2D2BA
    PATH2INPUTD=args[5] #Repository of results obtained by R2D2BA

    PATH3INPUTM=args[6] #file path to Dof cov file 1kb
    PATH3INPUTD=args[7] #file path to Dof cov file 1kb
    
    PATH2OUTPUT=args[8] #Path to a pdf output

}

#####################################################################################
#####Import fils#####
chr_info <- read.table("Resources/chr_nb.txt",header=F) #Chr#0#length#centromer
pdf(file = PATH2OUTPUT, onefile = TRUE, width=10, height=13)

pathtosegM<- PATH2INPUTM #Mother
pathtosegD <- PATH2INPUTD #Daughter

#####################################################################################
#first strain plot
###############
print("Reading : ")
listdirs<-list.dirs(pathtosegM)
table <- data.frame("RTG" = character(), "Pvalue" = numeric(),"PV" = numeric(),stringsAsFactors=FALSE)
for (t in 1:length(listdirs)){
    filespath<-listdirs[t]
    print (filespath)

    filesnameInvAll<-list.files(filespath, pattern = c("-InvarAllelels.txt"))
    filesnameVarAll<-list.files(filespath, pattern = c("-VarAlleles.txt" ))
    RTG <-basename(filespath) #name of the RTG according to the file name
    
    RTGP<- read.table(paste(filespath,"/",filesnameInvAll,sep=""),fill=TRUE,sep="\t")
    RTGVH <- read.table(paste(filespath,"/",filesnameVarAll,sep=""),fill=TRUE,sep="\t")
    newdata<-read.table(paste(PATH3INPUTM,sep=""),fill=TRUE,sep="\t")
    newdata[,1] <- suppressWarnings(as.numeric(as.character(newdata[,1])))
    DOCFile <- newdata[with(newdata,order(newdata$V1)),]
    
    ####################################################################################################################################################
    #plot invisible points to label Y axe
    ####################################################################################################################################################
    
    plot(x=chr_info[,2],y=chr_info[,1],pch=10, col = "white", xlab="Coordinate (kb)", ylab="Chromosomes", ylim=rev(c(1,17)),xlim=c(0,1550000),cex.lab=1, cex.main=1,frame=FALSE, yaxt="n", xaxt="n", cex.axis=2)
    redefinition=colors()[1:657]
    palette(redefinition)
    
    #####label axes
    axis(side=2, at=c(1:16), labels=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"), tick=FALSE, line=-1, pos=NA, cex.axis=1,las=1)
    axis(side=2, at=seq(0.65,16.15,0.5), labels=c("a","b","a","b","a","b","a","b","a","b","a","b","a","b","a","b","a","b","a","b","a","b","a","b","a","b","a","b","a","b","a","b"), tick=FALSE, line=-2, pos=NA, cex.axis=0.65,las=1)
    axis(side=1, at=c(0,250000,500000,750000,1000000,1250000,1500000), labels=c(0,250,500,750,1000,1250,1500), tick=TRUE, line=-2, pos=NA, cex.axis=1.6,tcl=0.5)
    #####legend
    legend(900000,7.75, c("Centromere."), cex= 0.7, col=c(153),pch=21,bty="n",pt.cex=1,pt.lwd=1);
    legend(900000,8, c("Invariable markers. "), cex= 0.7, col=c(314),pch= 124,bty="n",pt.cex=1,pt.lwd=1);
    legend(900000,8.75, c("Other varying markers ( allelic ratio != 0 and !=1)"), cex=0.7, bty="n",col=c(155),pch= 124,pt.cex=1,pt.lwd=1);
    legend(900000,8.25, c("LOH (allelic ratio =0)"), cex= 0.7, col=c(554),pch= 124,bty="n",pt.cex=1,pt.lwd=1);
    legend(900000,8.5, c("LOH (allelic ratio =1)"), cex= 0.7, col=c(490),pch= 124,bty="n",pt.cex=1,pt.lwd=1);
    legend(900000,9.25, c("Normalized Depth of Cov. <0.5"), cex= 0.7, col=c(630),pch= 20,bty="n",pt.cex=1,pt.lwd=1);
    legend(900000,9.5, c("Normalized Depth of Cov. >1.5"), cex= 0.7, col=c(448),pch= 20,bty="n",pt.cex=1,pt.lwd=1);
    mtext("")
    #####plot black line to define chromosomes
    #####Chromosome length should vary from one strain to another
    for (i in c(-0.3, 0.1)){
        z= chr_info[,3]
        #print(z[1])
        points(x=c(0,z[1]), y=c(1+i,1+i), type="l"); points(x=c(0,z[2]), y=c(2+i,2+i), type="l"); points(x=c(0,z[3]), y=c(3+i,3+i), type="l"); points(x=c(0,z[4]), y=c(4+i,4+i), type="l");
        points(x=c(0,z[5]), y=c(5+i,5+i), type="l"); points(x=c(0,z[6]), y=c(6+i,6+i), type="l"); points(x=c(0,z[7]), y=c(7+i,7+i), type="l"); points(x=c(0,z[8]), y=c(8+i,8+i), type="l");
        points(x=c(0,z[9]), y=c(9+i,9+i), type="l"); points(x=c(0,z[10]), y=c(10+i,10+i), type="l"); points(x=c(0,z[11]), y=c(11+i,11+i), type="l"); points(x=c(0,z[12]), y=c(12+i,12+i),
        type="l"); points(x=c(0,z[13]), y=c(13+i,13+i), type="l"); points(x=c(0,z[14]), y=c(14+i,14+i), type="l"); points(x=c(0,z[15]), y=c(15+i,15+i), type="l"); points(x=c(0,z[16]), y=c(16+i,16+i), type="l")}


    ####Invariable RTG alleles
    points(x=subset(RTGP$V2, RTGP$V10<=0.8),y=(subset(RTGP$V1, RTGP$V10<=0.8)-0.3),pch=124, col = 314, lwd = 1, cex=0.8) #gray
    #####Variable RTG alleles
    points(x=subset(RTGVH$V2, RTGVH$V10==0.0),y=(subset(RTGVH$V1, RTGVH$V10==0.0)-0.3),pch=124, col = 554, lwd = 1, cex=0.8)#rouge
    points(x=subset(RTGVH$V2, RTGVH$V10>0.9),y=(subset(RTGVH$V1, RTGVH$V10>0.9)-0.3),pch=124, col = 490, lwd = 1, cex=0.8) #bleu
    points(x=subset(RTGVH$V2, 0.29<=RTGVH$V10 & RTGVH$V10<=0.67),y=(subset(RTGVH$V1, 0.29<=RTGVH$V10 & RTGVH$V10<=0.67)-0.3),pch=124, col = 155, lwd = 1, cex=0.8)#Noir
    #######DOC
    points(x=subset(DOCFile$V4, DOCFile$V2>1.5 & DOCFile$V1!="ref|NC"), y=(as.numeric(subset(DOCFile$V1, DOCFile$V2>1.5 & DOCFile$V1!="ref|NC"))-0.1),pch=20, col =448,lwd = 1, cex=0.8)
    points(x=subset(DOCFile$V4, DOCFile$V2<0.7 & DOCFile$V1!="ref|NC"), y=(as.numeric(subset(DOCFile$V1, DOCFile$V2<0.7 & DOCFile$V1!="ref|NC"))-0.1),pch=20, col =630, lwd = 1, cex=0.8)
    #####plot centromere
    points(x=chr_info[,4],y=(chr_info[,1]-0.3),pch=21, col = 153, cex=1.5,lwd=1)
    #####close
}
###########
#Second strain plot
###########

listdirs<-list.dirs(pathtosegD)
table <- data.frame("RTG" = character(), "Pvalue" = numeric(),"PV" = numeric(),stringsAsFactors=FALSE)
for (t in 1:length(listdirs)){
    filespath<-listdirs[t]
    print (filespath)
    
    filesnameInvAll<-list.files(filespath, pattern = c("-InvarAllelels.txt"))
    filesnameVarAll<-list.files(filespath, pattern = c("-VarAlleles.txt" ))
    RTG <-basename(filespath) #name of the RTG according to the file name
    
    RTGP<- read.table(paste(filespath,"/",filesnameInvAll,sep=""),fill=TRUE,sep="\t")
    RTGVH <- read.table(paste(filespath,"/",filesnameVarAll,sep=""),fill=TRUE,sep="\t")
    newdata<-read.table(paste(PATH3INPUTD,sep=""),fill=TRUE,sep="\t")
    newdata[,1] <- suppressWarnings(as.numeric(as.character(newdata[,1])))
    DOCFile <- newdata[with(newdata,order(newdata$V1)),]
    ####################################################################################################################################################
    #plot invisible points to label Y axe
    ####################################################################################################################################################

    ####Invariable RTG alleles
    points(x=subset(RTGP$V2, RTGP$V10<=0.8),y=(subset(RTGP$V1, RTGP$V10<=0.8)+0.1),pch=124, col = 314, lwd = 1, cex=0.8)
    
    #####Variable RTG alleles
    points(x=subset(RTGVH$V2, RTGVH$V10==0.0),y=(subset(RTGVH$V1, RTGVH$V10==0.0)+0.1),pch=124, col = 554, lwd = 1, cex=0.8)#rouge
    points(x=subset(RTGVH$V2, RTGVH$V10>0.9),y=(subset(RTGVH$V1, RTGVH$V10>0.9)+0.1),pch=124, col = 490, lwd = 1, cex=0.8) #bleu
    points(x=subset(RTGVH$V2, 0.29<=RTGVH$V10 & RTGVH$V10<=0.67),y=(subset(RTGVH$V1, 0.29<=RTGVH$V10 & RTGVH$V10<=0.67)+0.1),pch=124, col = 155, lwd = 1, cex=0.8)#Noir
    
    
    #######DOC
    points(x=subset(DOCFile$V4, DOCFile$V2>1.5 & DOCFile$V1!="ref|NC"), y=(as.numeric(subset(DOCFile$V1, DOCFile$V2>1.5 & DOCFile$V1!="ref|NC"))+0.3),pch=20, col =448 , lwd = 1, cex=0.8)
    points(x=subset(DOCFile$V4, DOCFile$V2<0.7 & DOCFile$V1!="ref|NC"), y=(as.numeric(subset(DOCFile$V1, DOCFile$V2<0.7 & DOCFile$V1!="ref|NC"))+0.3),pch=20, col =630, lwd = 1, cex=0.8)
    
    #####plot centromere
    points(x=chr_info[,4],y=(chr_info[,1]+0.1),pch=21, col = 153, cex=1.5,lwd=1)
    #####close
}

print(" PDF plot generated in ")
print(PATH2OUTPUT)

garbage <- dev.off()






