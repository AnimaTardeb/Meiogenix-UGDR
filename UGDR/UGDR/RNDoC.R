#!/usr/local/bin/R
#####################################################################################
#####################################################################################
# -*- coding: utf8 -*-
# This file is part of UGDR-BA
# This script plot normalized depth of coverage between a cell and its reference.
# This script is part of the UGDR-BA pipeline
# Copyright (C) 2016 Bedrat Amina - Meiogenix - Institut Curie.
#RUN: Python NDoC.py -h
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
#####################################################################################
# cat script.R | /usr/local/bin/R --slave --args filein figout.pdf
#####################################################################################

#Input file : column 1 => chr name, column 2=> normelized DOC , column 3 => strain name column 4 => lokation X kb

#this script is used by NDoC.py automaticaly. however you can run it separatly

#The scrpt run with two arguments : the Xkbnormalised file and a name of a pdf output

#####################################################################################
#Start

args = commandArgs()
PATH2INPUT=args[4]
PATH2OUTPUT=args[5]


couvperkb<-read.table(paste(PATH2INPUT,sep=""),header=F)
pdf(file = PATH2OUTPUT, onefile = TRUE, width=13, height=9)
#Get the name according to the file's name
NAME <-basename(PATH2INPUT)
echant<-strsplit(NAME,"_")[[1]][1]

#####################################################################################

#oma A vector of the form c(bottom, left, top, right) giving the size of the outer margins in lines of text.
par(oma=c(.5,1.5,3,1.5))

ChrFile<-subset(couvperkb, couvperkb$V1!="chrmt")
ChrNameL=list(c(0,0), c("1","A"),c("2","B"), c("3","C"), c("4","D"),c("5", "E"),c("6", "Ff"), c("7","G"), c("8","H"), c("9","I") , c("10","J"), c("11","K"), c("12","L"), c("13","M"), c("14","N"), c("15","O"),c("16","P"))
ChrName<-c(levels(factor(ChrFile$V1)))
XVar=0
data<-data.frame()
dataifo <- data.frame()
for(i in c(1:length(ChrName))){
	for (j in c(1:length(ChrNameL))){
		if(any(ChrName[i]  == ChrNameL[[j]])){
            
			chr<-c(as.numeric(ChrNameL[[j]][1]))
			xifo<-c(subset(ChrFile$V2, ChrFile$V3=="Sk" & ChrFile$V1!="chrmt"& ChrFile$V1==ChrNameL[[j]][1]))
			xA<-c(subset(ChrFile$V2, ChrFile$V3=="Sc" & ChrFile$V1!="chrmt"& ChrFile$V1==ChrNameL[[j]][1]))
            
			newXVar=max(length(xifo), length(xA))+XVar
			yifo<-c(newXVar)
			length(xifo) <- length(xA) <-max(c(length(chr),length(xifo), length(xA),length(yifo)))
            
			data<-cbind(chr, xifo, xA, yifo)
			dataifo<-rbind(dataifo, data)
			XVar= newXVar
        }
    }
}
#####################################################################################
#Plot mean Depth of coverage for all the chromosomes in one figure

plot(x=dataifo$xifo, xlab="position (x1 kb)",ylab="", main=paste(echant, ": Normelized covrage mean", sep=""), type='l',lwd=0.5, ylim= c(-1,3), col="red")
points(x=dataifo$xA,type='l',lwd=0.5, ylim= c(-1,3), col="blue")

#Adding vertical line to delimit chromosomes and a legend

abline(v=dataifo$y, col="black", lty=1, lwd=0.5)
text(dataifo$y,3,dataifo$chr, col = "gray60", adj = c(1, -.1))
legend("bottomleft", legend=c("other","S288c"),pch=c(5,5), col = c('red','blue'), lwd=c(2.5,2.5),bg='white')

#####################################################################################
#####################################################################################
#Plot mean Depth of coverage for all the chromosomes separatly

par(oma=c(.5,1.5,3,0.5))
par(mfrow=c(4,4))
for(i in c(1:16)){
    ChrFile=subset(couvperkb,couvperkb$V1==i)
    matplot(x=subset(ChrFile$V2, ChrFile$V3=="Sc" ),xlab="position (x1 kb)",ylab="", main=paste(echant, ": Chr. ",i, sep=""), xlim= c(0, length(subset(ChrFile$V2, ChrFile$V3=="Sc"))), ylim= c(-1,3), pch=5, col="blue")
	points(x=subset(ChrFile$V2, ChrFile$V3=="Sk" ),xlab="position (x1 kb)",ylab="",xlim= c(0, length(subset(ChrFile$V2, ChrFile$V3=="Sk"))), ylim= c(-1,3),  pch=5,col="red")
    legend("topright", legend=c("other","S288c"),pch=c(5,5), col = c('red','blue'), lwd=c(.5,.5),bg='white',cex = 0.65)
}

dev.off()


#END

