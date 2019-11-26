#!/usr/bin/python 
#####################################################################################
#####################################################################################
# -*- coding: utf8 -*-
# This file is part of R2D2-BA
# This script extract informations from VCF files to give a human readable 
# files resuming alleles variations, and difining region of recombinaition
# Copyright (C) 2016 Bedrat Amina - Meiogenix - Institut Curie
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
import re, os
############################################################################################################    
#Lis les fichiers VCF et les transforme en dictionnaire avec toutes les information necessaire     
############################################################################################################
def readVCFFile(Filein):
    Liste_SEQ, Dic_File = [] , {}
    for line in Filein:
        seq_comment  =  re.compile(r"^\#") 
        iterator , Liste_l, i =  seq_comment.finditer(line) , [], 0
        try:   
            while seq_comment :
                n  =  iterator.next()
                pass
        except StopIteration: pass
        if not seq_comment.search(line):
            header      =   line.split("\n")[0]
            fields      =   header.split("\t")
            filter_info =   fields[7].split(";")
            dico        =   {}
            for I in range(len(filter_info)-1):
                info=filter_info[I].split("=")
                if len(info)>=2:
                    dico[info[0]]=info[1]
                else:
                    dico[info[0]]="100000"
                    
            ao  =   dico["AO"].split(",")
            type=   dico["TYPE"].split(",")
            if len(ao)==1:
                AO      =       ao[0]
                DP      =       dico["DP"]
                RO      =       dico["RO"]
                Types   =       type[0]#mettre le tt dans un dico pour une cle de nom du chromosome
                if fields[0] in  Dic_File.keys():
                    try:
                        ratio_value =   "%.3f" % round(float(AO)/int(DP), 3)
                    except ZeroDivisionError:
                        print "############################################################################"
                        print("ERROR: DP is equal to 0.")
                        print "############################################################################"
                        os._exit(1)
                    Liste_l.extend([fields[1],fields[5],Types, AO, DP,RO,fields[3], fields[4],ratio_value])
                    Dic_File[fields[0]].append(Liste_l)
                else:
                    Liste_SEQ  =   []
                    Dic_File[fields[0]] =   Liste_SEQ
                    try:
                        ratio_value =   "%.3f" % round(float(AO)/int(DP), 3)
                    except ZeroDivisionError:
                        print "############################################################################"
                        print "ERROR: DP is equal to 0."
                        print "############################################################################"
                        os._exit(1)

                    Liste_l.extend([fields[1],fields[5],Types, AO, DP, RO,fields[3], fields[4],ratio_value])
                    Dic_File[fields[0]].append(Liste_l)
            else:
                pass
            i = i+1
    return Dic_File

############################################################################
# Read the alleles stocked in a dictionary and write them in one file
############################################################################ 
"""
def distribution_alleles(Dic_alleles, file_alleles):
    liste = []
    try:
        for cle, val in Dic_alleles.items():
            for j in range(len (val)):
                liste_l = [cle, val[j][2], val[j][0], val[j][1], val[j][4], val[j][8], val[j][5],val[j][6], val[j][3],val[j][7]]
                liste.append(liste_l)
        print val[j][2]," : ",len(liste) #le nombre d allele par type
        #title = "chr\ttype\tPos\tQual\tCov\tRatio\tRO\tref\tAO\tAlt\n"
        #file.write(title)
        for i in range (len(liste)):
            for j in range (len (liste[i])):
                file_alleles.write(str(liste[i][j]))
                file_alleles.write("\t")
            file_alleles.write("\n")
    except:
        print "############################################################################"
        print "ERROR: Make sure your 'chr_name.txt' and 'VCF' files are correct"
        print "############################################################################"
        
        os._exit(1)
"""

def distribution_alleles(Dic_alleles, file_alleles):
    liste = []
    for cle, val in Dic_alleles.items():
        for j in range(len (val)):
            liste_l = [cle, val[j][2], val[j][0], val[j][1], val[j][4], val[j][8], val[j][5],val[j][6], val[j][3],val[j][7]]
            liste.append(liste_l)
    print val[j][2]," : ",len(liste) #le nombre d allele par type
    #title = "chr\ttype\tPos\tQual\tCov\tRatio\tRO\tref\tAO\tAlt\n"
    #file.write(title)
    for i in range (len(liste)):
        for j in range (len (liste[i])):
            file_alleles.write(str(liste[i][j]))
            file_alleles.write("\t")
        file_alleles.write("\n")
##################################################################
# Extract the alleles present in the repeated regions
##################################################################

def read_repeatedR(liste,filein):
    filein = open(filein,'r')
    dic = {}
    
    for lines in filein:
        List, chr, start, end = [], lines.split(" ")[0],lines.split(" ")[1],lines.split(" ")[2]
        List.extend([start,end])
        for i in range (len(liste)):
            if chr == liste[i]:
                chrname = i+1
                if chrname not in dic.keys():
                    GL = []
                    GL.append(List)
                    dic[chrname] = GL
                else:
                    GL.append(List)
                    dic[chrname] = GL
    filein.close()
    return dic
##################################################################
#Read Chromosome name
##################################################################
def read_chrname(filein):
    chromosome_name=[]
    LnomC=[]
    file = open(filein,"r")
    for lines in file:
        ngs_name, chr_name= lines.split("\t")[0], lines.split("\t")[1]
        LnomC.append(ngs_name)
        chromosome_name.append(chr_name)
    return chromosome_name, LnomC

##################################################################
# Write all the alleles informations in one file 
##################################################################         
def writefile(fileout, chro, posi, ratio, ro, ref, ao, alt, qal, chrortg, posirtg, ratiortg , ror, refr, aor, altr,qalr,couv, couvr, color, Lname):
    #title = "Chr\tPosition\tRatio\tRO\tRef-All\tAO\tAlt-All\tQual \tChr\tPosition\tRatio\tRO\tRef-All\tAO\tAlt-all\tQual\tColor\n"
    LinetoWrite =  str(chro)+"\t"+ str(posi) +"\t"+ str(ratio) +"\t"+ str(ro)+str(ref)+"\t"+ str(ao)+str(alt)+"\t"+str(couv)+str("X")+"\t" +str(qal)+"\t"+ str(chrortg) +"\t"+ str(posirtg) +"\t"+ str(ratiortg)+"\t"+ str(ror)+str(refr)+"\t"+ str(aor)+str(altr)+"\t"+str(couvr)+str("X")+"\t"+str(qalr)+"\n"
    fileout.write(LinetoWrite)

##################################################################
#Write the
##################################################################
def writeRecReg(liste, Infile):
    for i in range(len(liste)):
        for j in range(1,len(liste[i])):
            LinetoWrite = str(liste[i][0])+"\t"+str(liste[i][j][0])+"\t"+str(liste[i][j][1][0])+"\t"+str(liste[i][j][1][1])+"\n"
            Infile.write(LinetoWrite)





























