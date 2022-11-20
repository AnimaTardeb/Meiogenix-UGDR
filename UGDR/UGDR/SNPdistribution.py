#!/usr/bin/python
#####################################################################################
#####################################################################################
# -*- coding: utf8 -*-
# This file is part of R2D2-BA
# This script extracts informations from VCF files to give a human readable
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
import re
import os
import ReWr
import RR
import Ratio_variation
########################################################################
#Each allele is clustered in a specific dictionary 
#The names are changed into an adequate type for the first function
#but they are hardcoded 
########################################################################

def Extract_Alleles(Dic_File,liste_chr_name, qual, AO,DicRepReg, ratio):
    Dic_SNP, Dic_DEL, Dic_Com,Dic_Ins, Dic_Mnp   =   {},{},{},{},{}
    DicSNP, DicDEL, DicCom,DicIns, DicMnp   =   {},{},{},{},{}
    for cle , val in Dic_File.items():
        for i in range (len(liste_chr_name)):
            snp,dele,com,ins,mnp=[],[],[],[],[]
            
            if cle==liste_chr_name[i]: #on oublie la mitochondrie
                for j in range(len (val)):
                    #Multiple if's means your code would go and check all the if conditions, 
                    #where as in case of elif, if one if condition satisfies it would not check other conditions..
                    if float(val[j][1])>qual and val[j][2]=='snp' and float(val[j][3])>AO:
                        snp.append(val[j])
                    if float(val[j][1])>qual and val[j][2]=='del' and float(val[j][3])>AO:
                        dele.append(val[j])
                    if float(val[j][1])>qual and val[j][2]=='complex' and float(val[j][3])>AO:
                        com.append(val[j])
                    if float(val[j][1])>qual and val[j][2]=='mnp' and float(val[j][3])>AO:
                        mnp.append(val[j])
                    if float(val[j][1])>qual and val[j][2]=='ins' and float(val[j][3])>AO:
                        ins.append(val[j])  
                          
                DicSNP[i+1]=snp
                DicDEL[i+1]=dele
                DicCom[i+1]=com
                DicIns[i+1]=ins
                DicMnp[i+1]=mnp
    #remove repeated alleles, or alleles higher than a certain ratio
    #print len(DicRepReg)
    Dic_SNP=RMRalleles(DicSNP, DicRepReg,ratio)
    
    return Dic_SNP,DicDEL,DicCom,DicIns,DicMnp


##########################
# Eliminate the repeated snp
##########################
def RMRalleles(DicP, DicRepReg, Ratio):
        dic={}

        for cle , val in DicP.items():
            for C, V in DicRepReg.items():
                if  cle==C:
                    liste=[]
                    i=0
                    #for i in range(len(val)):
                    while i <len(val):
                        flag=False
                        for j in range (len(V)):
                            if  int(V[j][0]) < int(val[i][0]) < int(V[j][1]):
                                #print "cledicp", cle,"==",C, V[j][0],"< ",val[i][0], "<", V[j][1]
                                flag=True
                        if flag ==False and val[i] not in liste :
                            liste.append(val[i])
                        i+=1
                    dic[cle]=liste
        return dic

"""
    def RMRalleles(DicP, DicRepReg, Ratio):
    dic={}
        if len(DicRepReg)!=0:
            for cle , val in DicP.items():
                for C, V in DicRepReg.items():
                    if  cle==C:
                        liste=[]
                        i=0
                        #for i in range(len(val)):
                        while i <len(val):
                            flag=False
                            for j in range (len(V)):
                                if  int(V[j][0]) < int(val[i][0]) < int(V[j][1]) or float(val[i][8])>=Ratio:
                                    print "cledicp", cle,"==",C, V[j][0],"< ",val[i][0], "<", V[j][1]
                                    flag=True
                            if flag ==False and val[i] not in liste :
                                liste.append(val[i])
                                i+=1
                        dic[cle]=liste
                return dic
                else:
                print "############################################################################"
                print "You file of Repited Regions is empty "
                print "############################################################################"
    return DicP
"""
##################################################################
#function that extend listes
##################################################################

def extend_liste(dic_key, liste_chr, vhvp_liste, Liste):
    if dic_key in vhvp_liste and liste_chr!=[] :
        vhvp_liste.append(liste_chr)
    else:
        if vhvp_liste!=[]and liste_chr!=[]:
            Liste.append(vhvp_liste)
            vhvp_liste=[]
            vhvp_liste.extend([dic_key,liste_chr])
            
        elif liste_chr!=[]:
            vhvp_liste=[]
            vhvp_liste.extend([dic_key,liste_chr])
            
    return vhvp_liste
                                
##################################################################
##################################################################   

def Compute_ratio(dicP, dicRECb, fileout1, fileout2, Q, AO, Lname): #n =3 et k=4 nbr de copie de chromosome
    ListeP, ListeRECb= [],[]
    for cleP, valP in sorted(dicP.items()):
        for cleRECb , valRECb in sorted(dicRECb.items()):
            if cleP==cleRECb:
                LP1, LR2,i=[],[],0 #LP1=>parental-liste, LR2=>RECb-liste
                while  i <len(valP):
                    j=0
                    while j < len(valRECb):
                        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        #if cas of LOH
                        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        """if valP[i][0]!=valRECb[j][0] and int(valRECb[j-1][0])<int(valP[i][0])<int(valRECb[j][0]):
                            if  0.29<float(valP[i][-1])<=0.37:
                                lch1=[float(valP[i][-1]), float(valRECb[j][-1])]
                                lch2=[int(valP[i][0]), lch1]
                                ReWr.writefile(fileout2, cleP,valP[i][0], float(valP[i][-1]),valP[i][5],valP[i][6],valP[i][3],valP[i][7],valP[i][1],"", valP[i][0],"0.0","", "","","", "", valP[i][4],"",461, Lname)
                                LP1=extend_liste(cleP, lch1,LP1, ListeP)#Liste parental
                                LR2=extend_liste(cleRECb, lch2, LR2, ListeRECb)#Liste VH
                        """
                        if j!=0 and int(valRECb[j-1][0])<int(valP[i][0])<int(valRECb[j][0]):
                            if  0.45<=float(valP[i][-1])<=0.55:
                                lch1=[float(valP[i][-1]), float(valRECb[j][-1])]
                                lch2=[int(valP[i][0]), lch1]
                                ReWr.writefile(fileout2, cleP,valP[i][0], float(valP[i][-1]),valP[i][5],valP[i][6],valP[i][3],valP[i][7],valP[i][1],"", valP[i][0],"0.0","", "","","", "", valP[i][4],"",461, Lname)
                                LP1=extend_liste(cleP, lch1,LP1, ListeP)#Liste parental
                                LR2=extend_liste(cleRECb, lch2, LR2, ListeRECb)#Liste VH
                        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        #if two alleles are at the same position and have quality and num of reads supporting the mutation up to

                        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        if valP[i][0]==valRECb[j][0] and float(valRECb[j][1])>Q and float(valRECb[j][3])>AO and float(valP[i][1])>Q and float(valP[i][3])>AO: #ecriture des fileout1 & 2

                            lch1, lch2= Ratio_variation.ratio_variation(cleP, float(valP[i][-1]), float(valRECb[j][-1]),valP[i][0],valRECb[j][0],cleRECb,valP[i][5],valP[i][6],valP[i][3],valP[i][7],valP[i][1],valRECb[j][5],valRECb[j][6],valRECb[j][3],valRECb[j][7],valRECb[j][1], valP[i][4],valRECb[j][4],fileout1, fileout2,Lname)
                            LP1=extend_liste(cleP, lch1,LP1, ListeP)
                            LR2=extend_liste(cleRECb, lch2, LR2, ListeRECb)
                               
                        j=j+1
                

                    ##++++++++++++++++++++++++++++++++++++++++++++++
                    # fin de chromosome
                    ##++++++++++++++++++++++++++++++++++++++++++++++
                    if valP[i][0]> valRECb[-1][0] and int(valP[i][0])>int(valRECb[-1][0]):
                        #print valP[i][0], "=>",int(valRECb[-1][0])

                        if  0.45<=float(valP[i][-1])<=0.55 or 0.26<=float(valP[i][-1])<=0.38 :
                            lch1=[float(valP[i][-1]), float(valRECb[-1][-1])]
                            lch2=[int(valP[i][0]), lch1]
                            ReWr.writefile(fileout2, cleP,valP[i][0], float(valP[i][-1]),valP[i][5],valP[i][6],valP[i][3],valP[i][7],valP[i][1],"", valP[i][0],"0.0","TT>", "","","", "", valP[i][4],"",461, Lname)
                            LP1=extend_liste(cleP, lch1,LP1, ListeP)#Liste parental
                            LR2=extend_liste(cleRECb, lch2, LR2, ListeRECb)#Liste VH
                    ##++++++++++++++++++++++++++++++++++++++++++++++
                    # fin de chromosome
                    ##++++++++++++++++++++++++++++++++++++++++++++++
                    if int(valP[i][0])<int(valRECb[0][0]):
                        #print valP[i][0], "=>",int(valRECb[-1][0])
                        if  0.45<=float(valP[i][-1])<=0.55 or 0.24<=float(valP[i][-1])<=0.38:
                            lch1=[float(valP[i][-1]), float(valRECb[-1][-1])]
                            lch2=[int(valP[i][0]), lch1]
                            ReWr.writefile(fileout2, cleP,valP[i][0], float(valP[i][-1]),valP[i][5],valP[i][6],valP[i][3],valP[i][7],valP[i][1],"", valP[i][0],"0.0","rr<", "","","", "", valP[i][4],"",461, Lname)
                            LP1=extend_liste(cleP, lch1,LP1, ListeP)#Liste parental
                            LR2=extend_liste(cleRECb, lch2, LR2, ListeRECb)#Liste VH
 
                    i=i+1          
                ListeP.append(LP1)
                ListeRECb.append(LR2)
                
    return ListeP, ListeRECb
        
##############################################################################
## execution methodes
##############################################################################

def execute(filein1, filein2,AllelePfile,AlleleRECbfile, fileout1, fileout2, fileout3, fileout4, QUAL):
    #--------------------------------------------------------
    #Read file of chromosome name
    #--------------------------------------------------------
    filechrnamePath = os.path.dirname(os.path.realpath(__file__)) + "/Resources/chr_name.txt"
    chromosome_name,LChrName=ReWr.read_chrname(filechrnamePath)
    #--------------------------------------------------------
    #read file of repeated seq in the yeast
    #--------------------------------------------------------
    fileRepRegionsPath = os.path.dirname(os.path.realpath(__file__)) + "/Resources/100nt/mregions_100_annot_2011.bed"
    DicRepReg= ReWr.read_repeatedR(chromosome_name, fileRepRegionsPath)
    #--------------------------------------------------------
    #Read and convert the parental and RECb VCF files
    #--------------------------------------------------------
    DicParent=ReWr.readVCFFile(filein1)
    DicRECb=ReWr.readVCFFile(filein2)
    #--------------------------------------------------------
    #Extract Variant and Unvariant alleles > Qual and >AO
    #--------------------------------------------------------
    DicSNPP,DicDELP,DicComP,DicInsP,DicMnpP    =   Extract_Alleles(DicParent,LChrName, QUAL, 20,DicRepReg, 1)
    DicSNP,DicDEL,DicCom,DicIns,DicMnp  =   Extract_Alleles(DicRECb,LChrName, QUAL, 20,DicRepReg, 1)
    
    #--------------------------------------------------------
    print "\nReference file is proceeded. Alleles are summarized in -ReferencelAllele.txt file "
    #--------------------------------------------------------
    ReWr.distribution_alleles(DicSNPP, AllelePfile)
    ReWr.distribution_alleles(DicDELP,AllelePfile)
    ReWr.distribution_alleles(DicComP,AllelePfile)
    ReWr.distribution_alleles(DicInsP,AllelePfile)
    ReWr.distribution_alleles(DicMnpP,AllelePfile)
    #--------------------------------------------------------
    print "_______________________________________________________________\n"
    print "Recombined file (REC) is proceeded. Alleles are summarized in -RECAllele.txt file\n"
    #--------------------------------------------------------
    ReWr.distribution_alleles(DicSNP,AlleleRECbfile)
    ReWr.distribution_alleles(DicDEL,AlleleRECbfile)
    ReWr.distribution_alleles(DicCom,AlleleRECbfile)
    ReWr.distribution_alleles(DicIns,AlleleRECbfile)
    ReWr.distribution_alleles(DicMnp,AlleleRECbfile)
    #--------------------------------------------------------
    print "_______________________________________________________________\n"
    print " Computing Ratio with Qual >= ", QUAL, "and AO >= 20 (you can modify the Qual and AO in SNPdistribution.py lines 234, 235 and 259).\n"
    #--------------------------------------------------------
    #print DicSNP
    ListeP, ListeVH=Compute_ratio(DicSNPP, DicSNP, fileout1, fileout2, QUAL,20, LChrName)

    print "_______________________________________________________________\n"
    print " Extracting Recombined regions  "
    ###########
    #at least 6 alleles in one chromosome
    ###########
    ListeVHssREP=RR.RR(ListeVH, ListeP,6,fileout3, DicRepReg)
    #print len(ListeVHssREP)
    ReWr.writeRecReg(ListeVHssREP,fileout4 )

        
        
        
        
