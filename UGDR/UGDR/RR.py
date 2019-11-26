#!/usr/bin/python
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



def removeM(Dic, cle, position):
    flag =False
    for key, values in Dic.items():
        if key==cle:
            for i in range (len (values)):
                if  int(values[i][0])<int(position)<int(values[i][1]):
                    flag=True
                
    return flag

def RR(ListeRTG, ListeP, N, fileout, dicM):
    LRTG=[]
    NL, seq=[], ""
    for k in range (len(ListeRTG)):
        if len(ListeRTG[k])!=0:
            LRTG.append([ListeRTG[k][0]])
            for W in range(1,len(ListeRTG[k])):
                    flag=removeM(dicM, ListeRTG[k][0], ListeRTG[k][W][0])
                    if flag==False:
                        LRTG[k].extend([ListeRTG[k][W]])
        else:
            LRTG.append([ListeRTG[k]])       #print LRTG

    for elements in range(len(LRTG)):
        if len(LRTG[elements])-1>=N:
            NL.append(LRTG[elements])

    for I in range (len(NL)):
        for i in range(len(ListeP)):
            if len(ListeP[i])!=0 and len(NL[I])!=0:
                if NL[I][0]==ListeP[i][0]:
                    K,L, l ,yy=1,[], [], 1
                    if NL[I][0]!=9:
                        for J in range(1,len(NL[I])-1):
                            if NL[I][J][1][0]!= NL[I][J][1][1]:
                                score    = 0
                                for j in range(K, len(ListeP[i])-1):
                                    x=0
                                    if NL[I][J][0] <= ListeP[i][j][0] <= NL[I][J+1][0]:
                                        if  0.26<=ListeP[i][j][1][0]<=0.38 and 0.26<=ListeP[i][j][1][1]<=0.38:  
                                            x=1
                                        elif 0.58<=ListeP[i][j][1][0]<=0.72 and 0.58<=ListeP[i][j][1][1]<=0.72:
                                            x=2 
                                        score +=x
                                        K=j
                                if score <= int(33):
                                    L.insert(NL[I][J][0],NL[I][J+1][0])
                                    yy+=1
                                    if J+1== len(NL[I])-1:
                                        l=list(set(L))
                                        l.sort()
                                        if l!=0 and len(l)>=5:
                                            fileout.write(str(NL[I][0]))
                                            fileout.write("\t")
                                            seq=str(l[0])+ "\t"+str(l[-1])+"\t"+str(int(l[-1])-int(l[0]))+"\t"+str(len(l))+"\t"+str(score)+"\n"
                                            fileout.write(seq)
                                        L=[]
                                        yy=1    
                                else:
                                    l=list(set(L))
                                    l.sort()
                                    if l!=0 and len(l)>=5:
                                        fileout.write(str(NL[I][0]))
                                        fileout.write("\t")
                                        seq=str(l[0])+ "\t"+str(l[-1])+"\t"+str(int(l[-1])-int(l[0]))+"\t"+str(len(l))+"\t"+str(score)+"\n"
                                        fileout.write(seq)
                                    L=[]
                                    yy=1
                    if NL[I][0]==9:
                        for J in range(1,len(NL[I])-1):
                            if NL[I][J][1][0]!= NL[I][J][1][1]:
                                score    = 0
                                for j in range(K, len(ListeP[i])-1):
                                    x=0
                                    if NL[I][J][0] <= ListeP[i][j][0] <= NL[I][J+1][0]:
                                        if  0.2<=ListeP[i][j][1][0]<=0.31 and 0.2<=ListeP[i][j][1][1]<=0.31:  
                                            x=0.5
                                        elif 0.44<=ListeP[i][j][1][0]<=0.55 and 0.44<=ListeP[i][j][1][1]<=0.55:
                                            x=1
                                        elif 0.69<=ListeP[i][j][1][0]<=0.79 and 0.69<=ListeP[i][j][1][1]<=0.79:
                                            x=2
                                        score +=x
                                        K=j
        
                                if score <= int(10):
                                    L.insert(NL[I][J][0],NL[I][J+1][0])
                                    yy+=1
                                    if J+1== len(NL[I])-1:
                                        l=list(set(L))
                                        l.sort()
                                        if l!=0 and len(l)>=5:
                                            fileout.write(str(NL[I][0]))
                                            fileout.write("\t")
                                            seq=str(l[0])+ "\t"+str(l[-1])+"\t"+str(int(l[-1])-int(l[0]))+"\t"+str(len(l))+"\t"+str(score)+"\n"
                                            fileout.write(seq)
                                        L=[]
                                        yy=1    
                                else:
                                    l=list(set(L))
                                    l.sort()
                                    if l!=0 and len(l)>=5:
                                        fileout.write(str(NL[I][0]))
                                        fileout.write("\t")
                                        seq=str(l[0])+ "\t"+str(l[-1])+"\t"+str(int(l[-1])-int(l[0]))+"\t"+str(len(l))+"\t"+str(score)+"\n"
                                        fileout.write(seq)
                                    L=[]
                                    yy=1

    return NL


