#!/usr/bin/python
#####################################################################################
#####################################################################################
# -*- coding: utf8 -*-
# This file is part of UGDR-BA
# This script calculate normalized depth of coverage between a cell and its reference.
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

import os, sys,optparse
import re
from operator import truediv
import numpy as np



def ReadDCOVFile(Filein, kb, fileout):
    Dic_chro, liste_covrange, liste_meancov =   dict(), [], []
    couverture, length_file =   0.0, 0
    for line in Filein:
        """
        ########################IMPORTANT#########################
        #compiled on the basis that each list start with "chr"
        ########################IMPORTANT#########################
        """
        num  =  re.compile(r"^['chr']")
        iterator =  num.finditer(line)
        try:
            while num:
                n           =   iterator.next()
                chroname    =   line.split(":")[0]
                cov         =   int(line.split("\t")[1])
                couverture +=   cov
                length_file+=   1
                if chroname in Dic_chro.keys():
                    liste_covrange.append(cov)
                    Dic_chro[chroname]  =   liste_covrange
                else:
                    liste_covrange      =   []
                    Dic_chro[chroname]  =   []
                    liste_covrange.append(cov)
                    Dic_chro[chroname]  =   liste_covrange
    
        except StopIteration: pass
        if not num.search(line):
            pass

    TTavrgCov=couverture/length_file
    print "Number of base :", length_file, "\n Total Mean of Coverage: ", TTavrgCov
    
    for key, value in Dic_chro.items():
        for i in range(1, len(value), kb ):
            liste_intermediaire=value[i:i+kb]
            mean    =   sum(liste_intermediaire, 0.0)/len(liste_intermediaire)
            Lline   =   str(key)+ "\t"+str(mean)+"\n"
            fileout.write(Lline)
            liste_meancov.append(mean)
    return Dic_chro, TTavrgCov

######################################################################
#------------------------------------------------------
# Normalize compared to the mean-voc of ALL the sample
#------------------------------------------------------
######################################################################
"""
    ########################IMPORTANT#########################
    # Two methods:
    # Mean_dofcov_bp: in the case you analyze all the bases
    # Mean_dofcov_kb: in the case you analyse bases per x Kb
    ########################IMPORTANT#########################
    """
def Mean_dofcov_kb(dicp, dicrtg, covmp, covmrtg, fileout, kb):
    oldsettings = np.seterr(divide='ignore')#waring rised due to division on zero

    newdic  =   {}
    """------------------------------------------------------
        #  Normalization
       ------------------------------------------------------
    """
    for keyp, valuep in dicp.items():
        for i in range(len(valuep)):
            valuep[i]   =   valuep[i]/covmp
    for key, value in dicrtg.items():
        for i in range(len(value)):
            value[i]    =   value[i]/covmrtg
    """------------------------------------------------------
        #  Writing to a file with "write_meanCov_kb" method
       ------------------------------------------------------
    """
    for j in range(len(dicp.values())):
        oldsettings = np.seterr(divide='ignore',invalid='ignore')#waring rised due to division on zero
        value_normalized    =   np.array(dicrtg.values()[j])/np.array(dicp.values()[j])
        newdic.update({dicp.keys()[j]:value_normalized})
    dictionnaire    =   write_meanCov_kb(newdic, fileout, kb)
    return dictionnaire

####################################################################

def Mean_dofcov_bp(dicp, dicrtg, covmp, covmrtg, fileout):
    oldsettings = np.seterr(divide='ignore')#waring rised due to division on zero
    newdic  =   {}
    """------------------------------------------------------
    #  Normalization
    ------------------------------------------------------
    """
    for keyp, valuep in dicp.items():
        for i in range(len(valuep)):
            valuep[i]   =   valuep[i]/covmp

    for key, value in dicrtg.items():
        for i in range(len(value)):
            value[i]    =   value[i]/covmrtg

    """------------------------------------------------------
        #  Writing to a file with "write_meanCov_bp" method
        ------------------------------------------------------
    """

    for j in range(len(dicp.values())):
        np.seterr(divide='ignore', invalid='ignore')
        value_normalized    =   np.array(dicrtg.values()[j])/np.array(dicp.values()[j])
        newdic.update({dicp.keys()[j]:value_normalized})
        write_meanCov_bp(newdic, fileout)
    return newdic

######################################################################
#------------------------------------------------------
# Normalize compared to the mean of Xbp.
######################################################################
# Convert Chromosome name according the the two lists given as an agruments
#----------------------------------------------------------------------------
################################################################################

def convertDic(Dic_parental,nameSC ,nameSK, S1, S2):
    for key, value in Dic_parental.items():
        if key in nameSC:
            ik=str(nameSC.index(key)+1)+ "_"+S1
            Dic_parental[ik]= value
            Dic_parental.pop(key)
        elif key in nameSK:
            ik=str(nameSK.index(key)+1)+ "_"+S2
            Dic_parental[ik]= value
            Dic_parental.pop(key)
        else:
            pass
    return Dic_parental
#-------------------------------------
# Writing the files
#-------------------------------------
"""
########################IMPORTANT#########################
# Two methods:
# Write_meanCov_bp: in the case you analyze all the bases
# Write_meanCov_kb: in the case you analyse bases per x Kb
########################IMPORTANT#########################
"""
def write_meanCov_bp(newdic, fileout):
    
    for key, value in newdic.items():
        keyi    =   key.split("_")
        for i in range(1, len(value)):
            Lline   =   str(keyi[0])+ "\t"+str(value[i])+"\t"+ str(keyi[1])+"\n"
            fileout.write(Lline)
#------------------------------------

def write_meanCov_kb(newdic, fileout, kb):
    newdic_norm =   {}
    
    for key, value in sorted(newdic.items()):
        keyi    =   key.split("_")
        liste   =   []
        for i in range(1, len(value)-1,kb):
            liste_intermediaire=value[i:i+kb]
            mean    =   sum(liste_intermediaire, 0.0)/len(liste_intermediaire)
            Lline   =   str(keyi[0])+ "\t"+str(mean)+"\t"+ str(keyi[1])+"\t"+ str(i)+ "\n"
            fileout.write(Lline)
            liste.append(mean)
        newdic_norm.update({key:liste})

    return newdic_norm

##################################################################
#Read chromosome fine names
##################################################################
def read_chrname(filein):
    chromosome_name =   []
    nameS1  =   []
    nameS2  =   []
    file    =   open(filein,"r")
    
    for lines in file:
        ngs_name, chr_name, ngs2_name = lines.split("\t")[0], lines.split("\t")[1], lines.split("\t")[3]
        S1, S2  =   lines.split("\t")[2], lines.split("\t")[5]
        nameS1.append(ngs_name)
        nameS2.append(ngs2_name)
        chromosome_name.append(chr_name)
    
    return chromosome_name, nameS1, nameS2, S1, S2.split("\n")[0]
######################################################################
######################################################################

def main(argv):
    #Help
    parser = optparse.OptionParser(usage='\033[1m' +"\npython %prog [options] Referenc_Depth_of_coverage_file [options] RTG_Depth_of_coverage_file [options] output_Repository \n"+ "\033[0;0m", version="%prog 1.0",description="Normalized Depth of Coverage (NDoC): To analyze alleles coverage and identify regions of recombination. the inputs are two depth of coverage files generated for UGDR. The per chromosome depth profile is generated at the end of the analysis.")
    
    parser.add_option("-i", "--ref_f", default=None, action="store_true",help="Depth of coverage file of the reference strain")
    parser.add_option("-j", "--RTG_f", default=None, action="store_true",help="Depth of coverage file of Recombined strain ")
    parser.add_option("-o", "--out_dir", default=False, action="store_true",help="Results repository")
    
    (opts, args) = parser.parse_args() 
    if len(args) != 3:
        print "\n"
        parser.print_help()
        exit(1)
    parfile     =   args[0]
    rtgfile     =   args[1]
    outputDir   =   args[2]
    return  parfile,rtgfile, outputDir
#================================================================
if __name__ == "__main__":
    parfile, RTGdoc, outputDir= main(sys.argv[1:])
    
    #================================================================
    # These spliting to get file name is specific for the file generated
    # by Glocal or the file name should be : nameoFinterest.anythink.anythink
    #================================================================
    
    rtgname =   RTGdoc.split("/")[-1]
    filename=   rtgname.split(".")[0]
    
    pname   =   parfile.split("/")[-1]
    Pfilename   =  pname.split(".")[0]
    
    OutDirname  =   "DepthOfCov-Results-"+filename
    OPF         =   os.listdir(outputDir)
    flag    =   False
    for dir in OPF:
        if dir  ==  OutDirname:
            flag    =   True
    if flag ==  True:
        print "#######################################################################"
        print "#    "+OutDirname+" directory  already exist                           "
        print "#    "+OutDirname+" Crushed                                            "
        print "#    New Results directory Created                                     "
        print "#######################################################################"
    else:
        os.makedirs(outputDir+"/"+OutDirname+"/", mode=0777)
        print "########################################################################"
        print "#        "+OutDirname+" directory Created                              #"
        print "########################################################################"
    #================================================================
    #================================================================
    try:
        Res1file    =    open (outputDir+"/"+OutDirname+"/"+Pfilename+".txt", "w")
        Res2file    =    open (outputDir+"/"+OutDirname+"/"+filename+".txt", "w")
        Res3file    =    open (outputDir +"/"+OutDirname+"/"+filename+"__normalized1kb.txt", "w")
        
        filein      =    open(parfile,"r")
        filein2     =    open(RTGdoc, "r")
        
        ##################################################################
        filechrnamePath = os.path.dirname(os.path.realpath(__file__)) + "/chr_name.txt"
        chromosome_name,nameSC, nameSK, S1, S2 = read_chrname(filechrnamePath)
        ##################################################################
        print "Reference DOC file: ",Pfilename
        Dic_parental_one, covparental   =   ReadDCOVFile(filein, 1000, Res1file)
        print "\nRecombined DOC file: ",filename
        Dic_RTG_one, covRTG             =   ReadDCOVFile(filein2, 1000, Res2file)

        #Compile chromosome name
        Dic_parental    =   convertDic(Dic_parental_one,nameSC ,nameSK, S1, S2)
        Dic_RTG         =   convertDic(Dic_RTG_one,nameSC ,nameSK, S1, S2)


        NEWDIC  = Mean_dofcov_kb(Dic_parental, Dic_RTG,covparental, covRTG , Res3file, 1000)
        Res3file.close()
        #Representation
        
        ResRep      =   outputDir +"/"+OutDirname+"/"+filename+"__normalized1kb.txt"
        ResRepPdf   =   outputDir +"/"+OutDirname+"/"+filename+".pdf"
        #sortie  =os.popen( "cat RNDoC.R | /usr/local/bin/R --slave --args "+ResRep+" "+ResRepPdf, "r").read()
        os.system("cat RNDoC.R | /usr/local/bin/R --slave --args "+ResRep+" "+ResRepPdf)

    except IOError:
        print "\n"
        print "######################################################################"
        print "#                        Incorrect File\n"
        print "######################################################################"
        sys.exit(0)
    #================================================================

    Res1file.close()
    Res2file.close()
    Res3file.close()

    filein.close()
    filein2.close()


######################################################################
######################################################################
