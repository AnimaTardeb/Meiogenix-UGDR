#!/usr/bin/python 
#####################################################################################
#####################################################################################
# -*- coding: utf8 -*-
# This file is part of UGDR-BA
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

import os, sys, time,optparse
import SNPdistribution

def main(argv):
    parser = optparse.OptionParser(usage='\033[1m' +"\npython %prog [options] REC_vcf_file [options] Reference_vcf_file [options] output_Repository \n"+"\033[0;0m", version="%prog 0.75",
        description='\033[1m' +"UGDR for analyzing alleles variation and identifing regions of recombination in yeast. This script requires two vcf files format to be compared and to extract the regions of recombination.")
    parser.add_option("-i", "--REC_rep", default=None, action="store_true",
        help="Folder of more than one recombinant to test ")
    parser.add_option("-I", "--REC_file", default=None, action="store_true",
        help="One recombinant (VCF file) to test (REC)")
    parser.add_option("-j", "--par_file", default=False, action="store_true",
        help="Reference (VCF) file")
    parser.add_option("-o", "--out_dir", default=False, action="store_true",
        help="Results folder ")
    
    parser.add_option("-c","--DPQUAL", dest ="num",  help="Filter applied on DP for 80X use 200"+"\033[0;0m" +"\n\n")

    (options, args) = parser.parse_args()
    #to check the -c argument
    
    if (options.num == None):
            print (parser.usage)
            exit(0)
    else:
            number = options.num

    if len(args) != 3:
        print "\n"
        parser.print_help()
        exit(1)   
    RECr      =   args[0]
    parfile     =   args[1]
    outputDir   =   args[2]
    RECfile     =   options.REC_file
    QUAL   =   number
    return  RECr, parfile, outputDir, RECfile, int(QUAL)
#================================================================
if __name__ == "__main__":
    
    RECrep, parfile, outputDir, RECfile, QUAL = main(sys.argv[1:])
    #================================================================
    # In the case you use repository of VCF files -i
    #================================================================    
    if RECfile==None:
        files = os.listdir(RECrep)
        for i in files :
            filename=i.split(".")[0] 
            OutDirname="Results-"+filename
            OPF= os.listdir(outputDir)
            flag=False
            for dir in OPF:
                if dir== OutDirname:
                    flag=True
            if flag==True:
                print "\n########################################################################"
                print "#    "+OutDirname+" directory  already exist                            "
                print "#    "+OutDirname+" Crushed                                             "
                print "#    New Results directory Created                                      "  
                print "########################################################################"
            else:
                os.makedirs(outputDir+"/"+OutDirname+"/", mode=0777)        
                print "########################################################################"
                print "#        "+OutDirname+" directory Created                              #" 
                print "########################################################################"

            startTime = time.time()   
            try:
                AlleleL1file    =     open (outputDir +"/"+OutDirname+"/"+filename+"-ReferenceAllele.txt", "w")
                AlleleRECfile   =     open (outputDir +"/"+OutDirname+"/"+filename+"-RECAlleles.txt", "w")
        
                Res1file    =    open (outputDir+"/"+OutDirname+"/"+filename+"-InvarAllelels.txt", "w")
                Res2file    =    open (outputDir+"/"+OutDirname+"/"+filename+"-VarAlleles.txt", "w")
                Res3file    =    open (outputDir +"/"+OutDirname+"/"+filename+"-RR.txt", "w")
                Res4file    =    open (outputDir +"/"+OutDirname+"/"+filename+"-NEWFILE.txt", "w")
                filein      =    open (parfile,"r")
                filein2     =    open (RECrep+"/"+i, "r")
        
                SNPdistribution.execute(filein, filein2,AlleleL1file,AlleleRECfile, Res1file, Res2file, Res3file, Res4file, QUAL)
            
                ResRep= outputDir +"/"+OutDirname+"/"
                print ResRep
                ResRepPdf=ResRep+filename+".pdf"
                sortie=os.popen("cat Alleles-rep-onefile.r | /usr/local/bin/R --slave --args "+ResRep+" " + ResRepPdf, "r").read()
                print sortie
            
            
            except IOError:
                print "\n"
                print "######################################################################"
                print "#                        Incorrect Path File"
                print "#    Please check your paths, these can solve half of the problems"
                print "######################################################################"
                sys.exit(0)     
                #================================================================ 
            fin=time.time()
            #print "Files created in ." ,fin-startTime, "secondes", i
            Res1file.close()
            Res2file.close()
            Res3file.close()
            Res4file.close()
            
            filein.close()
            filein2.close()
            CC=outputDir +"/"+OutDirname+"/"+filename+"-RR.txt"
            print "merde",CC
            os.remove(CC)

    #================================================================
    # Case of using one VCF file -I
    #================================================================
    else:
        fname=RECrep.split("/")[-1]
        filename=fname.split(".")[0]

        OutDirname="Results-"+filename
        OPF= os.listdir(outputDir)
        flag=False
        for dir in OPF:
            if dir== OutDirname:
                flag=True
        if flag==True:
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
        startTime = time.time()   
        try:
            AlleleL1file    =     open (outputDir +"/"+OutDirname+"/"+filename+"-ReferenceAllele.txt", "w")
            AlleleRECfile   =     open (outputDir +"/"+OutDirname+"/"+filename+"-RECAlleles.txt", "w")
            
            Res1file    =    open (outputDir+"/"+OutDirname+"/"+filename+"-InvarAllelels.txt", "w")
            Res2file    =    open (outputDir+"/"+OutDirname+"/"+filename+"-VarAlleles.txt", "w")
            Res3file    =    open (outputDir +"/"+OutDirname+"/"+filename+"-RR.txt", "w")
            Res4file    =    open (outputDir +"/"+OutDirname+"/"+filename+"-NEWFILE.txt", "w")
            filein      =    open(parfile,"r")
            filein2     =    open(RECrep, "r")

            SNPdistribution.execute(filein, filein2,AlleleL1file,AlleleRECfile, Res1file, Res2file, Res3file, Res4file, QUAL)
            Res1file.close()
            Res2file.close()
            Res3file.close()
            Res4file.close()
            #os.remove(outputDir +"/"+OutDirname+"/"+filename+"-RR.txt")
            ResRep= outputDir +"/"+OutDirname+"/"
            print ResRep
            ResRepPdf=ResRep+filename+".pdf"
            os.system("cat Alleles-rep-onefile.r | /usr/local/bin/R --slave --args "+ResRep+" " + ResRepPdf)
            CC=outputDir +"/"+OutDirname+"/"+filename+"-RR.txt"
            #os.remove(CC)
            os.remove(outputDir +"/"+OutDirname+"/"+filename+"-NEWFILE.txt")

        except IOError:
            print "\n"
            print "######################################################################"
            print "#                        Incorrect File\n"
            print "#    Please check your paths, these can solve half of the problems"
            print "######################################################################"
            sys.exit(0)     
            #================================================================ 
        filein.close() 
        fin=time.time()
        #print "Files created in ." ,fin-startTime, "secondes", i
        filein2.close()


   
