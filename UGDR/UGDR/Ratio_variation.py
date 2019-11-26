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

import ReWr
import RR

##################################################################
########################################################################

ONE_ONE=[0.9, 1]
ONE_THIRD=[0.29, 0.37]
TWO_THIRD=[0.62, 0.7]
A_QUARTER=[0.21, 0.29]
HALF=[0.45, 0.55]
THREE_QUARTERS=[0.71, 0.79]

"""
ONE_ONE=[0.92, 1]
ONE_THIRD=[0.38-0.05, 0.38+0.05]
TWO_THIRD=[0.61-0.05, 0.61+0.05]
A_QUARTER=[0.28-0.05, 0.28+0.05]
HALF=[0.5-0.05, 0.5+0.05]
THREE_QUARTERS=[0.71-0.05, 0.71+0.05]
"""
########################################################################
########################################################################


def ratio_variation(cle, val1, val2,p1,p2, cle2,ro,ref,ao,alt, qual,ror,refr,aor,altr, qualr, couv, couvr, fileout1, fileout2, Lname):#val et val2 sont ratio p1, P2 position
        lch1,lsnp1, lch2,lsnp2 =[],[],[],[]
        #Invariant alleles ratio
        if A_QUARTER[0]<=val1<A_QUARTER[1] and A_QUARTER[0]<=val2<A_QUARTER[1]:#
            lsnp1=[val1,val2]
            lch1=[int(p1),lsnp1]
            ReWr.writefile(fileout1, cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr,couv, couvr, 461, Lname)
        
        if ONE_THIRD[0]<val1<=ONE_THIRD[1] and ONE_THIRD[0]<val2<=ONE_THIRD[1]:#
            lsnp1=[val1,val2]
            lch1=[int(p1),lsnp1]
            ReWr.writefile(fileout1,cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr,couv, couvr,552, Lname)
        if HALF[0]<val1<=HALF[1] and HALF[0]<val2<=HALF[1]:#
            lsnp1=[val1,val2]
            lch1=[int(p1),lsnp1]
            ReWr.writefile(fileout1, cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr,couv, couvr, 461, Lname)
        
        if TWO_THIRD[0]<val1<=TWO_THIRD[1] and TWO_THIRD[0]<val2<=TWO_THIRD[1]:#
            lsnp1=[val1,val2]
            lch1=[int(p1),lsnp1]
            ReWr.writefile(fileout1,cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr,couv, couvr,552, Lname)
        if THREE_QUARTERS[0]<val1<=THREE_QUARTERS[1] and THREE_QUARTERS[0]<val2<=THREE_QUARTERS[1]:#
            lsnp1=[val1,val2]
            lch1=[int(p1),lsnp1]
            ReWr.writefile(fileout1, cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr,couv, couvr, 461, Lname)

        #variant alleles ratio  
        
        if A_QUARTER[0]<=val1<A_QUARTER[1] and HALF[0]<=val2<=HALF[1]: # or 0.58<=val1<=0.72 and 0.26<=val2<=0.38  or 0.58<=val1<=0.72 and 0.98<=val2:
            lsnp2=[val1,val2]
            lch2=[int(p1),lsnp2]
            ReWr.writefile(fileout2,cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr,couv, couvr,614, Lname)
        
        if ONE_THIRD[0]<val1<=ONE_THIRD[1] and  TWO_THIRD[0]<val2<=TWO_THIRD[1]: # or ONE_THIRD[0]<val1<=ONE_THIRD[1]and ONE_ONE[0]<=val2:
            lsnp2=[val1,val2]
            lch2=[int(p1),lsnp2]
            ReWr.writefile(fileout2, cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr, couv, couvr,614, Lname)

        if HALF[0]<=val1<=HALF[1] and ONE_ONE[0]<=val2 or HALF[0]<=val1<=HALF[1] and THREE_QUARTERS[0]<val2<=THREE_QUARTERS[1] or HALF[0]<=val1<=HALF[1] and A_QUARTER[0]<=val2<A_QUARTER[1] or HALF[0]<=val1<=HALF[1] and ONE_THIRD[0]<val2<=ONE_THIRD[1] or HALF[0]<=val1<=HALF[1] and TWO_THIRD[0]<val2<=TWO_THIRD[1]:
            lsnp2=[val1,val2]
            lch2=[int(p1),lsnp2]
            ReWr.writefile(fileout2, cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr, couv, couvr,614, Lname)
        """if ONE_ONE<[0]<=val1<=ONE_ONE[1] and  HALF[0]<=val2<=HALF[1] :
            lsnp2=[val1,val2]
            lch2=[int(p1),lsnp2]
            ReWr.writefile(fileout2, cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr, couv, couvr,614, Lname)
        
        """
        if TWO_THIRD[0]<val1<=TWO_THIRD[1] and ONE_ONE[0]<=val2 or TWO_THIRD[0]<val1<=TWO_THIRD[1] and ONE_THIRD[0]<val2<=ONE_THIRD[1]:
            lsnp2=[val1,val2]
            lch2=[int(p1),lsnp2]
            ReWr.writefile(fileout2, cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr, couv, couvr,614, Lname)

        if THREE_QUARTERS[0]<val1<=THREE_QUARTERS[1] and ONE_ONE[0]<=val2 or THREE_QUARTERS[0]<val1<=THREE_QUARTERS[1] and HALF[0]<=val2<=HALF[1]:
            lsnp2=[val1,val2]
            lch2=[int(p1),lsnp2]
            ReWr.writefile(fileout2, cle,p1, val1, ro,ref,ao,alt, qual, cle2, p2, val2,ror,refr,aor,altr, qualr, couv, couvr,614, Lname)
        
        else:
            pass


        return lch1, lch2

##################################################################
################################################################## 


