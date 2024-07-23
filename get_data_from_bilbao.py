#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 12:11:29 2024

@author: tony
"""


import requests
import pandas
import re
import numpy as np
import math
import time

##This script contains the following functions for the space groups 1-230 from Bilbao:
## get maximal subgroups
## get the list of generators (coordinate convention = Bilbao = a preferred convention in ITC = our LSM paper)
## get the table of Wyckoff positions

"""
Maximalsubgroups230 = []

for i in range(1,231):
    print(i)

    url = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-lxi?client=wp&way=up&path=@maxsub&gnum=%d"%i

    html = requests.get(url).text
    
    a = [m.start() for m in re.finditer('<table', html)][1]
    b = [m.start() for m in re.finditer('</table>', html)][1]
    
    #print(a,b)
    
    #print(html[a:b+8])
    
    tb0 = html[a:b+8]
    
    tb1 = re.sub('<.*?>', '', tb0)
    tb2 = tb1.split('N IT number HM symbol Index Type Transformations')[1]
    tb3 = tb2.split('show..')[:-1]
    
    msubgrp_for_i = []
    while(len(tb3) >=1 and tb3[-1][-2] == 'k'):
        tb3.remove(tb3[-1])
    for item in tb3:
        msubgrp_for_i.append((item.split()[1],int(item.split()[3])))
    print(msubgrp_for_i)
    Maximalsubgroups230.append(msubgrp_for_i)
"""  

"""
Gens230 = []

for no in range(1,231):
    print(no)

    url = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen?gnum=%d"%no

    html = requests.get(url).text
    
    a = [m.start() for m in re.finditer('<td align="center">&nbsp; &nbsp; ', html)]
    b = [m.start() for m in re.finditer('&nbsp; &nbsp; </td>', html)]
    
    gens_for_i = []
    
    for i in range(len(a)):
    
        #print(html[a[i]+33:b[i]])
    
        tb0 = html[a[i]+33:b[i]]

        gens_for_i.append(tb0)
        
    print(gens_for_i)

    Gens230.append(gens_for_i)
"""


Wyckoff230 = []

#For each group no, the Wyckoff table is such that, in each row of the table, 
#the first three columns are "Multiplicity" and "Wyckoff letter" and "Site symmetry", then followed by
#the fourth column which is "Coordinates".
#We first use "<table><tr><td>" to obtain the "Site symmetry" in each row, and then back track 
#the first three columns that are  between the first "<td align=center>" of the row and "<table><tr><td>".
for no in range(1,231):
    print(no)

    url = "https://www.cryst.ehu.es/cgi-bin/cryst/programs//nph-wp-list?gnum=%d&what=wpos"%no

    html = requests.get(url).text
    
    a0 = [m.start() for m in re.finditer('Wyckoff Positions of Group',html)][0]
    b0 = [m.start() for m in re.finditer('<h3>Wyckoff position and site symmetry group of',html)][0]
    
    
    html0 = html[a0:b0]
    
    a = [m.start() for m in re.finditer('<td align=center>',html0)]
    #a = [m.start() for m in re.finditer('<table><tr><td>', html0)]
    b = [m.start() for m in re.finditer('</td></tr></table>', html0)]
    
    wyckoff_for_no = []
    
    for i in range(len(b)): #b runs over the rows of the Wyckoff table
    
        #print(html[a[i]+33:b[i]])

        tb0 = html0[a[0]:b[i]]   #this is the scope of the row i
        wyckoff_for_i = []
        
        a1 = [m.start() for m in re.finditer('<td align=center>',tb0)]
        b1 = [m.start() for m in re.finditer('</td>',tb0)]
    
        #for j in range(len(a1)):
        wyckoff_for_i.append((tb0[a1[0]+17:b1[0]])+(tb0[a1[1]+17:b1[1]]))
        wyckoff_for_i.append((tb0[a1[2]+17:b1[2]]))
        
        aa = [mm.start() for mm in re.finditer('&standard=1">',tb0)]
        bb = [mm.start() for mm in re.finditer('</a></nobr>',tb0)]
        
        for j in range(len(aa)):
            wyckoff_for_i.append(tb0[aa[j]+13:bb[j]])

        wyckoff_for_no.append(wyckoff_for_i)
        
        
        #find the start of the next row, always counted in html0:
        a = [(b[i]+10) + m.start() for m in re.finditer('<td align=center>',html0[(b[i]+10):])]
        
        
    print(wyckoff_for_no)

    Wyckoff230.append(wyckoff_for_no)
    
    
    
    
    
    
    
    
    
    
    