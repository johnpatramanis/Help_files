import numpy as np
import math
import os
import argparse
import time
import re
import random
import sys
import seaborn
import scipy.special

DATA=[]
for chr in range(1,22):
    FILE=open('refined_output{}.ibd'.format(chr),'r')
    for line in FILE:
        line=line.strip().split()
        DATA.append(line)
    
#Here we count samples per pop
SAMPLES={}
IND=[]
for seg in DATA:
    IND.append(seg[0])
    IND.append(seg[2])

UNIQUEIND=list(set(IND))


for kappa in range(0,len(UNIQUEIND)):
    myind=(UNIQUEIND[kappa]).split('_')
    
    if len(myind)>2:
        UNIQUEIND[kappa]='_'.join(myind[0:2])
    if len(myind)<=2:
        UNIQUEIND[kappa]=myind[0]

for ZED in UNIQUEIND:
    if ZED in SAMPLES.keys():
        SAMPLES[ZED]+=1
    else:
        SAMPLES[ZED]=1
#print(SAMPLES)


###################keep only poplabels

for seg in range(0,len(DATA)):

    seg0=DATA[seg][0].split('_')
    seg2=DATA[seg][2].split('_')
    

    if len(seg0)>2:
        DATA[seg][0]='_'.join(seg0[0:2])
    if len(seg0)<=2:
        DATA[seg][0]=seg0[0]
        
        
    if len(seg2)>2:
        DATA[seg][2]='_'.join(seg2[0:2])
    if len(seg2)<=2:
        DATA[seg][2]=seg2[0]
        

        
#creation of overal dictionary
IBDS={}
for seg in DATA:
    if seg[0] in IBDS.keys():
        pass
    if seg[0] not in IBDS.keys():
        IBDS[seg[0]]={}


    if seg[2] in IBDS.keys():
        pass
    if seg[2] not in IBDS.keys():
        IBDS[seg[2]]={}



#creation ,appending of inner dictionaries
#prepei na diairw me pop1*pop2
for seg in DATA:
    if float(seg[-1])>=5:
        if seg[2] in IBDS[seg[0]].keys():
            IBDS[seg[0]][seg[2]]+=float(seg[-1])/(SAMPLES[seg[0]]*SAMPLES[seg[2]])
            
        if seg[2] not in IBDS[seg[0]].keys():
            IBDS[seg[0]][seg[2]]=float(seg[-1])/(SAMPLES[seg[0]]*SAMPLES[seg[2]])
    
    
    
    
    
        if seg[0] in IBDS[seg[2]].keys():
            IBDS[seg[2]][seg[0]]+=float(seg[-1])/(SAMPLES[seg[0]]*SAMPLES[seg[2]])
            
            
        if seg[0] not in IBDS[seg[2]].keys():
            IBDS[seg[2]][seg[0]]=float(seg[-1])/(SAMPLES[seg[0]]*SAMPLES[seg[2]])


CHOSENPOP='Crete'
pop=[]
for y in IBDS[CHOSENPOP]:
    pop.append([IBDS[CHOSENPOP][y],y])

for y in sorted(pop):
    print(y)