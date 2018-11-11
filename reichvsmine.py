import argparse
import gzip
import re
import random
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import plotly
import plotly.plotly as py
import plotly.figure_factory as ff
import concurrent.futures
from scipy.cluster import hierarchy

#python metaptux.py --vcf peloplineage.vcf

#########################Eisagwgh parametrwnnnn##########################################################################
parser = argparse.ArgumentParser()  
parser.add_argument('--vcf',nargs='+',type=str)
parser.add_argument('--vcf2',nargs='+',type=str)




args = parser.parse_args()
#########################################################################################################################
print('wow')

myfile=open(args.vcf[0])            
mygenotypes={}

for line in myfile:
    if line[0]!='#':
        line=line.strip().split()
        position=str(line[0])+':'+str(line[1])
        genotype=str(line[9][0:3])
        mygenotypes[position]=genotype

                
                
print(len(mygenotypes))



genoset=set(mygenotypes)
        
        
        
        
        
if args.vcf2:
    myfile2=open(args.vcf2[0])
    bothalt=0
    bothref=0
    reichref=0
    mineref=0
    heteroz=0
    same=0
    diff=0
    
    
    for line in myfile2:
        if line[0]!='#':
            line=line.strip().split()
            genotypes=[x for x in line[10:11]]
            position=str(line[0])+':'+str(line[1])
            if position in genoset:
                REF=mygenotypes[position]
                for GENO in range(0,len(genotypes)):
                    if REF==genotypes[GENO] and REF=='1/1':
                        bothalt+=1
                    if REF==genotypes[GENO] and REF=='0/0':
                        bothref+=1
                    if REF!=genotypes[GENO] and REF=='1/1':
                        reichref+=1
                    if REF!=genotypes[GENO] and REF=='0/0':
                        mineref+=1
                    if REF!=genotypes[GENO] and ( REF=='0/1' or REF=='1/0'):
                        heteroz+=1
    print('\ncorrectly alt:',bothalt,'\n correctly ref:',bothref,'\n incorrectly alt:',reichref,'\n incorrectly ref:',mineref,'\n heterozygous:',heteroz)            
    print(sum([bothalt,bothref,reichref,mineref,heteroz]))        
