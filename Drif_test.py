import msprime
import numpy as np
import numpy.linalg
import math
import os
import time
import re
import random
import sys
import os.path
from multiprocessing import Process,Manager
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
from bokeh.io import output_notebook, show
from bokeh.plotting import figure
from bokeh.models import TapTool, CustomJS, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.io import output_file
from bokeh import colors


start_time = time.time()


REPS=1


##########################################################################################################################################


N_Boet1i=100000
N_Boet2i=10
N_Boet3i=10

N_Foreign1=100
N_Foreign2=100

generation_time = 20

T=3000/generation_time
T_old=5000/generation_time
T_old_old=10000/generation_time



r_boet=0.001


N_Boet1= N_Boet1i/ math.exp(-r_boet * (T/generation_time))
N_Boet2= N_Boet2i/ math.exp(-r_boet * (T/generation_time))
N_Boet3= N_Boet3i/ math.exp(-r_boet * (T/generation_time))



###############################################################################################################################



population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_Boet1i,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Boet2i,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Boet3i,growth_rate=r_boet),
    msprime.PopulationConfiguration(initial_size=N_Foreign1,growth_rate=0),
    msprime.PopulationConfiguration(initial_size=N_Foreign2,growth_rate=0)

]



migration_matrix = [
[0,0.001,0.001,0.0,0.00],
[0.01,0,0.01,0.0,0.0],
[0.1,0.1,0,0,0],
[0.00,0.0,0,0,0.00],
[0.00,0.0,0.,0.0,0]
]

N1=50
N2=50
N3=50




POPS=[N1,N2,N3]
samples=[msprime.Sample(0,0)]*N1 + [msprime.Sample(0,50)]*N2 + [msprime.Sample(0,100)] *N3




demographic_events = [

msprime.PopulationParametersChange(time =T , growth_rate = 0 , population_id = 0),
msprime.PopulationParametersChange(time =T ,growth_rate=0 , population_id = 1),
msprime.PopulationParametersChange(time =T ,growth_rate=0 , population_id = 2),

msprime.PopulationParametersChange(time =T , initial_size = N_Boet1i , population_id = 0),
msprime.PopulationParametersChange(time =T , initial_size = N_Boet2i, population_id = 1),
msprime.PopulationParametersChange(time =T , initial_size = N_Boet3i, population_id = 2),

msprime.MigrationRateChange(time=T , rate=0, matrix_index=(0,3)), 
msprime.MigrationRateChange(time=T , rate=0, matrix_index=(3,0)),
msprime.MigrationRateChange(time=T , rate=0, matrix_index=(4,0)),
msprime.MigrationRateChange(time=T , rate=0, matrix_index=(0,4)),


msprime.MassMigration(time=T,source=2,destination=1,proportion = 1.0),
msprime.MassMigration(time=T,source=1,destination=0,proportion = 1.0),
msprime.MassMigration(time=T_old,source=3,destination=0,proportion = 1.0),
msprime.MassMigration(time=T_old_old,source=4,destination=0,proportion = 1.0)




]






##################################################################################################################################################


def SIMULATE(L,argument,samples,population_configurations,migration_matrix,demographic_events):
    j=int(argument)
    #recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(j))
    dd = msprime.simulate(samples=samples,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,mutation_rate=1e-8,recombination_rate=2e-8,
        demographic_events=demographic_events,length=150000000)
    outfile=open('ms_prime_{}'.format(j),'w')
    for var in dd.variants():
        L.append([int(j),var.index,var.position])
        for genotype in var.genotypes:
            outfile.write(str(genotype))
        outfile.write('\n')
    outfile.close()    
    wow=open('mynewvcf{}.vcf'.format(j),'w')
    dd.write_vcf(wow,2,str(j))
    wow.close()
    
    population_labels= ["Boet1_"]*int(N1/2) + ["Boet2_"]*int(N2/2) + ["Boet3_"]*int(N3/2)
    d=0
    newlabels=[]
    for i in range(0,len(population_labels)):
        newlabels.append(population_labels[i]+str(d))
        d+=1
        if i==len(population_labels)-2:
            newlabels.append(population_labels[i]+str(d))
            break
        if population_labels[i]!=population_labels[i+1]:
            d=0
    population_labels=newlabels
    wow=open('mynewvcf{}.vcf'.format(j))
    wowzers=open('myvcf{}.vcf'.format(j),'w')
    for line in wow:
        line=line.strip().split()
        if line[0]=='#CHROM':
            line[9:]=population_labels
        wowzers.write("\t".join(line))
        wowzers.write("\n")
    wow.close()
    
    return j,L
L=[]
SIMULATE(L,0,samples,population_configurations,migration_matrix,demographic_events)



elapsed_time_1 = time.time() - start_time        
    
print('Step 1 : {} '.format(elapsed_time_1/60))        



vcf=open('myvcf0.vcf')



def get_colors(n):
    colorz=[]
    for k in range(0,n):
        c1=round(random.random(), 3)
        c2=round(random.random(), 3)
        c3=round(random.random(), 3)
        colorz.append((c1,c2,c3))
    return colorz
#####IIIINCOOOMPLEEEETE




for line in vcf:
    if line[0]=='#' and line[1]!='#':
        population_labels=line.strip().split()[9:]
        true_labels=[]
        for k in population_labels:
            familyname=k.split('_')[0]
            true_labels.append(familyname)
        genotypes=[[] for x in range(0,len(population_labels))]
    if line[0]!='#':
        line=line.strip().split()
        counter=0
        for j in line[9:]:
            if j=='0|0':
                genotypes[counter].append(0)
            if j=='1|0' or j=='0|1':
                genotypes[counter].append(1)
            if j=='1|1':
                genotypes[counter].append(2)
            if j=='.|.':
                genotypes[counter].append(9)
            counter+=1

print(len(genotypes),len(genotypes[0]))
for wow in genotypes:
    print(len(wow))            
#TRANSFORM TO tSNE
X = np.asarray(genotypes)    
X_embedded = TSNE(n_components=2,learning_rate=100.0,n_iter=1000000,perplexity=50.0).fit_transform(X)    
print(X_embedded.shape)   



#PLOTING
plt.figure(figsize=(100, 60))

COLORPALLETE=get_colors(len(set(true_labels)))
COLORZ_TO_LABELS={}

uniquelabels=[x for x in set(true_labels)]
for j in range(0,len(uniquelabels)):
    COLORZ_TO_LABELS[uniquelabels[j]]=COLORPALLETE[j]

colors=[ COLORZ_TO_LABELS[x] for x in true_labels]
   
plt.scatter([x[0] for x in X_embedded],[x[1] for x in X_embedded],label=true_labels,c=colors)
plt.show()


pca = PCA(n_components=2).fit_transform(genotypes)
X=[]
Y=[]
for j in pca:
    X.append(j[0])
    Y.append(j[1])

plt.figure(figsize=(100, 60))

   
plt.scatter(X,Y,label=population_labels,c=colors)
plt.show()
        


