import msprime
import numpy as np
import math
import os


N_OG=1000
N_OUT=1000
N_AB=1000
N_A0=1000
N_B0=1000

r_A=0.000
r_B=0.000

generation_time = 25

T_split_OUT_AB=25000/generation_time
T_split_AB=12500/generation_time




N_A=N_A0 / math.exp(-r_A * T_split_AB)
N_B=N_B0 / math.exp(-r_B * T_split_AB)



population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_OUT),
    msprime.PopulationConfiguration(initial_size=N_A, growth_rate=r_A),
    msprime.PopulationConfiguration(initial_size=N_B, growth_rate=r_B)
]

#migration of AB to OUT
m_OUT_AB=0.00000

migration_matrix = [
    [0,0.00000,0.000000],
    [0.00000,0,0.00000],
    [0.000000,0.00000,0]]

samples=[msprime.Sample(0,0)]*100 + [msprime.Sample(1,0)]*100 + [msprime.Sample(2,0)] *100


demographic_events = [
    # A and B merge
    msprime.MassMigration(time=T_split_AB, source=2, destination=1, proportion=1.0),
    msprime.MigrationRateChange(time=T_split_AB, rate=0),
    msprime.MigrationRateChange(time=T_split_AB, rate=m_OUT_AB, matrix_index=(0, 1)),
    msprime.MigrationRateChange(time=T_split_AB, rate=m_OUT_AB, matrix_index=(1, 0)),
    msprime.PopulationParametersChange(time=T_split_AB, initial_size=N_AB, growth_rate=0, population_id=1),
    # Population AB merges into OUT
    msprime.MassMigration(time=T_split_OUT_AB, source=1, destination=0, proportion=1.0),
    # Size changes to N_A at T_AF
    msprime.PopulationParametersChange(time=T_split_OUT_AB, initial_size=N_OUT, population_id=0)
]

dd = msprime.simulate(samples=samples,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,mutation_rate=0.01,
    demographic_events=demographic_events,length=1.0, recombination_rate=2e-8,num_replicates=22)

j=0
for i in dd:
    wow=open('mynewvcf{}.vcf'.format(j),'w')
    i.write_vcf(wow,2,str(j))
    wow.close()
#    print(dd.samples())
    population_labels= ["africa"]*50 + ["asia"]*50 + ["europe"]*50
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
    wowzers.close()
    j+=1

os.system('rm mynewvcf*.vcf')
