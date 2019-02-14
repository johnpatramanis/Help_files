import msprime
import numpy as np
import math
import os
import time
import re


N_OG=1000
N_OUT=1000
N_AB=1000
N_A0=1000
N_B0=1000

r_A=0.000
r_B=0.000

generation_time = 25

T_split_OUT_AB=100000/generation_time
T_split_AB=50000/generation_time




N_A=N_A0 / math.exp(-r_A * T_split_AB)
N_B=N_B0 / math.exp(-r_B * T_split_AB)



population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_OUT),
    msprime.PopulationConfiguration(initial_size=N_A, growth_rate=r_A),
    msprime.PopulationConfiguration(initial_size=N_B, growth_rate=r_B)
]

#migration of AB to OUT
m_OUT_AB=0.001

migration_matrix = [
    [0,0.0001,0.0001],
    [0.0001,0,0.001],
    [0.0001,0.001,0]]

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


reps=1000
totalf3=[]
for REPS in range(0,reps):
    dd = msprime.simulate(samples=samples,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,mutation_rate=1e-8,
        demographic_events=demographic_events,length=10000.0, recombination_rate=2e-8,num_replicates=22)

    j=1
    variants=[]
    for i in dd:
        
        for var in i.variants():
            variants.append(var.position)
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
    os.system('bcftools concat -o total_chroms.vcf myvcf*.vcf')
    #os.system('rm myvcf*.vcf')




    VCF=open('total_chroms.vcf','r')
    newVCF=open('newtotal_chroms.vcf','w')

    snpcount=0
    for line in VCF:
        if line[0]!='#' and snpcount<len(variants):
            line=line.strip().split()
            line[2]='rs{}'.format(snpcount)
            line[1]=str(variants[snpcount])
            line.append('\n')
            line='\t'.join(line)
            snpcount+=1
        newVCF.write(line)

    VCF.close
    newVCF.close

    os.system('mv newtotal_chroms.vcf total_chroms.vcf')

    os.system('plink --vcf total_chroms.vcf --make-bed --out simulation')



    

    import os.path
    if os.path.isfile('simulation.bed'):
        simulationfile='simulation'
    else:
        simulationfile='simulation-temporary'
    
    os.system('plink --bfile {} --pca 10 --out pcaofsimulation'.format(simulationfile))
    
    
    ####################################### 3 Pop Test ######################################################################################
    parfile=open('parfile.txt','w')

    parfile.write('genotypename:    {}.bed\n'.format(simulationfile))
    parfile.write('snpname:         {}.bim\n'.format(simulationfile))
    parfile.write('indivname:       {}.fam\n'.format(simulationfile))
    parfile.write('outputformat:   PACKEDANCESTRYMAP\n')
    parfile.write('genotypeoutname: simulation.geno\n')
    parfile.write('snpoutname:      simulation.snp\n')
    parfile.write('indivoutname:    simulation.ind\n')
    parfile.write('pordercheck: NO')

    parfile.close()


    os.system('convertf -p parfile.txt')

    IND=open('simulation.ind','r')
    newIND=open('newsimulation.ind','w')

    for line in IND:
        line=line.strip().split()
        label=re.search(r'([a-z]+)([0-9]+):[a-z]+[0-9]+',line[0])
        pop=label.group(1)
        number=label.group(0)
        line[0]=str(number)
        line[2]=str(pop)
        newIND.write('\t'.join(line))
        newIND.write('\n')
        
    IND.close
    newIND.close()
        
    os.system('mv newsimulation.ind simulation.ind')

    Pop3=open('qp3Poplist','w')
    Pop3.write('africa europe asia')
    Pop3.close()

    Parfilepop=open('3popparfile','w')
    Parfilepop.write('SSS: allmap\n')
    Parfilepop.write('indivname:   simulation.ind\n')
    Parfilepop.write('snpname:     simulation.snp\n')
    Parfilepop.write('genotypename: simulation.geno\n')
    Parfilepop.write('popfilename: qp3Poplist\n')

    Parfilepop.close()


    os.system('ls')



    os.system('qp3Pop -p 3popparfile >f3stat')
    f3file=open('f3stat','r')
    for line in f3file:
        line=line.strip().split()
        #print(line)
        if line[0]=='result:':
            totalf3.append(float(line[4]))
    f3file.close()
    
    
    
    
    os.system('rm simulation.*')
    os.system('rm simulation-temporary.*')
#########################################################################################################################################
        




        
    for k in range(1,21):
    
    
        os.system('plink --vcf myvcf{}.vcf --make-bed --out simulation'.format(k))
        if os.path.isfile('simulation.bed'):
            simulationfile='simulation'
        else:
            simulationfile='simulation-temporary'
        
        
        
        ####################################### 3 Pop Test ######################################################################################
        parfile=open('parfile.txt','w')

        parfile.write('genotypename:    {}.bed\n'.format(simulationfile))
        parfile.write('snpname:         {}.bim\n'.format(simulationfile))
        parfile.write('indivname:       {}.fam\n'.format(simulationfile))
        parfile.write('outputformat:   PACKEDANCESTRYMAP\n')
        parfile.write('genotypeoutname: simulation.geno\n')
        parfile.write('snpoutname:      simulation.snp\n')
        parfile.write('indivoutname:    simulation.ind\n')
        parfile.write('pordercheck: NO')

        parfile.close()


        os.system('convertf -p parfile.txt')

        IND=open('simulation.ind','r')
        newIND=open('newsimulation.ind','w')

        for line in IND:
            line=line.strip().split()
            label=re.search(r'([a-z]+)([0-9]+):[a-z]+[0-9]+',line[0])
            pop=label.group(1)
            number=label.group(0)
            line[0]=str(number)
            line[2]=str(pop)
            newIND.write('\t'.join(line))
            newIND.write('\n')
            
        IND.close
        newIND.close()
            
        os.system('mv newsimulation.ind simulation.ind')

        Pop3=open('qp3Poplist','w')
        Pop3.write('locals metropolis apoikia')
        Pop3.close()

        Parfilepop=open('3popparfile','w')
        Parfilepop.write('SSS: allmap\n')
        Parfilepop.write('indivname:   simulation.ind\n')
        Parfilepop.write('snpname:     simulation.snp\n')
        Parfilepop.write('genotypename: simulation.geno\n')
        Parfilepop.write('popfilename: qp3Poplist\n')

        Parfilepop.close()


        #SNP=open('simulation.snp','r')
        #newSNP=open('newsimulation.snp','w')
        #snpcounter=0
        #for line in SNP:
        #    if snpcounter<len(variants):
        #        line=line.strip().split()
        #        line[0]='rs{}'.format(snpcounter)
        #        line[2]=str(variants[snpcounter])
        #        line.append('\n')
        #        line='\t'.join(line)
        #        snpcounter+=1
        #        newSNP.write(line)




        #SNP.close
        #newSNP.close

        os.system('mv newsimulation.snp simulation.snp')
        os.system('ls')



        os.system('qp3Pop -p 3popparfile >f3stat_{}'.format(REPS))
        
        
        f3file=open('f3stat_{}'.format(REPS),'r')
        
        for line in f3file:
            line=line.strip().split()
            #print(line)
            if line[0]=='result:':
                totalf3.append(float(line[4]))
        f3file.close()

    
    
    
    
    
    
    
print(totalf3[0])
print(totalf3[1:])


###############################################################################################################################################
