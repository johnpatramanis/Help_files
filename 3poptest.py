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


start_time = time.time()


reps=1
for REPS in range(0,reps):

    totalf3=[]
    
##############################################################################################################################################
#Simulation Parameters
    
    parametersfile=open('PARAMETERS_{}'.format(REPS),'w')
    
    N_locals=int(round(random.uniform(500.0,1000.0)))
    N_metropolis=int(round(random.uniform(500.0,1000.0)))
    
    generation_time = 20
    T_COLONIZATION=700/generation_time
    
    
    COLONIZER=random.randint(0,1)
    if COLONIZER==0:
        N_initial_colony=int(round(random.uniform(200.0,float(N_locals))))
        while N_initial_colony>N_metropolis:
            N_initial_colony=int(round(random.uniform(200.0,float(N_metropolis))))
    if COLONIZER==1:
        N_initial_colony=int(round(random.uniform(200.0,float(N_metropolis))))



    r_locals=10**(-1*random.uniform(1,4))
    r_metropolis=10**(-1*random.uniform(1,4))
    r_colony=10**(-1*random.uniform(1,4))
    
    while (N_initial_colony / (math.exp(-r_colony * T_COLONIZATION)) ) > N_metropolis:
        r_colony=10**(-1*random.uniform(1,4))
    
    N_finale_colony=N_initial_colony / (math.exp(-r_colony * T_COLONIZATION))
    print(N_locals,N_metropolis,N_initial_colony,N_finale_colony)
    ###############################################################################################################################
    


    population_configurations = [
        msprime.PopulationConfiguration(initial_size=N_locals,growth_rate=r_locals),
        msprime.PopulationConfiguration(initial_size=N_metropolis, growth_rate=r_metropolis),
        msprime.PopulationConfiguration(initial_size=N_finale_colony, growth_rate=r_colony)
    ]



    migration_matrix = [
        [0,10**(-1*random.uniform(1,4)),10**(-1*random.uniform(1,4))],
        [10**(-1*random.uniform(1,4)),0,10**(-1*random.uniform(1,4))],
        [10**(-1*random.uniform(1,4)),10**(-1*random.uniform(1,4)),0]]

    N1=20
    N2=20
    N3=20
    POPS=[N1,N2,N3]
    samples=[msprime.Sample(0,0)]*N1 + [msprime.Sample(1,0)]*N2 + [msprime.Sample(2,0)] *N3

    demographic_events = [
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(0, 2)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(2, 0)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(1, 2)),
    msprime.MigrationRateChange(time=T_COLONIZATION, rate=0, matrix_index=(2, 1)),
    msprime.PopulationParametersChange(time=T_COLONIZATION, initial_size=N_initial_colony, growth_rate=0, population_id=2),
    msprime.MassMigration(time=T_COLONIZATION, source=2, destination=COLONIZER, proportion=1.0),
    
    
    ]

    parametersfile.write('\t'.join([str(x) for x in [COLONIZER,N_locals,N_metropolis,N_initial_colony,N_finale_colony,r_locals,r_metropolis,r_colony]]))
    parametersfile.write('\n')
    parametersfile.write('\t'.join([str(x) for x in migration_matrix]))
    
    print(migration_matrix)
    
    




######################################################################################################################################################
#RUN the simulation and output genotypes in vcfs and ms format files, one for each chrom 

    
    def SIMULATE(L,argument,samples,population_configurations,migration_matrix,demographic_events):
        j=int(argument)
        recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(j))
        dd = msprime.simulate(samples=samples,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,mutation_rate=1e-8,
            demographic_events=demographic_events,recombination_map=recomb_map)
        outfile=open('ms_prime_{}'.format(j),'w')   
        for var in dd.variants():
            L.append([j,var.index,var.position])
            for genotype in var.genotypes:
                outfile.write(str(genotype))
            outfile.write('\n')
        outfile.close()    
        wow=open('mynewvcf{}.vcf'.format(j),'w')
        dd.write_vcf(wow,2,str(j))
        wow.close()
        
        population_labels= ["locals"]*int(N1/2) + ["metropolis"]*int(N2/2) + ["apoikia"]*int(N3/2)
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
    if __name__ == '__main__':
        with Manager() as manager:
            L=manager.list(L)
            processes=[]
            for loop in range(1,23):
                p=Process(target=SIMULATE,args=(L,loop,samples,population_configurations,migration_matrix,demographic_events,))
                processes.append(p)
                
                p.start()
        
                
            for p in processes:
                p.join()
            #print(len(L),'1')
            sys.stdout.flush()
            variants=sorted(list(L))


    variantinfo=['{}\t{}\t{}\n'.format(x[0],x[1],x[2])for x in variants]
    print(len(variants),len(variantinfo))

    variantinformation=open('variants_info.txt','w')
    variantinformation.write('CHROM\tVARIANT\tPOSITION\n')
    for loop in variantinfo:
        variantinformation.write(loop)
    
    variantinformation.close()

    elapsed_time_1 = time.time() - start_time        
        
    print('Step 1 : {} '.format(elapsed_time_1/60))        
    


#Arrange the vcf files into one, fix labels , bed file , transform to eigen, calculate pca and f stats
        
    os.system('rm mynewvcf*.vcf')
    os.system('bcftools concat -o total_chroms.vcf myvcf*.vcf')
    os.system('rm myvcf*.vcf')




    VCF=open('total_chroms.vcf','r')
    newVCF=open('newtotal_chroms.vcf','w')

    snpcount=0
    
    variants=sorted(variants)
    for line in VCF:
        if line[0]!='#' and snpcount<len(variants):
            line=line.strip().split()
            if len(line)<=2:
                continue
            line[2]='rs{}'.format(snpcount)
            line[1]=str(variants[snpcount][2])
            line.append('\n')
            line='\t'.join(line)
            snpcount+=1
        newVCF.write(line)

    VCF.close
    newVCF.close

    os.system('mv newtotal_chroms.vcf total_chroms.vcf')

    os.system('plink --vcf total_chroms.vcf --make-bed --out simulation')

    
############################################################## RUN PCA ######################################################################
    if os.path.isfile('simulation.bed'):
        simulationfile='simulation'
    else:
        simulationfile='simulation-temporary'
    
    os.system('plink --bfile {} --pca 10 --out pcaofsimulation'.format(simulationfile))


########################################### 3 Pop Test ######################################################################################
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
    
    
    for k in range(1,22):
        segments=[[0,100],[200,300],[400,500]]
        for j in segments:
            os.system('plink --vcf total_chroms.vcf --chr {} --from-kb {} --to-kb {}  --make-bed --out simulation'.format(k,segments[0],segments[1]))

            
        ############################################################## RUN PCA ######################################################################
            if os.path.isfile('simulation.bed'):
                simulationfile='simulation'
            else:
                simulationfile='simulation-temporary'
            
            os.system('plink --bfile {} --pca 10 --out pcaofsimulation'.format(simulationfile))


        ########################################### 3 Pop Test ######################################################################################
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
    
    
    
    
    
    os.system('rm simulation.*')
    os.system('rm simulation-temporary.*')
    os.system('rm ms_prime_*')
    for x in range(1,23):
        os.system('rm ms_{}'.format(x))

    elapsed_time_3=time.time() - start_time
    print('step 3 : {}'.format(elapsed_time_3/60))        
    
    
###############################################################################################################################################
print(totalf3)

