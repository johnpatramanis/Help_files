import msprime
import numpy as np
import math
import os
import re
import argparse


#python3 Sample_Analysis.py --file_list list.txt



parser = argparse.ArgumentParser()
parser.add_argument('--file_list',nargs=1,type=str)


args = parser.parse_args()







file1=open(str(args.file_list[0]))



files=[]
for line in file1:
    line=line.strip().split()
    files.append(line)
file1.close()

print(len(files))



#os.system('cp ancientdatabase.geno for_merge_ancient.geno')
#os.system('cp ancientdatabase.snp for_merge_ancient.snp')
#os.system('cp ancientdatabase.ind for_merge_ancient.ind')

#os.system('cp vdata.geno for_merge_modern.geno')
#os.system('cp vdata.snp for_merge_modern.snp')
#os.system('cp vdata.ind for_merge_modern.ind')





for rep in files:
    os.system('samtools mpileup -R -B -q30 -Q30 -l /home/kluser2/datasets/wholegenome/1200Ksnps.bed  -f /home/kluser1/aDNA/Reference_Genomes/bwakit/hs37.fa  /home/kluser2/datasets/wholegenome/{} > {}_ancient'.format(rep[0],rep[0]))
    os.system('wc -l {}_ancient'.format(rep[0]))
    os.system('pileupCaller --sampleNames {}  --samplePopName TEST -f /home/kluser2/datasets/wholegenome/ancientdatabase.snp  -o EigenStrat -e {}_ancient < {}_ancient'.format(rep[1],rep[0],rep[0]))

    os.system('mv {}_ancient.geno.txt {}_ancient.geno'.format(rep[0],rep[1]))
    os.system('mv {}_ancient.snp.txt {}_ancient.snp'.format(rep[0],rep[1]))
    os.system('mv {}_ancient.ind.txt {}_ancient.ind'.format(rep[0],rep[1]))
    

    os.system('samtools mpileup -R -B -q30 -Q30 -l /home/kluser2/datasets/wholegenome/snpsforcallingHOA.bed  -f /home/kluser1/aDNA/Reference_Genomes/bwakit/hs37.fa  /home/kluser2/datasets/wholegenome/{} > {}_HOA'.format(rep[0],rep[0]))
    os.system('wc -l {}_HOA'.format(rep[0]))
    os.system('pileupCaller --sampleNames {}  --samplePopName TEST -f /home/kluser2/datasets/wholegenome/vdata.snp  -o EigenStrat -e {}_HOA < {}_HOA'.format(rep[1],rep[0],rep[0]))

    os.system('mv {}_HOA.geno.txt {}_HOA.geno'.format(rep[0],rep[1]))
    os.system('mv {}_HOA.snp.txt {}_HOA.snp'.format(rep[0],rep[1]))
    os.system('mv {}_HOA.ind.txt {}_HOA.ind'.format(rep[0],rep[1]))
    
    
    

    supp1=open('vdata.snp','r')
    supp2=open('ancientdatabase.snp','r')

    snp1=open('{}_HOA.snp'.format(rep[1]),'r')
    snp2=open('{}_ancient.snp'.format(rep[1]),'r')

    newsnp1=open('{}_HOA_new.snp'.format(rep[1]),'w')
    newsnp2=open('{}_ancient_new.snp'.format(rep[1]),'w')



    di={}
    for j in supp1:
            line=j.strip().split()
            position=str(line[1])+':'+str(line[3])
            rs=line[0]
            alt=line[4]
            ref=line[5]
            di[position]=[rs,alt,ref]

    print(len(di))



    for f in snp1:
            line=f.strip().split()
            position=str(line[1])+':'+str(line[3])
            try:
                    line[0]=di[position][0]
                    if line[4]=='<NON_REF>':
                            if line[5]==di[position][1]:
                                    line[4]=di[position][2]
                            else:
                                    line[4]=di[position][1]
            except KeyError:
                    line[0]='Missing_rs'
            newsnp1.write("\t".join(line))
            newsnp1.write("\n")
    di={}
    
    
    for j in supp2:
            line=j.strip().split()
            position=str(line[1])+':'+str(line[3])
            rs=line[0]
            alt=line[4]
            ref=line[5]
            di[position]=[rs,alt,ref]

    print(len(di))

    for f in snp2:
            line=f.strip().split()
            position=str(line[1])+':'+str(line[3])
            try:
                    line[0]=di[position][0]
                    if line[4]=='<NON_REF>':
                            if line[5]==di[position][1]:
                                    line[4]=di[position][2]
                            else:
                                    line[4]=di[position][1]
            except KeyError:
                    line[0]='Missing_rs'
            newsnp2.write("\t".join(line))
            newsnp2.write("\n")
            
            
            
            
    supp1.close()
    supp2.close()

    snp1.close()
    snp2.close()

    newsnp1.close()
    newsnp2.close()

    os.system('mv {}_HOA_new.snp {}_HOA.snp'.format(rep[1],rep[1]))
    os.system('mv {}_ancient_new.snp {}_ancient.snp'.format(rep[1],rep[1]))

    
    
    
    parfileancient=open('parfile_ancient','w')
    parfilemodern=open('parfile_HOA','w')
    
    parfilemodern.write('geno1: vdata.geno\n')
    parfilemodern.write('snp1: vdata.snp\n')
    parfilemodern.write('ind1: vdata.ind\n')
    parfilemodern.write('geno2: {}_HOA.geno\n'.format(rep[1]))
    parfilemodern.write('snp2: {}_HOA.snp\n'.format(rep[1]))
    parfilemodern.write('ind2: {}_HOA.ind\n'.format(rep[1]))
    parfilemodern.write('genooutfilename: {}_HOA_merged.geno\n'.format(rep[1]))
    parfilemodern.write('snpoutfilename: {}_HOA_merged.snp\n'.format(rep[1]))
    parfilemodern.write('indoutfilename: {}_HOA_merged.ind'.format(rep[1]))
    
    
    
    
    
    
    parfileancient.write('geno1: ancientdatabase.geno\n')
    parfileancient.write('snp1: ancientdatabase.snp\n')
    parfileancient.write('ind1: ancientdatabase.ind\n')
    parfileancient.write('geno2: {}_ancient.geno\n'.format(rep[1]))
    parfileancient.write('snp2: {}_ancient.snp\n'.format(rep[1]))
    parfileancient.write('ind2: {}_ancient.ind\n'.format(rep[1]))
    parfileancient.write('genooutfilename: {}_ancient_merged.geno\n'.format(rep[1]))
    parfileancient.write('snpoutfilename: {}_ancient_merged.snp\n'.format(rep[1]))
    parfileancient.write('indoutfilename: {}_ancient_merged.ind\n'.format(rep[1]))
    
    parfilemodern.close()
    parfileancient.close()
    
    
    
    os.system('mergeit -p parfile_ancient')
    os.system('mergeit -p parfile_HOA')
    
    os.system('rm parfile_HOA')
    os.system('rm parfile_ancient')
    
    
    
    convert_ancient=open('convert_ancient','w')
    convert_modern=open('convert_HOA','w')
    
    convert_ancient.write('genotypename: {}_ancient_merged.geno\n'.format(rep[1]))
    convert_ancient.write('snpname: {}_ancient_merged.snp\n'.format(rep[1]))
    convert_ancient.write('indivname: {}_ancient_merged.ind\n'.format(rep[1]))
    convert_ancient.write('outputformat: PACKEDPED\n')
    convert_ancient.write('genotypeoutname: {}_ancient_merged.bed\n'.format(rep[1]))
    convert_ancient.write('snpoutname: {}_ancient_merged.bim\n'.format(rep[1]))
    convert_ancient.write('indivoutname: {}_ancient_merged.fam'.format(rep[1]))

    
    convert_modern.write('genotypename: {}_HOA_merged.geno\n'.format(rep[1]))
    convert_modern.write('snpname: {}_HOA_merged.snp\n'.format(rep[1]))
    convert_modern.write('indivname: {}_HOA_merged.ind\n'.format(rep[1]))
    convert_modern.write('outputformat: PACKEDPED\n')
    convert_modern.write('genotypeoutname: {}_HOA_merged.bed\n'.format(rep[1]))
    convert_modern.write('snpoutname: {}_HOA_merged.bim\n'.format(rep[1]))
    convert_modern.write('indivoutname: {}_HOA_merged.fam'.format(rep[1]))

    convert_modern.close()
    convert_ancient.close()
    
    os.system('convertf -p convert_ancient')
    os.system('convertf -p convert_HOA')
    
    os.system('rm convert_ancient')
    os.system('rm convert_HOA')