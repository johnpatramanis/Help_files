import itertools
import numpy as np

#############################################################################################################################################

def PAIRWISE(individs_tupple):
    individs_list=list(individs_tupple)
    ind1=individs_list[0]
    ind2=individs_list[1]
    DIFF=0
    for l in range(0,len(ind1)):
        if ind1[l]!=ind2[l]:
            DIFF+=1
    PAIRWISEDIFF=DIFF/len(ind1)
    return PAIRWISEDIFF

def AVG_PAIRWISE(pop1,pop2):
    MEAN=[]
    combinations=list(itertools.product(pop1, pop2))
    for comb in combinations:
        pairwise=PAIRWISE(comb)
        MEAN.append(pairwise)
    MEAN=np.mean(MEAN)
    return MEAN


def F2(pop1,pop2):
    myF2=AVG_PAIRWISE(pop1,pop2)-((AVG_PAIRWISE(pop1,pop1)+AVG_PAIRWISE(pop2,pop2))/2)
    return myF2


def F3(pop1,pop2,popX):
    myF3=(F2(popX,pop1)+F2(popX,pop2)-F2(pop1,pop2))/2
    return myF3


#############################################################################################################################################


FILE=open('test.ms','r')
npop=3
Ns=[10,10,10]
set=0
datasets=[]
for line in FILE:
    if line[0]=='/':
        set+=1
        datasets.append([])
    if line[0]=='1' or line[0]=='0':
        datasets[(set-1)].append(line.strip())
        
print(datasets)


ALLF3=[]

for set in datasets:
    populations=[[]]
    counter=0
    popcount=0
    for ind in set:
        populations[popcount].append(ind)
        counter+=1
        if counter==Ns[popcount]:
            counter=0
            popcount+=1
            populations.append([])
    populations=populations[:-1]
    combinations=list(itertools.product(populations[0], populations[1]))
    print(F3(populations[1],populations[0],populations[2]))
    ALLF3.append(F3(populations[1],populations[0],populations[2]))
F3TOTAL=np.mean(ALLF3)
print(F3TOTAL)