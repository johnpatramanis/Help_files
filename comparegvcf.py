import re
import argparse

parser = argparse.ArgumentParser()  
parser.add_argument('--gvcf',nargs=1,type=str)
parser.add_argument('--bim',nargs=1,type=str)

args = parser.parse_args()




file1=open(str(args.gvcf[0]))
file2=open(str(args.bim[0]))
newfilename='new{}'.format(file1)
newfile=open('newgvcf.gvcf','w')


data2=[]
for f in file2:
	line=f.split()
	data2.append(str(line[0])+":"+str(line[3]))


setdata2=set(data2)
print(setdata2)


counts=0
for line in file1:
    if line[0]!='#':
        line=f.strip().split()
        position=str(line[0])+':'+str(line[1])
        if position in setdata2:
            newfile.write(line)
            newfile.write('\n')
            counts+=1
            print(counts)
        else:
            pass
        
    else:
        newfile.write(line)
        newfile.write('\n')





