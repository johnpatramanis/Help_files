import re
import argparse

parser = argparse.ArgumentParser()  
parser.add_argument('--file1',nargs=1,type=str)
parser.add_argument('--file2',nargs=1,type=str)

args = parser.parse_args()




file1=open(str(args.file1[0]))
file2=open(str(args.file2[0]))
newfilename='combined.txt'
newfile=open(str(newfilename),'w')

data1=[]
for f in file1:
	line=f.split()
	data1.append(line[1])

data2=[]
for f in file2:
	line=f.split()
	data2.append(str(line[0])+":"+str(line[3]))

data1=set(data1)
data2=set(data2)

data3=data1.intersection(data2)


print(len(data3))

for l in data3:
	newfile.write(str(l))
	newfile.write('\n')
