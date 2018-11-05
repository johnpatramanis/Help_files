import re
import argparse

parser = argparse.ArgumentParser()  
parser.add_argument('--file1',nargs=1,type=str)
parser.add_argument('--file2',nargs=1,type=str)

args = parser.parse_args()


mypops=['europe','ancient','near_east','north_africa']

file1=open(str(args.file1[0]))
file2=open(str(args.file2[0]))

newfile=open('mykeepfam.txt','w')

data1=[]
for f in file1:
	line=f.strip().split(" ")
	data1.append(line)

data2=[]
for f in file2:
	line=f.strip().split(" ")
	data2.append(line)

data3=[]
for x in data1:
	for y in data2:
		if x[0]==y[0]:
			if (y[1] in mypops or y[2] in mypops) and (y[0] not in data3):
				data3.append(y[0])


for l in data3:
	newfile.write(str(l))
	newfile.write('\n')