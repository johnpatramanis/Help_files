import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file',nargs=1,type=str)
parser.add_argument('--supp',nargs=1,type=str)
parser.add_argument('--out',nargs=1,type=str)

args = parser.parse_args()




file1=open(str(args.file[0]))
supp=open(str(args.supp[0]))
file2=open(str(args.out[0]),'w')

di={}
for j in supp:
	line=j.strip().split()
	position=str(line[0])+':'+str(line[3])
	rs=line[1]
	alt=line[4]
	ref=line[5]
	di[position]=[rs,alt,ref]
	
print(len(di))



for f in file1:
        line=f.strip().split()
	position=str(line[0])+':'+str(line[3])
        try:
		line[1]=di[position][0]
		if line[4]=='<NON_REF>':
			if line[5]==di[position][1]:
				line[4]=di[position][2]
			else:
				line[4]=di[position][1]
	except KeyError:
		line[1]='Missing_rs'
        file2.write("\t".join(line))
        file2.write("\n")




