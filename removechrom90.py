import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file',nargs=1,type=str)
parser.add_argument('--out',nargs=1,type=str)

args = parser.parse_args()




file1=open(str(args.file[0]))
file2=open(str(args.out[0]),'w')

mychroms=[x for x in range(0,24)]
mychroms.append('MT')

for f in file1:
	if f[0]!='#':
        	line=f.strip().split()
        	if str(line[0])not in mychroms:
			line[0]=str(0)
	file2.write("\t".join(line))
	file2.write("\n")
