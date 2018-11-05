import re
import argparse

parser = argparse.ArgumentParser()  
parser.add_argument('--file',nargs=1,type=str)
parser.add_argument('--out',nargs=1,type=str)

args = parser.parse_args()




file1=open(str(args.file[0]))
file2=open(str(args.out[0]),'w')



for f in file1:
	line=f.strip().split()
	line[1]=str(line[0])+':'+str(line[3])
	file2.write("\t".join(line))
	file2.write("\n")


