import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file',nargs=1,type=str)

args = parser.parse_args()

print(str(args.file[0]))


file=open(str(args.file[0]))
newfilename='new'+str(str(args.file[0]))
newfile=open(str(newfilename),'w')


for f in file:
        line=f.strip().split()
	if (line[4]=='A'  and line[5]=='C' ) or (line[4]=='C' and line[5]=='A') or (line[4]=='G'  and line[5]=='C' ) or (line[4]=='C'  and line[5]=='G' ) or (line[4]=='G' and line[5]=='T') or (line[4]=='T' and line[5]=='G') or (line[4]=='A' and line[5]=='T') or (line[4]=='T' and line[5]=='A')  :
        	line=" ".join(line)
        	newfile.write(line)
		newfile.write('\n')


