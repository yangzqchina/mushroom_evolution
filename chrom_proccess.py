# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 09:18:42 2016

@author: zqyang
"""
import sys

infile=open(sys.argv[1])
map_lines=[]
chrom_dict={}
chroms=[]
i=0
outfile=open(sys.argv[2],'w')
for line in infile:
    lines=line.rstrip().split('\t')
    chrom=lines[1].split(':')[0]
    #lines.append(chrom)
    if not chrom in chroms:
        chroms.append(chrom)
        i=i+1
        chrom_dict.setdefault(chrom,str(i))
    map_lines.append(lines)
    outfile.write(chrom_dict[chrom]+'\t'+'\t'.join(lines[1:])+'\n')
infile.close()
outfile.close()
