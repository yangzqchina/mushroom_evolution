# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 14:41:51 2017

@author: zqyang
"""

import pandas as pd
import argparse

Args=argparse.ArgumentParser(description="Get stats information of bam file.\n\
 It will print sample reads mapped_reads mapped_reate genome_coverage mean_depth")
Args.add_argument('-d','--depth',action='store',default=None,dest='depthfile',\
help='Input a depth file of samtools depth(-a)')
Args.add_argument('-s','--stats',action='store',default=None,dest='statsfile',\
help='Input a stats file.( samtools flagstats)')
Args.add_argument('-l','--label',action='store',default=None,dest='label',help='Sample name')
args=Args.parse_args()

sum_base=40162560
file1=pd.read_table(args.depthfile,header=None)
file1.columns=['chr','pos','frq']
cov=file1.loc[file1.frq!=0,:].shape[0]/float(sum_base)
mean_depth=file1.frq.sum()/float(sum_base)
file1=pd.read_table(args.statsfile,header=None)
sum_reads=int(file1.loc[0,0].split(' ')[0])
mapped_reads=int(file1.loc[4,0].split(' ')[0])
mapped_ratio=mapped_reads/float(sum_reads)
print('%s\t%i\t%i\t%f\t%f\t%f'%(args.label,sum_reads,mapped_reads,mapped_ratio,cov,mean_depth))
