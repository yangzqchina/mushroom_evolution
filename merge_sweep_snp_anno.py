# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 09:54:44 2017

@author: zqyang
"""

import pandas as pd
import os,sys
import argparse

ref_anno=pd.read_table('/public/home/zqyang/pangenome/ref_seq/genes.gtf',header=None)
ref_anno=ref_anno.loc[ref_anno.loc[:,2]!='gene',[0,2,3,4,8]]
os.chdir('/public/home/zqyang/pangenome/analysis/gwas')
ref_anno.columns=['chr','region','start','end','INFO']
k_sweeps=os.listdir('Q_K_sweep_snp_A_vs_B/')

anno_info=[]
for file_i in k_sweeps:
    trait=file_i.replace('.txt','').replace('snp_','')
    k_sweep=pd.read_table(os.path.join(sys.argv[1],file_i))
    k_sweep=k_sweep.loc[k_sweep.SWEEP==True,:]
    k_sweep.loc[:,'Trait']=trait

    for i in k_sweep.index:
        snp_info=list(k_sweep.loc[i,:])
        snp_chr,pos=k_sweep.loc[i,['Chromosome','Position']]
        can_region=ref_anno.loc[(ref_anno.chr==snp_chr)&(ref_anno.start<=pos)&(ref_anno.end>=pos),:]
        for j in can_region.index:
            gene_info=list(can_region.loc[j,:])
            anno_info.append(snp_info+gene_info)
results=pd.DataFrame(anno_info)
results.columns=list(k_sweep.columns)+['chr','region','start','end','INFO']
results.to_csv(sys.argv[2]+'_gwas_info.csv',index=False)
