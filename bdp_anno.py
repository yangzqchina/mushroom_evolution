# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 22:01:03 2016

@author: zqyang
"""
import numpy as np
import pandas as pd
#gff_file='/public/home/zqyang/pangenome/analysis/anno/unmapped_genome.gff3'
#bdpfile='/public/home/zqyang/pangenome/analysis/gwas/bdp.txt'
#gff_dat=pd.read_table(gff_file,header=None)
#region_info=['scaffold','marker','start','end']
#bdp=pd.read_table(bdpfile,sep=' ').loc[:,region_info]
class gff(object):
    def __init__(self,gff_file,funcfile):
        self.filename=gff_file
        self.funcfile=funcfile
    def get_info(self):
        gene_info=pd.read_table(self.filename,header=None)
        gene_mrna=gene_info.loc[gene_info.loc[:,2]=='mRNA',:]
        gene_id=[]
        for i in gene_info.index:
            line=gene_info.loc[i,8]
            gene_id.append(line.split(';')[0].replace('ID=',''))
        gene_info.loc[:,8]=gene_id
        gene_info=gene_info.loc[:,[0,2,3,4,8]]
        gene_mrna=gene_mrna.loc[:,[0,2,3,4,8]]
        gene_info.columns=['scaffold','region','start','end','geneid']
        gene_mrna.columns=['scaffold','region','start','end','geneid']
        gene_cds=gene_info.loc[gene_info.loc[:,'region']=='CDS',:]        
        mrna_dict={}
        mrnaid=[]
        geneid=[]
        for i in gene_mrna.index:
            line=gene_mrna.loc[i,'geneid'].split(';')
            mrna_i=line[0].replace('ID=','')
            gene_i=line[1].replace('Parent=','')
            mrna_dict.setdefault(gene_i,mrna_i)
            mrnaid.append(mrna_i)
            geneid.append(gene_i)
        gene_mrna.loc[:,'geneid']=geneid
        gene_mrna.loc[:,'mrnaid']=mrnaid
        gene_cds.loc[:,'geneid']=map(lambda x:x.replace('cds.',''),gene_cds.loc[:,'geneid'])
        gene=gene_info.loc[gene_info.loc[:,'region']=='gene',:]
        gene.loc[:,'cds_length']=gene_cds.groupby('geneid').apply(lambda x:\
        sum(x.end-x.start))[map(lambda x:mrna_dict[x],gene.geneid)].tolist()
        return(gene,gene_cds,gene_mrna)
    def merge_info(self):
        gene,gene_cds,gene_mrna=self.get_info()
        func=pd.read_table(self.funcfile)
        func.columns=['mrnaid']+func.columns.tolist()[1:]
        rna_all_info=pd.merge(left=gene_mrna,right=func,on='mrnaid')
        return(rna_all_info)

funcfile='/public/home/zqyang/pangenome/analysis/anno/function_annotation/function_annotation.txt'
gff_file='/public/home/zqyang/pangenome/analysis/anno/unmapped_genome.gff3'
bdpfile='/public/home/zqyang/pangenome/analysis/gwas/bdp.txt'

gff_dat=gff(gff_file,funcfile)
gff_info,gff_cds,mrna_info=gff_dat.get_info()
rna_info=gff_dat.merge_info()
sweep_region=[]

#bdpfile='/public/home/zqyang/pangenome/analysis/gwas/bdp.txt'
#gff_dat=pd.read_table(gff_file,header=None)
region_info=['scaffold','marker','start','end']
bdp=pd.read_table(bdpfile,sep=' ').loc[:,region_info]
header=list(set(rna_info.columns)-set(['scaffold','start','end']))
Filter=[]
for bdp_index in bdp.index:
    bdp_chr,bdp_start,bdp_end=bdp.loc[bdp_index,['scaffold','start','end']]
    can_rna=rna_info.loc[rna_info.scaffold==bdp_chr,:]
    for can_i in can_rna.index:
        rna_start,rna_end=can_rna.loc[can_i,['start','end']]
        if (bdp_start>rna_end)|(bdp_end<rna_start):
            Filter.append(list(bdp.loc[bdp_index,['scaffold','marker','start','end']])+\
            list(can_rna.loc[can_i,header])+['SAME SCAFFOLD'])
        else:
            Filter.append(list(bdp.loc[bdp_index,['scaffold','marker','start','end']])+\
            list(can_rna.loc[can_i,header])+['PASS'])
dat=pd.DataFrame(Filter)
dat.columns=['scaffold','marker','start','end']+header+['Filter']
dat=dat.loc[:,dat.columns[range(4)+range(16,19)+range(4,16)]]
dat.to_csv("/public/home/zqyang/pangenome/analysis/anno/bdp.anno.csv",index=False)
results=dat.loc[dat.Filter=='PASS',:]
results.to_csv('/public/home/zqyang/pangenome/analysis/anno/bdp_anno_pass.csv',index=False)