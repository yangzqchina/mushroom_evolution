# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 21:36:09 2016

@author: zqyang
"""

import os,sys
import pandas as pd
import numpy as np


def pan_gwas(workdir,input_file,files_name,windows_file,cdpfile,bdpfile,bin=500,step=250,maf=0.05,cutoff=0.05):
    os.chdir(workdir)
    depth_data=pd.read_table(input_file,header=None)
    samples=np.loadtxt(files_name,dtype=np.str).tolist()
    tags=[sample_i.split('/')[1].split('_')[0] for sample_i in samples]
    depth_data.columns=np.array(['scaffold','pos']+tags,dtype=np.str)
    scaffolds=list(set(depth_data.scaffold))
    
    results=[]
    results_index=["scaffold","marker","start","end"]+tags
    for scaffold_i in scaffolds:
        #scaffold_i=scaffolds[0]
        pos=depth_data.ix[depth_data.scaffold==scaffold_i,"pos"]
        #min_pos=np.min(pos)
        max_pos=np.max(pos)
        #windows=max_pos//step+1
        starts=range(0,max_pos,step)
        i=0
        for start_i in starts:
            marker=scaffold_i+"_"+str(i)
            end=min(start_i+bin,max_pos)
            values=np.mean(depth_data.ix[(depth_data.scaffold==scaffold_i)&(depth_data.pos<end)&(depth_data.pos>start_i),tags],axis=0).tolist()
            results.append([scaffold_i,marker,start_i,end]+values)
            i+=1
    win_data=pd.DataFrame(results)
    win_data.columns=results_index
    win_data.to_csv(windows_file,sep="\t")
    #cdp_data=win_data
    win_data.loc[:,tags]=win_data[tags]/win_data[tags].mean()
    bdp_data=win_data
    bdp_data[tags]=bdp_data[tags][bdp_data[tags]<=cutoff].fillna(1)
    bdp_data[tags]=bdp_data[tags][bdp_data[tags]>cutoff].fillna(0)
    marker_means=bdp_data[tags].mean(axis=1)
    select_marker=(marker_means<(1-cutoff))|(marker_means>cutoff)
    #filter_bdp=bdp_data
    bdp_data.loc[select_marker,:].to_csv(bdpfile, index=False)
    #filter_bdp.to_csv(bdpfile, index=False)
    win_data.loc[select_marker,:].to_csv(cdpfile,index=False)
    #filter_cdp.to_csv(cdpfile,index=False)

workdir="./"
#input_file=sys.argv[1]
input_file=sys.argv[1]
files_name="bam.list"
windows_file="windows_file.txt"
cdpfile="cdp.csv"
bdpfile="bdp.csv"
bin=500
step=250
maf=0.05
cutoff=0.1
pan_gwas(workdir,input_file,files_name,windows_file,cdpfile,bdpfile,bin=500,step=250,maf=0.05,cutoff=0.05)