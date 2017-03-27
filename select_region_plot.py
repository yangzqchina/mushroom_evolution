# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a plot fst_vs_pi script file.

Usage:

python filter_plot.py popA_vs_popB_Fst.windowed.weir.fst popA_5kb.pi.windowed.pi \
popB_5kb.pi.windowed.pi A B

output:
A_vs_B.site.csv
A_vs_B_fst_vs_pi.pdf

"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse

parses=argparse.ArgumentParser(description='This is a software to analysis Fst vs Pi')

parses.add_argument('--fst',dest='fst',action='store',help='A Fst file of vcftools',default=None)
parses.add_argument('--pi1',dest='pi1',action='store',help='A population\'s Pi file',default=None)
parses.add_argument('--pi2',dest='pi2',action='store',help='Another population\'s Pi file',default=None)
parses.add_argument('--l1',dest='label1',action='store',help='A label reprisents one population',default=None)
parses.add_argument('--l2',dest='label2',action='store',help='A label reprisents another population',default=None)
parses.add_argument('--region',dest='region',action='store',help='Which region of genome you are analysiing',default='')
args=parses.parse_args()
matplotlib.use('Agg') 
os.chdir('/public/home/zqyang/pangenome/analysis/gatk_all_results/population_computing')

fst=args.fst
pia=args.pi1
pib=args.pi2
labelA=args.label1
labelB=args.label2
region=args.region
outfile=labelA+'_vs_'+labelB
#fst="popA_vs_popB_Fst.windowed.weir.fst"
#pia="popA_5kb.pi.windowed.pi"
#pib="popB_5kb.pi.windowed.pi"
#outfile="popA_vs_popB"
fst=pd.read_table(fst)
pia=pd.read_table(pia)
pib=pd.read_table(pib)

pia.columns=pd.Index(list(pia.columns[:3])+['N_VARIANTS_'+labelA,'PI'+labelA])
pib.columns=pd.Index(list(pib.columns[:3])+['N_VARIANTS_'+labelB,'PI'+labelB])

dat=pd.merge(left=fst,right=pia,on=('CHROM', 'BIN_START', 'BIN_END'))
dat=pd.merge(left=dat,right=pib,on=('CHROM', 'BIN_START', 'BIN_END'))
dat.loc[:,"PI_Ratio"]=dat.loc[:,'PI'+labelA]/dat.loc[:,'PI'+labelB]

#compute cutoff values of pi ratio and Fst
cutoff_pi=dat.PI_Ratio.quantile(q=(0.05,0.9)).tolist()
cutoff_fst=dat.MEAN_FST.quantile(q=0.95)
dat.loc[:,"group"]="no_pass"
dat.loc[(dat.PI_Ratio>cutoff_pi[1])&(dat.MEAN_FST>cutoff_fst),"group"]=labelA
dat.loc[(dat.PI_Ratio<cutoff_pi[0])&(dat.MEAN_FST>cutoff_fst),"group"]=labelB
if not region=='':
    dat.to_csv("%s_in_%s_fst_vs_pi.site.csv"%(region,outfile),index=False)
else:
    dat.to_csv("%s_fst_vs_pi.site.csv"%(outfile),index=False)
#create a function to compute pi distribution and cumulative frequency distribution
#a bin-number is necessary
def stat_freq(dat,col_name="PI_Ratio",bins=1000):
    data=dat.loc[:,col_name]
    n=len(data)
    min_data=min(data)
    max_data=max(data)
    breaks=np.linspace(min_data,max_data,bins).tolist()
    results=[]
    for i in range(len(breaks)-1):
        if i==0:
            break_index=map(lambda x:(x>=breaks[i])&(x<=breaks[i+1]),data)
            freq_i=len(data[break_index])/float(n)*100
            cumsum_i=freq_i
            value_i=np.mean([breaks[i],breaks[i+1]])
            results.append([value_i,freq_i,cumsum_i])
        else:
            break_index=map(lambda x:(x>breaks[i])&(x<=breaks[i+1]),data)
            freq_i=len(data[break_index])/float(n)*100
            cumsum_i=freq_i+cumsum_i
            value_i=np.mean([breaks[i],breaks[i+1]])
            results.append([value_i,freq_i,cumsum_i])
    return(results)

dat_freq=np.array(stat_freq(dat,"PI_Ratio",1200),dtype=np.float)
select_area=(dat_freq[:,0]>cutoff_pi[0])&(dat_freq[:,0]<cutoff_pi[1])
pi_no=(dat_freq[:,0]<=cutoff_pi[0])|(dat_freq[:,0]>=cutoff_pi[1])

fst_freq=np.array(stat_freq(dat,"MEAN_FST",800),dtype=np.float)
fst_select_area=(fst_freq[:,0]>cutoff_fst)
fst_no=(fst_freq[:,0]<=cutoff_fst)
#fig1.get_yticklabels()[0].set_visible(False)

#plt.figure(figsize=(10,10))

fill_color1="#333366"
fill_color2="#999999"
cutline_color="#663300"
font = {'family' : 'serif',  
        'color'  : 'darkred',  
        'weight' : 'normal',  
        'size'   : 16,  
        } 

fig=plt.figure(figsize=(10,10))

fig0=fig.add_axes([0.8,0.8,0.2,0.2])
plt.xlim(0,1)
plt.ylim(0,1)
#text="Choice Standard:\n$"+str(cutoff_pi[0])+"<\pi_A/pi_B\<"+str(cutoff_pi[1])+"$"

#First, we should write a text on the topright of the picture
fig0.spines['right'].set_color('none')
fig0.spines['top'].set_color('none')
fig0.spines['left'].set_color('none')
fig0.spines['bottom'].set_color('none')
fig0.set_xticks([])
fig0.set_yticks([])
plt.text(s="Cumulative\nFrequency(%)",x=0.2,y=0.45,fontdict=font)

#Next, we'll plot a pi ratio's distribution
fig1=fig.add_axes([0,0.8,0.8,0.2])
#plt.plot(dat_freq[:,0],dat_freq[:,1],color='black',linewidth=0.5,ls='dashed')
plt.fill_between(dat_freq[pi_no,0],dat_freq[pi_no,1],y2=0,color=fill_color1,\
label="$Selected\ Region$")
plt.fill_between(dat_freq[select_area,0],dat_freq[select_area,1],y2=0,\
color=fill_color2,label="$Not\ Selected\ Region$")
plt.xlim((-2,25))
plt.ylabel(r'Frequency (%)',fontdict=font)
plt.axvline(cutoff_pi[0],ls='dashed',linewidth=1.5,color=cutline_color)
plt.axvline(cutoff_pi[1],ls='dashed',linewidth=1.5,color=cutline_color)
plt.legend(loc=4)
fig1_2=plt.twinx()
plt.plot(dat_freq[:,0],dat_freq[:,2],color="maroon",linewidth=2,label="$Cumulative\ Frequency$")
plt.xlim((-2,25))
plt.ylim((0,100))
plt.text(s="Choice Standard:\n$%f<\pi_{%s}/\pi_{%s}<%f$\n$Fst>%f$"%(cutoff_pi[0],\
labelA,labelB,cutoff_pi[1],cutoff_fst),x=7,y=30,fontdict=font)
#plt.ylabel(r'cumulative(%)',fontsize=15)
for i in fig1.get_xticklabels():
    i.set_visible(False)
fig1_2.get_yticklabels()[0].set_visible(False)
plt.legend(loc=1)

fig2=fig.add_axes([0,0,0.8,0.8])
p1=plt.scatter(x=dat.loc[dat.group==labelA,"PI_Ratio"],y=dat.loc[dat.group==labelA,\
"MEAN_FST"],c="#333366",marker='o',label="$\pi_{%s}<\pi_{%s}$"%(labelA,labelB),\
edgecolors='none',s=15)
p2=plt.scatter(x=dat.loc[dat.group==labelB,"PI_Ratio"],y=dat.loc[dat.group==labelB,\
"MEAN_FST"],c="#006633",marker='o',label="$\pi_{%s}<\pi_{%s}$"%(labelB,labelA),\
edgecolors='none',s=15)
p3=plt.scatter(x=dat.loc[dat.group=="no_pass","PI_Ratio"],y=dat.loc[dat.group=="no_pass",\
"MEAN_FST"],c=fill_color2,marker='o',label='$not\ be\ selected$',edgecolors='none',\
s=15,alpha=0.5)
plt.xlabel(r"$\pi\ ratio\  (\pi_{%s}/\pi_{%s})$"%(labelA,labelB),fontdict=font)
plt.ylabel(r'$Fst$',fontdict=font)
#fst_list=[cutoff_fst for i in range(dat.shape[0])]
#plt.axvline(cutoff_pi)
plt.axvline(cutoff_pi[0],ls='dashed',linewidth=1.5,color=cutline_color)
plt.axvline(cutoff_pi[1],ls='dashed',linewidth=1.5,color=cutline_color)
plt.axhline(cutoff_fst,ls='dashed',linewidth=1.5,color=cutline_color)
plt.xlim((-2,25))
plt.legend(fontsize=16)
yticks=fig2.get_yticklabels()
xticks=fig2.get_xticklabels()
fig2.get_yticklabels()[len(yticks)-2].set_visible(False)
fig2.get_xticklabels()[len(xticks)-1].set_visible(False)

fig3=fig.add_axes([0.8,0,0.2,0.8])
#plt.plot(fst_freq[:,1],fst_freq[:,0],color='black',linewidth=0.5,ls='dashed')
plt.fill_betweenx(fst_freq[fst_no,0],fst_freq[fst_no,1],0,color=fill_color2)
plt.fill_betweenx(fst_freq[fst_select_area,0],fst_freq[fst_select_area,1],0,color=fill_color1)
plt.xlabel(r'Frequency (%)',fontdict=font)
plt.axhline(cutoff_fst,ls='dashed',linewidth=1.5,color=cutline_color)
fig3.get_yticklabels()[0].set_visible(False)
fig3.get_yticklabels()[len(yticks)-2].set_visible(False)
#plt.ylim((-0.2,1))

fig3_2=plt.twiny()
plt.plot(fst_freq[:,2],fst_freq[:,0],color="maroon",linewidth=2)
#plt.xlim((-2,25))
#plt.xlabel(r'cumulative(%)',fontdict=font)
plt.xlim((0,100))
plt.ylim((-0.2,1))
yticks=fig3.get_yticklabels()
xticks=fig3.get_xticklabels()
fig3_2.get_xticklabels()[0].set_visible(False)

#fig=plt.figure(figsize=(10,10))

#plt.annotate()
if region!='':
    fig.savefig("%s_in_%s_fst_vs_pi.pdf"%(region,outfile),bbox_inches='tight')
else:
    fig.savefig(outfile+"_fst_vs_pi.pdf",bbox_inches='tight')
#fig3.get_xticklabels()[len(xticks)-1].set_visible(False)

