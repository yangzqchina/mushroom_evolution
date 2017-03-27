# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 19:24:32 2016

@author: zqyang
"""

import pandas as pd
import numpy as np
import argparse
from scipy import stats
from mne.stats import bonferroni_correction

parses=argparse.ArgumentParser(description='Use BDP file to make Chisq-test')
parses.add_argument('--bdp',dest='bdp',default=None,help='BDP File. First line\
 is column name.')
parses.add_argument('--pop1',dest='pop1',default=None,help='population1 File.\
One line one sample name')
parses.add_argument('--pop2',dest='pop2',default=None,help='population2 File.\
One line one sample name.')
parses.add_argument('--name1',dest='n1',default='pop1',type=str,help='the name of \
population1')
parses.add_argument('--name2',dest='n2',default='pop2',type=str,help='the name of \
population2')
parses.add_argument('-o','--output',dest='outfile',default='bdp',\
help='output file name. if you set results, it will produce a file named \
results.chisq.csv')

args=parses.parse_args()

def cdp_chisq(bdpfile,pop1,pop2,n1,n2):
    bdp=pd.read_table(bdpfile,sep=' ')
    pop1=pd.read_table(pop1,header=None).loc[:,0].tolist()
    n_pop1=len(pop1)
    pop2=pd.read_table(pop2,header=None).loc[:,0].tolist()
    n_pop2=len(pop2)
    bdp_4h=bdp.columns.tolist()[:4]
    #p1_dat=bdp.loc[:,bdp_4h]
    #p2_dat=bdp.loc[:,bdp_4h]
    p1=bdp.loc[:,pop1]
    p2=bdp.loc[:,pop2]
    p1[p1==1]=0
    p1[p1==2]=1
    p2[p2==1]=0
    p2[p2==2]=1
    results=bdp.loc[:,bdp_4h]
    results.loc[:,n1+'_yes']=p1.sum(axis=1)
    results.loc[:,n1+'_no']=n_pop1-p1.sum(axis=1)
    results.loc[:,n2+'_yes']=p2.sum(axis=1)
    results.loc[:,n2+'_no']=n_pop2-p2.sum(axis=1)
    chisq_results=[]
    fisher_results=[]
    for x in results.index:
        contingency=[results.loc[x,[n1+'_yes',n1+'_no']],\
            results.loc[x,[n2+'_yes',n2+'_no']]]
        try:
            chisq_results.append(stats.chi2_contingency(contingency)[:2])
            fisher_results.append(stats.fisher_exact(contingency))
        except:
            chisq_results.append(np.array([np.NaN,np.NaN]))
            fisher_results.append(np.array([np.NaN,np.NaN]))
    chisq_results=np.array(chisq_results)
    fisher_results=np.array(fisher_results)
    results.loc[:,'chisq']=chisq_results[:,0]
    results.loc[:,'chisq_p-value']=chisq_results[:,1]
    results.loc[:,'oddsratio']=fisher_results[:,0]
    results.loc[:,'fisher_p-value']=fisher_results[:,1]
    results.loc[:,'bonferroni_chisq']=bonferroni_correction(results.loc[:,'chisq_p-value'])[1]
    results.loc[:,'bonferroni_fisher']=bonferroni_correction(results.loc[:,'fisher_p-value'])[1]
    results.loc[results.loc[:,'bonferroni_chisq']>1,'bonferroni_chisq']=1
    results.loc[results.loc[:,'bonferroni_fisher']>1,'bonferroni_fisher']=1
    return(results)

bdpfile=args.bdp
pop1=args.pop1
pop2=args.pop2
n1=args.n1
n2=args.n2
results=cdp_chisq(bdpfile,pop1,pop2,n1,n2)    
results.to_csv(args.outfile+'.chisq.csv',index=False)    
