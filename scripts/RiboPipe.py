#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2021-02-15 13:23:59
LastEditors: Li Fajin
LastEditTime: 2022-05-28 12:21:37
Description: file content
'''

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os


'''
###################################
Date: 20220528.15:00
all polarity with Tranlation initiation and elongation

19636 trans x 49 conditionPairs
###################################
'''
ConditionPairs_of_initiation_and_elongation=pd.read_csv("../data/TranslationRelatedConsitionPairs.txt",sep="\t")
polarity_with_initiation_and_elongation=pd.read_csv("../data/total_merged_diff_polarity_202102.txt",sep="\t",index_col=0)
RPKM_with_initiation_and_elongation=pd.read_csv("../data/total_merged_cds_level_RPKM_202102.txt",sep="\t",index_col=0)
RPKM_with_initiation_and_elongation.apply(np.mean, axis=1).describe() ## trans describe
RPKM_with_initiation_and_elongation.apply(np.mean, axis=0).describe() ## sample describe

coverage_with_initiation_and_elongation=pd.read_csv("../data/total_merged_read_coverage_202102.txt",sep="\t",index_col=0)
coverage_with_initiation_and_elongation_sorted=coverage_with_initiation_and_elongation[sorted(coverage_with_initiation_and_elongation.columns)]
coverage_with_initiation_and_elongation_sorted.apply(np.mean, axis=1).describe() ## trans describe
coverage_with_initiation_and_elongation_sorted.apply(np.mean, axis=0).describe() ## sample describe


RPKM_with_initiation_and_elongation_filtered_gt0_in_each_samples=RPKM_with_initiation_and_elongation.loc[RPKM_with_initiation_and_elongation.apply(np.min,axis=1)>=1,:]
coverage_with_initiation_and_elongation_filtered_mean_of_each_trans_gt_pectMean=coverage_with_initiation_and_elongation_sorted.loc[coverage_with_initiation_and_elongation_sorted.apply(np.mean,axis=1)>=coverage_with_initiation_and_elongation_sorted.apply(np.mean,axis=1).describe()['mean'],:]
trans_in_common_associated_with_translation=set(RPKM_with_initiation_and_elongation_filtered_gt0_in_each_samples.index).intersection(set(coverage_with_initiation_and_elongation_filtered_mean_of_each_trans_gt_pectMean.index))
polarity_with_initiation_and_elongation_filtered=polarity_with_initiation_and_elongation.loc[trans_in_common_associated_with_translation,:]


polarity_with_initiation_and_elongation_filtered_vars=np.var(polarity_with_initiation_and_elongation_filtered,axis=1)
polarity_with_initiation_and_elongation_filtered_vars.describe()

trans_filtered_by_polarityVars=polarity_with_initiation_and_elongation_filtered.iloc[np.where(polarity_with_initiation_and_elongation_filtered_vars>polarity_with_initiation_and_elongation_filtered_vars.describe()[1]+polarity_with_initiation_and_elongation_filtered_vars.describe()[2])[0],:]
polarity_with_initiation_and_elongation_filtered_filtered_byPolarityVars=polarity_with_initiation_and_elongation_filtered.drop(trans_filtered_by_polarityVars.index)


plt.rc('font',weight='bold')
sns.set(font_scale=0.5)
fig=plt.figure()
col_c=dict(zip(ConditionPairs_of_initiation_and_elongation['TranslationType'].unique(), ['#2a93d4','#ffb5ba','#79bd9a']))
col_colors=ConditionPairs_of_initiation_and_elongation['TranslationType'].map(col_c)
ax=sns.clustermap(polarity_with_initiation_and_elongation_filtered_filtered_byPolarityVars,method ='ward',metric='euclidean',row_cluster=True,col_cluster=True,
                vmin=-0.5,vmax=0.5,center=0,cmap = 'RdBu_r',col_colors=col_colors.values,yticklabels=False,xticklabels=True)
ax.cax.set_visible(False)
plt.savefig("py20220528-polarity-with-Translation-initiation-and-elongation-no-labels.pdf")

