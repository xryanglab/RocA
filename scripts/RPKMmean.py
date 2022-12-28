#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2021-02-14 14:24:27
LastEditors: Li Fajin
LastEditTime: 2021-02-14 14:35:54
Description: file content
'''


import numpy as np
import pandas as pd
import sys
from optparse import OptionParser
from collections import defaultdict


def create_parser():
	'''argument parser'''
	usage="usage: python %prog [options]" +"\n" + __doc__+"\n"
	parser=OptionParser(usage=usage)
	parser.add_option("-i","--input",action="store",type="string",dest="RPKM",help="Input file in txt format.")
	parser.add_option("-o","--output",action="store",type="string",dest="output_prefix",help="Prefix of output files.[required]")
	parser.add_option('-g','--group',action="store",type="string",dest="group_name",help="Group name of each group separated by comma. e.g. 'si-control,si-eIF3e'")
	parser.add_option('-r','--replicate',action="store",type="string",dest="replicate_name",help="Replicate name of each group separated by comma. e.g. 'si_3e_1_80S,si_3e_2_80S__si_cttl_1_80S,si_ctrl_2_80S'")
	return parser

def calculate_mean_RPKM(data,groups,replicates,output_prefix):
	'''calculate the mean polarity of each group based on different replicates'''
	transcriptList=data.index.values
	labels_dict={}
	data_dict={}
	data_mean_dict=defaultdict(dict)
	for g,r in zip(groups,replicates):
		labels_dict[g]=r.strip().split(',')
	## separate different groups
	for g in groups:
		# data_dict[g]=data.loc[:,labels_dict[g]]
		data_dict[g]=data.reindex(columns=labels_dict[g])

	## calculate the mean polarity between different replicates
	for g in groups:
		for trans in transcriptList:
			index=np.where(~np.isnan(data_dict[g].loc[trans,:]))[0]
			if len(index) == len(labels_dict[g]):
				polarity=np.mean(data_dict[g].loc[trans,data_dict[g].columns[index]])
				data_mean_dict[g][trans]=polarity
			elif len(index) < len(labels_dict[g]) and len(index) > 1:
				polarity=np.mean(data_dict[g].loc[trans,data_dict[g].columns[index]])
				data_mean_dict[g][trans]=polarity
			elif len(index) == 1:
				polarity=np.mean(data_dict[g].loc[trans,data_dict[g].columns[index]])
				data_mean_dict[g][trans]=polarity
			elif len(index) == 0:
				continue
			else:
				raise KeyError("Key error, please check your data input!")

	## transform the dict to a python dataframe
	for g in groups:
		data_mean_dict[g]=pd.DataFrame(data_mean_dict[g],index=[g]).T

	## concatenate different data frame
	data_mean=pd.concat([v for v in data_mean_dict.values()],axis=1,sort=True)
	## write the mean density file
	data_mean.to_csv(output_prefix+"_mean_RPKM_dataframe.txt",sep="\t")


def main():
	'''main function'''
	parser=create_parser()
	(options,args)=parser.parse_args()
	(input,output_prefix,groups,replicates)=(options.RPKM,options.output_prefix,options.group_name.strip().split(','),options.replicate_name.strip().split('__'))
	data=pd.read_csv(input,sep="\t",index_col=0)
	samples=data.columns
	## calculate the mean polarity
	calculate_mean_RPKM(data,groups,replicates,output_prefix)

if __name__ == "__main__":
    main()