#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2021-02-22 10:15:28
LastEditors: Li Fajin
LastEditTime: 2021-02-22 11:23:17
Description: file content
'''



import sys
import numpy as np
import pandas as pd
from optparse import OptionParser



def create_parser_for_mean_ratio():
	'''argument parser'''
	usage="usage: python %prog --control control.txt --treat treat.txt --kmer kmer.txt -o outprefix  "
	parser=OptionParser(usage=usage)
	parser.add_option("-c","--control",action="store",type="string",dest="control",
		help="Control ribosome density file on polypurine sequence motifs.")
	parser.add_option("-t","--treat",action="store",type="string",dest="treat",
		help="treatment ribosome density file on polypurine sequence motifs.")
	parser.add_option("-o","--otput_prefix",action="store",type="string",dest="output_prefix",
		help="Prefix of output files.[required]")
	parser.add_option("--kmer",action="store",type="string",dest="kmer",help="kmer motifs to keep!")
	return parser

def parseKmer(kmer):
    kmers=pd.read_csv(kmer)
    kmers=kmers.values
    return kmers

def parseMetaReads(metaReads,kmer):
    metaReadsDict={}
    kmerDict={}
    kmers=parseKmer(kmer)
    with open(metaReads,'r') as f:
        for line in f:
            if line.strip()=="":
                continue
            trans=line.strip().split("\t")[0]
            motif=line.strip().split("\t")[1]
            start=line.strip().split("\t")[2]
            stop=line.strip().split("\t")[3]
            reads=np.array([i for i in line.strip().split("\t")[4:] if i.strip()!=""]).astype(np.float)
            key=trans+":"+motif+":"+str(start)+"-"+str(stop)
            metaReadsDict[key]=reads
            if motif in kmers:
                kmerDict[key]=reads
            else:
                pass
    print("There are " +str(len(metaReadsDict))+" motifs in "+str(metaReads))
    print("But, only "+str(len(kmerDict))+" motifs selected!")
    return metaReadsDict,kmerDict

def calculateDensityRatio(controlmetaReadsDict,treatmetaReadsDict,output):
    common=set(controlmetaReadsDict.keys()).intersection(set(treatmetaReadsDict.keys()))
    print("Finally, There are "+str(len(common))+" motif in common used!")
    ratio_list=[]
    for motif in common:
        control=controlmetaReadsDict[motif]
        treat=treatmetaReadsDict[motif]
        ratio=(treat+1)/(control+1)
        ratio_list.append(ratio)
    ratio_list=np.array(ratio_list)
    ratio_data=pd.DataFrame(sum(ratio_list)/len(ratio_list))
    ratio_data.to_csv(output,index=None)

def main():
    parser=create_parser_for_mean_ratio()
    (options,args)=parser.parse_args()
    (control,treat,output,kmer)=(options.control,options.treat,options.output_prefix,options.kmer)
    print("Start...")
    controlmetaReadsDict,controlKmerDict=parseMetaReads(control,kmer)
    treatmetaReadsDict,treatKmerDict=parseMetaReads(treat,kmer)
    calculateDensityRatio(controlKmerDict,treatKmerDict,output)
    print("Finish!")

if __name__=="__main__":
    main()