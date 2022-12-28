#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2021-02-20 15:34:03
LastEditors: Li Fajin
LastEditTime: 2021-11-04 11:43:39
Description: file content
'''


import sys
import numpy as np
import pandas as pd

control=sys.argv[1]
treat=sys.argv[2]
output=sys.argv[3]


def parseMetaReads(metaReads):
    metaReadsDict={}
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
    print("There are " +str(len(metaReadsDict))+" motifs in "+metaReads+" !")
    return metaReadsDict

def calculateDensityRatio(controlmetaReadsDict,treatmetaReadsDict,output):
    common=set(controlmetaReadsDict.keys()).intersection(set(treatmetaReadsDict.keys()))
    print("There are "+str(len(common))+" motif in common!")
    fout=open(output+".array","w")
    ratio_list=[]
    for motif in common:
        control=controlmetaReadsDict[motif]
        treat=treatmetaReadsDict[motif]
        ratio=(treat+1)/(control+1)
        ratio_list.append(ratio)
        fout.write("%s\t" %(motif))
        for i in range(len(ratio)):
            fout.write("%s\t" %(str(ratio[i])))
        fout.write("\n")
    fout.close()
    ratio_list=np.array(ratio_list)
    ratio_data=pd.DataFrame(sum(ratio_list)/len(ratio_list))
    ratio_data.to_csv(output,index=None)

def main():
    controlmetaReadsDict=parseMetaReads(control)
    treatmetaReadsDict=parseMetaReads(treat)
    calculateDensityRatio(controlmetaReadsDict,treatmetaReadsDict,output)

if __name__=="__main__":
    main()