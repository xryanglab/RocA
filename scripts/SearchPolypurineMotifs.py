#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2021-01-14 21:05:14
LastEditors: Li Fajin
LastEditTime: 2021-01-15 13:27:36
Description: search polypurine motifs
'''


import re
import sys
import numpy as np
import pandas as pd
from itertools import groupby,chain
from collections import defaultdict
from optparse import OptionParser

def create_parser_for_poly_purine():
	'''argument parser.'''
	usage="usage: python %prog [options]"
	parser=OptionParser(usage=usage)
	parser.add_option("-i","--input", action="store",type="string",default=None,dest="fastaFile",
			help="Input file(s) in fasta format. ")
	parser.add_option("-o","--otput_prefix",action="store",type="string",dest="output_prefix",
			help="Prefix of output files.[required]")
	parser.add_option('--kmer',action="store",type="int",dest="kmer",default=4,
			help="Length of kmer used for statistics. default=%default")
	parser.add_option("--base",action='store',type='string',dest='base',default='AG',help='Poly base. If AG, output poly-purine.default=%default')

	return parser

def fastaIter(transcriptFile):
	'''
	This function is used to get a dict of transcript sequence
	'''
	fastaDict={}
	f=open(transcriptFile,'r')
	faiter=(x[1] for x in groupby(f,lambda line: line.strip()[0]==">")) ## groupby returns a tuple (key, group)
	for header in faiter:
		geneName=header.__next__().strip(">").split(" ")[0].strip()
		seq=''.join(s.strip() for s in faiter.__next__())
		flag=0
		for nt in ['I','K','M','R','S','W','Y','B','D','H','V','N','X']:
			if nt in seq:
				flag+=1
				flag_nt=nt
		if flag != 0:
			print(geneName+" filtered"+"--"+"There is a ambiguous nucleotide",flag_nt,"in your sequence")
			continue
		fastaDict[geneName]=seq
	return fastaDict


def StatisticsPolyPurine(fastaFile,kmer,Base,output):
	'''
	statistics of  polypurine motifs.
	'''
	readSeqs=fastaIter(fastaFile)
	fout=open(output,'w')
	for read in readSeqs:
		read_seq=readSeqs[read]
		read_length=len(read_seq)
		for i in range(0,read_length): ## change the region as you wish
				motif=read_seq[i:(i+kmer)]
				if len(motif) < kmer:
					continue
				if all([base in list(Base) for base in list(motif)]):
					tmp=[read.strip(),read_length,motif,i,i+kmer]
					fout.write("\t".join(str(i) for i in tmp))
					fout.write("\n")
				else:
					continue


def parse_args_for_poly_purine():
	parsed=create_parser_for_poly_purine()
	(options,args)=parsed.parse_args()
	print("start...")
	StatisticsPolyPurine(options.fastaFile,options.kmer,options.base,options.output_prefix+"_poly"+options.base+"_"+str(options.kmer)+"_mer.txt")
	print("Finish!")



def main():
	"""main program"""
	parse_args_for_poly_purine()

if __name__ == "__main__":
		main()













