#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2021-01-15 13:22:06
LastEditors: Li Fajin
LastEditTime: 2021-01-15 14:05:49
Description: calculate motif score for polypurine motifs
'''


import sys
from itertools import groupby,chain
from collections import defaultdict
from optparse import OptionParser

def create_parser_for_poly_purine_score():
	'''argument parser.'''
	usage="usage: python %prog [options]"
	parser=OptionParser(usage=usage)
	parser.add_option("-i","--input", action="store",type="string",dest="fastaFile",
			help="Input file(s) in fasta format. ")
	parser.add_option("-m","--motifs", action="store",type="string",dest="motifFile",
			help="Collapsed motif file. ")
	parser.add_option("-o","--otput_prefix",action="store",type="string",dest="output_prefix",
			help="Prefix of output files.[required]")
	return parser

def GetAllBases(fastaFile):
    all_bases=0
    with open(fastaFile,'r') as f:
        for line in f:
            if line.strip()=="":
                continue
            if line.strip().startswith(">"):
                continue
            seq=line.strip()
            seq_len=len(seq)
            all_bases+=seq_len
    return all_bases

def GetScores(fastaFile,motifFile,output_prefix):
    all_bases=GetAllBases(fastaFile)
    with open(motifFile,'r') as fin,open(output_prefix+"_motif_scores.txt",'w') as fout:
        fout.write("motif"+"\t"+"length"+"\t"+"number"+"\t"+"all_bases"+"\t"+"Ref_all_bases"+"\t"+"score"+"\n")
        for line in fin:
            if line.strip()=="":
                continue
            kmer=line.strip().split("\t")[0]
            kmer_len=len(kmer)
            number=int(line.strip().split("\t")[1])
            score=number*kmer_len/all_bases
            fout.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(kmer,str(kmer_len),str(number),str(number*kmer_len),str(all_bases),str(score)))



def parse_args_for_poly_purine_score():
	parsed=create_parser_for_poly_purine_score()
	(options,args)=parsed.parse_args()
	print("start...")
	GetScores(options.fastaFile,options.motifFile,options.output_prefix)
	print("Finish!")



def main():
	"""main program"""
	parse_args_for_poly_purine_score()

if __name__ == "__main__":
		main()