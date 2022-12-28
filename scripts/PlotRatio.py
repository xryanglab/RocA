#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2021-02-20 15:39:33
LastEditors: Li Fajin
LastEditTime: 2021-02-20 20:31:32
Description: file content
'''

import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from optparse import OptionParser
from optparse import OptionParser

def create_parser_for_ratioPlot():
	'''argument parser'''
	usage="usage: python %prog -i file1.txt,file2.txt,file3.txt -o outprefix --format pdf -u 50 -d 50 -t file1,file2,file3 "
	parser=OptionParser(usage=usage)
	parser.add_option("-i","--input",action="store",type="string",dest="Ratio",
		help="ratio files. multiple files are separated by comma.")
	parser.add_option("-o","--otput_prefix",action="store",type="string",dest="output_prefix",
		help="Prefix of output files.[required]")
	parser.add_option("-f","--format",action="store",type="string",dest="format",default='pdf',
		help="Plot format of metaplot.[pdf/png/jpg/tiff/]. default=%default.[required]")
	parser.add_option("-u","--upstream",action="store",type="int",dest="upstream",default=100,
		help="upstream base number from first base of polypurine sequence motifs. Consistent with that in RiboDensityAroundPolyPurineMotifs.py.")
	parser.add_option("-d","--downstream",action="store",type="int",dest="downstream",default=100,
		help="downstream base number from 5'end of exon. Consistent with that in CalDepthOfFiveEnds.py")
	parser.add_option("-t","--legend",action="store",type="string",dest="legend",
		help="Legends separated by comma corresponding to --input")
	parser.add_option("--ymax",action="store",type="float",dest="ymax",default=None,help="The max of ylim. default=%default")
	parser.add_option("--ymin",action="store",type="float",dest="ymin",default=None,help="The min of ylim. default=%default")
	parser.add_option("--log2",action="store_true",dest="log2",help="log2 or not")
	parser.add_option("--axvline",action="store",type="float",default=None,dest="axvline",help="plot axvline")
	parser.add_option("--start",action="store",type="int",dest="start_position",default=3,help="The start position need to be averaged.default=%default")
	parser.add_option("--window",action="store",type="int",dest="window",default=7,help="The length of silde window. ddefault=%default")
	parser.add_option("--step",action="store",type='int',dest="step",default=1,help="The step length of slide window. default=%default")
	parser.add_option("--slide-window",action="store_true",dest="slideWindow",help="Using slide window to average the density.")

	return parser



def parseRatio(Ratio,legends):
	metaFiles=Ratio.strip().split(",")
	metaLegend=legends.strip().split(",")
	meanReadsDict={}

	for i in range(len(metaFiles)):
		tmpfile=metaFiles[i]
		tmplegend=metaLegend[i]
		tmpReadsList=pd.read_csv(tmpfile)
		meanReadsDict[tmplegend]=np.array(tmpReadsList.iloc[:,0])
	return meanReadsDict

def slide_window_average(meanReadsDict,upstream,downstream,start,window,step):
	''' Used for calculating mean density with a slide window'''
	winLen=upstream+downstream+1
	mean_data={}
	for k,v in meanReadsDict.items():
		tmp_data=np.zeros(winLen)
		tmp_data[0:int(start)]+=v[0:int(start)]
		tmp_data[-int(start):]+=v[-int(start):]
		for j in np.arange(start,winLen-start,step):
			tmp_data[j]+=np.mean(v[(j-int((window-1)/2)):(j+int((window-1)/2))])
		mean_data[k]=tmp_data
	return mean_data


def metaplotForReadsOnMotifs(meanReadsDict,Format,upstream,downstream,ymin,ymax,outprefix,axvline,log2):
	plt.rc('font',weight='bold')
	winLen=upstream+downstream+1
	text_font={"size":20,"family":"Arial","weight":"bold"}
	legend_font={"size":20,"family":"Arial","weight":"bold"}
	if len(meanReadsDict) <=8:
		colors=["b","orangered","green","c","m","y","k","w"]
	else:
		colors=sns.color_palette('husl',len(meanReadsDict))

	if not log2:
		pass
	else:
		for k in meanReadsDict.keys():
			meanReadsDict[k]=np.log2(meanReadsDict[k]+1)

	with PdfPages(outprefix + "_metaplots.pdf") as pdf:
		figs=[]
		i=0
		fig=plt.figure(figsize=(16,8))
		ax=fig.add_subplot(111)
		for k in meanReadsDict.keys():
			i+=1
			plt.plot(np.arange(0,winLen),meanReadsDict[k],color=colors[i],linewidth=2,label=k)

		ax.set_ylabel("Enrichment Ratio",fontdict=text_font)
		ax.set_xticks(np.arange(0,winLen,50))
		ax.set_xticklabels((np.arange(0,winLen,50)-upstream))
		ax.axvline(upstream,color="r",dashes=(3,2),alpha=0.5,label="x=0")
		if axvline:
			ax.axvline(upstream+axvline,color="b",dashes=(3,2),alpha=0.5,label="x="+str(axvline))
		else:
			pass
		ax.set_xlabel("Distance from 1st base of motifs (nt)",fontdict=text_font)
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_linewidth(2)
		ax.spines["left"].set_linewidth(2)
		ax.tick_params(which="both",width=2,labelsize=20)
		if not ymin and not ymax:
			pass
		elif not ymin and ymax:
			ax.set_ylim(0,ymax)
		elif ymin and not ymax:
			raise IOError("Please offer the ymax parameter as well!")
		elif ymin and ymax:
			ax.set_ylim(ymin,ymax)
		else:
			raise IOError("Please enter correct ymin and ymax parameters!")
		plt.legend(loc="best",prop=legend_font)
		plt.tight_layout()
		figs.append(fig)
		for fig in figs:
			pdf.savefig(fig)
		plt.close()


def main():
	parser=create_parser_for_ratioPlot()
	(options,args)=parser.parse_args()
	print("start...")
	meanReadsDict=parseRatio(options.Ratio,options.legend)
	if options.slideWindow:
		meanReadsDict_average=slide_window_average(meanReadsDict,options.upstream,options.downstream,options.start_position,options.window,options.step)
		metaplotForReadsOnMotifs(meanReadsDict_average,options.format,options.upstream,options.downstream,options.ymin,options.ymax,options.output_prefix+"_average",options.axvline,options.log2)
		metaplotForReadsOnMotifs(meanReadsDict,options.format,options.upstream,options.downstream,options.ymin,options.ymax,options.output_prefix,options.axvline,options.log2)
	else:
		metaplotForReadsOnMotifs(meanReadsDict,options.format,options.upstream,options.downstream,options.ymin,options.ymax,options.output_prefix,options.axvline,options.log2)

	print("Finish!")

if __name__=="__main__":
	main()