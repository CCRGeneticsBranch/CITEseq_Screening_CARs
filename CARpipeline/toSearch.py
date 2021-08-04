import re
import numpy as np
import regex
import sys 
from Bio.Seq import Seq
import pandas as pd 
from Bio import AlignIO
from Bio.Align import AlignInfo
import time
import os 
import gzip

samFile = sys.argv[1]

df_sam =  pd.read_csv(samFile,delimiter = "\t",header=None,error_bad_lines=False,skiprows=2,warn_bad_lines=False,usecols = [0,1,2,3,4,5,6,7,8,9],engine='python' )
df_sam.columns = ["read", "one", "mapTo", "three", "four", "five", "six", "seven", "eight", "seq"]


def reverse_complement(dna):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return ''.join([complement[base] for base in dna[::-1]])



n = -1
for x in df_sam.seq.str.split('\n', 1):
	n +=1

	col1 = df_sam.iloc[n]["read"]
	col2 = str(df_sam.iloc[n]["one"]) 
	col3 = str(df_sam.iloc[n]["mapTo"]) 
	col4 = str(df_sam.iloc[n]["three"]) 
	col5 = str(df_sam.iloc[n]["four"]) 
	col6 = str(df_sam.iloc[n]["five"]) 
	col7 = str(df_sam.iloc[n]["six"]) 
	col8 = str(df_sam.iloc[n]["seven"]) 
	col9 = str(df_sam.iloc[n]["eight"])
	#print(type(col1)) 
	
	
	capture2 = regex.search(r'(\D{16})(\D{10})(TTTCTTATATGGG){s<=1,d<=1}(\D*)(CCCATATAAGAAA){s<=1,d<=1}(\D{10})(\D{16})',x[0])
	if capture2 != None:
		print(col1 + "FWDinREV" + "\t" + col2  + "\t" + col3  + "\t" + col4  + "\t" + col5  + "\t" + col6  + "\t" + col7  + "\t" + col8  + "\t" + col9  + "\t" + capture2[1] + capture2[2] + capture2[3] + capture2[4])
		print(col1 + "REVinFWD" + "\t" + col2  + "\t" + col3  + "\t" + col4  + "\t" + col5  + "\t" + col6  + "\t" + col7  + "\t" + col8  + "\t" + col9  + "\t" + reverse_complement(capture2[7]) + reverse_complement(capture2[6]) + capture2[5] + capture2[4])

	capture4 = regex.search(r'(\D*)(CCCATATAAGAAA){s<=1,d<=1}(\D{10})(\D{16})(\D*)(\D{16})(\D{10})(TTTCTTATATGGG){s<=1,d<=1}(\D*)',x[0])
	if capture4 != None and capture2  == None:	
		print(col1 + "REVfirst" + "\t" + col2  + "\t" + col3  + "\t" + col4  + "\t" + col5  + "\t" + col6  + "\t" + col7  + "\t" + col8  + "\t" + col9  + "\t" + reverse_complement(capture4[4]) + reverse_complement(capture4[3]) + capture4[2] + capture4[1])
		print(col1 + "FWDsecond" + "\t" + col2  + "\t" + col3  + "\t" + col4  + "\t" + col5  + "\t" + col6  + "\t" + col7  + "\t" + col8  + "\t" + col9  + "\t" + capture4[6] + capture4[7] + capture4[8] + capture4[9])

	capture1 = regex.search(r'(\D{16})(\D{10})(TTTCTTATATGGG){s<=1,d<=1}(\D*)(\D{16})(\D{10})(TTTCTTATATGGG){s<=1,d<=1}(\D*)',x[0])
	if capture1 != None and capture2 == None and capture4 == None:
		print(col1 + "FWD" + "\t" + col2  + "\t" + col3  + "\t" + col4  + "\t" + col5  + "\t" + col6  + "\t" + col7  + "\t" + col8  + "\t" + col9  + "\t" + capture1[1] + capture1[2] + capture1[3] + capture1[4])
		print(col1 + "FWD2" + "\t" + col2  + "\t" + col3  + "\t" + col4  + "\t" + col5  + "\t" + col6  + "\t" + col7  + "\t" + col8  + "\t" + col9  + "\t" + capture1[5] + capture1[6] + capture1[7] +  capture1[8])


	capture3 = regex.search(r'(\D*)(CCCATATAAGAAA){s<=1,d<=1}(\D{10})(\D{16})(\D*)(CCCATATAAGAAA){s<=1,d<=1}(\D{10})(\D{16})',x[0])
	if capture3 != None and capture2  == None and capture1 == None and capture4 == None:
		print(col1 + "REV" + "\t" + col2  + "\t" + col3  + "\t" + col4  + "\t" + col5  + "\t" + col6  + "\t" + col7  + "\t" + col8  + "\t" + col9  + "\t" + reverse_complement(capture3[4]) + reverse_complement(capture3[3]) + capture3[2] + capture3[1])
		print(col1 + "REV2" + "\t" + col2  + "\t" + col3  + "\t" + col4  + "\t" + col5  + "\t" + col6  + "\t" + col7  + "\t" + col8  + "\t" + col9  + "\t" + reverse_complement(capture3[8]) + reverse_complement(capture3[7]) + capture3[6] + capture3[5])
		
	captureFwd = regex.search(r'(\D{16})(\D{10})(TTTCTTATATGGG){s<=1,d<=1}(\D*)',x[0]) 
	if captureFwd != None and capture1 == None and capture2 == None and capture3 == None and capture4 == None:
		print(col1 + "FWDonly" + "\t" + col2  + "\t" + col3  + "\t" + col4  + "\t" + col5  + "\t" + col6  + "\t" + col7  + "\t" + col8  + "\t" + col9  + "\t" + captureFwd[1] + captureFwd[2] + captureFwd[3] + captureFwd[4])


	revCapture = regex.search(r'(\D*)(CCCATATAAGAAA){e<=0}(\D{10})(\D{16})',x[0])
	if revCapture != None and capture1 == None and capture2 == None and capture3 == None and captureFwd == None and capture4 == None:
		print(col1 + "REVonly" + "\t" + col2  + "\t" + col3  + "\t" + col4  + "\t" + col5  + "\t" + col6  + "\t" + col7  + "\t" + col8  + "\t" + col9  + "\t" + reverse_complement(revCapture[4]) + reverse_complement(revCapture[3]) + revCapture[2] + revCapture[1])
