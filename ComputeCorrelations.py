### As the project is under development, you cannot access the code of the "ComputeCorelation" method.
### Please contact me to get the full version.

#@asaaidi 2020
# -*- coding: utf-8 -*-
import numpy as np
import math
from collections import defaultdict
import glob, os, subprocess
import argparse
import operator
import csv 
#import Prediction as PR

parser = argparse.ArgumentParser()                                               
parser.add_argument("--file", "-f", type=str, required=True)
#parser.add_argument("--profile", "-p", type=str, required=True)
parser.add_argument("--output", "-o", type=str, required=True)
parser.add_argument("--fasta", "-R", type=str, required=True)
parser.add_argument("--amplicon", "-a",  type=int, required=True)
args = parser.parse_args()

def ListBasePairsFromStruct(Struct):  # return dic={structure:[liste de pairs de base ],....}
    lista = []
    stack = []
    ListUmp =[]
    for i in range(len(Struct)):  # sequence length
	if  Struct[i] ==".":
	    ListUmp.append(i)
        if Struct[i] == '(':  # Opening base pair
            stack.append(i)
        elif Struct[i] == ')':  # Closing base pair
            k = stack.pop()
            lista.append((k, i))#0-based
	
    return lista,ListUmp

def GetPseudoknots(Struct):
    lista = []
    stack = []
    for i in range(len(Struct)):  # sequence length
        if Struct[i] == '<':  # Opening base pair
            stack.append(i)
        elif Struct[i] == '>':  # Closing base pair
            k = stack.pop()
            lista.append((k, i))#0-based
    return lista

def classification(mutationrate):
	R="Unkown"
	if mutationrate<0.02:
		R="L"
	if 0.02<=mutationrate<=0.04:
		R="M"
	if mutationrate>0.04:
		R="H"
	return R
	
	("P3a",(21,68,8))

if __name__=="__main__":
	#python2.7 ComputeCorrelationsMaximalhelicesElimination.py --file Incell/Modified/5S_INCELL1_5SFREQUENCY_COUNTex.txt  --output 2020-10-21 --fasta alignment_sequences/5S.fa --amplicon 20 --p /home/user/Documents/GeorgiaTech2020/Generate_normalized_prior_Input/2020-10-20/Profiling/5SProfiling.csv
	setofhelices=[]
	EmbeddedHelices==False
	MinLenHelix=0
	File=args.file #-f
	OutputFolder=args.output #-o
	Fastafile=args.fasta #-R
	Amplicon=args.amplicon #-a #20
	#Profile=args.profile #-p
	if not os.path.exists(OutputFolder):
    		os.makedirs(OutputFolder)
    	RNA=File.split("/")[-1].split(".")[0]
	if 1==1:
		
                outputIndexation=os.path.join(OutputFolder,RNA+"_MutationIndex.csv")
		#print outputIndexation
		######### Parse RNA sequence and structure
		fp=open(Fastafile)
		for i, line in enumerate(fp):
			if i==1:
				seq=line
			if i==2:
				SS=line
		fp.close()
                lenRNA=len(seq)

		try:
   			SS
		except NameError:# #structure not provided
			SS="".join(["na" for i in range(lenRNA)])	
		#parse MM UU MU UM file		
		MM = np.zeros((lenRNA,lenRNA)) #MM dict
		UU= np.zeros((lenRNA,lenRNA)) # uU dict
                UM= np.zeros((lenRNA,lenRNA))
		MU= np.zeros((lenRNA,lenRNA))
		N= np.zeros((lenRNA,lenRNA)) 
		Prior=np.zeros((lenRNA,lenRNA)) 
		ListPairs=ListBasePairsFromStruct(SS)[0]+ GetPseudoknots(SS)
		ListUnpaired=ListBasePairsFromStruct(SS)[1]
				

		with open(File) as f:  #i-1based 	 j-1based 	 UU 	 UM 	 MU 	 MM 	 #Reads 	 MM/#Reads  	 UU/#Reads 
		    next(f)
		    for line in f:
		       #(i,j,bb,kk,dd,aa,X,Y,Z)= line.split() # format with 9 columns 
		       (i,j,bb,kk,dd,aa)= line.split()
		       # change to 0-based
		       I=int(i)-1
		       J=int(j)-1
		       #print I,J,int(MM)
		       N[I,J]=int(bb+kk+dd+aa)
		       MM[I,J]= int(aa)
		       UU[I,J]= int(bb) #0-based
		       UM[I,J]=int(kk)
		       MU[I,J]= int(dd)

		with open (outputIndexation, 'w') as outfile:	
			outfile.write("%s,%s,%s,%s,%s\n"% ("position","nucleotide","structure","MutationRate","MutationRange"))				
			for i in range(Amplicon, lenRNA-Amplicon):
				if i in ListUnpaired:
					L="Up"
				else:
					if i+1 not in  ListUnpaired and i-1 not in ListUnpaired:
						L="P"
					else:
						L="HE"
				X=0
				if (UU[i,i]+MM[i,i])!=0:
					X=MM[i,i]/float(UU[i,i]+MM[i,i])
					
				outfile.write("%i,%s,%s,%.4f,%s\n"% (i+1,seq[i],L,X,classification(X)))

				####################### Compute Prior on the range [amplicon,lenRNA-amplicon]x[i+7 , lenRNA-amplicon]: 
				for j in range(i+7, lenRNA-Amplicon):
					Totalcount=UU[i,j]+MM[i,j]+UM[i,j]+MU[i,j]
					if (float(UU[i,i]+MM[i,i])!=0 and float(UU[j,j]+MM[j,j])!=0):
						Prior[i,j]=MM[i,i]/float(UU[i,i]+MM[i,i])*MM[j,j]/float(UU[j,j]+MM[j,j])*Totalcount
					
		outfile.close()
		
		Computecorrelations(MM,UU,Prior,ListPairs,lenRNA,File,MinLenHelix,setofhelices)
		
		print ("Prior and Relative differences have been successfully computed for "+(File.split("/")[-1]).split("_")[0]+ " RNA for the condition "+ (File.split("/")[-1]).split("_")[1] ) 
