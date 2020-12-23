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

def Computecorrelations(MM,UU,Prior,Pairs,lenRNA,Name,MinLenHelix,setofhelices):
	Norm= np.zeros((lenRNA,lenRNA)) 
	for i in range( len(MM)):
		for j in range(len(MM)):
			if Prior[i,j]!=0:
				Norm[i,j]=(MM[i,j]-Prior[i,j])/float(Prior[i,j])
	#set diagonals to Zero
        for i in range( len(MM)):
		MM[i,i]=np.nan
		Norm[i,i]=0
	
	Norm[Norm==0]=np.nan # set null normalized values to NA because Norm was initially defined as matrix of zeros
	######### Output MM count and Rltdiff
	with open (os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+"_Rawdata.csv"), 'w') as outfile: 
		outfile.write("%s,%s,%s,%s,%s\n"% ( "i","j","ncts","MM","Rltdiff" ))
		for i in range( Amplicon, lenRNA-Amplicon):
			for j in range( i+7,lenRNA-Amplicon):
				outfile.write("%i,%i,%s,%.4f,%.4f \n"% (i+1,j+1,str(seq[i]+seq[j]),round(MM[i,j],4),round(Norm[i,j],4)))
	outfile.close()
	######## output data for specific helices
	'''
	outputcorrelation=os.path.join(OutputFolder,(Name.split("/")[-1]).split(".")[0]+"_Norm_Featured_helices.csv")
	Lista=[]
	for Helix in setofhelices:
		if Helix[1][2]>MinLenHelix: 	
			listtuple=[]	
			i=Helix[1][0]-1
			j=Helix[1][1]-1
			l=Helix[1][2]
			for k in range(l):# neigborhood [i-1,i+1]x[j-1,j+1]
				for (n,t) in [(i+k,j-k),(i+k,j-k-1),(i+k,j-k+1),(i+k-1,j-k),(i+k-1,j-k-1),(i+k-1,j-k+1),(i+k+1,j-k),(i+k+1,j-k-1),(i+k+1,j-k+1)]:
					if j+1<lenRNA and math.isnan(Norm[n,t])!=True and (n,t,Norm[n,t])not in listtuple:  
						Lista.append((Helix[0],n,t,Norm[n,t],l))
						#Listadit.append(((Helix[0],l),n,t,Norm[n,t]))
						listtuple.append((n,t,Norm[n,t]))	
	
        # Get the #of common pairs between helices NA not included.
	dicto= defaultdict()
	dictobis= defaultdict() #oneadditional element that is the length of the helix
	CommonPairs=defaultdict()
	#SortedAverageHelix=defaultdict()
	#initialization
	for elem in Lista:
		dicto[elem[0]]=list()
		dictobis[elem[0]]=list()
		CommonPairs[elem[0]]=list()	
	#incrementation
	for elem in Lista:
		dicto[elem[0]].append((elem[1],elem[2],elem[3]))
		dictobis[elem[0]].append((elem[1],elem[2],elem[3],elem[4]))
	Averagehelix=defaultdict() 
	for key in dictobis:
		
		toaverage=[]
		## CC content and AA content 
		listAA=0
		listCC=0
		for elem in dictobis[key]:
			# verification if key=="eP19plus1BP":
			#	print elem, seq[elem[0]],seq[elem[1]]
			if seq[elem[0]].upper()=="C" and seq[elem[1]].upper()=="C":
				listCC+=1
			if seq[elem[0]].upper()=="A" and seq[elem[1]].upper()=="A":
				listAA+=1
			#print "loloo",elem[0],seq[elem[0]],seq[elem[1]]
			toaverage.append(elem[2])
		#tuple of average, nbr of observations, expected number of observations
		Averagehelix[key]=(round(np.mean(toaverage),3),len(toaverage),dictobis[key][0][3]*5+4,listCC,listAA)
	
	
	
	#1 ----------first filtering layer not enough data, the nbr of observations is strctely less than the half of the expected number of observations 
	for elem in Averagehelix.keys():
		
		if Averagehelix[elem][1]<Averagehelix[elem][2]/2:
			del Averagehelix[elem]
	#print "After first filtering", Averagehelix
	#print Averagehelix
	#2 ------------ look for common pairs and remove those that are dominated by other helices
	for key1 in Averagehelix:
		firstlist=[elem for elem in dicto[key1]]
		for key2 in Averagehelix:
			if key2!=key1:
				secondlist=[elem for elem in dicto[key2]]
				#get pairs in common btw helices
				
				if len(set(firstlist) & set(secondlist))!=0:
					CommonPairs[key1].append((key2,len(set(firstlist) & set(secondlist))))
	#print  set(firstlist) ,"\n",  set(secondlist),"\n", set(firstlist) & set(secondlist)
	#print CommonPairs	
	removed=[]
	for Helix in Averagehelix:
		for Helixrelated in CommonPairs[ Helix]:
			Helice=Helixrelated[0]	
			commonpairs=Helixrelated[1]
			# if dominated by another helix in average and most than half of observed values are in common, remove this helix
			#print Averagehelix[Helice][0],Averagehelix[Helix][0],commonpairs,Averagehelix[Helix][1],Averagehelix[Helice][0]<Averagehelix[Helix][0] , commonpairs>=Averagehelix[Helix][1]
			if Averagehelix[Helix][0]<Averagehelix[Helice][0] and 2*commonpairs>=Averagehelix[Helix][1]:
				
				removed.append(Helix)
	for Helix in list(set(removed)): # a trick to remove redundant helices
		if Averagehelix[Helix]:
			del Averagehelix[Helix]
			
	#sort the dictionary:
	SortedAverageHelix = Averagehelix.items()
	SortedAverageHelix.sort(key=lambda x:-x[1][0])
	#print  "initila,",SortedAverageHelix
	#print "Averagehelix \n",Averagehelix	
		
	with open (outputcorrelation, 'w') as outfile:
		outfile.write("%s,%s,%s,%s,%s,%s,%s,%s\n"% ( "Dataset","Helix","i","j","seqij","Rltdiff","ijdist","commonpairs" ))
		for elem in Lista:
			outfile.write("%s,%s,%s,%s,%s,%.4f,%i,%s \n"% ( (Name.split("/")[-1]).split("_")[1],elem[0],str(elem[1]+1),str(elem[2]+1),str(seq[elem[1]]+seq[elem[2]]),round(elem[3],4),elem[2]-elem[1]+1,CommonPairs[elem[0]]))
	'''
	Order=[]
	'''
	with open (outputcorrelation, 'w') as outfile:
		outfile.write("%s,%s,%s,%s,%s,%s\n"% ( "Helix","averageRLtdiff","nbrobservations","Length","nbrCCpairs","nbrAApairs" ))
		for elem in SortedAverageHelix:
			Order.append(elem[0])
			outfile.write("%s,%.4f,%i,%i,%i,%i \n"% (elem[0], elem[1][0] , elem[1][1], (elem[1][2]-4)/5, elem[1][3], elem[1][4]))
	outfile.close()
	'''
	return Order
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
	
def EmbeddedHelices(listofHelices): # Embedded helices are defined as follow: H2 is inside the space [H1+L1,H1-L1] the maximal distance is two nucleotides from each side.
	List=[]
	#listofHelices=[("Cx",(79,97,8)),("Cx2",(83,92,4)),("Cx3",(80,98,4))]
	for H1 in listofHelices:
		for H2 in listofHelices:
			#print  H1,H2,H2[1][0]>H1[1][0],H2[1][0]-(H1[1][0]+H1[1][2])<3,H1[1][1]-H1[1][2]-H2[1][1]<3
			if H2[1][0]>=H1[1][0]+ H1[1][2] and H2[1][1]<=H1[1][1]- H1[1][2]and  H2[1][0]-H1[1][0]-H1[1][2]<3 and H1[1][1]-H1[1][2]-H2[1][1]<3:
				# length of the new embeededhelix is the sum of the two helices lengths and the maximal gap value that is less than 3
				LengthEmbeddedhelix=H1[1][2]+H2[1][2]+np.amax([H2[1][0]-H1[1][0]-H1[1][2],H1[1][1]-H1[1][2]-H2[1][1]])
				List.append((str(H1[0]+H2[0]),(H1[1][0],H1[1][1],LengthEmbeddedhelix)))
	print "List embeeded", List
	return List

    
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
                lenRNA=len(seq)-1

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
