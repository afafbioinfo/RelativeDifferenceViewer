
# RelativeDifferenceViewer

RelativeDifferenceViewer is a program for viewing raw mutation data and the corresponding relative differences.
It outputs raw mttaion count and Relative difference heatmaps.

## Installing RelativeDifferenceViewer

 RelativeDifferenceViewer depends mainly on Python 2.7 and R.

## Executing RelativeDifferenceViewer:

Rscript --vanilla Distributions2D.R --RD1=RNA_Experiment1_RMRPFREQUENCY_COUNTex_Rawdata.csv  --MD1=RNA_Experiment1_RNAFREQUENCY_COUNTex_MutationIndex.csv  --RD2=RNA_Experiment2_RNAPFREQUENCY_COUNTex_Rawdata.csv  --MD2=RNA_Experiment2_RNAFREQUENCY_COUNTex_MutationIndex.csv

## Input files

Please follow the two main steps bellow to get the input file in the required format:

  - Get raw mutations
  run shapemapper then run  Ringmapper while replacing ringmapper.py file with the new edits.
  
Example of use:

./RingMapper-master/ringmapper.py output/ /Pipeline_Modified_parsed.mut  RNAname  --fasta  fastafile.fa    --untreated  output/Pipeline_Untreated_$1_parsed.mut  --window 1 --metric g

You will notice  2 additional  files that will appear  FREQUENCY_COUNTex (for DMS experiment) and  FREQUENCY_COUNTbg (for the control)).

Those are the two files that we will use for the computation with the fasta file.


  - To compute mutation rates and correlations: you can use the python script , for the arguments  please type on the terminal  python2.7 ComputeCorrelations.py -h 
  
Example of use:

ComputeCorrelations.py --file RNAname_Experiment1_RNAnameFREQUENCY_COUNTex.txt --output 2020-12-23 --fasta  RNAname.fa --amplicon 20

Note that the amplicon here was set to 20; this means that the sequence will be cut by 20 nucleotides from each side.

Two files will be generated  RNAname_Experiment1_RNAnameFREQUENCY_COUNTex_MutationIndex.csv, RNAname_Experiment1_RNAnameFREQUENCY_COUNTex_Rawdata.csv that you can find in 2020-12-23 folder.

  - As a last step, use the R script , please type  Rscript --vanilla  -h on the terminal to see the required arguments .
Note that here you are comparing two conditions, so you need to run step 2  for each condition (Cell free and In cell, 10mol and 0 mol ...).
Example of use:
Rscript --vanilla Distributions2D.R --RD1=rpsM_Cell-Free_RMRPFREQUENCY_COUNTex_Rawdata.csv  --MD1=rpsM_Cell-Free_rpsMFREQUENCY_COUNTex_MutationIndex.csv  --RD2=rpsM_INCELL1_rpsMPFREQUENCY_COUNTex_Rawdata.csv  --MD2=rpsM_INCELL1_rpsMFREQUENCY_COUNTex_MutationIndex.csv

## OUtputs:
4 plots: 
Grid.png, Mutation_rate.png, Raw_mutation_count.png and Rdiff.png

