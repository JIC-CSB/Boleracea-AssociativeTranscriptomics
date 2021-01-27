# B. oleracea Associative Transcriptomics

___

## Instructions

### Download
Download the files.

The SNP folder should contain:
- SNPGrapher.R
- ArabidopsisAnnotation.tsv
- GrapherDirections.txt

The GEM folder should contain:
- GEMGrapher.R
- GEMRegress.R
- ArabidopsisAnnotation.tsv
- BolRPKMTASSEL.txt
- GrapherDirections.txt
- Marker_to_At_Bol.csv

### SNP
To graph your TASSEL outputs, open SNPGrapher.R in RStudio.
Set your working directory to that which contains SNPGrapher.R.
Press source to run.
You will be prompted to select the model you used for GWAS analysis, either GLM or MLM.
You will then need to navigate to and select your ArabidopsisAnnotation.tsv file and your the GrapherDirections.txt file.
The script will prompt you to enter an output identifier, this can be anything, it will be the name of your output file. You will also be asked for the trait you wish to graph, this needs to be input exactly as it is written in your TASSEL file.
The script will then run and output a Manhattan plot with your chosen output identifier to your working directory.
- Should you wish to visualise a specific chromosome, select it in the prompt.
- Once plotted you can select specific points with your mouse to view the gene model name.
- Press Esc when you are finished and you will be asked if you want to save this chromosome specific plot.
- Repeat this for all chromosomes you wish to individually graph.



### GEM
To run GEMRegress.R to run GEM analysis on your data, open it in RStudio.
Set your working directory to that which contains the script. BolRPKMTASSEL.txt and Marker_to_At_Bol.csv should also be in here.
Press source to run the script.
You will be prompted to select your trait file.
The script will prompt you to enter an output identifier, this can be anything, it will be the name of your output file.
The script will then run and output your results.

To plot this data you will need to run the GEMGrapher.R script, open this in RStudio.
Set your working directory to that which contains the GEMGrapher.R.
Press source to run the script.
Select the GEMRegress results file you wish to graph.
The script will prompt you to enter an output identifier, this can be anything, it will be the name of your output file.
The script will then run and output a Manhattan plot, with your chosen output identifier, to your working directory.
- Should you wish to visualise a specific chromosome, select it in the prompt.
- Once plotted you can select specific points with your mouse to view the gene model name.
- Press Esc when you are finished and you will be asked if you want to save this chromosome specific plot.
- Repeat this for all chromosomes you wish to individually graph.
