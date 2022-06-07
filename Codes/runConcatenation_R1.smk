# Read 1 concatenation
# Pipeline for run concatenation - Concatenate 5 runs of fastq files together

# Command for dry run: snakemake --snakefile runConcatenation_R1.smk -n -p -j 16

# Command for cluster execution: snakemake --snakefile runConcatenation_R1.smk -j 16 --latency-wait 60

import os

# Define the name for output directory
concatDir = "runConcatenation"

# Create a directory for concatenated fastq files
os.system("mkdir -p {concatDir}".format(concatDir=concatDir))

# Extract sample ID from pair-end fastq files
run1Dir = "/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/RNAseq/laneConcatenation/CIA_Paw_Run1_Concat_20210521"
run2Dir = "/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/RNAseq/laneConcatenation/CIA_Paw_Run2_Concat_20210527"
run3Dir = "/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/RNAseq/laneConcatenation/CIA_Paw_Run3_Concat_20210529"
run4Dir = "/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/RNAseq/laneConcatenation/CIA_Paw_Run4_Concat_20210601"
run5Dir = "/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/RNAseq/laneConcatenation/CIA_Paw_Run5_Concat_20210605"

extractList = []

for fileName in os.listdir(run5Dir):
	extractList.append('_'.join(fileName.split("_")[0:2]))

# Get unique sample ID
ID = list(set(extractList))
print(ID)

rule all:
	input:
	  expand(concatDir + "/{sample}_R1_001.fastq.gz", sample = ID)
   
rule runCon_R1:
	input:
	  R1_run1 = run1Dir + "/{sample}_R1_001.fastq.gz",
	  R1_run2 = run2Dir + "/{sample}_R1_001.fastq.gz",
	  R1_run3 = run3Dir + "/{sample}_R1_001.fastq.gz",
	  R1_run4 = run4Dir + "/{sample}_R1_001.fastq.gz",
	  R1_run5 = run5Dir + "/{sample}_R1_001.fastq.gz"
	output:
	  concatDir + "/{sample}_R1_001.fastq.gz"
	message:"R1: Concatenate lanes on {wildcards.sample}"
	shell:
		"cat {input.R1_run1} {input.R1_run2} {input.R1_run3} {input.R1_run4} {input.R1_run5} > {output}"
