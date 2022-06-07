# Pipeline for lane concatenation

# Command for dry run: snakemake --snakefile laneConcatenation.smk -n -p --cluster qsub -j 12 --verbose
# Command for dag: snakemake --snakefile laneConcatenation.smk --dag | dot -Tpdf > dag.pdf
# Command to execute: snakemake --snakefile laneConcatenation.smk -j 48
# Command for cluster execution: snakemake --snakefile laneConcatenation.smk --cluster qsub -j 12 --latency-wait 30 --verbose
# Command for run report: snakemake --snakefile laneConcatenation.smk --report report.html

import os

# Define the name for output directory
concatDir = "CIA_Paw_Run5_Concat_20210605"

# Create a directory for concatenated fastq files
os.system("mkdir -p {concatDir}".format(concatDir=concatDir))

# Extract sample ID from pair-end fastq files
rawDataDir = "/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/RNAseq/illuminaOutput/CIA_Paw_Run5_20210605"

extractList = []

for fileName in os.listdir(rawDataDir):
	extractList.append('_'.join(fileName.split("_")[0:2]))

# Get unique sample ID
ID = list(set(extractList))
print(ID)

rule all:
	input:
	  expand(concatDir + "/{sample}_R1_001.fastq.gz", sample = ID),
	  expand(concatDir + "/{sample}_R2_001.fastq.gz", sample = ID)
   
rule laneCon_R1:
	input:
	  lane1 = rawDataDir + "/{sample}_L001_R1_001.fastq.gz",
	  lane2 = rawDataDir + "/{sample}_L002_R1_001.fastq.gz",
	  lane3 = rawDataDir + "/{sample}_L003_R1_001.fastq.gz",
	  lane4 = rawDataDir + "/{sample}_L004_R1_001.fastq.gz"
	output:
	  concatDir + "/{sample}_R1_001.fastq.gz"
	message:"R1: Concatenate lanes on {wildcards.sample}"
	shell:
		"cat {input.lane1} {input.lane2} {input.lane3} {input.lane4} > {output}"

rule laneCon_R2:
	input:
	  lane1 = rawDataDir + "/{sample}_L001_R2_001.fastq.gz",
	  lane2 = rawDataDir + "/{sample}_L002_R2_001.fastq.gz",
	  lane3 = rawDataDir + "/{sample}_L003_R2_001.fastq.gz",
	  lane4 = rawDataDir + "/{sample}_L004_R2_001.fastq.gz"
	output:
	  concatDir + "/{sample}_R2_001.fastq.gz"
	message:"R2: Concatenate lanes on {wildcards.sample}"
	shell:
		"cat {input.lane1} {input.lane2} {input.lane3} {input.lane4} > {output}"
