#RNAseq data insert size estimation
#Firstly, perform alignment using bwa through coding sequence(cds)
#cds provided from Ensembl CanFam3.1 was used

#CDS parameter below(CDS) should be assigned with your directory and file
CDS = "/Directory/Of/Ensembl/CanFam3.1/CDS/Canis_familiaris.CanFam3.1.cdna.all.fa"

#bwa mem
#cds,output_file,input_fastq files should be assigned.
#parameters with capital letters should be replaced with your values.
#input: cds, forward and reverse fastq files
#output: sam file

bwa mem CDS INPUT_FWD_FQ INPUT_REV_FQ > OUTPUT_SAM

#samtools stats
#output_file,input_sam files should be assigned.
#parameters with capital letters should be replaced with your values.
#input: sam file
#output: information of statistics from sam file

samtools stats INPUT_SAM | grep ^SN > OUTPUT_TXT
