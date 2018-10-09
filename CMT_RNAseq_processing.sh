#RNAseq processing pipeline
#Firstly,perform alignment using Tophat
#Next, assemble transcripts, estimates their abundances using Cufflinks

#two parameters below(GTF,BOWTIE2INDEX) should be assigned with your directory and files
GTF = "/Directory/Of/Ensembl/CanFam3.1/Annotation/genes.gtf"
BOWTIE2INDEX = "/Directory/Of/Ensembl/CanFam3.1/Sequence/Bowtie2Index/genome"


#tophat
#output_directory,gtf,bowtieindex,input_fastq files should be assigned.
#parameters with capital letters should be replaced with your values.

tophat -p 10 -o TOPHAT_OUTPUT_DIR -G GTF --library-type fr-firststrand BOWTIE2INDEX INPUT_FWD_FQ INPUT_REV_FQ


#cufflinks
#output_directory,input_bam files should be assigned.
#parameters with capital letters should be replaced with your values.

cufflinks -p 10 -G GTF --library-type fr-firststrand --output-dir CUFFLINKS_OUTPUT_DIR INPUT_BAM
