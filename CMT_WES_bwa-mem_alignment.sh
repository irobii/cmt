#!/bin/bash
#$ -cwd
#$ -S /bin/bash


# Basic argument should be assigned with your directory and files
PLATFORM
TYPE
PROJECT
REF=/Directory/Of/Ensembl/CanFam3.1/Sequence/genome.fa
DataPath
INPUT_FWD_FQ1
INPUT_FWD_FQ2
AlignPath


DIR=/Directory/Of/PROJECT
SAMPLE=${INPUT_FWD_FQ1%_1*}

gatk3=/Directory/Of/GATK_VERSION3/GenomeAnalysisTK.jar

# fixed annotation data path
dbSNP=/Directory/Of/Ensembl/CanFam3.1/Annotation/canis_familiaris.vcf.gz
SV=/Directory/Of/Ensembl/CanFam3.1/Annotation/canis_familiaris_structural_variations.vcf.gz # All structural variations
VAR=/Directory/Of/Ensembl/CanFam3.1/Annotation/canis_familiaris_incl_consequences.vcf.gz # All consequences of the variations on the Ensembl transcriptome, as called by the variation consequence pipeline
BOWTIE2INDEX = /Directory/Of/Ensembl/CanFam3.1/Sequence/Bowtie2Index/genome.fa

# 1. Alignment
# 1-1. Actual Alignment. -I option to use illumina 1.3+ quailities. For the latest version, we don't need -I option. 
bwa mem -t 10 -M $BOWTIE2INDEX $INPUT_FWD_FQ1 $INPUT_FWD_FQ2 > $AlignPath$SAMPLE.sam
gatk SortSam -SO coordinate -I $AlignPath$SAMPLE.sam -O $AlignPath$SAMPLE.bam -VALIDATION_STRINGENCY LENIENT

# 2. Add or replace read groups
gatk AddOrReplaceReadGroups -SO coordinate -I $AlignPath$SAMPLE.bam -O $AlignPath$SAMPLE.RGadded.bam \
-CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT \
-RGLB $PROJECT -RGPL $PLATFORM -RGPU $PLATFORM -RGSM $SAMPLE

# 3. Marking PCR duplicates
gatk MarkDuplicates -I $AlignPath$SAMPLE.RGadded.bam -O $AlignPath$SAMPLE.RGadded.marked.bam \
-METRICS_FILE metrics -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT \

# 4. [Optional] Local realignment around indels
#step1. To create a table of possible indels
java -Xmx16g -Djava.io.tmpdir=$DIR/tmp -jar $gatk3 \
-T RealignerTargetCreator \
-R $REF \
-o $AlignPath$SAMPLE.bam.list \
-I $AlignPath$SAMPLE.RGadded.marked.bam


# 5. [step2] To realign reads around indels targets
java -Xmx16g -Djava.io.tmpdir=$DIR/tmp -jar $gatk3 \
-I $AlignPath$SAMPLE.RGadded.marked.bam \
-R $REF \
-T IndelRealigner \
-targetIntervals $AlignPath$SAMPLE.bam.list \
-o $AlignPath$SAMPLE.RGadded.marked.realigned.bam


# 6. The mate information must be fixed
gatk FixMateInformation \
-INPUT $AlignPath$SAMPLE.RGadded.marked.realigned.bam \
-OUTPUT $AlignPath$SAMPLE.RGadded.marked.realigned.fixed.bam \
-SO coordinate \
-VALIDATION_STRINGENCY LENIENT \
-CREATE_INDEX true

# 7. Quality score recalibration. CountCovariates is not available in the new GATK 2.0 version. Instead, use 'BaseRecalibrator' followed by 'PrintReads'
# 7.1) BaseRecalibrator
gatk BaseRecalibrator \
-R $REF \
-I $AlignPath$SAMPLE.RGadded.marked.realigned.fixed.bam \
--known-sites $dbSNP \
--known-sites $SV \
--known-sites $VAR \
-O $AlignPath$SAMPLE.recal_data.grp

# 7.2) PrintReads
gatk ApplyBQSR \
-R $REF \
-I $AlignPath$SAMPLE.RGadded.marked.realigned.fixed.bam \
-bqsr $AlignPath$SAMPLE.recal_data.grp \
-O $AlignPath$SAMPLE.RGadded.marked.realigned.fixed.recal.bam

# 8. Remove samfile
rm $AlignPath$SAMPLE.sam
