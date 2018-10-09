#!/bin/bash
#$ -cwd
#$ -S /bin/bash

# basic argument
PROJECT=$1
REF=$2
AlignPath=$3
TUMOR=$4
# ID=${TUMOR%%.bam}
ID=${TUMOR%%-tumor.RGadded.marked.realigned.fixed.recal.bam}
DIR=/data/project/$PROJECT

if [ ! -d $DIR/05_filtering/OxoG/$ID ]
then
	mkdir -p $DIR/05_filtering/OxoG/$ID
fi
AnalysisPath=$DIR/05_filtering/OxoG/$ID/

gatk CollectSequencingArtifactMetrics -R $REF -I $AlignPath$ID"-tumor/"$TUMOR --FILE_EXTENSION ".txt" -O $AnalysisPath$ID.gatkOxoG

gatk FilterByOrientationBias -V 03.somatic/muTect2/WES3rd/$ID.mutect2.somatic.filter.vcf \
--artifact-modes 'G/T' \
-P $AnalysisPath$ID.gatkOxoG.pre_adapter_detail_metrics.txt -O $AnalysisPath$ID.mutect2.somatic.gatkOxoG_filtered.vcf
