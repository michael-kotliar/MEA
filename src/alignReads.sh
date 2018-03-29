#!/bin/bash

MEA_BIN=`dirname $0`
source $MEA_BIN/mea.config

##############################################################################
#############   Module 3: allele-specific alignment
##############################################################################

PARAM_VALID=1
PARAM_SINGLE_READS=1
if [ "$1" = "-s" ]; then
    if [ $# -eq 7 ]; then
        PARAM_FASTQ_FILE=$2
        PARAM_GENOME=$3
        PARAM_REFERENCE_GENOME=$4
        PARAM_STRAIN1=$5
        PARAM_STRAIN2=$6
        PARAM_BAM_PREFIX=$7
    else
        PARAM_VALID=0
    fi
elif [ "$1" = "-p" ]; then
    if [ $# -eq 8 ]; then
        PARAM_SINGLE_READS=0
        PARAM_FASTQ_FILE1=$2
        PARAM_FASTQ_FILE2=$3
        PARAM_GENOME=$4
        PARAM_REFERENCE_GENOME=$5
        PARAM_STRAIN1=$6
        PARAM_STRAIN2=$7
        PARAM_BAM_PREFIX=$8
    else
        PARAM_VALID=0
    fi
else
    PARAM_VALID=0
fi

if [ $PARAM_VALID = 0 ]; then
    echo "
Usage:
        mea alignReads <-s/-p> input_reads_1 [input_reads_2] genome_concat reference_genome strain1 strain2 outputPrefix

Options:
    -s              to align single-end reads (requires one input file)
    -p              to align paired-end reads (requires two input files)
                    
    input_reads_1
                    the 1st input reads file in fastq.
                    (fastq.gz or bam is supported when using BWA)
                    
    input_reads_2
                    (paired end) the 2nd input reads file in fastq.
                    (fastq.gz or bam is supported when using BWA)
                    
    genome_concat
                    path to the indexed reference for concatenated insilico genome.
                    for BWA, specifiy path to the fasta.
                    for Bowtie2 and Tophat2, specify path and basename of index files
                    for Bismark, specify genome folder, excluding <Bisulphite_Genome>
    reference_genome
                    path to the reference genome indices folder
    strain1         name of strain1
                    (e.g. hap1 or CASTEiJ)
                    
    strain2         name of strain2
                    (e.g. hap2 or C57BL6J)
                    
    outputPrefix    prefix for output files, including the full path, without an extension
                    (e.g. ./TSC_H3K36me3 )
"
exit 1
fi

#------------------------------------------------------------------------------------
# detect the allelic reads from the reads aligned to the concatenated genome
#------------------------------------------------------------------------------------
function detectAllelicConcatenated {
    
    printProgress "[detectAllelicConcatenated] Started"
    
    local PARAM_INPUT_SAM=$1
    local PARAM_STRAIN=$2
    local PARAM_REFFASTA=$3
    local PARAM_QUALITY=$4
    local PARAM_OUT_PREFIX=$5
    
    # output header first
    samtools view -SH "$PARAM_INPUT_SAM" \
        | awk -v ref="$PARAM_STRAIN" '($0 ~ ref) {print $0}' \
        | sed 's/'"$PARAM_STRAIN"'_//g' \
        > "$PARAM_OUT_PREFIX".sam
    
    # append reads
    samtools view -S $PARAM_INPUT_SAM \
        | awk -v ref="$PARAM_STRAIN" '(($3 ~ ref)&&($5'"$PARAM_QUALITY"')) {print $0}' \
        | sed 's/'"$PARAM_STRAIN"'_//g' \
        >> "$PARAM_OUT_PREFIX".sam
    
    # convert to bam
    samtools view -bt $PARAM_REFFASTA "$PARAM_OUT_PREFIX".sam > "$PARAM_OUT_PREFIX".unsorted.bam
    
    # sort by coordinates
    samtools sort "$PARAM_OUT_PREFIX".unsorted.bam "$PARAM_OUT_PREFIX"
    
    if [ -f "$PARAM_OUT_PREFIX".bam ]
    then
        printProgress "[detectAllelicConcatenated] Filtered BAM file created."
        # remove temp files
        printProgress "[detectAllelicConcatenated] Removing temporary SAM files."
        rm -f "$PARAM_OUT_PREFIX".sam
        rm -f "$PARAM_OUT_PREFIX".unsorted.bam
    else
        printProgress "[detectAllelicConcatenated] ERROR: Filtered BAM file was not created."
    fi
    printProgress "[detectAllelicConcatenated] Done"
}


#------------------------------------------------------------------------------------
# alignReads
#------------------------------------------------------------------------------------

if [ $PARAM_SINGLE_READS = 1 ]; then
    #align single-end reads to concatenated insilico genome
    printProgress "align to concatenated insilico genome"
     STAR --runMode alignReads --genomeDir "$PARAM_GENOME" "$MEA_STAR_ALN_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE"
    mv Aligned.out.sam "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
    mv Log.out "$PARAM_BAM_PREFIX"_STAR_RunParameters.tsv
    mv Log.final.out "$PARAM_BAM_PREFIX"_STAR_AlignmentSummary.tsv
    printProgress "align to reference genome"
    STAR --runMode alignReads --genomeDir "$PARAM_REFERENCE_GENOME" "$MEA_STAR_ALN_TOTAL_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE"
    mv Aligned.out.sam "$PARAM_BAM_PREFIX"_total.sam
    samtools view -bShu "$PARAM_BAM_PREFIX"_total.sam | samtools sort - "$PARAM_BAM_PREFIX"_total
    samtools index "$PARAM_BAM_PREFIX"_total.bam
    mv Log.out "$PARAM_BAM_PREFIX"_total_STAR_referenceRunParameters.tsv
    mv Log.final.out "$PARAM_BAM_PREFIX"_total_STAR_referenceAlignmentSummary.tsv
    printProgress "detecting allelic reads"
    detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
    detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
    rm SJ.out.tab Log.progress.out
else #[ $PARAM_SINGLE_READS = 0 ]
    #align paired-end reads to concatenated insilico genome
    printProgress "align to the insilico concatenated genome"
    STAR --runMode alignReads --genomeDir "$PARAM_GENOME" "$MEA_STAR_ALN_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE1" "$PARAM_FASTQ_FILE2"
    mv Aligned.out.sam "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
    mv Log.out "$PARAM_BAM_PREFIX"_STAR_RunParameters.tsv
    mv Log.final.out "$PARAM_BAM_PREFIX"_STAR_AlignmentSummary.tsv
    printProgress "align to the reference genome"
    STAR --runMode alignReads --genomeDir "$PARAM_REFERENCE_GENOME" "MEA_STAR_ALN_TOTAL_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE1" "$PARAM_FASTQ_FILE2"
    mv Aligned.out.sam "$PARAM_BAM_PREFIX"_total.sam
    samtools view -bShu "$PARAM_BAM_PREFIX"_total.sam | samtools sort - "$PARAM_BAM_PREFIX"_total
    samtools index "$PARAM_BAM_PREFIX"_total.bam
    mv Log.out "$PARAM_BAM_PREFIX"_total_STAR_referenceRunParameters.tsv
    mv Log.final.out "$PARAM_BAM_PREFIX"_total_STAR_referenceAlignmentSummary.tsv
    printProgress "detecting allelic reads"
    detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
    detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
    rm SJ.out.tab Log.progress.out
fi
