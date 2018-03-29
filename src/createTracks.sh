#!/bin/bash

MEA_BIN=`dirname $0`
source $MEA_BIN/mea.config

if test $# -ne 8
then
    echo "
Usage:   
         mea createTracks <-s/-p> bamPrefix chrLengthFile strain1 strain2 strain1.refmap strain2.refmap outputDIR
         
Options:
         -s to create tracks for the single-end aligned reads
         -p to create tracks for the paired-end aligned reads

         bamPrefix      prefix used for the output of alignReads command
         chrLengthFile  path to the reference genome chrLengthFile
         strain1        name of strain1 (e.g. hap1)
         strain2        name of strain2 (e.g. hap2)
         genome1.refmap path to the refmap file created for insilico genome 1
         genome1.refmap path to the refmap file created for insilico genome 2
         outputDIR      output directory (where to create track files)
         
Output:
         outputDIR/outputPrefix_strain1.bedGraph
         outputDIR/outputPrefix_strain1.bw        read profiles for strain1 projected to reference genome
         
         outputDIR/outputPrefix_strain2.bedGraph
         outputDIR/outputPrefix_strain2.bw        read profiles for strain2 projected to reference genome
         
         outputDIR/outputPrefix_strain1.wig.gz
         outputDIR/outputPrefix_strain2.wig.gz    unprojected read profiles for strain1 and strain2
"
exit 1
fi

##############################################################################
#############   Module 4: projection to reference genome
##############################################################################

###converts bam to wig using bedtools
function BAM2WIGbedtools {
    local PARAM_INPUT_PREFIX=$1
    local PARAM_OUTPUT_DIR=$2

    local VAR_q=$MEA_BAM2WIG_PARAM_MIN_QUALITY     # min read quality [0]
    local VAR_F=$MEA_BAM2WIG_PARAM_FILTERING_FLAG  # filtering flag [0]
    local VAR_x=$MEA_BAM2WIG_PARAM_SE_EXTENSION    # average fragment length used for fixed length of the read extension [0]. used for ChIP-seq (SET) only
    local VAR_INPUT_BASENAME=`basename $PARAM_INPUT_PREFIX`
    
    printProgress "[bam2wig of $VAR_INPUT_BASENAME] started"
    samtools view -bh -F "$VAR_F" -q "$VAR_q" "$PARAM_INPUT_PREFIX".bam > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bam
    
    if [ "$VAR_OPTION" = "-s" ]; then
        bedtools genomecov -bg -split -ibam "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bam -g "$PARAM_CHROM_SIZES" \
        > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q"_unsorted.bedGraph
    elif [ "$VAR_OPTION" = "-p" ]; then
        bedtools genomecov -bg -split -ibam "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bam -g "$PARAM_CHROM_SIZES" \
        > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q"_unsorted.bedGraph
    else
        echo "Invalid option $VAR_OPTION"
        exit 1
    fi
    sort -k1,1 -k2,2n "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q"_unsorted.bedGraph > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bedGraph
    bedGraphToBigWig "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bedGraph "$PARAM_CHROM_SIZES" "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bw #output this file for viz
    
    awk 'BEGIN {
        print "track type=wiggle_0"
    }
    NF == 4 {
        print "fixedStep chrom="$1" start="$2+1" step=1 span=1"
        for(i = 0; i < $3-$2; i++) {
            print $4
        }
    }' "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bedGraph > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".wig
    bgzip -c "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".wig > "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".wig.gz
    mv "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bedGraph "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".bedGraph
    mv "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bw "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".bw
    rm "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".wig
}

### projects a wig profile to reference genome
function projectToReferenceGenome {
    local PARAM_WIG_FILE=$1
    local PARAM_REFMAP_FILE=$2
    local PARAM_BEDGRAPH_FILE=$3
    
    printProgress "Started projectToReferenceGenome"
    
    $ALEA_JAR project\
        --input-wig=$PARAM_WIG_FILE\
        --input-refmap=$PARAM_REFMAP_FILE\
        --output-bedgraph=$PARAM_BEDGRAPH_FILE
    
    printProgress "Finished projectToReferenceGenome"
}

VAR_OPTION=$1
shift

PARAM_BAM_PREFIX=$1
PARAM_CHROM_SIZES=$2
PARAM_STRAIN1=$3
PARAM_STRAIN2=$4
PARAM_REFMAP_FILE1=$5
PARAM_REFMAP_FILE2=$6
PARAM_OUTPUT_DIR=$7

meaCreateDir "$PARAM_OUTPUT_DIR"

VAR_OUTPUT_BASENAME=`basename $PARAM_BAM_PREFIX`
VAR_OUTPUT_PREFIX1="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"
VAR_OUTPUT_PREFIX2="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"

BAM2WIGbedtools "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1" "$PARAM_OUTPUT_DIR" "$PARAM_CHROM_SIZES"
BAM2WIGbedtools "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2" "$PARAM_OUTPUT_DIR" "$PARAM_CHROM_SIZES"
mv "$VAR_OUTPUT_PREFIX1".wig.gz "$VAR_OUTPUT_PREFIX1"_preProject.wig.gz
mv "$VAR_OUTPUT_PREFIX2".wig.gz "$VAR_OUTPUT_PREFIX2"_preProject.wig.gz
mv "$VAR_OUTPUT_PREFIX1".bedGraph "$VAR_OUTPUT_PREFIX1"_preProject.bedGraph
mv "$VAR_OUTPUT_PREFIX2".bedGraph "$VAR_OUTPUT_PREFIX2"_preProject.bedGraph
mv "$VAR_OUTPUT_PREFIX1".bw "$VAR_OUTPUT_PREFIX1"_preProject.bw
mv "$VAR_OUTPUT_PREFIX2".bw "$VAR_OUTPUT_PREFIX2"_preProject.bw
mv "$VAR_OUTPUT_PREFIX1".wig "$VAR_OUTPUT_PREFIX1"_preProject.wig
mv "$VAR_OUTPUT_PREFIX2".wig "$VAR_OUTPUT_PREFIX2"_preProject.wig

projectToReferenceGenome "$VAR_OUTPUT_PREFIX1"_preProject.wig.gz "$PARAM_REFMAP_FILE1" "$VAR_OUTPUT_PREFIX1".bedGraph
bedGraphToBigWig "$VAR_OUTPUT_PREFIX1".bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX1".bw

projectToReferenceGenome "$VAR_OUTPUT_PREFIX2"_preProject.wig.gz "$PARAM_REFMAP_FILE2" "$VAR_OUTPUT_PREFIX2".bedGraph
bedGraphToBigWig "$VAR_OUTPUT_PREFIX2".bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX2".bw

BAM2WIGbedtools "$PARAM_BAM_PREFIX"_total "$PARAM_OUTPUT_DIR" "$PARAM_CHROM_SIZES"
printProgress "[bam2wig] finished successfully"


### generate normalized bigwigs
### and place them in UCSC track hub

VAR_q=$MEA_BAM2WIG_PARAM_MIN_QUALITY     # min read quality [0]
VAR_F=$MEA_BAM2WIG_PARAM_FILTERING_FLAG  # filtering flag [0]
VAR_OUTPUT_TOTAL="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_total_F"$VAR_F"_q"$VAR_q"
printProgress "[TrackHub Generate] started"

meaCreateDir "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"
touch "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt "$PARAM_OUTPUT_DIR"/Track_Hub/hub.txt "$PARAM_OUTPUT_DIR"/Track_Hub/genomes.txt

RPM_SCALING_FACTOR_tmp=$(samtools view -c "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".bam)
RPM_SCALING_FACTOR=$(echo "scale=25;1000000/$RPM_SCALING_FACTOR_tmp" | bc)
echo "scaling factor: $RPM_SCALING_FACTOR"

#ref
awk -F "\t" -v RPM_SCALE="$RPM_SCALING_FACTOR" 'BEGIN {OFS="\t"; print $0} {
    $4 *= RPM_SCALE
    print $0;}' "$VAR_OUTPUT_PREFIX1".bedGraph > "$VAR_OUTPUT_PREFIX1"_tmp1.bedGraph
grep -v "type" "$VAR_OUTPUT_PREFIX1"_tmp1.bedGraph > "$VAR_OUTPUT_PREFIX1"_tmp2.bedGraph
sort -k1,1 -k2,2n "$VAR_OUTPUT_PREFIX1"_tmp2.bedGraph > "$VAR_OUTPUT_PREFIX1"_tmp3.bedGraph
bedGraphToBigWig "$VAR_OUTPUT_PREFIX1"_tmp3.bedGraph "$PARAM_CHROM_SIZES" "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_RPM.bw
rm "$VAR_OUTPUT_PREFIX1"_tmp1.bedGraph "$VAR_OUTPUT_PREFIX1"_tmp2.bedGraph "$VAR_OUTPUT_PREFIX1"_tmp3.bedGraph
#alt
awk -F "\t" -v RPM_SCALE="$RPM_SCALING_FACTOR" 'BEGIN {OFS="\t"; print $0} {
    $4 *= RPM_SCALE
    print $0;}' "$VAR_OUTPUT_PREFIX2".bedGraph > "$VAR_OUTPUT_PREFIX2"_tmp1.bedGraph
grep -v "type" "$VAR_OUTPUT_PREFIX2"_tmp1.bedGraph > "$VAR_OUTPUT_PREFIX2"_tmp2.bedGraph
sort -k1,1 -k2,2n "$VAR_OUTPUT_PREFIX2"_tmp2.bedGraph > "$VAR_OUTPUT_PREFIX2"_tmp3.bedGraph
bedGraphToBigWig "$VAR_OUTPUT_PREFIX2"_tmp3.bedGraph "$PARAM_CHROM_SIZES" "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_RPM.bw
rm "$VAR_OUTPUT_PREFIX2"_tmp1.bedGraph "$VAR_OUTPUT_PREFIX2"_tmp2.bedGraph "$VAR_OUTPUT_PREFIX2"_tmp3.bedGraph
#total
mv "$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_total.bedGraph "$VAR_OUTPUT_TOTAL".bedGraph
awk -F "\t" -v RPM_SCALE="$RPM_SCALING_FACTOR" 'BEGIN {OFS="\t"; print $0} {
    $4 *= RPM_SCALE
print $0;}' "$VAR_OUTPUT_TOTAL".bedGraph > "$VAR_OUTPUT_TOTAL"_RPM.bedGraph
bedGraphToBigWig "$VAR_OUTPUT_TOTAL"_RPM.bedGraph "$PARAM_CHROM_SIZES" "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_total_F"$VAR_F"_q"$VAR_q"_RPM.bw

### code adapted from Aaron Bogutz, Louis Lefebvre lab (UBC)
printf "genome "$MEA_BUILD"\ntrackDb "$MEA_BUILD"/trackDb.txt" > "$PARAM_OUTPUT_DIR"/Track_Hub/genomes.txt
printf "hub <name>\nshortLabel <short name>\nlongLabel <Hub to display data at UCSC>\ngenomesFile genomes.txt\nemail <email>" > "$PARAM_OUTPUT_DIR"/Track_Hub/hub.txt
printf "track %s\ncontainer multiWig\nshortLabel %s\nlongLabel %s\ntype bigWig\nvisibility full\nmaxHeightPixels 150:60:32\nconfigurable on\nautoScale on\naggregate transparentOverlay\nshowSubtrackColorOnUi on\npriority 1.0\n\n" $VAR_OUTPUT_BASENAME $VAR_OUTPUT_BASENAME $VAR_OUTPUT_BASENAME | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 215,215,215\n\taltColor 225,225,225\n\tautoScale on\n\n" $VAR_OUTPUT_BASENAME"_total" $VAR_OUTPUT_BASENAME $VAR_OUTPUT_BASENAME"_total" $VAR_OUTPUT_BASENAME"_total" $VAR_OUTPUT_BASENAME"_total.bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 215,215,215\n\taltColor 225,225,225\n\tautoScale on\n\n" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1 $VAR_OUTPUT_BASENAME $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1 $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1 $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1".bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 215,215,215\n\taltColor 225,225,225\n\tautoScale on\n\n" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2 $VAR_OUTPUT_BASENAME $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2 $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2 $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2".bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
hubCheck -noTracks "$PARAM_OUTPUT_DIR"/Track_Hub/hub.txt

printProgress "[TrackHub generate] finished successfully"


