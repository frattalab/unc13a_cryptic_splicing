#!/bin/bash


FOLDER=$1
FILES=$FOLDER*.bam


for filename in $FILES; do
 # get the basename of the file
  b=$(basename $filename)
  # cut off eeverythign after a period to get the sample
  sample=$(echo $b | cut -d"." -f1)
  # run bcftools mpileup and do some awk magic to skip the header lines
  # if you want to run you'll need to change the location of your fasta obvs
  bcfoutput="$( bcftools mpileup --no-BAQ -Q 0 -r chr19:17642430-17642430 -f /SAN/vyplab/vyplab_reference_genomes/sequence/human/gencode/GRCh38.primary_assembly.genome.fa $filename | awk 'BEGIN{FS=OFS="\t"} f{print $8} /^#CHROM/{f=1}' | awk -F';' '{print $1 " "$2}')"

  depth=$(echo $bcfoutput | cut -d" " -f1 | cut -d"=" -f2)

  reads=$(echo $bcfoutput | cut -d" " -f2 | cut -d"=" -f2)
  IFS=, read var1 var2 var3 var4 var5 <<< $reads

  ((creads=var1 + var2))
  ((greads=var3 + var4))

  echo "$sample $creads $greads $depth"

done > $FOLDER"unc13.snp.out.tsv"
awk 'BEGIN { print "sample  ref_reads alt_reads total"} { print }' $FOLDER"unc13.snp.out.tsv"   > /tmp/out && mv /tmp/out $FOLDER"unc13.snp.out.tsv"
