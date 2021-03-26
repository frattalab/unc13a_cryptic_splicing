for i in bowtie_aligned_2/*.sam
do
samtools sort $i > $i\.bam
samtools index $i\.bam
done