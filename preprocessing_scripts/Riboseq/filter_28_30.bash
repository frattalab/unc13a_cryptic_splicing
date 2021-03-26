for i in *out.bam 
do 
samtools view -h $i | awk 'length($10) >= 28 || length($10) <= 30 || $1 ~ /^@/' | samtools view -bS - > $i\filtered_28_30.bam
done

