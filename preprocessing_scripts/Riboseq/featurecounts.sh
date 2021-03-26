for i in *30.bam
do
featureCounts -t CDS -a /camp/lab/ulej/home/users/farawar/genomes/hs/annotation/gencode.v29.annotation.gtf -o $i\featureCounts.txt $i
done
