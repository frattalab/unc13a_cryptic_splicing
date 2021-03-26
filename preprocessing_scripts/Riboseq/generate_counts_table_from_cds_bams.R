library(tidyverse)

filenames <- Sys.glob("A:/home/users/wilkino/ALL_RIBOSEQ/sarah_hill/full_riboseq_snakemake/snakemake_output/results/genome_dedup/*ounts.txt.gz")

for(filename in filenames){
  sample <- str_split(str_split(filename, "/", simplify = T)[11], "\\.", simplify = T)[1]
  
  all_content = readLines(filename)
  skip = all_content[-1]
  this_counts <- read_tsv(skip)
  
  this_counts <- this_counts[,c(1,7)]  # get only gene_id and the number of counts
  colnames(this_counts) <- c("gene_id", sample) 
  
  if(filename == filenames[1]){
    all_counts <- this_counts
  } else {
    all_counts <- full_join(all_counts, this_counts, by = "gene_id")
  }
}

write_csv(all_counts, "C:/Users/ogw/Google Drive/UCL PhD/Year 2/sarah hill/cds_feature_counts.csv")
