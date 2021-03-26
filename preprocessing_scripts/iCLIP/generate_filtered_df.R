library(Rsamtools)
library(tidyverse)

# this script uses the new bowtie2 aligned ones with the larger penalties

setwd("A:/home/users/wilkino/Unc13a/minigene_clip/bowtie_aligned_2/")

files <- Sys.glob("*.bam")

for(file in files){

  sample <- str_sub(file, 1,2)
  this_df <- data.frame(scanBam(file)) %>%
    mutate(sample = sample)

  if(file == files[1]){
    all_df <- this_df
  } else {
    all_df <- bind_rows(all_df, this_df)
  }
}

###### extra filtering #########

filtered_df <- all_df %>%
  filter(str_length(cigar) == 3) %>%  # ensures no gaps (can still be mismatches - cigar doesn't include this info)
  filter(qwidth > 25) %>%  # reasonably stringent to ensure matches are good
  filter(strand == "+") %>%  # unecessary as used "norc" option
  select(pos, sample) %>%  # 1-based coordinate of start
  group_by(sample, pos) %>%
  mutate(xls = n()) %>%
  unique() %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(condition = str_sub(sample, 1, 1)) %>%
  mutate(normalised_counts = 1000*xls/sum(xls))

write_csv(filtered_df, "C:/Users/ogw/Google Drive/UCL PhD/Year 2/unc13a/iCLIP/final_data/bowtie2_aligned_2_filtered_df.csv")

