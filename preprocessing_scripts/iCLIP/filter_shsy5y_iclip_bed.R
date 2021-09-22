library(tidyverse)
library(patchwork)

get_bedgraph <- function(filename){
  bg <- read_tsv(filename, col_names = c("chr", "start", "end", "score"))
  return(bg)
}

bg1 <- get_bedgraph("shsy_1_full.bedgraph")
bg2 <- get_bedgraph("shsy_2_full.bedgraph")

combined <- full_join(bg1 %>% rename(score1 = score), 
                      bg2 %>% rename(score2 = score),
                      by=c("chr", "start"))

combined[is.na(combined)] <- 0

combined2 <- combined %>%
  mutate(total =score1+score2 ) %>%
  filter((chr == "chr19" & start > 17599327 & ((end.x < 17690344 & end.y == 0) | 
                                                 (end.y < 17690344 & end.x == 0))) | 
           chr == "chr9" & start > 35160008 & ((end.x < 35407338 & end.y == 0) | 
                                                 (end.y < 35407338 & end.x == 0)))

write_csv(combined2, "shsy5y_tdp43_iclip_unc13a_b_only.csv")