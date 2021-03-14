library(tidyverse)
library(Rsamtools)
library(DescTools)
library(clustringr)
library(stringdist)

# The purpose of this script is to find the number of unique UMIs  (i.e. unique cDNAs) in the 
# targeted RNA sequencing of UNC13A, using a random hexamer RT primer.
#
# This script takes .bam files aligned to the human (HG38) genome. It filters for alignments near
# the UNC13A cryptic exon. It works out how many unique UMIs there are, and whether these 
# cDNAs contain the risk or non-risk alleles of rs12973192
#
# The script extracts the UMIs associated with each read from the read name (these were already moved
# from the read header during demultiplexing with Ultraplex), and checks each read aligns to the cryptic.
#
# Because these libraries were very low complexity/over-amplified (especially for the reads which contained the 
# cryptic exon) there is a high risk of erroneous false-positive UMIs due to sequencing/PCR errors.
# To minimise this, the code first ranks UMIs by their abundance (because the most abundant UMIs
# are unlikely to be derived from sequencing errors). It then works down the list of UMIs, 
# ensuring that each one is sufficiently different from the UMIs which have been added to the 
# list of "true" UMIs. It does this by finding the minimum string distance of each UMI to 
# the UMIs currently in the list of "true" UMIs. It requires that for a UMI sequence which 
# has a string distance of n bases, it is at least 1/(10^n) times the abundance of 
# the most similar UMI. It also requires that each UMI is at least 0.01x times as abundant
# as the most abundant, and is present in at least 50 reads.
#
# This script differs from the specific RT script in that it does not detect correctly spliced transcripts.
# It also doesn't have a lower limit on the number of reads containing each UMI (though 
# each must still be at least 0.01x the abundance of the most abundant).
#
# Note that although the forward reads in the fastq files for these samples correspond to the forward
# strand of the RNA, and thus the MINUS strand of the genome, the bam file format automatically reverse
# complements the reads, thus "C" is still healthy, and "G" is still risk.


get_sample <- function(string){
    return(str_split(string,":")[[1]][2])
}

find_unique_umis <- function(df, distance_multiplier = 20, break_ratio = 100){
    # df should have two columns - "strings" which are the UMIs, and "n" which
    # is the copy number of each
    #
    # higher values of distance multiplier increase ability to detect low abundance
    # UMIs, at the risk of false positives. 100 works well for 1% error rate,
    # 250 is good for 0.2%
    #
    # for longer UMIs, there's increased risk of false positives for some reason
    # so use a lower value
    #
    # break ratio - stop searching once UMIs are very low abundance

    df_remaining <- df
    df_added <- data.frame(strings = c(), n = c())

    biggest <- max(df$n)

    while(T){
        # Find the one with the largest number
        largest <- df_remaining$strings[which.max(df_remaining$n)]
        this_biggest <- max(df_remaining$n)


        if(this_biggest < biggest/break_ratio){
            break
        }

        this_n <- df_remaining$n[which.max(df_remaining$n)]

        # find its distance from others
        if(nrow(df_added) > 0){

            df_added <- df_added %>%
                mutate(distances = stringdist::stringdist(strings, largest))
            min_distance <- min(df_added$distances[which(df_added$distances > 0)])

            # find how many there are of these
            already_added_sum <- sum(df_added$n[which(df_added$distances==min_distance)]) # PERHAPS THIS SHOULD BE ALL OF THEM? AND MAYBE SHOULD BE MAX

        } else {
            min_distance <- 8 # set to large so defo included
            already_added_sum <- 0
        }

        ratio <- already_added_sum/this_n

        if(ratio < distance_multiplier^min_distance & !(str_detect(largest, "N"))){
            # then add
            df_added <- bind_rows(df_added, df_remaining %>%
                                      filter(strings == largest))
        }

        # always remove regardless
        df_remaining <- df_remaining[which(!df_remaining$strings == largest),]

        if(nrow(df_remaining)%%100==0){
            print(nrow(df_remaining))
        }

        if(nrow(df_remaining) == 0){
            break
        }
    }

    if(nrow(df_added) > 0){
        return(df_added %>% select(-distances))
    } else {
        return(df_added)
    }
}

##### RUN

setwd("C:/Users/ogw/Google Drive/UCL PhD/Year 2/unc13a/analyse_targeted_RNASEQ/tg3/")

sample_names = read_csv("sample_names.csv")  # which barcode was used for each sample

samples = data.frame(n = 1:14) %>% mutate(sample_name = paste0("A:/home/users/wilkino/Unc13a/tg3/ultraplex_demux_tg3_s",
                                                               n, ".fastq.gzAligned.sortedByCoord.out.bam"))

for(sample in samples$n){
    
    bam_name = paste0("A:/home/users/wilkino/Unc13a/tg3/ultraplex_demux_tg3_s",
                      sample, ".fastq.gzAligned.sortedByCoord.out.bam")

    message(sample)

    bam_df <- data.frame(scanBam(bam_name))

    df2 <- bam_df %>%
        mutate(umi = str_sub(qname, str_length(qname)-7, str_length(qname))) %>%
        filter(rname == "chr19", pos == 17641536) %>% # corresponds to primer GGTCACGAAGTGGAACAGG
        mutate(read_type = ifelse(str_detect(cigar, "21M857N46M"), "cryptic", "not_cryptic")) # should have 857 nt gap if cryptic
    
    if(nrow(df2) > 20){
    
        umi_df <- df2 %>%
            filter(read_type == "cryptic") %>%
            select(umi) %>%
            group_by(umi) %>%
            mutate(umi_n = n()) %>% ungroup() %>%
            select(umi, umi_n) %>% unique() %>%
            filter(!str_detect(umi, "N")) %>%
            arrange(desc(umi_n)) %>%
            select(strings = umi, n = umi_n)
        
        # find genuine UMIs
        selected_umis <- find_unique_umis(umi_df, distance_multiplier = 10)
    
        # Now filter the bam_df for these umis
        df3 <- df2 %>%
            filter(umi %in% selected_umis$strings)
    
        # worked out which base is in the crucial position
        # in the forward strand, healthy = c
        df4 <- df3 %>%
            mutate(nt = str_sub(seq, 38, 38))  # would expect rs12973192 to be 38th nucleotide of aligned sequence 
    
        min_acceptable_fract <- 0.95
    
        nt_list <- c()
        fracs <- c()
    
        for(this_umi in selected_umis$strings){
            this_df <- filter(df4, umi == this_umi) %>%
                group_by(nt) %>%
                mutate(nt_n = n()) %>%
                ungroup() %>%
                mutate(nt_frac = nt_n/n()) %>%
                select(nt, nt_frac) %>%
                unique() %>%
                arrange(desc(nt_frac))
            
            if(max(this_df$nt_frac) >= min_acceptable_fract){
                this_nt <- this_df$nt[1]
            } else {
                this_nt <- "dunno"
            }
            
            nt_list <- c(nt_list, this_nt)
    
            fracs <- c(fracs, max(this_df$nt_frac))
    
        }
    
        final_df <- data.frame(umis = selected_umis$strings, nt = nt_list) %>%
            group_by(nt) %>%
            mutate(n = n()) %>%
            mutate(risk = ifelse(nt == "G", "Risk", "Non-Risk")) %>%  # in bam file, reads are reverse complemented, so correspond to forward genomic
            ungroup() %>%
            select(nt, risk, n) %>%
            unique() %>%
            mutate(sample_name = sample)
    
        if(sample == samples$n[1]){
            all_df <- final_df
        } else {
            all_df <- bind_rows(all_df, final_df)
        }
    } else {
        print(sample)
    }
}


to_join <- samples %>%
    select(sample_name = n)

to_join2 <- bind_rows(to_join %>% mutate(nt = "G", risk = "Risk"),
                      to_join %>% mutate(nt = "C", risk = "Non-Risk"))

just_C_or_G_tmp <- all_df %>%
    full_join(to_join2, by = c("nt", "sample_name", "risk")) %>%
    filter(nt == "G"|nt=="C")  %>%
    left_join(sample_names, by = "sample_name") %>%
    filter(!(actual_sample_name %in% c("P17_07","P35_07","P64_11","P47_11"))) # remove controls

just_C_or_G_tmp[is.na(just_C_or_G_tmp)] <- 0

just_C_or_G <- just_C_or_G_tmp %>%
    group_by(sample_name) %>%
    mutate(total = sum(n)) %>%
    mutate(n_healthy = ifelse(nt=="C", n, 0)) %>% # C is healthy in forward genomic strand
    mutate(p_value = pbinom(max(n_healthy), total, 0.5))

# Write final data out
write_csv(just_C_or_G, "values_for_plot.csv")
