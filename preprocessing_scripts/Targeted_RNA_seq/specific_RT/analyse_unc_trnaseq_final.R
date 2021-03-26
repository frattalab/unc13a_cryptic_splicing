library(tidyverse)
library(Rsamtools)
library(DescTools)
library(clustringr)
library(stringdist)

# The purpose of this script is to find the number of unique UMIs  (i.e. unique cDNAs) in the 
# targeted RNA sequencing of UNC13A, using a specific RT primer.
#
# This script takes .bam files aligned to the human (HG38) genome. It filters for alignments near
# the UNC13A cryptic exon. It works out how many unique UMIs there are, and whether these 
# cDNAs contain the risk or non-risk alleles of rs12973192
#
# The script first identifies whether reads contain the cryptic exon, or are "correctly" spliced. 
# It then extracts the UMIs associated with each read from the read name (these were already moved
# from the read header during demultiplexing with Ultraplex)
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
# This script differs from the random hexamer script in that it detects whether reads are cryptic
# or correctly spliced by looking at the CIGAR string.

get_sample <- function(string){
    # Detects sample name from the file name
    return(str_split(string,":")[[1]][2])
}

find_unique_umis <- function(df, distance_multiplier = 20, break_ratio = 100,
                             min_n = 50){
    # df should have two columns - "strings" which are the UMIs, and "n" which
    # is the copy number of each
    #
    # higher values of distance multiplier increase ability to detect low abundance
    # UMIs, at the risk of false positives. 100 works well for 1% error rate,
    # 250 is good for 0.2%

    df_remaining <- df
    df_added <- data.frame(strings = c(), n = c())
    biggest <- max(df$n)

    while(T){
        
        # Find the one with the largest number
        largest <- df_remaining$strings[which.max(df_remaining$n)] # the umi sequence
        this_biggest <- max(df_remaining$n)


        if(this_biggest < biggest/break_ratio | this_biggest < min_n){
            break
        }

        this_n <- df_remaining$n[which.max(df_remaining$n)]

        # find its distance from others
        if(nrow(df_added) > 0){
            
            df_added <- df_added %>%
                mutate(distances = stringdist::stringdist(strings, largest))
            
            min_distance <- min(df_added$distances[which(df_added$distances > 0)])

            # find how many there are of these
            already_added_sum <- sum(df_added$n[which(df_added$distances==min_distance)])

        } else {  # if this is the first one
            min_distance <- 8 # set to large so definitely included
            already_added_sum <- 0
        }

        ratio <- already_added_sum/this_n

        if(ratio < distance_multiplier^min_distance & !(str_detect(largest, "N"))){
            # then add
            df_added <- bind_rows(df_added, df_remaining %>%
                                      filter(strings == largest))
        }

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

# a csv containing the names of the samples
samples <- read_csv("//data.thecrick.org/lab-ulej/home/users/wilkino/Unc13a/Targeted_RNAseq/barcodes2.csv",
                   col_names = "sample_with_barcode")

# the directory containing the bam files
dir <- "C:/Users/ogw/Downloads/bams/ultraplex_demux_"

samples$sample <- unlist(lapply(samples$sample_with_barcode, get_sample))

for(sample in samples$sample){

    message(sample)

    bam_name <- paste0(dir,
                        sample, ".fastq.gzAligned.out.sam.bam")

    bam_df <- data.frame(scanBam(bam_name))

    df2 <- bam_df %>%
        mutate(umi = str_sub(qname, str_length(qname)-7, str_length(qname))) %>%
        filter(rname == "chr19", pos == 17641520) %>%  # position that primer "CATTGTTCTGCACGTCGGTC" starts from in HG38 genome
        mutate(read_type = ifelse(str_detect(cigar, "37M857N48M"), "cryptic", "not_cryptic"))  # cryptic leaves gap of 857 nucleotides
    
    # filter just for cryptic reads
    umi_df <- df2 %>%
        filter(read_type == "cryptic") %>%
        select(umi) %>%
        group_by(umi) %>%
        mutate(umi_n = n()) %>% 
        ungroup() %>%
        select(umi, umi_n) %>% 
        unique() %>%
        filter(!str_detect(umi, "N")) %>%
        arrange(desc(umi_n)) %>%
        select(strings = umi, n = umi_n)

    # Find UMIs which appear genuine
    selected_umis <- find_unique_umis(umi_df, distance_multiplier = 10)
    
    if(nrow(selected_umis) > 0){

        # Now filter the bam_df for these umis
        df3 <- df2 %>%
            filter(umi %in% selected_umis$strings)

        # worked out which SNP variant we have
        
        # in the forward strand, healthy = C
        df4 <- df3 %>%
            mutate(nt = str_sub(seq, 54, 54))  # rs12973192 will be 54th base in read
    
        min_acceptable_fract <- 0.9 # remove ambiguous bases
    
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
                this_nt <- this_df$nt[1]  # most abundant because this_dt is ranked by abundance
            } else {
                this_nt <- "dunno"
            }
    
            nt_list <- c(nt_list, this_nt)
            fracs <- c(fracs, max(this_df$nt_frac))
        }
    
        final_df <- data.frame(umis = selected_umis$strings, nt = nt_list) %>%
            group_by(nt) %>%
            mutate(n = n()) %>%
            mutate(risk = ifelse(nt == "G", "Risk", "Non-Risk")) %>%  # in forward genomic strand, G is risk
            ungroup() %>%
            select(nt, risk, n) %>%
            unique() %>%
            mutate(sample_name = sample)
    
        if(sample == samples$sample[1]){
            all_df <- final_df
        } else {
            all_df <- bind_rows(all_df, final_df)
        }
    } else { #if no umis
        # do nothing
    }

}

# Make a dataframe to add values to
to_join <- samples %>%
    select(sample_name = sample)
to_join2 <- bind_rows(to_join %>% mutate(nt = "G", risk = "Risk"),
                      to_join %>% mutate(nt = "C", risk = "Non-Risk"))


just_C_or_G_tmp <- all_df %>%
    full_join(to_join2, by = c("nt", "sample_name", "risk")) %>%
    filter(nt == "G"|nt=="C")  %>%
    filter(!(sample_name == "P17_07"|sample_name =="P35_07"|sample_name=="P64_11"|sample_name == "P47_11")) # remove controls

just_C_or_G_tmp[is.na(just_C_or_G_tmp)] <- 0

just_C_or_G <- just_C_or_G_tmp %>%
    group_by(sample_name) %>%
    mutate(total = sum(n)) %>%
    mutate(n_healthy = ifelse(nt=="C", n, 0)) %>%
    mutate(p_value = pbinom(max(n_healthy), total, 0.5))

# Write final data out

write_csv(just_C_or_G, "for_plot.csv")