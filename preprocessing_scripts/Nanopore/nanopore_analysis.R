library(tidyverse)
library(Rsamtools)

find_coverage_density <- function(vec, q_start, q_end, 
                                  buffer = 3, l = 100000){
    # The purpose of this function is to detect 
    
    q_start <- q_start + buffer
    q_end <- q_end - buffer
    if(q_end < q_start){
        q_start <- round((q_start+q_end)/2)
        q_end <- q_start
    }
    
    counts_this_region = sum(vec[q_start:q_end])
    density <- counts_this_region/(q_end-q_start+1)
    
    return(density)
}

make_coverage_vector <- function(cigar, start, vector_l=1000000,
                                 find_pos = F){
    # The purpose of this function is to return a vector of 0s and 1s where
    # if the read is aligned at the position it's a 1, otherwise it's a 0
    
    # First split the cigar string into its components
    # S = soft clip, M = match (though could contain mismatches), 
    # I = insert in record, # D = deletion in record, N = skipped region 
    # from the reference, H = hard clip
    cig_split <- unlist(strsplit(cigar, "(?<=S|M|I|D|N|H)", perl=TRUE))
    
    # generate a vector of 0s
    v <- rep(0, vector_l)
    
    # initialise the first position
    pos <- start
    
    for(cig in cig_split){
        # Determine which type of alignment this section is
        type <- str_sub(cig, str_length(cig), str_length(cig))
        
        # Find the length of this section
        cig_length <- as.numeric(str_sub(cig, 1, str_length(cig)-1))
        
        
        # If it's aligned (note that this may include mismatches) then it's 1
        if(type == "M"){
            v[pos:(pos+cig_length-1)] <- 1
        } 
        
        # Move the position along the reference by the cig length, unless 
        # this section is an insert or clipping
        if(!type %in% c("I", "S", "H")){
            pos <- pos + cig_length
        }
    }
    
    if(find_pos){
        return(pos)
    } else{
        return(v)
    }
    
}

### State the positions of interest for the reference

# The intron retention event that is TDP-43 regulated
ir_start <- 61106 
ir_end <- 62405

# the cryptic exon (just the short one)  CHECK THIS
cryptic_start <- 47804
cryptic_end <- cryptic_start+127

# The exon upstream of the cryptic
upstream_exon_start <- 47459  # not actually start - the amplicon start
upstream_exon_end <- 47500

# The exon downstream of the IR
downstream_exon_start <- 62405
downstream_exon_end <- downstream_exon_start+20 # not actually end - the amplicon end


setwd("C:/Users/ogw/Google Drive/UCL PhD/Year 2/unc13a/nanopore/")
samples <- c("05", "06", "07", "08", "09", "10", "11", "12")

frac <- 1  # optional subsampling to improve speed - set to 1 to ignore

for(sample in samples){
    print(sample)
    
    # Read in the bam file and filter for high quality & long reads
    bam_df <- data.frame(scanBam(paste0("bams/barcode", sample, ".bam"))) %>%
        filter(strand %in% c("+", "-")) %>%
        filter(!is.na(str_length(cigar))) %>%
        filter(mapq > 50) %>%
        filter(qwidth > 500) %>%
        group_by(qname) %>%
        arrange(desc(mapq)) %>%
        distinct(qname, .keep_all = T) %>%
        mutate(barcode = paste0("bc",sample)) %>%
        ungroup() %>%
        sample_frac(frac)
    
    
    # Initialise columns for alignment density
    bam_df$ir_d <- 0 # The whole retained intron
    bam_df$cryptic_d <- 0 # The whole (short) cryptic
    bam_df$upstream_d<- 0 # (Part of the) exon upstream of the cryptic
    bam_df$downstream_d<- 0 # (Part of the) exon downstream of the IR
    bam_df$start_ir_d<- 0 # The start of the retained intron (5' trancript-wise)
    bam_df$end_cryptic_d<- 0 # The end of the cryptic (3' transcript-wise)
    
    
    # Iterate through each record in the filtered bam and analyse
    for(i in 1:nrow(bam_df)){
        if(i%%100 == 0){
            print(paste(i, "of", nrow(bam_df)))
        }
        
        # Generate a coverage vector
        this_v <- make_coverage_vector(bam_df$cigar[i],
                          bam_df$pos[i])
        
        # Find the coverage density for each region (specified above)
        bam_df$ir_d[i] <- find_coverage_density( 
                              this_v,
                              ir_start,
                              ir_end)
        
        bam_df$cryptic_d[i] <- find_coverage_density(
                                    this_v,
                                    cryptic_start,
                                    cryptic_end)
        
        bam_df$upstream_d[i] <- find_coverage_density(
                                   this_v,
                                   upstream_exon_start,
                                   upstream_exon_end)
        
        bam_df$downstream_d[i] <- find_coverage_density( 
                                  this_v,
                                  downstream_exon_start,
                                  downstream_exon_end)
        
        bam_df$start_ir_d[i] <- find_coverage_density( 
                                 this_v,
                                  ir_start,
                                  ir_start+50) # first 50 nucleotides of IR
        
        bam_df$end_cryptic_d[i] <- find_coverage_density(
                                       this_v,
                                       cryptic_end-20,
                                       cryptic_end) # last 20 nucleotides of IR
    }
    
    full_reads <- bam_df %>%
        filter(upstream_d > 0.7 | end_cryptic_d > 0.9) %>% # some -ve strand reads may terminate within cryptic
        filter(downstream_d > 0.5 | start_ir_d > 0.3) %>% # some +ve strand reads might terminate in intron - should still include these
        mutate(classification = ifelse(end_cryptic_d < 0.5 & start_ir_d < 0.1, "Neither",
                                       ifelse(end_cryptic_d > 0.8 & start_ir_d < 0.1, "Cryptic",
                                              ifelse(end_cryptic_d < 0.1 & start_ir_d > 0.3, "IR",
                                                     ifelse(end_cryptic_d > 0.8 & start_ir_d > 0.3, "both", NA)))))
    
    write_csv(full_reads, paste0("processed_again/barcode", sample, ".csv"))
    
}


# read in processed files

csvs <- Sys.glob("processed_again/*.csv")

for(csv in csvs){
    if(csv == csvs[1]){
        combined <- read_csv(csv)
    } else {
        combined <- bind_rows(combined, read_csv(csv))
    }
}

sample_names <- data.frame(barcode = c("bc05", "bc06", "bc07", "bc08", "bc09", "bc10", 
                                       "bc11", "bc12"),
                           name = c("SHSY_NT1", "SHSY_Dox1", 
                                    "SHSY_NT2", "SHSY_Dox2",
                                    "P56", "P45", "P63", "P28"))


df <- combined %>%
    filter(!is.na(classification)) %>%
    dplyr::select(barcode, classification) %>%
    group_by(barcode, classification) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(fraction = n/n()) %>%
    ungroup() %>%
    complete(barcode, nesting(classification), fill = list(fraction = 0)) %>%
    left_join(sample_names) %>%
    unique()

p1 <- ggplot(df, aes(x = name, y = 100*fraction, fill = classification)) +
    geom_bar(stat="identity", position="dodge") +
    ggeasy::easy_rotate_x_labels() +
    ylab("Percentage of reads") +
    xlab("") +
    ggtitle("Read distribution")
    

p2 <- ggplot(df %>% filter(classification != "Neither"), aes(x = name, y = 100*fraction, fill = classification)) +
    geom_bar(stat="identity", position="dodge") +
    ggeasy::easy_rotate_x_labels() +
    ylab("Percentage of reads") +
    ggtitle("'Neither' removed") +
    xlab("")

library(patchwork)

p1+p2

ggsave("combined_plot_nanopore.png")

