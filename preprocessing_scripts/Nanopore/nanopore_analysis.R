library(tidyverse)
library(Rsamtools)
library(patchwork)


find_coverage_density <- function(vec, q_start, q_end, 
                                  buffer = 3){
    # The purpose of this function is to calculate the fraction of the region
    # to which the record is aligned. Eg if this is a 100 nt region and 80 nt
    # of the reference are aligned then it would return 0.8 .
    
    # Buffer is a small amount that is trimmed from the 5' and 3' end to reduce
    # noise at the ends of regions.
    
    # Update start and ends with the buffer
    q_start <- q_start + buffer
    q_end <- q_end - buffer
    
    # if, after applying the buffer, the size is now very small, set the size 
    # to 1
    if(q_end < q_start){
        q_start <- round((q_start+q_end)/2)
        q_end <- q_start
    }
    
    # Find the number of aligned bases in this region
    counts_this_region = sum(vec[q_start:q_end])
    # Divide by the length of the region
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

# the cryptic exon (just the short one)  TODO check this
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
    
    # Classify each read based on whether it contains the cryptic or IR
    # Due to the relatively high error rate of nanopore this is performed 
    # in a permissive way (in both directions)
    full_reads <- bam_df %>%
        filter(upstream_d > 0.7 | end_cryptic_d > 0.9) %>% # some -ve strand reads may terminate within cryptic
        filter(downstream_d > 0.5 | start_ir_d > 0.3) %>% # some +ve strand reads might terminate in intron - should still include these
        mutate(classification = ifelse(end_cryptic_d < 0.5 & start_ir_d < 0.1, "Neither",
                                       ifelse(end_cryptic_d > 0.8 & start_ir_d < 0.1, "Cryptic",
                                              ifelse(end_cryptic_d < 0.1 & start_ir_d > 0.3, "IR",
                                                     ifelse(end_cryptic_d > 0.8 & start_ir_d > 0.3, "both", NA)))))
    
    # write processed reads to csv
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


# make a dataframe linking barcodes to the sample names
sample_names <- data.frame(barcode = c("bc05", "bc06", "bc07", "bc08", "bc09", "bc10", 
                                       "bc11", "bc12"),
                           name = c("SHSY_NT1", "SHSY_Dox1", 
                                    "SHSY_NT2", "SHSY_Dox2",
                                    "FTD1", "FTD2", "FTD3", "FTD4"))

# Find the number of reads of each type for each sample
df <- combined %>%
    filter(!is.na(classification)) %>%  # remove reads that didn't classify clearly
    dplyr::select(barcode, classification) %>%
    group_by(barcode, classification) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    group_by(barcode) %>%
    mutate(fraction = n/n()) %>%
    ungroup() %>%
    complete(barcode, nesting(classification), fill = list(fraction = 0)) %>%
    left_join(sample_names) %>%
    unique() %>%
    mutate(class2 = ifelse(classification == "both", "CE+IR", classification)) %>%
    mutate(type = ifelse(str_detect(name, "SHSY"), "SHSY5Y", "Patient")) %>%
    ungroup()

write_csv(df, "combined_nanopore_dataframe.csv")


p2_patient <- ggplot(df %>% filter(classification != "Neither",
                                   type=="Patient"), aes(x = name, y = 100*fraction, fill = class2)) +
    geom_bar(stat="identity", position="dodge") +
    ggpubr::theme_pubr() +
    ggeasy::easy_rotate_x_labels(side = "right") +
    ylab("Percentage of reads") +
    ggtitle("Patient RNA") +
    xlab("") +
    ggeasy::easy_add_legend_title("")

p2_shsy <- ggplot(df %>% filter(classification != "Neither",
                                type=="SHSY5Y"), aes(x = name, y = 100*fraction, fill = class2)) +
    geom_bar(stat="identity", position="dodge") +
    ggpubr::theme_pubr() +
    ggeasy::easy_rotate_x_labels(side = "right") +
    ylab("Percentage of reads") +
    ggtitle("SH-SY5Y") +
    xlab("") +
    ggeasy::easy_add_legend_title("") 

p2_shsy | p2_patient
ggsave("nanopore_neither_removed.svg", height = 15, width = 20, units="cm")

# without neither removed

p2_patient <- ggplot(df %>% filter(classification != "",
                                   type=="Patient"), aes(x = name, y = 100*fraction, fill = class2,
                                                         colour = class2)) +
    geom_bar(stat="identity", position="dodge") +
    ggpubr::theme_pubr() +
    ggeasy::easy_rotate_x_labels(side = "right") +
    ylab("Percentage of reads") +
    ggtitle("Patient RNA") +
    xlab("") +
    ggeasy::easy_add_legend_title("")
p2_patient

p2_shsy <- ggplot(df %>% filter(classification != "",
                                type=="SHSY5Y"), aes(x = name, y = 100*fraction, fill = class2,
                                                     colour = class2)) +
    geom_bar(stat="identity", position="dodge") +
    ggpubr::theme_pubr() +
    ggeasy::easy_rotate_x_labels(side = "right") +
    ylab("Percentage of reads") +
    ggtitle("SH-SY5Y") +
    xlab("") +
    ggeasy::easy_add_legend_title("") 
p2_shsy | p2_patient
ggsave("nanopore_neither_not_removed.svg", height = 15, width = 20, units="cm")



