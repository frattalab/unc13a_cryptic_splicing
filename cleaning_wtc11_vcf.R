library(data.table)
library(VariantAnnotation)
library(tidyverse)
#get the genomic coords of the revelvant ALS/FTD genes 
als_ftd_risk_genes = c("AARS2", "ALS2", "ALS3", "ALS7", "ANG", "APBB2", "APOE", "APP", 
                       "ATP13A2", "ATXN2", "C19orf12", "C6orf10", "C9orf72", "CAPZA1", 
                       "CCNF", "CDH22", "CFAP410", "CHCHD10", "CHMP2B", "CHRNA3,4,B4", 
                       "CHRNA4", "CNTN6", "CRYM", "CSF1R", "CTSF", "CX3CR1", "DAO", 
                       "DCTN1", "ELP3", "EPB41L1", "ERBB4", "EWSR1", "FIG4", "FLNC", 
                       "FUS", "GFRA2", "GLE1", "GRN", "GSN", "HNRNPA1", "HNRNPA2B1", 
                       "ITM2B", "KIF5A", "LMNB1", "LOC101929163", "LRRK2", "MAPT", "MATR3", 
                       "MPO", "NEFH", "NEK1", "NIPA1", "NOTCH3", "NPC2", "OPTN", "PFN1", 
                       "PLD3", "PON1-3", "PRKAR1B", "PRNP", "PRPH", "PSEN1", "PSEN2", 
                       "RNPLA6", "SETX", "SIGMAR1", "SOD1", "SORT1", "SPG11", "SQSTM1", 
                       "SRR", "SS18L1", "TAF15", "TARDBP", "TBK1", "TIA1", "TMEM106B", 
                       "TNIP1", "TREM2", "TRPM2", "TSC1", "TUBA4A", "UBQLN2", "UNC13A", 
                       "VAPB", "VCP")

gene_locations = annotables::grch38 %>% filter(symbol %in% als_ftd_risk_genes) %>% dplyr::select(symbol,chr,start,end) %>% filter(
    chr %in% as.character(1:22)
) %>% mutate(chr = paste0("chr",chr))

wtc11_vcf = "/Users/annaleigh/Documents/GitHub/unc13a_cryptic_splicing/data/wtc11_vcf.gz"

gene_gr = makeGRangesFromDataFrame(gene_locations,keep.extra.columns = T)

tab <- TabixFile(wtc11_vcf)
vcf_rng <- VariantAnnotation::readVcf(tab, "hg38", param=gene_gr)

library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb = keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene)
coding <- predictCoding(vcf_rng, txdb, seqSource=Hsapiens)

coding_mutations = annotatr::annotate_regions(coding,gene_gr) %>% as.data.table()

coding_mutations = coding_mutations %>% 
    filter(CONSEQUENCE != "synonymous") %>% 
    dplyr::select(seqnames,start,end,annot.symbol,REFAA,PROTEINLOC,VARAA,strand) %>% 
    makeGRangesFromDataFrame(,keep.extra.columns = T,ignore.strand = F)

names(coding_mutations) = paste(seqnames(coding_mutations),start(coding_mutations),sep = ":")
library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

seqlevelsStyle(coding_mutations) = "UCSC"  # necessary
cur19 = liftOver(coding_mutations, ch)
class(cur19)

cur19 = unlist(cur19)
genome(cur19) = "hg19"

cur19$orig = names(cur19)
cur19 = cur19 %>% as.data.table()

als_od = fread("/Users/annaleigh/Documents/GitHub/unc13a_cryptic_splicing/data/als_od_mutations.csv")

als_od %>% 
    mutate(seqnames = paste0("chr",Chromosome)) %>% 
    mutate(pos2 = as.numeric(`Position hg19`)) %>% 
    left_join(cur19, by = c("seqnames" = "seqnames","pos2" = "start")) %>% 
    filter(!is.na(orig))
