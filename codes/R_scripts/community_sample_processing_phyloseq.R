#################################################
## Initial community process using phyloseq    ##
#################################################
## Author: Fan Yang
## Required: 1) otu table with taxa in rows; 2) taxonomy table with otu label matching those in the otu table
## This script runs in interminal using "Rscript" command and system inputs
## This script outputs phyloseq objects (1 with unmodified information, the other with taxa sums >= 5 across all samples) in "RDS" format that can be directly loaded into R for downstream analysis. It will also out put histograms for sequencing depth and sequencing coverage.
##
## to use: navigate to the directory you wish R objects and histograms will be saved 
## in terminal, type command `Rscrpt /PATH/TO/community_sample_processing_phyloseq.R </PATH/TO/otu_table.txt> </PATH/TO/taxa_table.txt>
## for example: 
## `Rscript ~/Documents/repos/R_code/rscripts/community_sample_processing_phyloseq.R ../R/cdhit_otu_table_wide.txt ../R/cdhit_taxa_table_w_repseq.txt`


library(phyloseq)
library(vegan)

## create three directories in current location to store different types of outputs:
dir.create("RDS_objects")
dir.create("Figures")
dir.create("Text_outputs")

args <- commandArgs(TRUE)
## input arguments:
### 1: otu table with taxa in rows
### 2: taxonomy table, taxa are also in rows

otu <- read.delim(args[1], row.names = 1) #read in otu table
# check otu table dimension
print("the dimension of the OTU table is: ")
dim(otu)

tax <- read.delim(args[2], row.names=1) #read in taxa table
# check taxa table dimension
print("the dimension of the taxa table is: ")
dim(tax)

data.phy<-phyloseq(otu_table(as.matrix(otu),taxa_are_rows=T), tax_table(as.matrix(tax)))
# check phyloseq object size
print("the phyloseq object contains: ")
data.phy
print("saving raw sequence phyloseq object to current directory ... ")
saveRDS(data.phy, "RDS_objects/raw_sequence_phyloseq.RDS")

## get rid of any taxa that summing up to be less than 5 across all samples
data.taxmin5.phy<-prune_taxa(taxa_sums(data.phy) >= 5, data.phy)
data.taxmin5.phy<-prune_samples(sample_sums(data.taxmin5.phy) > 0, data.taxmin5.phy) #get rid of samples with all 0's, if it exists
print("saving phyloseq object with taxa sum >= 5 across all samples to current directory ... ")
saveRDS(data.taxmin5.phy, "RDS_objects/taxsum_min5_sequence_phyloseq.RDS")

## save the histogram of sample sequencing depth distribution
print("Generating histogram on sample sequencing depth to current directory ... ")
pdf("Figures/data.taxmin5.sequencing_depth_hist.pdf")
hist(sample_sums(data.taxmin5.phy), breaks = nsamples(data.taxmin5.phy)/2)
dev.off()

## calculate Good's coverage for each sample  
si <- data.frame(sample_sums(data.taxmin5.phy)) #get sample sums
names(si) <- "sample_sums"
totu<-t(data.frame(otu_table(data.phy))) #transpose subsetted otu table so that samples are in row

si$n1 <- rowSums(totu == 1) #counts of singletons in each sample
si$C <- 1-(si$n1 / si$sample_sums) #calculate Good's coverage
si$SAMPLES <- row.names(si)
## save the sample coverage information to file
print("Saving the sample coverage information as a text file in current directory ... ")
write.table(si, "text_outputs/data.taxmin5.sample_goods_coverage.txt", sep = "\t", quote = F, row.names = F)

## plot Good's coverage histogram
print("Generating histogram on Good's estimated coverage to current directory ... ")
pdf("Figures/data.taxmin5.goods_coverage_hist.pdf")
hist(si$C, breaks = nrow(si)/2)
dev.off()

 
