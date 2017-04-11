# Preprocessing community data
 
## Objectives:  
1. Create phyloseq objects using raw OTU table and taxonomy table.   
2. Eliminate ambiguous rare OTUs that summed up to less than 5 across all samples.  
3. Generate a histogram to evaluate sample sequencing depths.   
4. Calculate the Good's Estimate of Coverage for each sample. 
5. Generate a histogram to evaluate sample coverage from sequencing.   
6. Save any output that could be used downstream directly.  

### 1. This step does all of the above and generates 3 different directories. In terminal, type:   
```bash
cd ~/pit_foaming_manuscript  
Rscript codes/R_scripts/community_sample_processing_phyloseq.R raw_data/otu_table.txt raw_data/taxa_table.txt
```

### 2. Sample trimming in `R` (from outputs generated from step 1):  
```R
library(ggplot2)

# import the sample coverage data:
C <- read.delim("Text_outputs/data.taxmin5.sample_goods_coverage.txt")

# quick plot to visual the sequencing depth and coverage:
ggplot(C, aes(x=sample_sums, y=C)) +
	ylim(0, 1) +
	geom_point() +
	geom_vline(xintercept=10000)
```

We can see that a sequencing depth greater and equal to 10000 good quality reads yielded a minimal coverage of 97% and the community coverages plateued. Therefore, we will remove samples that had less than 10000 sequences.    

```R   
library(phyloseq) 

# import phyloseq object from step 1:  
data.phy <- readRDS("RDS_objects/taxsum_min5_sequence_phyloseq.RDS")
```

    ##> data.phy
    ##phyloseq-class experiment-level object
    ##otu_table()   OTU Table:         [ 8338 taxa and 547 samples ]
    ##tax_table()   Taxonomy Table:    [ 8338 taxa by 6 taxonomic ranks ]

```R 
# removing sample with total sequences less than 10000. Also remove the taxa that are all 0 across samples:   
data.min10k <- prune_samples(sample_sums(data.phy) >= 10000, data.phy)
```  

    ##> data.min10k
    ##phyloseq-class experiment-level object
    ##otu_table()   OTU Table:         [ 8338 taxa and 503 samples ]
    ##tax_table()   Taxonomy Table:    [ 8338 taxa by 6 taxonomic ranks ]

```R  
data.min10k <- prune_taxa(taxa_sums(data.min10k) > 0, data.min10k)
```     
    ##> data.min10k
    ##phyloseq-class experiment-level object
    ##otu_table()   OTU Table:         [ 8336 taxa and 503 samples ]
    ##tax_table()   Taxonomy Table:    [ 8336 taxa by 6 taxonomic ranks ]


### 3. Incorporating sample metadata (samples without any metadata will be removed):    
```R
si <- read.delim("raw_data/sample_meta_data.txt")
```    

    ## > dim(si)
    ##[1] 488 112

```R
# add row names to sample metadata
row.names(si) <- si$SAMPLES
# add metadata to phyloseq object
sample_data(data.min10k) <- si
# because some samples were removed, remove OTUs that summed up 0 across all samples
data.min10k <- prune_taxa(taxa_sums(data.min10k) > 0, data.min10k)
```

    ##> data.min10k
    ##phyloseq-class experiment-level object
    ##otu_table()   OTU Table:         [ 8328 taxa and 488 samples ]
    ##sample_data() Sample Data:       [ 488 samples by 112 sample variables ]
    ##tax_table()   Taxonomy Table:    [ 8328 taxa by 6 taxonomic ranks ]


### note to Fan ###
save.image("sorting_things_out.R")
