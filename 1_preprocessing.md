# Preprocessing community data

## Objectives:  
1. Create phyloseq objects using raw OTU table and taxonomy table.   
2. Eliminate ambiguous rare OTUs that summed up to less than 5 across all samples.  
3. Generate a histogram to evaluate sample sequencing depths.   
4. Calculate the Good's Estimate of Coverage for each sample. 
5. Generate a histogram to evaluate sample coverage from sequencing.   
6. Save any output that could be used downstream directly.  

### 1. This step does all of the above and generates 3 different directories. In terminal, type:   
```
cd ~/pit_foaming_manuscript  
Rscript codes/R_scripts/community_sample_processing_phyloseq.R raw_data/otu_table.txt raw_data/taxa_table.txt
```


