library(phyloseq)
input.env<-read.delim("FM_16S_env_data.txt", sep ="\t", header=T, fill=T, stringsAsFactors=F)
input.tax<-read.delim("taxa_table.txt", sep="\t", header=T)
input.otu<-read.delim("../master_otus/cdhit_otu_table_wide.txt", sep="\t", header=T, row.names=1)
input.tax$repseq<-NULL
rownames(input.tax)<-input.tax$OTUS
head(input.tax)
input.tax$OTUS<-NULL
head(input.tax)
data0.phy<-phyloseq(otu_table(as.matrix(input.otu), taxa_are_rows=T), tax_table(as.matrix(input.tax)))
data0.phy

test<-data.frame(do.call('rbind', strsplit(as.character(input.env$ID_Miseq_ID), '-', fixed=T)))
test$X4<-gsub("(?<![0-9])0+", "", test$X3, perl=T)
test$X5<-paste("S_", test$X4, sep="") 
test$X5[547]
test$X5[547]<-"S_552"
input.env$ID_Miseq_ID[547]
input.env$ID_Miseq_ID[547]<-"FM-16s-552"
colnames(test)[5]<-"SAMPLES"
input.si<-cbind(test$SAMPLES, input.env)
colnames(input.si)[1]
colnames(input.si)[1]<-"SAMPLES"
write.table(input.si, "sample_information.txt", sep="\t",row.names=F, quote=F)

colnames(input.otu)[1:3]
colnames(input.otu)[1:3]<-c("S_97", "S_98", "S_99")
out<-cbind(row.names(input.otu), input.otu)
write.table(out, "otu_table.txt", sep="\t", row.names=F, quote=F)

out<-cbind(row.names(input.tax), input.tax)
write.table(out, "taxa_table.txt", sep="\t", row.names=F, quote=F)

data0.phy<-phyloseq(otu_table(as.matrix(input.otu), taxa_are_rows=T), tax_table(as.matrix(input.tax)))
row.names(input.si)<-input.si$SAMPLES
sample_data(data0.phy)<-input.si
data0.phy

names(input.si)
unique(input.si$Foaming.Status)
input.si[input.si==""]<-NA
unique(input.si$Description)
data1.phy<-subset_samples(data0.phy, !Foaming.Status %in% "DO NOT USE - no feed data")
data1.phy
data1.phy<-prune_samples(sample_sums(data1.phy)>=10000, data1.phy)
data1.phy
summary(sample_sums(data0.phy))
input.si$Long_Sample_ID
sample_data(data1.phy)$Long_Sample_ID
otu<-data.frame(otu_table(data1.phy))
si<-data.frame(sample_data(data1.phy))
totu<-data.frame(t(otu))

library(vegan)
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
names(si)
trt<-data.frame(si$Foaming.Status, si$Description, si[, 46:50])

library(ggplot2)
library(RColorBrewer)
 for (i in names(trt)){
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        colorCount = length(unique(trt[, i]))
        getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
        pdf(paste("16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS(data.mds, trt[, i], colors)
       print(p)
        dev.off()
print(adonis(data.trans~trt[,i]))
}

unique(trt$Description)
trt<-trt[rowSums(is.na(trt))==0, ]
dim(trt)
data2.phy<-data1.phy
data2.phy
sample_data(data2.phy)<-trt
data2.phy
otu<-data.frame(otu_table(data2.phy))
si<-data.frame(sample_data(data2.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
for (i in names(si)){
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        colorCount = length(unique(si[, i]))
        getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
        pdf(paste("16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS(data.mds, si[, i], colors)
        print(p)
        dev.off()
        print(adonis(data.trans~si[,i]))
for (i in names(si)){
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        colorCount = length(unique(si[, i]))
        getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
        pdf(paste("16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS(data.mds, si[, i], colors)
        print(p)
        dev.off()
        print(adonis(data.trans~si[,i]))
} 
getwd()
savehistory("all_16s_preprocess.R")
