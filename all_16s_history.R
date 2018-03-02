NMDS<-data.frame(MDS1,MDS2,Treatment)
NMDS.mean=aggregate(NMDS[,1:2],list(group=Treatment),mean)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  df_ell <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
  }
X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment),size=3,alpha=0.75) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=2, linetype=5)+theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
X1    
}
for (i in names(trt)){
library(ggplot2)
library(RColorBrewer)
for (i in names(trt)){
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
ls()
library(phyloseq)
library(ggplot2)
library(DESeq2)
library(plyr)
library(RColorBrewer)
ls()
data2.phy
test<-prune_taxa(taxa_sums(data2.phy)==0, data2.phy)
test
otu_table(test)[1:5, 1:5]
head(rowSums(otu_table(test)[, 1:5]
)
)
test<-prune_taxa(taxa_sums(data2.phy)!=0, data2.phy)
test
head(taxa_sums(test))
data2.phy<-prune_taxa(taxa_sums(data2.phy)!=0, data2.phy)
data2.phy
dim(trt)
trt[]<-lapply(trt, factor)
str(trt)
unique(trt$Description)
sample_data(data2.phy)<-trt
data2.phy
source("~/code/R/ggplot_nmds_eclipse.R")
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
library(vegan)
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
}
unique(si$FK_Pump.Interval)
names(si)
data1.phy
data0.phy
test<-prune_taxa(taxa_sums(data1.phy)!=0, data1.phy)
test
data1.phy<-prune_taxa(taxa_sums(data1.phy)!=0, data1.phy)
data1.phy
otu<-data.frame(otu_table(data1.phy))
si<-data.frame(sample_data(data1.phy))
totu<-data.frame(t(otu))
getwd()
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
head(trt)
data1.py
data1.phy
namesa(sample_data(data1.phy)
)
names(sample_data(data1.phy)a
names(sample_data(data1.phy))
str(sample_data(data1.phy)
)
names(sample_data(data1.phy))
si_1<-data.frame(data1.phy)
si_1<-data.frame(sample_data(data1.phy))
si_1[1:5, 1:5]
si_1[1:5, 1:10]
si_1$X..of.Max
head(si_1$Name)
names(si_1)
head(si_1$Long_Sample_ID)
test<-data.frame(do.call('rbind',strsplit(as.character(si_1$Long_Sample_ID), "-", fixed=T))
)
head(test)
tail(test)
unique(test$X5)
unique(test$X1)
unique(test$X2)
unique(test$X3)
dim(test)
save.image("all_16s_physeq.RData")
savehistory("all_16s_history.R")
ls()
data0.phy
library(phyloseq)
library(ggplot2)
library(vegan)
library(plyr)
library(RColorBrewer)
data0.phy
data1.phy
data2.phy
head(test)
dim(test)
dim(si_1)
unique(test$X5)
count(test, c(X5))
count(test, c('X5'))
count(test, c('X3'))
unique(test$X2)
count(test, c('X2'))
count(test, c('X1'))
dim(test)
dim(si_1)
head(test)
names(si_1)
head(si_1$FK_Pump.Interval)
ls()
head(si_1$Classification.Number)
head(si_1$Name)
si_1<-cbind(test, si_1)
head(si_1[, 1:12]
)
tail(si_1[, 1:12])
test<-data.frame(si_1$X3)
head(test)
test<-substring(si_1$X3, seq(1, nchar(x),2), seq(2, nchar(x),2))
test<-substring(si_1$X3, seq(1, nchar(si_1$X3),2), seq(2, nchar(si_1$X3),2))
str(si_1$X3)
test<-substring(as.character(si_1$X3), seq(1, nchar(as.character(si_1$X3)),2), seq(2, nchar(as.character(si_1$X3)),2))
lapply()
lapply
x<-"051713"
substr(x, seq(1, nchar(x), 2), seq(2, nchar(x),2))
substr(x, seq(1, nchar(x), 2), seq(2, nchar(x),2), seq(3, nchar(x),2))
substring(x, seq(1, nchar(x), 2), seq(2, nchar(x),2))
x<-si_1$X3
str(x)
x<-as.character(x)
str(x)
x<-si_1$X3
x<-as.character(as.numeric(x))
str(x)
x<-si_1$X3
x[]<-lapply(x, as.character)
head(x)
head(si_1$X3)
str(x)
save.image("all_16s_physeq.RData")
savehistory("all_16s_history.R")
ls()
library(phyloseq)
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(dplyr)
library(plyr)
getwd()
data0.phy
data1.phy
data2.phy
dim(si_1)
names(si_1)
si_1[1:5, 1:5]
head(x)
test<-data.frame(do.call('rbind',substring(x, seq(1, nchar(x), 2), seq(2, nchar(x),2)))
)
lapply
test<-lapply(x, substring(x, seq(1, nchar(x), 2), seq(2, nchar(x),2)))
test<-sapply(x, function(x) substring(x, seq(1, nchar(x), 2), seq(2, nchar(x),2)))
test<-sapply(as.character(x), function(x) substring(x, seq(1, nchar(x), 2), seq(2, nchar(x),2)))
head(test)
test<-lapply(as.character(x), function(x) substring(x, seq(1, nchar(x), 2), seq(2, nchar(x),2)))
head(test)
test<-sapply(as.character(x), function(x) substring(x, seq(1, nchar(x), 2), seq(2, nchar(x),2)))
res<-sapply(as.character(x), function(x) substring(x, seq(1, nchar(x), 2), seq(2, nchar(x),2)))
test<-data.frame(t(res))
test<-data.frame(t(res), x)
test<-data.frame(t(res), names(res))
dim(res)
res['4', ]<-colnames(res)
res<-rbin(res, colnames(res))
res<-rbind(res, colnames(res))
dim(res)
res[1:5, ]
res[, 1:5]
dim(res_
dim(res)
colnames(res)<-[1:493]
colnames(res)<-c(1:493)
res[, 1:5]
test<-data.frame(t(res))
head(test)
head(si_1$X3)
names(test)<-c("month", "date", "year")
head(si_1$X3)
head(test)
colnames(test)[4]<-"full"
head(test)
si_1<-data.frame(test[, 1:3], si_1)
si_1[1:5, 1:10]
tail(si_1)
tail(si_1[, 1:10])
colnames(si_1)[4:8]<-c("state", "some_number", "full_date", "id", "layer"
)
tail(si_1[, 1:10])
names(si_1)
si<-data.frame(si_1[, 9:12], si_1[, 1:8], si_1[, 13], si[, 51:58], si[, 14:50], si[, 59:109])
si<-data.frame(si_1[, 9:12], si_1[, 1:8], si_1[, 13], si_1[, 51:58], si_1[, 14:50], si_1[, 59:109])
names(si_1)
names(si)
colnames(si)[13]<-"Foaming.Status"
names(si)
head(si[, 1:21]
)
si$date<-NULL
head(si[, 1:21])
length(unique(si$id))
count(si, c(id, layer))
count(si, c("id", "layer"))
write.table(si, "fixed_inputs/sample_information_donotuseRemoved_rearranged.txt", sep="\t", row.names=F, quote=F)
dim(si)
str(head(si[, 1:21]))
head(si[, 1:20])
str(si[, 1:20])
test[]<-lapply(si[, 1:20], factor)
head(test)
test<-lapply(si[, 1:20], factor)
str(test)
si[, 1:20]<-lapply(si[, 1:20], factor)
str(si[, 1:20])
str(si[, 1:21])
sample_data(data1.phy)<-si
si[1:5, 1:20]
si_1[1:5, 1:20]
rownames(si)<-si$SAMPLES
si_1[1:5, 1:20]
si[1:5, 1:20]
data1.phy
sample_data(data1.phy)[1:5, 1:10]
sample_data(data1.phy)<-si
data1.phy
unique(sample_data(data1.phy)$Description)
si[, is.na(si$Description)]
subset(si, is.na(Description))
data2.phy<-subset_samples(data1.phy, is.na(Description))
data2.phy
data2.phy<-subset_samples(data1.phy, !is.na(Description))
source("~/code/R/ggplot_nmds_eclipse.R")
otu<-data.frame(otu_table(data1.phy))
si<-data.frame(sample_data(data1.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
dim(totu)
taxa_sums(data1.phy)
test<-subset_taxa(taxa_sums(data1.phy)==0, data1.phy)
data1.phy
subse
test<-subset_taxa(data1.phy, taxa_sums(data1.phy)==0)
dim(si)
names(si)
for (i in names(si[, 5:20])){
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
names$si
names(si)
si$full_date<-NULL
names(si)
for (i in names(si[, 5:19)){
for (i in names(si[, 5:19])){
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
for (i in names(si[, 10:12])){
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
for (i in names(si[, 15:19])){
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
for (i in names(si[, 16:19])){
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
names(si)
source("~/code/R/ordisurf_extraction.R")
source("~/code/R/ggplot_nmds_ordisurf.R")
getwd()
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(data1.si[, 19:107])){
         print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
         sf<-ordi.sf(data.mds, data1.si[, i])
         pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
for (i in names(data1.si[, 19:107])){
         print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
         sf<-ordi.sf(data.mds, si[, i])
         pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
         p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}
for (i in names(si[, 19:107])){
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
         sf<-ordi.sf(data.mds, si[, i])
         pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
         p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}
si[, 19]<-lapply(as.character(si[, 19], numberic)
)
si[, 19]<-lapply(as.character(si[, 19]), numberic)
si[, 19]<-lapply(as.character(si[, 19]), int)
head(si[, 19]
)
for (i in names(si[, 120:107])){
for (i in names(si[, 20:107])){
       print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
       sf<-ordi.sf(data.mds, si[, i])
       pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))      p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
for (i in names(si[, 20:107])){
       print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
       sf<-ordi.sf(data.mds, si[, i])
       pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
      p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}
head(si$Viscosity)
str(si$Viscosity)
length(unique(si$Viscosity)
)
for (i in names(si[, 47:107])){
       print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
       sf<-ordi.sf(data.mds, si[, i])
       pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
      p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
  print(p)
dev.off()
}
str(si[, 46:107])
test<-as.numeric(si$Viscosity)
head(test)
head(si$Viscosityq)
head(si$Viscosity)
test<-as.numeric(si[, 19:107])
test<-sapply(si[, 19:107], as.numeric)
str(test)
head(test)
test<-data.frame(sapply(si[, 19:107], as.numeric))
str(test)
si[, 19:107]<-data.frame(sapply(si[, 19:107], as.numeric))
str(si)
for (i in names(si[, 19:107])){
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
        print(p)
        dev.off()
}
length(unique(si$Select.Fish.Meal))
(unique(si$Select.Fish.Meal))
for (i in names(si[, 83:107])){
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
        print(p)
        dev.off()
}
length(unique(si$Tryptophan))
(unique(si$Tryptophan))
length(unique(si[, 83:107]))
sapply(si[, 83:107], function(x) length(unique(x)))
sapply(si[, 1:107], function(x) length(unique(x)))
sapply(si[, 83:107], function(x) unique(x))
for (i in names(si[, 83:107])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
savehistory("all_16s_history.R")
       sf<-ordi.sf(data.mds, si[, i])
       pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
      p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}
head(si$Viscosity)
str(si$Viscosity)
length(unique(si$Viscosity)
)
for (i in names(si[, 47:107])){
       print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
       sf<-ordi.sf(data.mds, si[, i])
       pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
      p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
  print(p)
dev.off()
}
str(si[, 46:107])
test<-as.numeric(si$Viscosity)
head(test)
head(si$Viscosityq)
head(si$Viscosity)
test<-as.numeric(si[, 19:107])
test<-sapply(si[, 19:107], as.numeric)
str(test)
head(test)
test<-data.frame(sapply(si[, 19:107], as.numeric))
str(test)
si[, 19:107]<-data.frame(sapply(si[, 19:107], as.numeric))
str(si)
for (i in names(si[, 19:107])){
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
        print(p)
        dev.off()
}
length(unique(si$Select.Fish.Meal))
(unique(si$Select.Fish.Meal))
for (i in names(si[, 83:107])){
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
        print(p)
        dev.off()
}
length(unique(si$Tryptophan))
(unique(si$Tryptophan))
length(unique(si[, 83:107]))
sapply(si[, 83:107], function(x) length(unique(x)))
sapply(si[, 1:107], function(x) length(unique(x)))
sapply(si[, 83:107], function(x) unique(x))
for (i in names(si[, 83:107])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
savehistory("all_16s_history.R")
save.image("all_16s_physeq.RData")
names(si)
for (i in names(si[, 5:18]){
for (i in names(si[, 5:18])){
print(i)}
for (i in names(si[, 5:18])){
print(head(si[, i])}
for (i in names(si[, 5:18])){
print(head(si[, i])
}
for (i in names(si[, 5:18])){
print(head(si[, i]))
}
test<-data1.phy
test<-subset_samples(data1.phy, is.na(Description))
test
test<-subset_samples(data1.phy, !is.na(Description))
test<-prune_taxa(data1.phy, taxa_sums(data1.phy)>0)
test<-prune_taxa(taxa_sums(test)>0, test)
test
data1.phy
for (i in names(si[, 5:18])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        working.phy<-subset_samples(data1.phy, !is.na(i))
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
        colorCount = length(unique(si1[, i]))
        getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
        pdf(paste("16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS(data.mds, si1[, i], colors)
        print(p)
        dev.off()
        print(adonis(data.trans~si1[,i]))
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
test<-subset_samples(data1.phy, !is.na(Description))
test
data1.phy
test<-prune_taxa(taxa_sums(test)>0, test)
test
otu<-data.frame(otu_table(test))
si<-data.frame(sample_data(test))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
colorCount = length(unique(si$Description))
        getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
pdf('test.pdf')
ggplot.NMDS(data.mds, si$Description, colors)
dev.off()
head(si$Name)
unique(si$Name)
data2.phy<-subset_samples(data1.phy, grepl("CF|CNF", Name))
data2.phy
data2.phy<-prune_samples(taxa_sums(data2.phy)>0, data2.phy)
data2.phy<-prune_taxa(taxa_sums(data2.phy)>0, data2.phy)
data2.phy
otu<-data.frame(otu_table(data2.phy))
si<-data.frame(sample_data(data2.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
names(si)
for (i in names(si[, 21:107])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
unique(si$Name)
colorCount = length(unique(si$Name))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si[, 21:107])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Name, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
savehistory("all_16s_history2.R")
save.image("all_16s_physeq.RData")
getwd()
ls()
colorCount = length(unique(si$ID_Feedmill))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si[, 20:107])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$ID_Feedmill, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
savehistory("all_16s_history2.R")
ls()
library(phyloseq)
library(plyr)
data0.phy
data1.phy
head(tax_table(data1.phy))
data1.melt<-psmelt(data1.phy)
dim(data1.melt)
data1.melt[1:5, 1:20]
head(data1.melt$genus)
savehistory("all_16s_history2)
savehistory("all_16s_history2")
library(phyloseq)
library(vegan)
library(plyr)
library(ggplot2)
library(RColorBrewer)
ls()
data0.phy
data1.
data1.phy
sample_data(data1.phy)[1:5, 1:20]
names(sample_data(data1.phy))
)
names(sample_data(data1.phy)[, 1:21])
names(sample_data(data1.phy)[, 21:107])
str(sample_data(data1.phy)[, 21:107])
test<-lapply(sample_data(data1.phy)[, 21:107], as.numeric)
str(test)
str(sample_data(data1.phy)[, 21:107])
str(test)
dim(si)
names(si)
si<-data.frame(sample_data(data1.phy))
names(si)
str(si[, 1:21])
str(si[, 21:107])
si[. 21:107]<-lapply(si[, 21:107], as.numeric)
si[, 21:107]<-lapply(si[, 21:107], as.numeric)
str(si[, 1:21])
str(si[, 21:107])
sample_data(data1.phy)<-si
data1.phy
str(sample_data(data1.phy)[, 1:23])
str(sample_data(data1.phy)[, 21:107])
source("~/code/R/ggplot_nmds_eclipse.R")
source("~/code/R/ggplot_nmds_ordisurf.R")
source("~/code/R/ordisurf_extraction.R")
str(si[, 1:21])
for (i in names(si[, 5:20])){
sink("16s_nmds_perfactor_adonis.txt", append = TRUE)
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        working.phy<-subset_samples(data1.phy, !is.na(i))
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
        colorCount = length(unique(si1[, i]))
        getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
        pdf(paste("16s_nmds_perfactor_",i,".pdf", sep=""))
        p<-ggplot.NMDS(data.mds, si1[, i], colors)
        print(p)
        dev.off()
        print(adonis(data.trans~si1[,i]))
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
        sink()
}
is.na(NA)
is.na("NA")
is.na(si$Description)
si$Description
heads(si$Description, n=50)
head(si$Description, n=50)
is.na(si$Description)
working.phy<-subset_samples(data1.phy, !is.na(Description))
working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
colorCount = length(unique(si1$Description))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
pdf("16s_nmds_perfactor_Description.pdf", length=12, width=15)
pdf("16s_nmds_perfactor_Description.pdf", height=12, width=15)
p<-ggplot.NMDS(data.mds, si1$Description, colors)
p
dev.off()
working.phy<-subset_samples(data1.phy, !is.na(Name))
working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
colorCount = length(unique(si1$Name))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
pdf("16s_nmds_perfactor_Name.pdf", length=12, width=15)
pdf("16s_nmds_perfactor_Name.pdf", height=12, width=15)
ggplot.NMDS(data.mds, si1$Name, colors)
dev.off()
working.phy<-subset_samples(data1.phy, !is.na(FK_Pump.Interval))
working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
colorCount = length(unique(si1$FK_Pump.Interval))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
pdf("16s_nmds_perfactor_FK_Pump.Interval.pdf", height=12, width=15)
ggplot.NMDS(data.mds, si1$FK_Pump.Interval, colors)
data1.phy
data2.phy
sample_data(data2.phy)$Foaming.Status
sample_data(data2.phy)$Name
sample_data(data1.phy)$Foaming.Status
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
data1.phy
dim(si)
si[1:5, 20:25]
for (i in names(si[, 20:108])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_byFS_493Samples",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$ID_Feedmill, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
     working.phy<-subset_samples(data1.phy, !is.na(Foaming.Status))
working.phy
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
colorCount = length(unique(si1$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si1[, 20:108])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si1[, i])
        pdf(paste("ordi_16s_nmds_byFS_493Samples",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
si1$FK_Pump.Interval<-as.numeric(si1$FK_Pump.Interval)
str(si1[, 1:5])
str(si1[, 19:25])
str(si[, 19:25])
for (i in names(si1[, 20:108])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si1[, i])
        pdf(paste("ordi_16s_nmds_byFS_493Samples",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
is.na(si$Tylan)
(si$Tylan)
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
data1.phy
dim(si)
str(si[, 1:22]
)
dim(si1)
str(si1[, 1:22])
sample_data(data1.phy)<-si1
si<-data.frame(sample_data(data1.phy))
str(si[, 1:22])
for (i in names(si[, 20:108])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        working.phy<-subset_samples(data1.phy, !is.na(i))
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
        sf<-ordi.sf(data.mds, si1[, i])
        pdf(paste("ordi_16s_nmds_byFS_493Samples_perfactor_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si1$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
savehistory("all_16s_history2.R")
library(phyloseq)
librara(ggplot2)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(plyr)
source("~/code/R/ggplot_nmds_eclipse.R")
source("~/code/R/ggplot_nmds_ordisurf.R")
source("~/code/R/ordisurf_extraction.R")
data1.phy
data2.phy
str(sample_data(data2.phy)[, 1:21]
)
si<-data.frame(sample_data(data2.phy))
dim(si)
str(si[, 1:20])
si[, 20:108]<-lapply(si[, 20:108], as.numeric)
str(si[, 1:25])
sample_data(data2.phy)<-si
str(sample_data(data2.phy)[, 1:21])
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si[, 20:108])){
sink("ordi_16s_nmds_byFS_cfVcnf_perfactor_envfit.txt", append = TRUE)
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        working.phy<-subset_samples(data2.phy, !is.na(i))
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si2<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
        sf<-ordi.sf(data.mds, si2[, i])
        pdf(paste("ordi_16s_nmds_byFS_cfVcnf_perfactor_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si2$Foaming.Status, colors, sf)
         print(p)
         dev.off()
         ef<-envfit(data.mds, data.frame(si2[, i]), permu=999, na.rm=T)
            print(ef)
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
savehistory("all_16s_history2.R")
save.image("all_16s_physeq.RData")
data2.phy
dim(si)
unique(si$Treated.Status)
count(si, c("Treated.Status")
)
data1.phy
test<-data.frame(sample_data(data1.phy))
count(test, c("Treated.Status")
)
myfunc <- function(v1) {
  deparse(substitute(v1))
}
myfunc(test$Treated.Status)
for (i in names(si[, 1:5])){
treatment<-myfunc(si[,i])
print(treatment)
}
for (i in names(si[, 1:5])){
i<-si[, i]
print(names(i))
print(head(i))
}
for (i in names(si[, 1:5])){
colnames(si[, i])<-si[, i]
print(names(i)
for (i in names(si[, 1:5])){
colnames(si[, i])<-si[, i]
for (i in names(si[, 1:5])){
print(si[, i])
}
for (i in names(si[, 1:2])){
print(si[, i])
}
for (i in names(si[, 1:2])){
print(i)
}
source("~/code/R/ggplot_nmds_ordisurf_2treatments.R")
pdf("test.pdf")
ggplot.NMDS.ordisurf(data.mds, si2$Foaming.Status, si2$Treated.Status, colors, sf)
dev.off()
getwd()
removeSource(ggplot.NMDS.ordisurf)
source("~/code/R/ggplot_nmds_ordisurf.R")
ls()
data3.phy<-subset_samples(data2.phy, Treated.Status=="0")
data3.phy
data2.phy
data2.phy<-prune_taxa(taxa_sums(data2.phy)>0, data2.phy)
data.ph
data2.phy
data3.phy<-prune_taxa(taxa_sums(data2.phy)>0, data3.phy)
data3.phy
data3.phy<-subset_samples(data2.phy, Treated.Status=="0")
data3.phy<-prune_taxa(taxa_sums(data3.phy)>0, data3.phy)
data3.phy
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si[, 20:108])){
sink("ordi_16s_nmds_byFS_cfVcnf_0trted_perfactor_envfit.txt", append = TRUE)
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        working.phy<-subset_samples(data3.phy, !is.na(i))
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
        si3<-data.frame(sample_data(working.phy))
        totu<-data.frame(t(otu))
        data.trans<-decostand(totu, "total")
        data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
        sf<-ordi.sf(data.mds, si3[, i])
        pdf(paste("ordi_16s_nmds_byFS_cfVcnf_0trted_perfactor_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si3$Foaming.Status, colors, sf)
        print(p)
        dev.off()
        ef<-envfit(data.mds, data.frame(si3[, i]), permu=999, na.rm=T)
        print(ef)
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
ls()
savehistory("all_16s_history2.R")
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Name, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
savehistory("all_16s_history2.R")
save.image("all_16s_physeq.RData")
getwd()
ls()
colorCount = length(unique(si$ID_Feedmill))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si[, 20:107])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$ID_Feedmill, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
savehistory("all_16s_history2.R")
ls()
library(phyloseq)
library(plyr)
data0.phy
data1.phy
head(tax_table(data1.phy))
data1.melt<-psmelt(data1.phy)
dim(data1.melt)
data1.melt[1:5, 1:20]
head(data1.melt$genus)
savehistory("all_16s_history2)
savehistory("all_16s_history2")
library(phyloseq)
library(vegan)
library(plyr)
library(ggplot2)
library(RColorBrewer)
ls()
data0.phy
data1.
data1.phy
sample_data(data1.phy)[1:5, 1:20]
names(sample_data(data1.phy))
)
names(sample_data(data1.phy)[, 1:21])
names(sample_data(data1.phy)[, 21:107])
str(sample_data(data1.phy)[, 21:107])
test<-lapply(sample_data(data1.phy)[, 21:107], as.numeric)
str(test)
str(sample_data(data1.phy)[, 21:107])
str(test)
dim(si)
names(si)
si<-data.frame(sample_data(data1.phy))
names(si)
str(si[, 1:21])
str(si[, 21:107])
si[. 21:107]<-lapply(si[, 21:107], as.numeric)
si[, 21:107]<-lapply(si[, 21:107], as.numeric)
str(si[, 1:21])
str(si[, 21:107])
sample_data(data1.phy)<-si
data1.phy
str(sample_data(data1.phy)[, 1:23])
str(sample_data(data1.phy)[, 21:107])
source("~/code/R/ggplot_nmds_eclipse.R")
source("~/code/R/ggplot_nmds_ordisurf.R")
source("~/code/R/ordisurf_extraction.R")
str(si[, 1:21])
for (i in names(si[, 5:20])){
sink("16s_nmds_perfactor_adonis.txt", append = TRUE)
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        working.phy<-subset_samples(data1.phy, !is.na(i))
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
        colorCount = length(unique(si1[, i]))
        getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
        pdf(paste("16s_nmds_perfactor_",i,".pdf", sep=""))
        p<-ggplot.NMDS(data.mds, si1[, i], colors)
        print(p)
        dev.off()
        print(adonis(data.trans~si1[,i]))
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
        sink()
}
is.na(NA)
is.na("NA")
is.na(si$Description)
si$Description
heads(si$Description, n=50)
head(si$Description, n=50)
is.na(si$Description)
working.phy<-subset_samples(data1.phy, !is.na(Description))
working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
colorCount = length(unique(si1$Description))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
pdf("16s_nmds_perfactor_Description.pdf", length=12, width=15)
pdf("16s_nmds_perfactor_Description.pdf", height=12, width=15)
p<-ggplot.NMDS(data.mds, si1$Description, colors)
p
dev.off()
working.phy<-subset_samples(data1.phy, !is.na(Name))
working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
colorCount = length(unique(si1$Name))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
pdf("16s_nmds_perfactor_Name.pdf", length=12, width=15)
pdf("16s_nmds_perfactor_Name.pdf", height=12, width=15)
ggplot.NMDS(data.mds, si1$Name, colors)
dev.off()
working.phy<-subset_samples(data1.phy, !is.na(FK_Pump.Interval))
working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
colorCount = length(unique(si1$FK_Pump.Interval))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
        colors = getPalette(colorCount)
pdf("16s_nmds_perfactor_FK_Pump.Interval.pdf", height=12, width=15)
ggplot.NMDS(data.mds, si1$FK_Pump.Interval, colors)
data1.phy
data2.phy
sample_data(data2.phy)$Foaming.Status
sample_data(data2.phy)$Name
sample_data(data1.phy)$Foaming.Status
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
data1.phy
dim(si)
si[1:5, 20:25]
for (i in names(si[, 20:108])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si[, i])
        pdf(paste("ordi_16s_nmds_byFS_493Samples",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$ID_Feedmill, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
     working.phy<-subset_samples(data1.phy, !is.na(Foaming.Status))
working.phy
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
colorCount = length(unique(si1$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si1[, 20:108])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si1[, i])
        pdf(paste("ordi_16s_nmds_byFS_493Samples",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
si1$FK_Pump.Interval<-as.numeric(si1$FK_Pump.Interval)
str(si1[, 1:5])
str(si1[, 19:25])
str(si[, 19:25])
for (i in names(si1[, 20:108])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si1[, i])
        pdf(paste("ordi_16s_nmds_byFS_493Samples",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
is.na(si$Tylan)
(si$Tylan)
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
data1.phy
dim(si)
str(si[, 1:22]
)
dim(si1)
str(si1[, 1:22])
sample_data(data1.phy)<-si1
si<-data.frame(sample_data(data1.phy))
str(si[, 1:22])
for (i in names(si[, 20:108])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        working.phy<-subset_samples(data1.phy, !is.na(i))
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si1<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
        sf<-ordi.sf(data.mds, si1[, i])
        pdf(paste("ordi_16s_nmds_byFS_493Samples_perfactor_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si1$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
savehistory("all_16s_history2.R")
library(phyloseq)
librara(ggplot2)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(plyr)
source("~/code/R/ggplot_nmds_eclipse.R")
source("~/code/R/ggplot_nmds_ordisurf.R")
source("~/code/R/ordisurf_extraction.R")
data1.phy
data2.phy
str(sample_data(data2.phy)[, 1:21]
)
si<-data.frame(sample_data(data2.phy))
dim(si)
str(si[, 1:20])
si[, 20:108]<-lapply(si[, 20:108], as.numeric)
str(si[, 1:25])
sample_data(data2.phy)<-si
str(sample_data(data2.phy)[, 1:21])
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si[, 20:108])){
sink("ordi_16s_nmds_byFS_cfVcnf_perfactor_envfit.txt", append = TRUE)
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        working.phy<-subset_samples(data2.phy, !is.na(i))
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si2<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
        sf<-ordi.sf(data.mds, si2[, i])
        pdf(paste("ordi_16s_nmds_byFS_cfVcnf_perfactor_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si2$Foaming.Status, colors, sf)
         print(p)
         dev.off()
         ef<-envfit(data.mds, data.frame(si2[, i]), permu=999, na.rm=T)
            print(ef)
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
savehistory("all_16s_history2.R")
save.image("all_16s_physeq.RData")
data2.phy
dim(si)
unique(si$Treated.Status)
count(si, c("Treated.Status")
)
data1.phy
test<-data.frame(sample_data(data1.phy))
count(test, c("Treated.Status")
)
myfunc <- function(v1) {
  deparse(substitute(v1))
}
myfunc(test$Treated.Status)
for (i in names(si[, 1:5])){
treatment<-myfunc(si[,i])
print(treatment)
}
for (i in names(si[, 1:5])){
i<-si[, i]
print(names(i))
print(head(i))
}
for (i in names(si[, 1:5])){
colnames(si[, i])<-si[, i]
print(names(i)
for (i in names(si[, 1:5])){
colnames(si[, i])<-si[, i]
for (i in names(si[, 1:5])){
print(si[, i])
}
for (i in names(si[, 1:2])){
print(si[, i])
}
for (i in names(si[, 1:2])){
print(i)
}
source("~/code/R/ggplot_nmds_ordisurf_2treatments.R")
pdf("test.pdf")
ggplot.NMDS.ordisurf(data.mds, si2$Foaming.Status, si2$Treated.Status, colors, sf)
dev.off()
getwd()
removeSource(ggplot.NMDS.ordisurf)
source("~/code/R/ggplot_nmds_ordisurf.R")
ls()
data3.phy<-subset_samples(data2.phy, Treated.Status=="0")
data3.phy
data2.phy
data2.phy<-prune_taxa(taxa_sums(data2.phy)>0, data2.phy)
data.ph
data2.phy
data3.phy<-prune_taxa(taxa_sums(data2.phy)>0, data3.phy)
data3.phy
data3.phy<-subset_samples(data2.phy, Treated.Status=="0")
data3.phy<-prune_taxa(taxa_sums(data3.phy)>0, data3.phy)
data3.phy
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si[, 20:108])){
sink("ordi_16s_nmds_byFS_cfVcnf_0trted_perfactor_envfit.txt", append = TRUE)
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        working.phy<-subset_samples(data3.phy, !is.na(i))
        working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
        si3<-data.frame(sample_data(working.phy))
        totu<-data.frame(t(otu))
        data.trans<-decostand(totu, "total")
        data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
        sf<-ordi.sf(data.mds, si3[, i])
        pdf(paste("ordi_16s_nmds_byFS_cfVcnf_0trted_perfactor_",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si3$Foaming.Status, colors, sf)
        print(p)
        dev.off()
        ef<-envfit(data.mds, data.frame(si3[, i]), permu=999, na.rm=T)
        print(ef)
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
ls()
savehistory("all_16s_history2.R")
ls()
data1.phy
library(phyloseq)
data1.phy
data2.phy
sample_data(data2.phy)$Treated.Status
unique(sample_data(data2.phy)$Foaming.Status)
unique(sample_data(data2.phy)$Description)
library(vegan)
data(BCI)
dim(BCI)
head(BCI[, 1:5])
simp <- diversity(BCI, "simpson")
head(simp)
length(simp)
dim(totu)
head(totu[, 1:5])
data1.phy
otu<-data.frame(otu_table(data1.phy))
totu<-t(otu)
head(totu[, 1:5])
dim(totu)
plot_richness
data1.rich<-estimate_richness(data1.phy, measures=c("Observed", "Chao1", "Shannon", "InvSimpson"))
dim(data1.rich)
head(data1.rich)
head(sample_data(data1.phy)[, 1:5])
data1.phy
si<-data.frame(sample_data(data1.phy))
head(si[, 1:5])
data1.rich$SAMPLES<-row.names(data1.rich)
head(data1.rich)
si1<-merge(si, data1.rich, "SAMPLES")
dim(si1)
head(si[, 1:5])
dim(si)
head(si[, 109:113])
head(si1[, 109:113])
data1.phy
sample_data(data1.phy)<-si1
head(si1[, 109:113])
row.names(si1)<-si1$SAMPLES
head(si1[, 1:5])
sample_data(data1.phy)<-si1
data1.phy
write.table(si1, "fixed_inputs/sample_information_donotuseRemoved_rearranged_w_alpha_indices.txt", quote=F, sep="\t", row.names=F)
source("~/code/R/ggplot_nmds_ordisurf.R")
source("~/code/R/ordisurf_extraction.R")
otu<-data.frame(otu_table(data1.phy))
totu<-t(otu)
head(totu[, 1:5])
si1<-data.frame(sample_data(data1.phy))
dim(si1)
dim(totu)
colnames(si1[, 108:113])
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
for (i in names(data1.si[, 109:113])){
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        ef<-envfit(data.mds, data.frame(data1.si[, i]), permu=999, na.rm=T)
        print(ef)
}
for (i in names(si1[, 109:113])){
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        ef<-envfit(data.mds, data.frame(si1[, i]), permu=999, na.rm=T)
        print(ef)
}
colorCount = length(unique(si1$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
library(RColorBrewer)
colorCount = length(unique(si1$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si1[, 109:113])){
        tryCatch({
        print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
        sf<-ordi.sf(data.mds, si1[, i])
        pdf(paste("ordi_16s_nmds_byFS_493Samples",i,".pdf", sep=""))
        p<-ggplot.NMDS.ordisurf(data.mds, si1$Foaming.Status, colors, sf)
         print(p)
         dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
getwd()
ls()
data3.phy
data2.phy
sample_sums(data1.phy)
sample_names(sample_sums(data1.phy)<1000)
sample_names(data1.phy, sample_sums(data1.phy)<1000)
sample_names(sample_sums(data1.phy)<1000, data1.phy)
sample_names()
sample_names
totu[rowSums(totu]<10000,]
totu[rowSums(totu)<10000,]
totu[rowSums(totu)<10000]
length(totu[rowSums(totu)<10000,])
totu[,rowSums(totu)<10000]
length(totu[,rowSums(totu)<10000])
totu[!rowSums(totu)>=10000, ]
head(totu[, 1:5])
rowSums(totu)
rowSums(totu)<10000
row.names(totu)[rowSums(totu)<10000, ]
rownames(totu)[rowSums(totu)<10000, ]
test<-prune_samples(sample_sums(data1.phy)>=10000, data1.phy)
test
data1.phy
data2.phy
unique(sample_data(data2.phy)$Description)
si<-data.frame(sample_data(data2.phy))
dim(si)
si2<-merge(si, data1.rich, "SAMPLES")
dim(si2)
row.names(si2)<-si2$SAMPLES
head(si2[, 1:5])
sample_data(data2.phy)<-si2
data2.phy
otu<-data.frame(otu_table(data2.phy))
totu<-t(otu)
dim(totu)
si2<-data.frame(sample_data(data2.phy))
dim(si2)
data.trans<-decostand(totu, "total")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
colorCount = length(unique(si2$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
for (i in names(si2[, 109:113])){
        tryCatch({
                print(c("Processing treatment:", i, "!!!!!!!!!!!!!!!!"))
                sf<-ordi.sf(data.mds, si2[, i])
            pdf(paste("ordi_16s_nmds_byFS_cfVcnf_",i,".pdf", sep=""))
                p<-ggplot.NMDS.ordisurf(data.mds, si2$Foaming.Status, colors, sf)
                print(p)
            dev.off()
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
pdf("cfVcnf_alpha_indices.pdf", width=15, height=8)
plot_richness(data2.phy, "Foaming.Status", "Name", measures=c("Observed", "Chao1", "Shannon", "InvSimpson")
)
dev.off()
pdf("all493_alpha_indices.pdf", width=15, height=8)
plot_richness(data1.phy, "Foaming.Status", "Name", measures=c("Observed", "Chao1", "Shannon", "InvSimpson"))
dev.off()
pdf("cfVcnf_alpha_indices.pdf", width=15, height=8)
p<-plot_richness(data2.phy, "Foaming.Status", "Name", measures=c("Observed", "Chao1", "Shannon", "InvSimpson"))
p+geom_boxplot(data=p$data, aes(x=Foaming.Status, y=value, color=NULL), alpha=0.1)
dev.off()
pdf("all493_alpha_indices.pdf", width=15, height=8)
p<-plot_richness(data1.phy, "Foaming.Status", "Name", measures=c("Observed", "Chao1", "Shannon", "InvSimpson"))
p+geom_boxplot(data=p$data, aes(x=Foaming.Status, y=value, color=NULL), alpha=0.1)
dev.off()
savehistory("all_16s_history3.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(vega)
library(vegan)
library(ggplot2)
library(RColorBrewer)
ls()
rm(BCI)
data2.phy
dim(si)
dim(totu)
decorana
decorana(totu)
hist(totu)
head(totu[, 1:5])
hist(colSums(totu))
hist(colSums(totu), breaks=5)
hist(totu[, 1])
hist(t(colSums(totu))
)
hist(rowSums(totu))
hist(rowSums(totu), breaks=5)
hist(rowSums(totu), breaks=1)
hist(rowSums(totu), breaks=2)
hist(rowSums(totu), breaks=100)
hist(rowSums(totu), breaks=1000)
hist(rowSums(totu), breaks=50)
hist(rowSums(totu), breaks=100)
hist(log(rowSums(totu), breaks=100))
hist(log(rowSums(totu)), breaks=100)
hist(log(rowSums(totu)), breaks=50)
hist(log(rowSums(totu)), breaks=20)
hist(log(rowSums(totu)), breaks=10)
hist(sqrt(rowSums(totu)), breaks=10)
hist(sqrt(rowSums(totu)), breaks=50)
hist(log(rowSums(totu)), breaks=50)
source("~/Documents/repos/code/R/ggplot_cca_ordisurf_2factors.R")
dim(data.trans)
dim(totu)
hist(rowSums(data.trans), break=50)
hist(rowSums(data.trans), breaks=50)
hist(rowSums(data.trans))
data.trans[, 1:5])
str(data.trans[, 1:5])
hist(data.frame(rowSums(data.trans)))
hist(data.frame(rowSums(data.trans)))
hist(log(rowSums(totu)), breaks=50)
test<-data.trans*10^5
head(test[, 1:5])
class(test)
test<-data.frame(test)
str(head(test[, 1:5]))
hist(rowSums(test))
head(test[, 1:5])
hist(test$OTU_0)
test$sum<-rowSums(test)
head(test$sum)
 data.trans<-decostand(totu, "hellinger")
head(data.trans[, 1:5])
hist(rowSums(data.trans))
hist(rowSums(data.trans), break= 50)
hist(rowSums(data.trans), breaks= 50)
hist(data.trans[, 1], breaks= 50)
hist(data.trans[, 2], breaks= 50)
hist(data.trans[, 3], breaks= 50)
hist(data.trans[, 100], breaks= 50)
hist(data.trans[, 300], breaks= 50)
ls()
ggplot.cca.ordisurf.2f
data.trans<-decostand(totu, "total")
data.cca<-capscale(data.trans ~ si$Foaming.Status, dist="bray")
scrs <- data.frame(scores(CCA, display="species"))
scrs <- data.frame(scores(data.cca, display="species"))
dim(scrs)
scrs[1:5, ]
plot(data.cca)
max(scrs)
min(scrs)
test<-subset(scrs, CAP1>0.5 | CAP1<-0.5)
test<-subset(scrs, CAP1>0.5 | CAP1< -0.5)
dim(test)
test
ls9)
ls()
tax<-data.frame(tax_table(data1.phy))
dim(tax)
cca_species<-merge(test, tax, row.names)
head(tax)
cca_species<-merge(test, tax, by.x=row.names(test), by.y=row.names(tax))
cca_species<-merge(test, tax, by="row.names")
dim(cca_species)
cca_species
test<-totu[row.names(totu)[cca_species$Row.names],]
dim(otu)
dim(totuu)
dim(totu)
test<-otu[row.names(otu)[cca_species$Row.names],]
dim(test)
head(test[ , 1:5])
head(row.names(otu))
test<-otu[cca_species$Row.names,]
head(row.names(otu))
head(row.names(test))
dim(test)
test[, 1:5]
rowSums(test>0)
CAP1<-scrs$CAP1
MDS1<-scrs$MDS1
cca.df<-data.frame(CAP1,MDS1)
dim(cca.df(
dim(cca.df)
head(cca.df)
scrs <- data.frame(scores(data.cca, display="sites"))
CAP1<-scrs$CAP1
MDS1<-scrs$MDS1
cca.df<-data.frame(CAP1,MDS1,DF_2F)
cca.df<-data.frame(CAP1,MDS1)
head(cca.df)
cca.df$type<-"sites"
cca.df$type<-"Samples"
ls()
cca_species
test<-data.frame(cca_species[, c("CAP1", "MDS1", "genus")])
head(test)
colnames(test)[3]<-"type"
cca.df<-rbind(test, cca.df)
dim(cca.df)
head(cca.df)
head(cca.df, n=20)
cca.df<-data.frame(CAP1,MDS1)
dim(cca.df)
cca.df$type<-"Samples"
cca.df<-data.frame(CAP1,MDS1, si(, c("Foaming.Status", "Name"))
)
cca.df<-data.frame(CAP1,MDS1, si[, c("Foaming.Status", "Name")])
head(cca.df)
head(scrs)
cca.df$type<-"Samples"
test
colnames(test)[3]<-"Foaming.Status"
test$Name<-"Species"
head(test)
head(cca.df)
cca.df$type<-NULL
head(cca.df)
cca.df.w.species<-rbind(test, cca.df)
dim(cca.df.w.species)
head(cca.df.w.species)
colnames(cca.df.w.species)[3]<-"FoamingStatus.and.Species"
head(cca.df.w.species)
X1<-ggplot() +
colorCount = length(unique(cca.df.w.species$FoamingStatus.and.Species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
X1<-ggplot() +
geom_poit(data=cca.df.w.species, aes(x=CAP1, y=MDS1, fill=FoamingStatus.and.Species, shape=Name), size=3, alpha=0.75) +
scale_shape_manual(values=c(21:30)) +
scale_fill_manual(values=c(COLORS), guide = guide_legend(override.aes = list(shape = 23)))+
theme_bw()
X1<-ggplot() +
geom_point(data=cca.df.w.species, aes(x=CAP1, y=MDS1, fill=FoamingStatus.and.Species, shape=Name), size=3, alpha=0.75) +
scale_shape_manual(values=c(21:30)) +
scale_fill_manual(values=c(COLORS), guide = guide_legend(override.aes = list(shape = 23)))+
theme_bw()
X1
X1<-ggplot() +
geom_point(data=cca.df.w.species, aes(x=CAP1, y=MDS1, fill=FoamingStatus.and.Species, shape=Name), size=3, alpha=0.75) +
scale_shape_manual(values=c(21:30)) +
scale_fill_manual(values=c(COLORS), guide = guide_legend(override.aes = list(shape = 23)))+
X1<-ggplot() +
geom_point(data=cca.df.w.species, aes(x=CAP1, y=MDS1, fill=FoamingStatus.and.Species, shape=Name), size=3, alpha=0.75) +
scale_fill_manual(values=c(colors), guide = guide_legend(override.aes = list(shape = 23)))+
theme_bw()
X1
X1<-ggplot() +
geom_point(data=cca.df.w.species, aes(x=CAP1, y=MDS1, fill=FoamingStatus.and.Species, shape=Name), size=3, alpha=0.75) +
scale_shape_manual(values=c(21:30)) +
scale_fill_manual(values=c(colors), guide = guide_legend(override.aes = list(shape = 23)))+
theme_bw()
X1
savehistory("all_16s_history4.R")
library(phyloseq)
library(ggplot2)
library(vegan)
library(RColorBrewer)
pdf("cca_sample_boundary_species.pdf")
X1
dev.off()
ls()
dim(data1.rich)
si<-data.frame(sample_data(data2.phy))
head(si[, 110:113])
test<-data.frame(si[, c("Foaming.Status","Observed", "Shannon", "InvSimpson")]) 
head(test)
library(plyr)
test.re<-ddply(test, .(Foaming.Status), t.test.plyr, "sds")
t.test.plyr <- function(x, var, mean=0 ){
  y <- rep(NA,10)
  y[6] <- nrow(x)[1]              # count observations
  if(nrow(x) < 2) return(y)       # exits if too less observations
  res <- t.test(x[var], mu=mean)  # doing the test
  y[1] <- res$statistic           # extract values of interest
  y[2] <- res$p.value      
  y[3] <- res$estimate     
  y[4] <- res$conf.int[1]  
  y[5] <- res$conf.int[2]  
  y[7] <- res$parameter    
  y[8] <- res$method       
  y[9] <- res$alternative  
  y[10] <- res$null.value   
  names(y) <- c("statistic","p.value","estimate","conf.int1", "conf.int2", "nobs","dof","method","alternative","null.value")
  y 
}
test.re<-ddply(test, .(Foaming.Status), t.test.plyr, "sds")
test.re<-ddply(test, .(Foaming.Status), t.test.plyr, "Observed")
test.re
test.re<-ddply(test, .(Foaming.Status), summarise, pval=t.test(Observed)$p.value)
test.re
summary(test)
test.re<-ddply(test, .(), summarise, pval=t.test(Foaming.Status ~ Observed)$p.value)
test.re<-ddply(test, .(Foaming.Status), summarise, pval=t.test(Foaming.Status ~ Observed)$p.value)
test.re<-ddply(test, .(Foaming.Status), summarise, pval=t.test(Observed)$p.value)
test.re
test0<-subset(test, Foaming.Status==0)
test1<-subset(test, Foaming.Status==1)
head(test0)
head(test1)
wilcox.test(test0$Observed, test1$Observed)
wilcox.test(test0$Shannon, test1$Shannon)
wilcox.test(test0$InvSimpson, test1$InvSimpson)
hist(test0)
hist(log(test0))
hist(log(test0[, 2:4]))
wilcox.test(log(test0$Observed), log(test1$Observed))
hist(test1)
wilcox.test(log(test0$Observed), test1$Observed)
hist(test0)
savehistory("all_16s_history4)
savehistory("all_16s_history4.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(vega)
library(vegan)
library(ggplot2)
library(RColorBrewer)
ls()
rm(BCI)
data2.phy
dim(si)
dim(totu)
decorana
decorana(totu)
hist(totu)
head(totu[, 1:5])
hist(colSums(totu))
hist(colSums(totu), breaks=5)
hist(totu[, 1])
hist(t(colSums(totu))
)
hist(rowSums(totu))
hist(rowSums(totu), breaks=5)
hist(rowSums(totu), breaks=1)
hist(rowSums(totu), breaks=2)
hist(rowSums(totu), breaks=100)
hist(rowSums(totu), breaks=1000)
hist(rowSums(totu), breaks=50)
hist(rowSums(totu), breaks=100)
hist(log(rowSums(totu), breaks=100))
hist(log(rowSums(totu)), breaks=100)
hist(log(rowSums(totu)), breaks=50)
hist(log(rowSums(totu)), breaks=20)
hist(log(rowSums(totu)), breaks=10)
hist(sqrt(rowSums(totu)), breaks=10)
hist(sqrt(rowSums(totu)), breaks=50)
hist(log(rowSums(totu)), breaks=50)
source("~/Documents/repos/code/R/ggplot_cca_ordisurf_2factors.R")
dim(data.trans)
dim(totu)
hist(rowSums(data.trans), break=50)
hist(rowSums(data.trans), breaks=50)
hist(rowSums(data.trans))
data.trans[, 1:5])
str(data.trans[, 1:5])
hist(data.frame(rowSums(data.trans)))
hist(data.frame(rowSums(data.trans)))
hist(log(rowSums(totu)), breaks=50)
test<-data.trans*10^5
head(test[, 1:5])
class(test)
test<-data.frame(test)
str(head(test[, 1:5]))
hist(rowSums(test))
head(test[, 1:5])
hist(test$OTU_0)
test$sum<-rowSums(test)
head(test$sum)
 data.trans<-decostand(totu, "hellinger")
head(data.trans[, 1:5])
hist(rowSums(data.trans))
hist(rowSums(data.trans), break= 50)
hist(rowSums(data.trans), breaks= 50)
hist(data.trans[, 1], breaks= 50)
hist(data.trans[, 2], breaks= 50)
hist(data.trans[, 3], breaks= 50)
hist(data.trans[, 100], breaks= 50)
hist(data.trans[, 300], breaks= 50)
ls()
ggplot.cca.ordisurf.2f
data.trans<-decostand(totu, "total")
data.cca<-capscale(data.trans ~ si$Foaming.Status, dist="bray")
scrs <- data.frame(scores(CCA, display="species"))
scrs <- data.frame(scores(data.cca, display="species"))
dim(scrs)
scrs[1:5, ]
plot(data.cca)
max(scrs)
min(scrs)
test<-subset(scrs, CAP1>0.5 | CAP1<-0.5)
test<-subset(scrs, CAP1>0.5 | CAP1< -0.5)
dim(test)
test
ls9)
ls()
tax<-data.frame(tax_table(data1.phy))
dim(tax)
cca_species<-merge(test, tax, row.names)
head(tax)
cca_species<-merge(test, tax, by.x=row.names(test), by.y=row.names(tax))
cca_species<-merge(test, tax, by="row.names")
dim(cca_species)
cca_species
test<-totu[row.names(totu)[cca_species$Row.names],]
dim(otu)
dim(totuu)
dim(totu)
test<-otu[row.names(otu)[cca_species$Row.names],]
dim(test)
head(test[ , 1:5])
head(row.names(otu))
test<-otu[cca_species$Row.names,]
head(row.names(otu))
head(row.names(test))
dim(test)
test[, 1:5]
rowSums(test>0)
CAP1<-scrs$CAP1
MDS1<-scrs$MDS1
cca.df<-data.frame(CAP1,MDS1)
dim(cca.df(
dim(cca.df)
head(cca.df)
scrs <- data.frame(scores(data.cca, display="sites"))
CAP1<-scrs$CAP1
MDS1<-scrs$MDS1
cca.df<-data.frame(CAP1,MDS1,DF_2F)
cca.df<-data.frame(CAP1,MDS1)
head(cca.df)
cca.df$type<-"sites"
cca.df$type<-"Samples"
ls()
cca_species
test<-data.frame(cca_species[, c("CAP1", "MDS1", "genus")])
head(test)
colnames(test)[3]<-"type"
cca.df<-rbind(test, cca.df)
dim(cca.df)
head(cca.df)
head(cca.df, n=20)
cca.df<-data.frame(CAP1,MDS1)
dim(cca.df)
cca.df$type<-"Samples"
cca.df<-data.frame(CAP1,MDS1, si(, c("Foaming.Status", "Name"))
)
cca.df<-data.frame(CAP1,MDS1, si[, c("Foaming.Status", "Name")])
head(cca.df)
head(scrs)
cca.df$type<-"Samples"
test
colnames(test)[3]<-"Foaming.Status"
test$Name<-"Species"
head(test)
head(cca.df)
cca.df$type<-NULL
head(cca.df)
cca.df.w.species<-rbind(test, cca.df)
dim(cca.df.w.species)
head(cca.df.w.species)
colnames(cca.df.w.species)[3]<-"FoamingStatus.and.Species"
head(cca.df.w.species)
X1<-ggplot() +
colorCount = length(unique(cca.df.w.species$FoamingStatus.and.Species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
X1<-ggplot() +
geom_poit(data=cca.df.w.species, aes(x=CAP1, y=MDS1, fill=FoamingStatus.and.Species, shape=Name), size=3, alpha=0.75) +
scale_shape_manual(values=c(21:30)) +
scale_fill_manual(values=c(COLORS), guide = guide_legend(override.aes = list(shape = 23)))+
theme_bw()
X1<-ggplot() +
geom_point(data=cca.df.w.species, aes(x=CAP1, y=MDS1, fill=FoamingStatus.and.Species, shape=Name), size=3, alpha=0.75) +
scale_shape_manual(values=c(21:30)) +
scale_fill_manual(values=c(COLORS), guide = guide_legend(override.aes = list(shape = 23)))+
theme_bw()
X1
X1<-ggplot() +
geom_point(data=cca.df.w.species, aes(x=CAP1, y=MDS1, fill=FoamingStatus.and.Species, shape=Name), size=3, alpha=0.75) +
scale_shape_manual(values=c(21:30)) +
scale_fill_manual(values=c(COLORS), guide = guide_legend(override.aes = list(shape = 23)))+
X1<-ggplot() +
geom_point(data=cca.df.w.species, aes(x=CAP1, y=MDS1, fill=FoamingStatus.and.Species, shape=Name), size=3, alpha=0.75) +
scale_fill_manual(values=c(colors), guide = guide_legend(override.aes = list(shape = 23)))+
theme_bw()
X1
X1<-ggplot() +
geom_point(data=cca.df.w.species, aes(x=CAP1, y=MDS1, fill=FoamingStatus.and.Species, shape=Name), size=3, alpha=0.75) +
scale_shape_manual(values=c(21:30)) +
scale_fill_manual(values=c(colors), guide = guide_legend(override.aes = list(shape = 23)))+
theme_bw()
X1
savehistory("all_16s_history4.R")
library(phyloseq)
library(ggplot2)
library(vegan)
library(RColorBrewer)
pdf("cca_sample_boundary_species.pdf")
X1
dev.off()
ls()
dim(data1.rich)
si<-data.frame(sample_data(data2.phy))
head(si[, 110:113])
test<-data.frame(si[, c("Foaming.Status","Observed", "Shannon", "InvSimpson")]) 
head(test)
library(plyr)
test.re<-ddply(test, .(Foaming.Status), t.test.plyr, "sds")
t.test.plyr <- function(x, var, mean=0 ){
  y <- rep(NA,10)
  y[6] <- nrow(x)[1]              # count observations
  if(nrow(x) < 2) return(y)       # exits if too less observations
  res <- t.test(x[var], mu=mean)  # doing the test
  y[1] <- res$statistic           # extract values of interest
  y[2] <- res$p.value      
  y[3] <- res$estimate     
  y[4] <- res$conf.int[1]  
  y[5] <- res$conf.int[2]  
  y[7] <- res$parameter    
  y[8] <- res$method       
  y[9] <- res$alternative  
  y[10] <- res$null.value   
  names(y) <- c("statistic","p.value","estimate","conf.int1", "conf.int2", "nobs","dof","method","alternative","null.value")
  y 
}
test.re<-ddply(test, .(Foaming.Status), t.test.plyr, "sds")
test.re<-ddply(test, .(Foaming.Status), t.test.plyr, "Observed")
test.re
test.re<-ddply(test, .(Foaming.Status), summarise, pval=t.test(Observed)$p.value)
test.re
summary(test)
test.re<-ddply(test, .(), summarise, pval=t.test(Foaming.Status ~ Observed)$p.value)
test.re<-ddply(test, .(Foaming.Status), summarise, pval=t.test(Foaming.Status ~ Observed)$p.value)
test.re<-ddply(test, .(Foaming.Status), summarise, pval=t.test(Observed)$p.value)
test.re
test0<-subset(test, Foaming.Status==0)
test1<-subset(test, Foaming.Status==1)
head(test0)
head(test1)
wilcox.test(test0$Observed, test1$Observed)
wilcox.test(test0$Shannon, test1$Shannon)
wilcox.test(test0$InvSimpson, test1$InvSimpson)
hist(test0)
hist(log(test0))
hist(log(test0[, 2:4]))
wilcox.test(log(test0$Observed), log(test1$Observed))
hist(test1)
wilcox.test(log(test0$Observed), test1$Observed)
hist(test0)
savehistory("all_16s_history4)
savehistory("all_16s_history4.R")
library(phyloseq)
library(ggplot)
library(ggplot2)
dim(si)
dim(totu)
names(si)
names(tax_table(data2.phy))
colnames(tax_table(data2.phy))
 cytophaga<-subset_taxa(data2.phy, genus = "Cytophaga")
dim(cytophaga)
cytophaga
 cytophaga<-prune_taxa(data2.phy, genus == "Cytophaga")
 cytophaga<-subset_taxa(data2.phy, genus == "Cytophaga")
cytophaga
 cytophaga<-prune_samples(cytophaga, sample_sums(cytophaga)>0)
prune_samples
 cytophaga<-prune_samples(sample_sums(cytophaga)>0, cytophaga)
cytophaga
head(otu_table(cytophaga))
data3.phy
 cytophaga<-prune_taxa(taxa_sums(cytophaga)>0, cytophaga)
cytophaga
data1.phy
data2.phy
data1.relative.phy<-data1.phy
data1.relative.phy<-transform_sample_counts(data1.relative.phy, function(x) x/sum(x))
head(sample_sums(data1.relative.phy))
data2.relative.phy<-subset_samples(data1.relative.phy, grepl("CF|CNF", Name))
data2.relative.phy
data2.phy
data2.relative.phy<-prune_taxa(taxa_sums(data2.relative.phy)>0, data2.relative.phy)
data2.relative.phy
cytophaga<-subset_samples(data2.relative.phy, genus == "Cytophaga")
cytophaga<-subset_taxa(data2.relative.phy, genus == "Cytophaga")
cytophaga
cytophaga<-prune_samples(sample_sums(cytophaga)>0, cytophaga)
cytophaga
unique(sample_data(cytophaga)$Foaming.Status)
library(plyr)
count(sample_data(cytophaga)$Foaming.Status)
cytophaga.nf<-subset_samples(cytophaga, Foaming.Status == "0")
cytophaga.f<-subset_samples(cytophaga, Foaming.Status == "1")
cytophaga.nf
cytophaga.nf<-prune_taxa(taxa_sums(cytophaga.nf)>0, cytophaga.nf)
cytophaga.f<-prune_taxa(taxa_sums(cytophaga.f)>0, cytophaga.f)
cytophaga.nf
cytophaga.f
test<-data.frame(otu_table(cytophaga.nf))
dim(test)
test
data.frame(otu_table(cytophaga.f))
test<-sample_sums(cytophaga.nf)
dim(test)
test
test<-data.frame(sample_sums(cytophaga.nf))
head(test)
test$SAMPLES<-row.names(test)
head(test)
colnames(test)[1]<-"Cytophaga.NF"
head(test)
dim(si)
names(si)
test<-merge(test, si[, c("SAMPLES", "Surface.Tension")], "SAMPLES")
head(test)
lenth(unique(test$Surface.Tension))
length(unique(test$Surface.Tension))
(unique(test$Surface.Tension))
test<-sample_sums(cytophaga)
head(test)
test<-data.frame(sample_sums(cytophaga))
head(test)
test$SAMPLES<-row.names(test)
test<-merge(test, si[, c("SAMPLES", "Surface.Tension")], "SAMPLES")
dim(test)
head(test)
length(is.na(test$Surface.Tension))
length(test$Surface.Tension == "NA")
length(unique(test$Surface.Tension))
test$Surface.Tension==NA
test$Surface.Tension=="NA"
test<-subset(test, test$Surface.Tension != "NA")
test
plot(x=test[, 2], y=test[, 3])
lm(test[,2]~test[,3])
anova(lm(test[,2]~test[,3]))
test<-merge(test, si[, c("SAMPLES", "Foaming.Status")], "SAMPLES")
test
ggplot(test, aes(x=sample_sums.cytophaga., y=Surface.Tension, color=Foaming.Status, group=Foaming.Status))
ggplot(test, aes(x=sample_sums.cytophaga., y=Surface.Tension, color=Foaming.Status, group=Foaming.Status))+geom_line()
ggplot(test, aes(x=sample_sums.cytophaga., y=Surface.Tension, color=Foaming.Status, group=Foaming.Status))+geom_point()
save.image("all_16s_physeq.RData")
ruminococcus2<-subset_taxa(data2.relative.phy, genus=="Ruminococcus2")
ruminococcus2
ruminococcus2<-prune_samples(sample_sums(ruminococcus2)>0, ruminococcus2)
ruminococcus2
test<-data.frame(sample_sums(ruminococcus2)
)
head(test)
test$SAMPLES<-row.names(test)
test<-merge(test, si[, c("SAMPLES", "Foaming.Status", "Acetic.Acid")], "SAMPLES")
head(test)
length(unique(test$Acetic.Acid))
test<-subset(test, test$Acetic.Acid != "NA")
dim(test)
test
ggplot(test, aes(x=sample_sums.ruminococcus2., y=Acetic.Acid, color=Foaming.Status, group=Foaming.Status))+geom_point()
ddply(test, .(Foaming.Status), lm, formula = Acetic.Acid ~ sample_sums.ruminococcus2.)
dlply(test, .(Foaming.Status), lm, formula = Acetic.Acid ~ sample_sums.ruminococcus2.)
mods<-dlply(test, .(Foaming.Status), lm, formula = Acetic.Acid ~ sample_sums.ruminococcus2.)
ldply(mods, coef)
ldply(mods, anova)
length(unique(si$LCFA))
length(unique(si$Oleic.Acid))
length(unique(si$SCFA))
length(unique(si$Acetic.Acid))
for (i in colnames(si)){
print(i)
print(length(unique(si[, i])))
}
cytophaga
test<-data.frame(sample_sums(cytophaga))
head(test)
test$SAMPLE<-row.names(test)
head(test)
test<-merge(test, si[, c("SAMPLES", "Foaming.Status", "Foam.Stability")], "SAMPLES")
colnames(test)[2]<-"SAMPLES")
colnames(test)[2]<-"SAMPLES"
test<-merge(test, si[, c("SAMPLES", "Foaming.Status", "Foam.Stability")], "SAMPLES")
head(test)
ggplot(test, aes_string(x=test[, 2], y=test[, 4], color=test[, 3], group=test[, 3]))+geom_point()
ggplot(test, aes_string(x=colnames(testt)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_line()
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()
afipia<-subset_taxa(data2.relative.phy, genus == "Afipia")
afipia
afipia<-prune_samples(sample_sums(afipia), afipia)
afipia<-prune_samples(sample_sums(afipia)>0, afipia)
afipia
test<-data.frame(sample_sums(afipia))
test$SAMPLES<-row.names(test)
head(test)
test<-merge(test, si[, c("SAMPLES", "Foaming.Status", "Foam.Stability")], "SAMPLES")
head(test)
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()
ugamma<-subset_taxa(data2.relative.phy, genus == "unclassified_Gammaproteobacteria")
ugamma<-prune_samples(sample_sums(ugamma)>0, ugamma)
ugamma
test<-data.frame(sample_sums(ugamma))
test$SAMPLES<-row.names(test)
test<-merge(test, si[, c("SAMPLES", "Foaming.Status", "Foam.Stability")], "SAMPLES")
head(test)
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()
head(otu_table(ugamma)[, 1:5])
otu_table(ugamma)[, 1:5]
max(taxa_sums(ugamma))
prune_taxa(taxa_sums(ugamma)==max(taxa_sums(ugamma)), ugamma)
otu_table(prune_taxa(taxa_sums(ugamma)==max(taxa_sums(ugamma)), ugamma))[, 1:5]
head(test)
test<-data.frame(otu_table(ugamma)["OTU_4157",])
dim(test)
test<-t(data.frame(otu_table(ugamma)["OTU_4157",]))
head(test)
test$SAMPLES<-row.names(test)
test<-data.frame(t(data.frame(otu_table(ugamma)["OTU_4157",])))
head(test)
test$SAMPLES<-row.names(test)
head(test)
test<-merge(test, si[, c("SAMPLES", "Foaming.Status", "Foam.Stability")], "SAMPLES")
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()
anova(lm(test$OTU_4157 ~ test$Foam.Stability))
anova(lmer(test$OTU_4157 ~ test$Foam.Stability))
library(lme4)
anova(lmer(test$OTU_4157 ~ test$Foam.Stability + ))
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()+theme_bw()
save.image("all_16s_physeq.RData")
test1<-data.frame(sample_sums(ugamma))
rm(test1)
pdf("ugamma_otu4157_to_foam_stability.pdf")
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()+theme_bw()
dev.off()
test<-data.frame(sample_sums(ugamma))
test$SAMPLES<-row.names(test)
test<-merge(test, si[, c("SAMPLES", "Foaming.Status", "Foam.Stability")], "SAMPLES")
pdf("ugamma_sum_to_foam_stability.pdf")
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()+theme_bw()
dev.off()
anova(lm(test$sample_sums.ugamma. ~ test$Foam.Stability))
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()+theme_bw()
savehistory("all_16s_history5.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
ls()
dim(si)
names(si)
for (i in colnames(si[, 20:113])){
l = length(unique(si[, i]))
print (paste(i, l, sep=": "))
}
for (i in colnames(si[, 20:113])){
l = length(unique(si[, i]))
if (l > 40){
print (paste(i, l , sep = ": "))
}}
for (i in colnames(si[, 20:113])){
l = length(unique(si[, i]))
if (l > 40){
l.m<-character()
for (i in colnames(si[, 20:113])){
l = length(unique(si[, i]))
if (l > 40){
l.m<-append(l.m, i)
}
print(l.m)
}
l.m<-character()
for (i in colnames(si[, 20:113])){
l = length(unique(si[, i]))
if (l > 40){
l.m<-append(l.m, i)
}
}
l.m
l.m<-l.m[-c(42:46)]
l.m
str(l.m)
write(l.m, "list_of_measurements_more_than_40values.txt", quote=F)
write(l.m, "list_of_measurements_more_than_40values.txt")
test<-data.frame(si[, c("Foam.Stability", "Foaming.Capacity", "Foaming.Status")])
head(test)
ggplot(test, aes_string(x=colnames(test)[2], y=colnames(test)[4], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()+the
ggplot(test, aes_string(x=colnames(test)[1], y=colnames(test)[2], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()+the
ggplot(test, aes_string(x=colnames(test)[1], y=colnames(test)[2], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()+theme_bw()
library(plyr)
count(test$Foaming.Status)
test<-data.frame(si[, c("Foam.Stability", "Manure.Temperature", "Foaming.Status")])
ggplot(test, aes_string(x=colnames(test)[1], y=colnames(test)[2], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()+theme_bw()
mods<-dlply(test, .(Foaming.Status), lm, formula = Foam.Stabiligy ~ Manure.Temperature)
mods<-dlply(test, .(Foaming.Status), lm, formula = Foam.Stability ~ Manure.Temperature)
ldply(mods, anova)
ggplot(test, aes_string(x=colnames(test)[1], y=colnames(test)[2], color=colnames(test)[3], group=colnames(test)[3]))+geom_point()+theme_bw()+geom_smooth(method='lm',formula=y~x)
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load(
"all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
ls()
data1.phy
data2.phy
data2.relative.phy
names(si)
head(si[, c("month", "years", "full_data")])
head(si[, c("month", "years", "full_date")])
head(si[, c("month", "year", "full_date")])
head(si[, c("month", "year", "full_date", "Long_Sample_ID")])
head(si[, c("month", "year", "full_date", "Long_Sample_ID", "id")])
si[, c("month", "year", "full_date", "Long_Sample_ID", "id")][id == "F498", ]
si[, c("month", "year", "full_date", "Long_Sample_ID", "id")][si$id == "F498", ]
ls()
data2.relative.phy
ls
head(otu_table(data2.relative.phy)[, 1:5])
dim(totu)
head(totu[, 1:5])
library(ineq)
s97.lcolc<-Lc(totu["S_97", ])
head(s97.lcolc)
s97.lcdf<-data.frame(L=rev(1-s97.lcolc$L), p=s97.lcolc$p, Uprob=c(1:length(s97.lcolc$L)/length(s97.lcolc$L)))
dim(s97.lcdf)
head(s97.lcdf)
index  <- which(lcdf$L >= 0.499 & lcdf$L <= 0.501)[1]
index  <- which(s97.lcdf$L >= 0.799 & s97.lcdf$L <= 0.801)[1]
index
ypos<-s97.lcdf$L[index]
yposs<-c(0, ypos)
xpos<-index/length(s97.lcdf$L)
xposs<-c(0, xpos)
xposs
ypositions <- data.frame(x = xposs, y = c(ypos,ypos))
ypositions
xpositions <- data.frame(x = c(xpos,xpos), y = yposs)
p <- ggplot(s97.lcdf, aes(x = Uprob, y = L)) + geom_line(colour = hcl(h=15, l=65, c=100)) + geom_line(aes(x = p, y = p))
p
head(si[, 1"5]
head(si[, 1:5]
)
head(si[, c("SAMPLES", "Foaming.Status")]
)
s98.lcolc<-Lc(totu["S_97", ])
s98.lcdf<-data.frame(L=rev(1-s98.lcolc$L), p=s98.lcolc$p, Uprob=c(1:length(s98.lcolc$L)/length(s98.lcolc$L)))
p <- ggplot(s98.lcdf, aes(x = Uprob, y = L)) + geom_line(colour = hcl(h=15, l=65, c=100)) + geom_line(aes(x = p, y = p))
p
p <- ggplot(s98.lcdf, aes(x = Uprob, y = L)) + geom_line() + geom_line(aes(x = p, y = p))
p
p <- ggplot(s98.lcdf, aes(x = Uprob, y = L)) + geom_line(colour = red)) + geom_line(aes(x = p, y = p))
p <- ggplot(s98.lcdf, aes(x = Uprob, y = L)) + geom_line(colour = "red")) + geom_line(aes(x = p, y = p))
p <- ggplot(s98.lcdf, aes(x = Uprob, y = L)) + geom_line(colour = "red") + geom_line(aes(x = p, y = p))
p
p <- ggplot(s98.lcdf, aes(x = Uprob, y = L)) + geom_line(colour = hcl(h=15, l=65, c=100)) + geom_line(aes(x = p, y = p))
p
head(s98.lcdf)
s98.lcdf$fs<-"nf"
s97.lcdf$fs<-"f"
test<-rbind(s98.lcdf, s97.lcdf)
p <- ggplot(test, aes(x = Uprob, y = L, color=fs)) + geom_line(colour = fs) + geom_line(aes(x = p, y = p))
head(test)
p <- ggplot(test, aes(x = Uprob, y = L, color=fs)) + geom_line(aes(x = p, y = p))
p
p <- ggplot(test, aes(x = Uprob, y = L)) + geom_line(colour = fs) + geom_line(aes(x = p, y = p, color=fs))
p <- ggplot(test, aes(x = p, y = p, color=fs))  + geom_line(aes(x = Uprob, y = L))
p
p <- ggplot(test, aes(x = Uprob, y = L)) + geom_line(aes(x = p, y = p, color=fs))
p
min(taxa_sums(data2.phy))
test<-prune_taxa(taxa_sums(data2.relative.phy)>=0.01, data2.relative.phy)
test
min(sample_sums(test))
max(sample_sums(test))
s97.lcolc<-Lc(totu["S_97", ])
s98.lcolc<-Lc(totu["S_98", ])
s98.lcdf<-data.frame(L=rev(1-s98.lcolc$L), p=s98.lcolc$p, Uprob=c(1:length(s98.lcolc$L)/length(s98.lcolc$L)))
s97.lcdf<-data.frame(L=rev(1-s97.lcolc$L), p=s98.lcolc$p, Uprob=c(1:length(s97.lcolc$L)/length(s97.lcolc$L)))
s98.lcdf$fs<-"nf"
s97.lcdf$fs<-"f"
p <- ggplot(test, aes(x = Uprob, y = L)) + geom_line(colour = blk) + geom_line(aes(x = p, y = p, color=fs))
test<-rbind(s98.lcdf, s97.lcdf)
p <- ggplot(test, aes(x = Uprob, y = L)) + geom_line(colour = blk) + geom_line(aes(x = p, y = p, color=fs))
p <- ggplot(test, aes(x = Uprob, y = L)) + geom_line(colour = "black") + geom_line(aes(x = p, y = p, color=fs))
p
p <- ggplot(test, aes(x = Uprob, y = L, color=fs)) + geom_line(colour = "black") + geom_line(aes(x = p, y = p))
p
 test<-prune_taxa(taxa_sums(data2.relative.phy)>=0.1, data2.relative.phy)
 test<-prune_taxa(taxa_sums(data2.relative.phy)>=0.01, data2.relative.phy)
totu<-t(data.frame(otu_table(test)))
dim(totu)
totu[1:5, 1:5]
 s97.lcolc<-Lc(totu["S_97", ])
 s98.lcolc<-Lc(totu["S_98", ])
s98.lcdf<-data.frame(L=rev(1-s98.lcolc$L), p=s98.lcolc$p, Uprob=c(1:length(s98.lcolc$L)/length(s98.lcolc$L)))
s97.lcdf<-data.frame(L=rev(1-s97.lcolc$L), p=s98.lcolc$p, Uprob=c(1:length(s97.lcolc$L)/length(s97.lcolc$L)))
s98.lcdf$fs<-"nf"
s97.lcdf$fs<-"f"
test<-rbind(s98.lcdf, s97.lcdf)
p <- ggplot(test, aes(x = Uprob, y = L, color=fs)) + geom_line(aes(x = p, y = p))
p
p <- ggplot(test, aes(x = Uprob, y = L, color=fs)) + geom_line(colour = "black") + geom_line(aes(x = p, y = p))
p
p <- ggplot(test, aes(x = Uprob, y = L, color=fs)) + geom_line(colour = fs) + geom_line(aes(x = p, y = p))
p <- ggplot(test, aes(x = Uprob, y = L, color=fs)) + geom_line() + geom_line(aes(x = p, y = p))
p
p <- ggplot(test, aes(x = Uprob, y = L, color=fs)) + geom_line() + geom_line(aes(x = p, y = p, color="black"))
p
p <- ggplot(test, aes(x = Uprob, y = L, color=fs)) + geom_line() + geom_line(aes(x = p, y = p))
p
head(rownames(totu)F)
head(rownames(totu))
for (i in rownames(totu[1:5, ])){
print (head(totu["i",]))
}
for (i in rownames(totu[1:5, ])){
print (head(totu[i,]))
}
head(si[, "name"])
head(si[, "Name"])
final_lc<-data.frame()
for (i in rownames(totu)){
lcolc<-Lc(totu[i, ])
lcdf<-data.frame(SAMPLES=i, L=rev(1-lcolc$L), p=lcolc$p, Uprob=c(1:length(lcolc$L)/length(lcolc$L)))
lcdf.merge<-merge(lcdf, si[, c("SAMPLES", "Foaming.Status", "Name"), "SAMPLES"]
final_lc<-rbind(final_lc, lcdf.merged)
print(paste("finished sample ", i, sep=""))
}
final_lc<-data.frame()
for (i in rownames(totu)){
lcolc<-Lc(totu[i, ])
lcdf<-data.frame(SAMPLES=i, L=rev(1-lcolc$L), p=lcolc$p, Uprob=c(1:length(lcolc$L)/length(lcolc$L)))
lcdf.merge<-merge(lcdf, si[, c("SAMPLES", "Foaming.Status", "Name")], "SAMPLES")
final_lc<-rbind(final_lc, lcdf.merged)
print(paste("finished sample ", i, sep=""))
}
final_lc<-data.frame()
for (i in rownames(totu)){
lcolc<-Lc(totu[i, ])
lcdf<-data.frame(SAMPLES=i, L=rev(1-lcolc$L), p=lcolc$p, Uprob=c(1:length(lcolc$L)/length(lcolc$L)))
lcdf.merge<-merge(lcdf, si[, c("SAMPLES", "Foaming.Status", "Name")], "SAMPLES")
final_lc<-rbind(final_lc, lcdf.merge)
print(paste("finished sample ", i, sep=""))
}
dim(final_lc)
head(final_lc)
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line() + geom_line(aes(x = p, y = p))
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.2) + geom_line(aes(x = p, y = p))
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.02) + geom_line(aes(x = p, y = p))
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.01) + geom_line(aes(x = p, y = p))
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.1) + geom_line(aes(x = p, y = p))
p
library(grid)
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.1) + geom_line(aes(x = p, y = p))+facet_grid(~Name)
p
library(plyr)
dim(s97.lcdf)
dim(s98.lcdf)
test<-prune_taxa(taxa_sums(data2.relative.phy)>=0.1, data2.relative.phy)
test
min(sample_sums(test))
max(sample_sums(test))
test<-prune_taxa(taxa_sums(data2.relative.phy)>=0.05, data2.relative.phy)
max(sample_sums(test))
min(sample_sums(test))
test
totu<-t(data.frame(otu_table(test)))
final_lc<-data.frame()
for (i in rownames(totu)){
lcolc<-Lc(totu[i, ])
lcdf<-data.frame(SAMPLES=i, index=c(1:length(lcolc$L)), L=rev(1-lcolc$L), p=lcolc$p, Uprob=c(1:length(lcolc$L)/length(lcolc$L)))
lcdf.merge<-merge(lcdf, si[, c("SAMPLES", "Foaming.Status", "Name")], "SAMPLES")
final_lc<-rbind(final_lc, lcdf.merge)
print(paste("finished sample ", i, sep=""))
}
head(final_lc)
length(unique(final_lc$index))
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.1) + geom_line(aes(x = p, y = p))
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.2) + geom_line(aes(x = p, y = p))
p
test1<-ddply(final_lc, .(index, Name), summarize, avgL=mean(L), avgP=mean(p), avgUprob=mean(Uprob))
head(test1)
p <- ggplot(test1, aes(x = avgUprob, y = avgL, color=Name)) + geom_line(alpha=0.2) + geom_line(aes(x = avgP, y = avgP))
p
data2.phy
min(taxa_sums(data2.phy))
hist(taxa_sums(data2.phy))
max(taxa_sums(data2.phy))
test<-prune_taxa(taxa_sums(data2.phy)>=20, data2.phy)
test
min(sample_sums(test)
)
totu<-t(data.frame(otu_table(test)))
final_lc<-data.frame()
for (i in rownames(totu)){
lcolc<-Lc(totu[i, ])
lcdf<-data.frame(SAMPLES=i, index=c(1:length(lcolc$L)), L=rev(1-lcolc$L), p=lcolc$p, Uprob=c(1:length(lcolc$L)/length(lcolc$L)))
lcdf.merge<-merge(lcdf, si[, c("SAMPLES", "Foaming.Status", "Name")], "SAMPLES")
final_lc<-rbind(final_lc, lcdf.merge)
print(paste("finished sample ", i, sep=""))
}
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.2) + geom_line(aes(x = p, y = p))
p
head(final_lc)
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.2) + geom_line(aes(x = p, y = p)) + theme_bw() + scale_color
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.2) + geom_line(aes(x = p, y = p)) + theme_bw() + scale_color_manual(values=c(
"black", "red"))
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(alpha=0.2) + geom_line(aes(x = p, y = p)) + theme_bw() + scale_color_manual(values=c("black", "red"), alpha=0.2)
savehistory("temp.R")
load("all_16s_physeq.RData")
ls()
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(aes(x = p, y = p)) + theme_bw() 
library(ggplot2)
library(grid)
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name)) + geom_line(aes(x = p, y = p)) + theme_bw() 
p
p <- ggplot(final_lc, aes(x = p, y = p, color=Name))+geom_line()
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name))+geom_line(alpha=0.2)
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name))+geom_line(alpha=0.1)
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name))+geom_line(alpha=0.1)+geom_line(aes(x = p, y = p, color="black"))
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name))+geom_line(alpha=0.1)+geom_line(aes(x = p, y = p), color="black")
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name))+geom_line(alpha=0.1)+geom_line(aes(x = p, y = p), color="black", style="dash") + theme_bw()
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name))+geom_line(size=2, alpha=0.1)+
geom_line(aes(x = p, y = p), color="black", style="dash") + 
theme(aspect.ratio=1)+theme_bw()+
theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+
theme(legend.title=element_text(size=15),legend.text=element_text(size=15))+
xlab("Percentage of OTU") + ylab("Percentage of Sequence Abundance")
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name))+geom_line(size=2, alpha=0.1)+
geom_line(aes(x = p, y = p), color="black", linetype="dash") + 
theme(aspect.ratio=1)+theme_bw()+
theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+
theme(legend.title=element_text(size=15),legend.text=element_text(size=15))+
xlab("Percentage of OTU") + ylab("Percentage of Sequence Abundance")
p
p <- ggplot(final_lc, aes(x = Uprob, y = L, color=Name))+geom_line(size=2, alpha=0.1)+
geom_line(aes(x = p, y = p), color="black") + 
theme(aspect.ratio=1)+theme_bw()+
theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+
theme(legend.title=element_text(size=15),legend.text=element_text(size=15))+
xlab("Percentage of OTU") + ylab("Percentage of Sequence Abundance")
p
dim(test)
library(phyloseq)
ls()
head(final_lc)
test<-final_lc[final_lc$Uprob==1, ]
dim(test)
head(test)
test<-final_lc[final_lc$SAMPLES=="S_97", ]
tail(test)
savehistory("temp.R")
load("all_16s_physeq.RData")
ls()
library(phyloseq)
data1.rich
head(data1.rich)
data1.phy
totu<-t(data.frame(otu_table(data1.phy))
)
dim(totu)
head(totu[, 1:5])
tail(data1.rich)
tail(totu[, 1:5])
library(vegan)
data1.rich$Pielou<-diversity(totu)/log(specnumber(totu))
head(data1.rich)
dim(sample_data(data1.phy))
sample_data(data1.phy)[1:5, 110:113]
sample_data(data1.phy)$Pielou<-data1.rich$Pielou
sample_data(data1.phy)[1:5, 110:114]
data2.ph
data2.phy
test<-prune_samples(grepl("CF|CNF", Name), data1.phy)
data2.phy<-subset_samples(data1.phy, grepl("CF|CNF", Name))
data2.phy<-prune_taxa(taxa_sums(data2.phy)>0, data2.phy)
data2.phy
si<-data.frame(sample_data(data2.phy))
head(si[, 1:5])
working.phy<-subset_samples(data2.phy, !is.na(Peilou))
working.phy<-subset_samples(data2.phy, !is.na(Pielou))
working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.cca<-capscale(data.trans ~ si$Pielou, dist="bray")
source("~/Documents/repos/code/R/ordisurf_extraction.R")
source("~/Documents/repos/code/R/ggplot_cca_ordisurf_2factors.R")
sf<-ordi.sf(data.cca, si[, i])
sf<-ordi.sf(data.cca, si$Pielou)
sf<-ordi.sf(data.cca, si$Foaming.Status)
data.cca<-capscale(data.trans ~ si$Pielou, dist="bray")
sf<-ordi.sf(data.cca, si$Pielou)
pdf("all_sample_nmds/env_trim_cca_cfVcnf/ordi_byFS/ordi_16s_cca_2f_byFoaming.Status_perfactor_Pielou.pdf")
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
library(RColorBrewer)
library(ggplot2)
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
p<-ggplot.cca.ordisurf.2f(data.cca, si[,c("Name", "Foaming.Status"), drop=F], colors, sf, Pielou)
p<-ggplot.cca.ordisurf.2f(data.cca, si[,c("Name", "Foaming.Status"), drop=F], colors, sf, si$Pielou)
dev.off()
pdf("all_sample_nmds/env_trim_cca_cfVcnf/ordi_byFS/ordi_16s_cca_2f_byFoaming.Status_perfactor_Pielou.pdf")
p<-ggplot.cca.ordisurf.2f(data.cca, si[,c("Name", "Foaming.Status"), drop=F], colors, sf, "Pielou")
dev.off()
plot(data.cca)
p<-ggplot.cca.ordisurf.2f(data.cca, si[,c("Name", "Foaming.Status"), drop=F], colors, sf, "Pielou")
pdf("all_sample_nmds/env_trim_cca_cfVcnf/ordi_byFS/ordi_16s_cca_2f_byFoaming.Status_perfactor_Pielou.pdf")
p
dev.off()
working.phy
data2.
data2.phy
data.cca<-capscale(data.trans ~ si$Foaming.Status, dist="bray")
data.cca
plot(data.cca)
sf<-ordi.sf(data.cca, si$Pielou)
pdf("all_sample_nmds/env_trim_cca_cfVcnf/ordi_byFS/ordi_16s_cca_2f_byFoaming.Status_perfactor_Pielou.pdf")
p<-ggplot.cca.ordisurf.2f(data.cca, si[,c("Name", "Foaming.Status"), drop=F], colors, sf, "Pielou")
print(p)
dev.off()
savehistory("temp.R")
load("all_16s_physeq.RData")
ls()
library(phyloseq)
data1.rich
head(data1.rich)
data1.phy
totu<-t(data.frame(otu_table(data1.phy))
)
dim(totu)
head(totu[, 1:5])
tail(data1.rich)
tail(totu[, 1:5])
library(vegan)
data1.rich$Pielou<-diversity(totu)/log(specnumber(totu))
head(data1.rich)
dim(sample_data(data1.phy))
sample_data(data1.phy)[1:5, 110:113]
sample_data(data1.phy)$Pielou<-data1.rich$Pielou
sample_data(data1.phy)[1:5, 110:114]
data2.ph
data2.phy
test<-prune_samples(grepl("CF|CNF", Name), data1.phy)
data2.phy<-subset_samples(data1.phy, grepl("CF|CNF", Name))
data2.phy<-prune_taxa(taxa_sums(data2.phy)>0, data2.phy)
data2.phy
si<-data.frame(sample_data(data2.phy))
head(si[, 1:5])
working.phy<-subset_samples(data2.phy, !is.na(Peilou))
working.phy<-subset_samples(data2.phy, !is.na(Pielou))
working.phy<-prune_taxa(taxa_sums(working.phy)>0, working.phy)
        otu<-data.frame(otu_table(working.phy))
si<-data.frame(sample_data(working.phy))
totu<-data.frame(t(otu))
data.trans<-decostand(totu, "total")
data.cca<-capscale(data.trans ~ si$Pielou, dist="bray")
source("~/Documents/repos/code/R/ordisurf_extraction.R")
source("~/Documents/repos/code/R/ggplot_cca_ordisurf_2factors.R")
sf<-ordi.sf(data.cca, si[, i])
sf<-ordi.sf(data.cca, si$Pielou)
sf<-ordi.sf(data.cca, si$Foaming.Status)
data.cca<-capscale(data.trans ~ si$Pielou, dist="bray")
sf<-ordi.sf(data.cca, si$Pielou)
pdf("all_sample_nmds/env_trim_cca_cfVcnf/ordi_byFS/ordi_16s_cca_2f_byFoaming.Status_perfactor_Pielou.pdf")
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
library(RColorBrewer)
library(ggplot2)
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
p<-ggplot.cca.ordisurf.2f(data.cca, si[,c("Name", "Foaming.Status"), drop=F], colors, sf, Pielou)
p<-ggplot.cca.ordisurf.2f(data.cca, si[,c("Name", "Foaming.Status"), drop=F], colors, sf, si$Pielou)
dev.off()
pdf("all_sample_nmds/env_trim_cca_cfVcnf/ordi_byFS/ordi_16s_cca_2f_byFoaming.Status_perfactor_Pielou.pdf")
p<-ggplot.cca.ordisurf.2f(data.cca, si[,c("Name", "Foaming.Status"), drop=F], colors, sf, "Pielou")
dev.off()
plot(data.cca)
p<-ggplot.cca.ordisurf.2f(data.cca, si[,c("Name", "Foaming.Status"), drop=F], colors, sf, "Pielou")
pdf("all_sample_nmds/env_trim_cca_cfVcnf/ordi_byFS/ordi_16s_cca_2f_byFoaming.Status_perfactor_Pielou.pdf")
p
dev.off()
working.phy
data2.
data2.phy
data.cca<-capscale(data.trans ~ si$Foaming.Status, dist="bray")
data.cca
plot(data.cca)
sf<-ordi.sf(data.cca, si$Pielou)
pdf("all_sample_nmds/env_trim_cca_cfVcnf/ordi_byFS/ordi_16s_cca_2f_byFoaming.Status_perfactor_Pielou.pdf")
p<-ggplot.cca.ordisurf.2f(data.cca, si[,c("Name", "Foaming.Status"), drop=F], colors, sf, "Pielou")
print(p)
dev.off()
savehistory("temp.R")
load("all_16s_physeq.RData")
ls()
library(phyloseq)
data1.phy
data2.phy
head(sample_data(data2.phy)[, 1:5])
head(sample_data(data1.phy)[, 1:5])
head(sample_data(data1.phy)[, 1:15])
unique(sample_data(data2.phy)$Treated.Status)
non_trt.phy<-prune_samples(Treated.Status=="0", data2.phy)
non_trt.phy<-subset_samples(data2.phy, Treated.Status=="0")
non_trt.phy
min(taxa_sums(non_trt.phy)
)
non_trt.phy<-prune_taxa(taxa_sums(non_trt.phy)>0, non_trt.phy)
non_trt.phy
library(plyr)
sample_names(non_trt.phy)
names(sample_data(non_trt.phy))
count(sample_data(non_trt.phy), "FK_Integrator")
plot_richness(non_trt.phy)
si<-data.frame(sample_data(non_trt.phy))
dim(si)
library(ggplot2)
si<-si[order(si$month, si$year), ]
head(si[, 1:15])
si<-si[order(si$year, si$month), ]
head(si[, 1:15])
library(grid)
 p<-ggplot(si, aes(x=month, y=Observed)) + geom_boxplot()+ facet_grid(~Foaming.Status)
p
 p<-ggplot(si, aes(group=month, y=Observed)) + geom_boxplot()+ facet_grid(~Foaming.Status)
p
 p<-ggplot(si, aes(month, Observed)) + geom_boxplot()+ facet_grid(~Foaming.Status)
p
 p<-ggplot(si, aes(month, Observed)) + geom_boxplot()
p
head(si[, c("full_date", "Observed")]
)
 p<-ggplot(si, aes(month, Observed, group=month)) + geom_boxplot()
p
 p<-ggplot(si, aes(month, Observed, group=month)) + geom_boxplot()+facet_grid(Foaming.Status~.)
p
si$group<-paste(si$month, si$year, sep="_")
 p<-ggplot(si, aes(group, Observed, group=group)) + geom_boxplot()+facet_grid(Foaming.Status~.)
p
str(si$group)
si$group<-factor(si$group, levels=levels(si$group))
str(si$group)
si$group<-paste(si$month, si$year, sep="_")
si$group<-factor(si$group, levels=unique(si$group))
str(si$group)
 p<-ggplot(si, aes(group, Observed, group=group)) + geom_boxplot()+facet_grid(Foaming.Status~.)
p
 p<-ggplot(si, aes(group, Pielou, group=group)) + geom_boxplot()+facet_grid(Foaming.Status~.)
p
 p<-ggplot(si, aes(group, Pielou, group=group, fill=Foaming.Status)) + geom_boxplot()
p
 p<-ggplot(si, aes(group, Pielou, group=group, color=Foaming.Status)) + geom_boxplot()
p
 p<-ggplot(si, aes(group, Pielou, group=group)) + geom_boxplot(aes(color=Foaming.Status)
)
p
 p<-ggplot(si, aes(group, Pielou)) + geom_boxplot(aes(color=Foaming.Status))
p
unique(si$Foaming.Status)
 p<-ggplot(si, aes(group, Pielou)) + geom_boxplot(aes(color=as.factor(Foaming.Status)))
p
 p<-ggplot(si, aes(group, Pielou)) + geom_ribbon(aes(ymin=min(Pielou), ymax=max(Pielou), alpha=0.2) + geom_line()
 p<-ggplot(si, aes(group, Pielou)) + geom_ribbon(aes(ymin=min(Pielou), ymax=max(Pielou)), alpha=0.2) + geom_line()
p
 p<-ggplot()+ geom_blank(si, aes(group, Pielou)) + geom_ribbon(aes(ymin=min(Pielou), ymax=max(Pielou)), alpha=0.2) + geom_line()
 p<-ggplot(si, aes(x=month, Pielou)) + geom_ribbon(aes(ymin=min(Pielou), ymax=max(Pielou)), alpha=0.2) + geom_line()
p
 p<-ggplot(si, aes(x=month, Pielou)) + geom_point(aes(ymin=min(Pielou), ymax=max(Pielou)), alpha=0.2)
p
si$index<-seq(1, length(si[,1]))
head(si$index)
 p<-ggplot(si, aes(index, Pielou)) + geom_ribbon(aes(ymin=min(Pielou), ymax=max(Pielou)), alpha=0.2) + geom_line()
p
si$index<-NULL
head(si[, 1:15])
library(plyr)
test<-ddply(si, .(Foaming.Status, group), avg=mean(Pielou), ste=sd(Pielou)/sqrt(length(Pielou))
)
dim(test)
head(test)
test<-ddply(si, .(Foaming.Status, group), summarize, avg=mean(Pielou), ste=sd(Pielou)/sqrt(length(Pielou))
)
head(test)
test$index<-seq(1, length(test[, 1]))
head(test)
ggplot(test, aes(index, avg)) + geom_ribbon(aes(ymin=min(avg-ste), ymax=max(avg+ste)), alpha=0.2) + geom_line(color=Foaming.Status)
ggplot(test, aes(index, avg)) + geom_ribbon(aes(ymin=min(avg-ste), ymax=max(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))
test$Foaming.Status<-factor(test$Foaming.Status, levels=c("0", "1"))
head(test)
ggplot(test, aes(index, avg)) + geom_ribbon(aes(ymin=min(avg-ste), ymax=max(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))
ggplot(test, aes(group, avg)) + geom_ribbon(aes(ymin=min(avg-ste), ymax=max(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))
test$index<-NULL
unique(test$group)
month_index<-data.frame(unique(test$group))
head(month_index)
names(month_index)<-"group"
month_index$index<-sep(1, length(month_index$group))
month_index$index<-seq(1, length(month_index$group))
head(month_index)
test<-merge(test, month_index, "group")
dim(test)
ggplot(test, aes(index, avg)) + geom_ribbon(aes(ymin=min(avg-ste), ymax=max(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))
ggplot(test, aes(index, avg)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))
ggplot(test, aes(index, avg, color=Foaming.Status)) + stat_smooth(method="loess", span=0.1, se=T, aes(fill=Foaming.Status), alpha=0.2) + teme_bw()
ggplot(test, aes(index, avg, color=Foaming.Status)) + stat_smooth(method="loess", span=0.1, se=T, aes(fill=Foaming.Status), alpha=0.2) + theme_bw()
ggplot(test, aes(index, avg, color=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))
ggplot(test, aes(index, avg)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), color=Foaming.Status, alpha=0.2) + geom_line(aes(color=Foaming.Status))
ggplot(test, aes(index, avg)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste), color=Foaming.Status), alpha=0.2) + geom_line(aes(color=Foaming.Status))
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))+scale_x_discrete(breaks=1:12, labels=unique(test$group))
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))+scale_x_discrete(breaks=1:13, labels=unique(test$group))
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))+scale_x_discrete(breaks=1:13, labels=unique(test$group))+xlab(NULL)
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))+scale_x_discrete(labels=group)
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))+scale_x_discrete(labels=test$group)
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))+scale_x_continuous(breaks=1:13, labels=test$group)
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))+scale_x_continuous(breaks=1:13, labels=unique(test$group))
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))+scale_x_continuous(breaks=1:13, labels=unique(test$group))+theme_bw()+scale_color_discrete(value=c("red", "black")+scale_fill_discrete(value="red", "black")
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))+scale_x_continuous(breaks=1:13, labels=unique(test$group))+theme_bw()+scale_color_discrete(value=c("red", "black"))+scale_fill_discrete(value=c("red", "black"))
ggplot(test, aes(index, avg, fill=Foaming.Status)) + geom_ribbon(aes(ymin=(avg-ste), ymax=(avg+ste)), alpha=0.2) + geom_line(aes(color=Foaming.Status))+scale_x_continuous(breaks=1:13, labels=unique(test$group))+theme_bw()
ddply(si, .(group), function(x) wilcox.test(x[Foaming.Status="0",], x[Foaming.Status=="1", ])
)
ddply(si, .(group), function(x) wilcox.test(x[Foaming.Status=0,], x[Foaming.Status==1, ])
)
head(si)
ddply(si, .(group), function(x) wilcox.test(x[Foaming.Status=0,], x[Foaming.Status==1, ])
)
ddply(si, .(group), function(x) wilcox.test(x[Foaming.Status==0,], x[Foaming.Status==1, ])
)
ddply(si, .(group), function(x) wilcox.test(x[si$Foaming.Status==0,], x[si$Foaming.Status==1, ])
)
ddply(si, .(group), function(x) wilcox.test(x[si$Foaming.Status=="0",], x[si$Foaming.Status=="1", ]))
ddply(idata.frame(si[, c("group", "Pielou", "Foaming.Status")]), .(group), function(x) wilcox.test(x[si$Foaming.Status=="0",], x[si$Foaming.Status=="1", ]))
ddply(idata.frame(si[, c("group", "Pielou", "Foaming.Status")]), .(group), function(x) wilcox.test(x[Foaming.Status=="0",], x[Foaming.Status=="1", ]))
ddply(idata.frame(si[, c("group", "Pielou", "Foaming.Status")]), .(group), function(x) wilcox.test(x[Foaming.Status=="0",Pielou], x[Foaming.Status=="1", Pielou]))
ddply(idata.frame(si[, c("group", "Pielou", "Foaming.Status")]), .(group), function(x) wilcox.test(x[si$Foaming.Status=="0",si$Pielou], x[si$Foaming.Status=="1", si$Pielou]))
test<-structure(list(Category = c("A", "C", 
"B", "C", "D", "E", 
"C", "A", "F", "B", 
"E", "C", "C", "A", 
"C", "A", "B", "H", 
"I", "A"), Type = c("POST", "POST", 
"POST", "POST", "PRE", "POST", "POST", "PRE", "POST", 
"POST", "POST", "POST", "POST", "PRE", "PRE", "POST", 
"POST", "POST", "POST", "POST"), Value = c(1560638113, 
1283621, 561329742, 2727503, 938032, 4233577690, 0, 4209749646, 
111467236, 174667894, 1071501854, 720499, 2195611, 1117814707, 
1181525, 1493315101, 253416809, 327012982, 538595522, 3023339026
)), .Names = c("Category", "Type", "Value"), row.names = c(21406L, 
123351L, 59875L, 45186L, 126720L, 94153L, 48067L, 159371L, 54303L, 
63318L, 104100L, 58162L, 41945L, 159794L, 57757L, 178622L, 83812L, 
130655L, 30860L, 24513L), class = "data.frame")
head(test)
ddply(idata.frame(data), .(Category), 
    function(x) wilcox.test(x[Type == "PRE",], x[Type == "POST",])
)
ddply(idata.frame(si[, c("group", "Pielou", "Foaming.Status")]), .(group), function(x) {w<-wilcox.test(subset(x, Foaming.Status=="0", select=Pielou)[,1], subset(x, Foaming.Status="1", select=Pielou)[,1])) return(w[c(1,3)])})
test<-data.frame(si[, c("group", "Pielou", "Foaming.Status")])
head(test)
test<-table(test$group)
head(test)
library(reshape)
library(reshape2)
test<-data.frame(si[, c("group", "Pielou", "Foaming.Status")])
test1<-dcast(test, value.va=Foaming.Status)
test1<-dcast(test, id.vars="Foaming.Status")
test1<-dcast(test, group+Pielou~Foaming.Status)
head(test1)
test1<-dcast(test, group~Foaming.Status)
head(test1)
test1<-dcast(test, row.names(test)+group~Foaming.Status)
head(test1)
test1<-dcast(test, row.names(test)+group~Foaming.Status, value.var=Pielou)
test1<-dcast(test, row.names(test)+group~Foaming.Status, value.var="Pielou")
head(test1)
wx <- function(d){
 w <- wilcox.test(
  # First vector (x)
    subset(d, Type == "PRE", select = Value )[,1], 
    subset(d, Type == "POST", select = Value )[,1]
      )
  # c(1,3) returns the Stat and the P-value (tweak that if you want something else)
  return(w[c(1,3)])
  }
ddply(test1, .(goup), function(x) wilcox.test(x[, c("0")], x[, c("1")]))
ddply(test1, .(rgoup), function(x) wilcox.test(x[, c("0")], x[, c("1")]))
ddply(test1, .(group), function(x) wilcox.test(x[, c("0")], x[, c("1")]))
colnames(test1)[3:4]<-c("non-foaming", "foaming")
ddply(test1, .(group), function(x) wilcox.test(x[, c("non-foaming")], x[, c("foaming")]))
ddply(test1, .(group), function(x) wilcox.test(x[, c("non-foaming")], x[, c("foaming")], na.rm=T))
str(test1)
wx <- function(d){
 w <- wilcox.test(
  # First vector (x)
    subset(d, Foaming.Status == "0", select = Pielou  )[,1], 
    subset(d, Foaming.Status == "1", select = Pielou  )[,1]
      )
  # c(1,3) returns the Stat and the P-value (tweak that if you want something else)
  return(w[c(1,3)])
  }
ddply(test, .(group), wx)
head(test)
wx(test1[, 3:4])
wx(test)
ddply(test, .(group), summarize, wx)
ddply(test, .(group), .fun= wx)
wx <- function(d){
 w <- wilcox.test(
  # First vector (x)
    subset(d, Foaming.Status == "0", select = Pielou  )[,1], 
    subset(d, Foaming.Status == "1", select = Pielou  )[,1],
    na.rm=T  )
  # c(1,3) returns the Stat and the P-value (tweak that if you want something else)
  return(w[c(1,3)])
  }
ddply(test, .(group), summarize, wx)
ddply(test, .(group), .function= wx)
subset(test, Foaming.Status == "0", select = Pielou  )[,1]
ddply(test, .(group), .fun= wx)
ddply(test, .(group), wx)
subset(test, Foaming.Status == "0", select = Pielou  )
subset(test, Foaming.Status == "0", select = Pielou  )[,1]
subset(test, Foaming.Status == "1", select = Pielou  )[,1]
head(test)
count(test, "group")
head(test1)
for (month in unique(test1$group)){
temp<-subset(test1, group==month)
print (dim(temp))
}
for (month in unique(test1$group)){
temp<-subset(test1, group==month)
p<-wilcox.test(temp$non-foaming, temp$foaming, na.rm=T)
print(month)
print(p)}
for (month in unique(test1$group)){
temp<-subset(test1, group==month)
print(month)
print(head(temp))
}
names(test1)
names(test1)<-c("sample", "group", "nonfoaming", "foaming")
for (month in unique(test1$group)){
temp<-subset(test1, group==month)
p<-wilcox.test(temp$nonfoaming, temp$foaming)
print(month)
print(p)}
ddply(test, .(group), .fun=wx)
p<-ggplot(si, aes(group, Observed, group=group)) + geom_boxplot()+facet_grid(Foaming.Status~.)
p
p<-ggplot(si, aes(group, Observed, group=group, color=as.factor(Foaming.Status)) + geom_boxplot()
)
p<-ggplot(si, aes(group, Observed, group=group, color=as.factor(Foaming.Status))) + geom_boxplot()
p
p<-ggplot(si, aes(group, Observed, group=group)) + geom_boxplot(aes(color=as.factor(Foaming.Status)))
p
p<-ggplot(si, aes(group, Observed)) + geom_boxplot(aes(color=as.factor(Foaming.Status)))
p
p<-ggplot(si, aes(group, Shannon)) + geom_boxplot(aes(color=as.factor(Foaming.Status)))
p
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
ls()
library(phyloseq)
data<-readRDS("RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/all_valid_samples_min_taxasums_5_min_seq_10k_physeq.RDS")
data
library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)
si<-data.frame(sample_data(data))
dim(si)
count(si, "id")
count(si, c("id", "Foaming.Status")
)
physeq<-subset_samples(data, id=="SP")
physeq
physeq<-prune_taxa(taxa_sums(physeq)>0, physeq)
physeq
physeq<-tax_glom(physeq, "genus")
physeq<-subset_taxa(physeq, domain!="Archaea" & domain!="unclassified_Root")
otu<-data.frame(otu_table(physeq))
si<-data.frame(sample_data(physeq))
tax<-data.frame(tax_table(physeq))
row.names(otu)<-tax$genus
totu<-data.frame(t(otu))
dataset<-merge(si, totu, by.x="SAMPLES", by.y="row.names")
print(dim(dataset))
head(dataset[, 1:9])
head(dataset[, 1:20])
dataset<-merge(si, totu, by.x="SAMPLES", by.y="row.names")
print(dim(dataset))
results_sp<-rcorr(as.matrix(dataset[,-c(1:20)]),type="spearman")
results_sp<-rcorr(as.matrix(dataset[,-c(1:20)]),type="spearman", na.omit=T)
str(dataset[, 21:30])
str(dataset[, 21:410])
str(dataset[, 400:410])
str(dataset[, 100: 117])
str(dataset[, 107: 116])
results_sp<-rcorr(as.matrix(dataset[,-c(1:20, 107:116)]),type="spearman", na.omit=T)
results_sp<-rcorr(as.matrix(dataset[,-c(1:20, 107:116)]),type="spearman")
dim(results_sp)
results_sp
results_hd<-hoeffd(as.matrix(dataset[,-c(1:20, 107:116)])
)
 rhos<-results_sp$r
        sp_ps<-results_sp$P
        ds<-results_hd$D
        ds_ps<-results_hd$P
       sp_melt<-na.omit(melt(sp_ps))
        ds_melt<-na.omit(melt(ds_ps))
        #creating a qvalue (adjusted pvalue) based on FDR
        sp_melt$spearman_qval<-fdrtool(sp_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        ds_melt$hoeffding_qval<-fdrtool(ds_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        #making column names more relevant
        names(sp_melt)[3]<-"spearman_pval"
        names(ds_melt)[3]<-"hoeffding_pval"
        # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
        sp_sub<-subset(sp_melt, spearman_qval < 0.05)
        ds_sub<-subset(ds_melt, hoeffding_qval < 0.05)
        # now melting the rhos, note the similarity between ps_melt and rhos_melt
        rhos_melt<-na.omit(melt(rhos))
        ds_melt<-na.omit(melt(ds))
        names(rhos_melt)[3]<-"rho"
        names(ds_melt)[3]<-"D"
       sp_merged<-merge(sp_sub,rhos_melt,by=c("Var1","Var2"))
        ds_merged<-merge(ds_sub, ds_melt,by=c("Var1","Var2"))
        merged<-merge(sp_merged, ds_merged, by=c("Var1", "Var2")
)
head(merge)
head(merged)
dim(sp_merged)
dim(ds_merged)
tail(merged)
length(unique(merged$rho))
      sp_melt$spearman_qval<-p.adjust(sp_melt$value, "fdr")
        ds_melt$hoeffding_qval<-p.adjust(ds_melt$value, "fdr")
#       sp_melt$spearman_qval<-fdrtool(sp_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
#       ds_melt$hoeffding_qval<-fdrtool(ds_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        #making column names more relevant
        names(sp_melt)[3]<-"spearman_pval"
        names(ds_melt)[3]<-"hoeffding_pval"
        # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
        sp_sub<-subset(sp_melt, spearman_qval < 0.05) 
        ds_sub<-subset(ds_melt, hoeffding_qval < 0.05)
        # now melting the rhos, note the similarity between ps_melt and rhos_melt
        rhos_melt<-na.omit(melt(rhos))
        ds_melt<-na.omit(melt(ds))
head(sp_me)
head(sp_melt)
       sp_melt<-na.omit(melt(sp_ps))
        ds_melt<-na.omit(melt(ds_ps))
head(sp_melt)
        sp_melt$spearman_qval<-p.adjust(sp_melt$value, "fdr")
        ds_melt$hoeffding_qval<-p.adjust(ds_melt$value, "fdr")
        names(sp_melt)[3]<-"spearman_pval"
        names(ds_melt)[3]<-"hoeffding_pval"
        # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
        sp_sub<-subset(sp_melt, spearman_qval < 0.05)
        ds_sub<-subset(ds_melt, hoeffding_qval < 0.05)
        # now melting the rhos, note the similarity between ps_melt and rhos_melt
        rhos_melt<-na.omit(melt(rhos))
        ds_melt<-na.omit(melt(ds))
        names(rhos_melt)[3]<-"rho"
        names(ds_melt)[3]<-"D"
        #merging together and remove negative rhos
        sp_merged<-merge(sp_sub,rhos_melt,by=c("Var1","Var2"))
        ds_merged<-merge(ds_sub, ds_melt,by=c("Var1","Var2"))
        merged<-merge(sp_merged, ds_merged, by=c("Var1", "Var2"))
head(merged)
head(sp_sub)
max(sp_sub$spearman_qval)
for (i in si$id){
x=i
print(x)}
for (i in unique(si$id)){
x=i
print(x)}
for (i in unique(data.frame(sample_data(data))$id)){
print(i)}
ls()
head(merged)
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
ls()
combined_cc<-read.delim("test/all_valid_samples_min_taxasums_5_min_seq_10k_physeq_combined_barn_cc_results.txt", sep="\t", header=T)
dim(combined_cc)
head(combined_cc)
min(combined_cc$D)
length(combined_cc$D[combined_cc$D<0.5])
length(combined_cc$D[combined_cc$D<0.6])
h<-hist(combined_cc$rho)
library(ggplot2)
ggplot(combined_cc, aes(rho)+geom_histogram(bin=100)
ggplot(combined_cc, aes(rho))+geom_histogram(bin=100)
ggplot(combined_cc, aes(rho))+geom_histogram()
ggplot(combined_cc, aes(rho))+geom_histogram(binwidth=100)
ggplot(combined_cc, aes(rho))+geom_histogram(binwidth=10)
ggplot(combined_cc, aes(rho))+geom_histogram()
ggplot(combined_cc, aes(rho))+geom_histogram(bins=100)
ggplot(combined_cc, aes(D))+geom_histogram(bins=100)
plot(combined_cc$rho ~ combined_cc$D)
plot(abs(combined_cc$rho) ~ combined_cc$D)
plot(combined_cc$D ~ abs(combined_cc$rho) )
ls()
D0.6<-subset(combined_cc, D>=0.6)
dim(D0.6)
dim(combined_cc)
6960*100/28798
ggplot(combined_cc, aes(D))+geom_histogram(bins=100)
min(D0.6$rho)
min(abs(D0.6$rho))
D0.7<-subset(combined_cc, D>=0.7)
dim(D0.7)
4656*100/28798
plot(combined_cc$D ~ abs(combined_cc$rho) )
plot(D0.7$D ~ abs(D0.7$rho) )
ls()
meta<-read.delim("foaming_status_cc/meta_w_genus_information.txt", sep="\t", header=T)
head(meta)
head(combined_cc)
temp<-merge(final_results, meta[, c("genus", "domain")], by.x="Var1", by.y="genus")
temp<-merge(combined_cc, meta[, c("genus", "domain")], by.x="Var1", by.y="genus")
head(temp)
temp<-merge(temp, meta[, c("genus", "domain")], by.x="Var2", by.y="genus")
head(temp)
bac.bac<-subset(temp, domain.x=="Bacteria" & domain.y=="Bacteria")[, 1:6]
head(bac.bac)
bac.bac<-subset(temp, domain.x=="Bacteria" & domain.y=="Bacteria")[, 1:8]
head(bac.bac)
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
test<-subset_taxa(phyloseq, taxa_names="Garciella")
test<-subset_taxa(data, taxa_names="Garciella")
test
data
test<-subset_taxa(data, genus="Garciella")
test
tax_table(data)
head(tax_table(data))
test<-subset_taxa(data, genus="Garciella")
test
data
dim(tax)
head(tax)
test<-prune_taxa(genus="Garciella", data)
test<-prune_taxa(data, genus=="Garciella")
data
head(tax_table(data))
test<-subset_taxa(data, genus=="Garciella")
test
test<-prune_samples(sample_sums(test)>0, test)
test
tax_table(test)
otu_table(test)[, 1:5]
otu_table(test)[, 10:30]
test<-tax_glom(test, taxrank="genus")
test
garciella<-t(data.frame(otu_table(test)))
head(garciella)
si<-data.frame(sample_data(test))
names(si)
ls9)
ls()
garciella<-merge(garciella, si[, c("SAMPLES", "id", "Crude.Fiber")], by.x="row.names", by.y="SAMPLES")
head(garciella)
plot(garciella$Crude.Fiber ~ garciella$OTU_3904)
test<-subset_taxa(data, kingdom=="Archaea")
test<-subset_taxa(data, domain=="Archaea")
test
test<-prune_samples(sample_sums(test)>0, test)
test
tax_table(test)
head(otu_table(test))
dim(si)
names(si)
head(si$FK_DNA.ID)
si$SAMPLES[si$FK_DNA.ID == "1973"0
si$SAMPLES[si$FK_DNA.ID == "1973")
si$SAMPLES[si$FK_DNA.ID == "1973"]
head(si$Long_Sample_ID)]
head(si$Long_Sample_ID)
si$SAMPLES[si$Long_Sample_ID == "IA-2-080813-F104-1B"]
head(si$layer)
head(si$fullmonth)
head(si$myear)
mcra_si<-read.delim("mcra_si.txt", header=T, sep="\t")
dim(mcra_si)
head(mcra_si)
mcra_si<-read.delim("mcra_si.txt", header=T, sep="\t")
head(mcra_si)
head(si$Long_Sample_ID)
si$to_match<-paste(si$state, si$some_number, si$full_date, si$id, sep="_")
head(si$to_match)
dim(mcra_si)
mcra_si<-merge(mcra_si, si, "to_match")
dim(mcra_si)
mcra_si<-read.delim("mcra_si.txt", header=T, sep="\t")
test<-data.frame(do.call('rbind', strsplit(as.character(si$Long_Sample_ID), "-", fixed=T)))
head(test)
si$to_match<-paste(test$X1, test$X2, test$X3, test$X4, sep="-")
mcra_si<-merge(mcra_si, si, "to_match")
dim(mcra_si)
unique(si$layer)
count(si, "layer")
library(plyr)
count(si, "layer")
si_no_2b<-subset(si, layer != "2B")
dim(si_no_2b)
mcra_si<-read.delim("mcra_si.txt", header=T, sep="\t")
mcra_si<-merge(mcra_si, si_no_2b, "to_match")
dim(mcra_si)
mcra_si<-read.delim("mcra_si.txt", header=T, sep="\t")
mcra_si<-merge(mcra_si, si_no_2b, "to_match", all.x=T)
dim(mcra_si)
head(mcra_si)
test<-mcra_si[is.na(mcra_si$Long_Sample_ID), ]
test
test$to_match
si$Long_Sample_ID[si$to_match %in% test$to_match]
si$Long_Sample_ID[si$id == "F458"]
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)
ls()
data
data1
data1.phy
physeq<-tax_glom(data, "phylum")
physeq<-subset_taxa(physeq, domain!="Archaea" & domain!="unclassified_Root")
physeq
min(sample_sums(physeq))
otu<-data.frame(otu_table(physeq_sub))
        si<-data.frame(sample_data(physeq_sub))
        tax<-data.frame(tax_table(physeq_sub))
        row.names(otu)<-tax$genus
otu<-data.frame(otu_table(physeq)
)
si<-data.frame(sample_data(physeq))
tax<-data.frame(tax_table(physeq))
row.names(otu)<-tax$genus
        totu<-data.frame(t(otu))
row.names(otu)<-tax$phylum
        totu<-data.frame(t(otu))
head(totu[, 1:5])
colSums(totu)
totu.melt<-melt(totu, id.var="row.names")
totu.melt<-melt(totu, id="row.names")
totu.melt<-melt(totu, row.names)
totu.melt<-melt(totu)
head(totu.melt)
totu$SAMPLES<-row.names(totu)
head(totu)
totu.melt<-melt(totu, id="SAMPLES")
head(totu.melt)
colnames(totu.melt)[2:3]<-c("phylum", "count")
names(si)
test<-merge(totu.melt, si[ , c("SAMLES", "DDGS")], "SAMPLES")
test<-merge(totu.melt, si[ , c("SAMPLES", "DDGS")], "SAMPLES")
head(test)
library(ggplo2)
library(ggplot2)
ggplot(test, aes(x=DDGS, y=count, color=phylum)) + geom_points()
ggplot(test, aes(x=DDGS, y=count, color=phylum)) + geom_point()
ggplot(test, aes(x=DDGS, y=count, color=phylum)) + geom_point()
test<-merge(totu.melt, si[ , c("SAMPLES", "Crude.Fiber")], "SAMPLES")
ggplot(test, aes(x=DDGS, y=count, color=phylum)) + geom_point()
ggplot(test, aes(x=Crude.Fiber, y=count, color=phylum)) + geom_point()
test<-merge(totu.melt, si[ , c("SAMPLES", "TDF")], "SAMPLES")
ggplot(test, aes(x=TDF, y=count, color=phylum)) + geom_point()
test<-merge(totu.melt, si[ , c("SAMPLES", "Foaming.Capacity")], "SAMPLES")
ggplot(test, aes(x=Foaming.Capacity, y=count, color=phylum)) + geom_point()
data_1e8<-transform_sample_counts(data, function(x) 1e8*x/sum(x))
head(sample_sums(data_1e8))
physeq<-tax_glom(data_1e8, "phylum")
physeq<-subset_taxa(physeq, domain!="Archaea" & domain!="unclassified_Root")
min(sample_sums(physeq))
otu<-data.frame(otu_table(physeq))
si<-data.frame(sample_data(physeq))
tax<-data.frame(tax_table(physeq))
row.names(otu)<-tax$phylum
totu<-data.frame(t(otu))
head(totu)
totu$SAMPLES<-row.names(totu)
totu.melt<-melt(totu, id="SAMPLES")
head(totu.melt)
names(totu.melt)[2:3]<-c("phylum", "counts")
head(totu.melt)
test<-merge(totu.melt, si[ , c("SAMPLES", "Foaming.Capacity")], "SAMPLES")
ggplot(test, aes(x=Foaming.Capacity, y=count, color=phylum)) + geom_point()
ggplot(test, aes(x=Foaming.Capacity, y=counts, color=phylum)) + geom_point()
ggplot(test, aes(x=Foaming.Capacity, y=counts, color=phylum)) + geom_point()+facet_grid(phylum~.)
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
ls()
rm(data_1e8)
data_1e5<-readRDS("RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/all_valid_samples_min_taxasums_5_min_seq_10k_1e5standardized_physeq.RDS")
library(phyloseq)
data1_1e5
data_1e5
test<-subset_taxa(data_1e5, genus == "unclassified_Ruminococcaceae")
test
min(sample_sums(test))
head(otu_table(test)[, 1:5])
plot_heatmap(test, sample.label="Name")
names(sample_data(test))
plot_heatmap(sort(taxa_sums(gpt),TRUE), sample.label="Name")
plot_heatmap(sort(taxa_sums(test),TRUE), sample.label="Name")
test <- prune_taxa(names(sort(taxa_sums(test),TRUE), test)
)
test <- prune_taxa(names(sort(taxa_sums(test),TRUE)[1:676], test)
)
test <- prune_taxa(names(sort(taxa_sums(test),TRUE)), test)
plot_heatmap(test, sample.label="Name")
plot_heatmap(test, sample.label="myear")
totu<-t(data.frame(otu_table(test))
)
si<-data.frame(sample_data(test))
head(totu[, 1:5])
test1<-merge(totu, si, by.x="row.names", by.y="SAMPLES")
head(test1[, 1:5])
plot(test1$OTU_1001 ~ test1$DDGS)
plot(test1$OTU_10064 ~ test1$DDGS)
barns_si<-read.delim("foaming_status_cc/spearman_hoeffding_cc/barn_foaming_rate.txt", sep="\t", header=T)
test1<-merge(test1, barns_si[, c("id", "category")], "id")
dim(test1)
plot(test1$OTU_10064 ~ test1$category)
plot(lm(test1$OTU_10064 ~ test1$category))
plot(lm(log(test1$OTU_10064) ~ test1$category))
plot(lm(log(test1$OTU_10064+1) ~ test1$category))
dim(si)
si<-merge(si, barns_si[, c("id", "category")], "id")
names(si)
dim(si)
sample_data(test)<-si
head(si[, 1:4])
row.names(si)<-"SAMPLES"
dim(si)
test
head(row.names(si))
class(si)
rownames(si)<-"SAMPLES"
rownames(si)<-si$SAMPLES
sample_data(test)<-si
TEST
test
test1 <- prune_taxa(names(sort(taxa_sums(test),TRUE)[1:20]), test)
plot_heatmap(test1, sample.label="category")
head(tax_table(test))
test.psmelt<-psmelt(otu_table(test, taxa_are_rows=T), sample_data(test))
test.psmelt<-psmelt(phyloseq(otu_table(test, taxa_are_rows=T), sample_data(test)))
head(test.psmelt[, 1:5])
library(plyr)
test1<-ddply(test.psmelt, .(OTU, category), summarize, abun=sum(Abundance))
head(test1)
test1<-ddply(test.psmelt, .(OTU, category), summarize, abun=mean(Abundance))
head(test1)
dim(test1)
ggplot(test1, aes(category, OTU)) + geom_raster(aes(fill=density))
library(ggplot2)
ggplot(test1, aes(category, OTU)) + geom_raster(aes(fill=density))
ggplot(test1, aes(category, OTU)) + geom_raster(aes(fill=abun))
test1 <- prune_taxa(names(sort(taxa_sums(test),TRUE)[1:20]), test)
test1
min(sample_sums(test1))
plot_heatmap(test1, "category")
test.psmelt<-psmelt(phyloseq(otu_table(test1, taxa_are_rows=T), sample_data(test1)))
head(test.psmetl[, 1:4])
head(test.psmelt[, 1:4])
test1<-ddply(test.psmelt, .(OTU, category), summarize, abun=mean(Abundance))
ggplot(test1, aes(category, OTU)) + geom_raster(aes(fill=abun))
test1$category<-factor(test1$category, levels=c("T0", "T0-10", "T10-50", "T50-75", "T75-100", "T100"))
test<-test.psmelt[test.psmelt$category=="T100", ]
head(test[, 1:5])
test$OTU<-reorder(test$OTU, -test$Abundance)
str(test$OTU)
test1$OTU<-factor(test1$OTU, levels=levels(test$OTU))
ggplot(test1, aes(category, OTU)) + geom_raster(aes(fill=abun))
data
data_1e5
ls()
si<-data.frame(sample_data(data))
test
test1
si<-data.frame(sample_data(data))
si<-merge(si, barns_si[, c("id", "category")], "id")
row.names(si)<-si$SAMPLES
sample_data(data)<-si
data
rumino<-subset_taxa(data_1e5, genus == "unclassified_Ruminococcaceae")
rumino
test<-data.frame(taxa_sums(rumino))
head(test)
test<-data.frame(sample_sums(rumino))
head(test)
test<-merge(test, si, by.x="row.names", by.y="SAMPLES")
plot(test$sample_sums.rumino.~test$category)
anova(lm(test$sample_sums.rumino.~test$category))
plot(lm(test$sample_sums.rumino.~test$category))
plot(lm(log(test$sample_sums.rumino.+1)~test$category))
anova(lm(log(test$sample_sums.rumino.+1)~test$category))
plot(test$sample_sums.rumino.~test$DDGS)
plot(test$sample_sums.rumino.~test$TDF)
plot(test$sample_sums.rumino.~test$Crude.Fiber)
plot(test$sample_sums.rumino.~ test$Adj.MPR)
plot(test$sample_sums.rumino.~ test$Acetic.Acid)
t0_t100<-subset_samples(data, category=="T0" & category=="T100")
data
t0_t100<-subset_samples(data, category %in% c("T0", "T100"))
t0_t100
min(taxa_sums(t0_t100))
t0_t100<-prune_taxa(taxa_sums(t0_t100)>0, t0_t100)
t0_t100
t0_t100.dds<-phyloseq_to_deseq2(t0_t100, ~ category)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(t0_t100.dds), 1, gm_mean)
library(plyr)
geoMeans = apply(counts(t0_t100.dds), 1, gm_mean)
library(deseq2)
library(DESeq2)
geoMeans = apply(counts(t0_t100.dds), 1, gm_mean)
t0_t100.dds<-estimateSizeFactors(t0_t100.dds, geoMeans=geoMeans)
t0_t100.dds<-DESeq(t0_t100.dds, fitType="local")
res<-results(t0_t100.dds)
res = res[order(res$padj, na.last=NA), ]
head(res)
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = merge(as(sigtab, "data.frame"), as(tax_table(t0_t100)[rownames(t0_t100), ], "matrix"))
head(sigtab)
sigtab<-as.data.frame(sigtab)
head(sigtab)
tax<-data.frame(tax_table(t0_t100))
head(tax)
sigtab<-merge(sigtab, tax, by="row.names")
head(sigtab)
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
head(posigtab)
dim(posigtab)
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
dim(x)
x
x = tapply(posigtab$log2FoldChange, posigtab$phylum, function(x) max(x))
x
posigtab$genus
unique(posigtab$genus)
plotMA(res)
t100_genera_names<-unique(posigtab$genus)
t100_genera_names
write(t100_genera_names, "positive_sig_2fold_abunndant_t100vst0.txt", quote=F)
write(t100_genera_names, "positive_sig_2fold_abunndant_t100vst0.txt")
t100_genera_names<-droplevels.factor(t100_genera_names)
str(t100_genera_names)
head(t100_genera_names)
t100_genera_names<-droplevels(t100_genera_names)
head(t100_genera_names)
saveRDS(t100_genera_names, "positive_sig_2fold_abunndant_t100vst0.RDS")
t100_genera_names
ls()
data
dim(si)
head(si$category)
head(barns_si)
barns_si$group<-barns_si$category
barns_si$group<-gsub("T10-50", "NCF", barns_si$group)
head(barns_si)
barns_si$group<-gsub("T0", "NCF", barns_si$group)
barns_si$group<-gsub("T0-10", "NCF", barns_si$group)
barns_si$group<-gsub("50-75", "NCF", barns_si$group)
barns_si$group<-gsub("75-100", "NCF", barns_si$group)
head(barns_si)
barns_si$group<-gsub("TNCF", "NCF", barns_si$group)
head(barns_si)
unique(barns_si$group)
barns_si$group<-gsub("NCF-10", "NCF", barns_si$group)
unique(barns_si$group)
barns_si$group<-gsub("T100", "CF", barns_si$group)
unique(barns_si$group)
si<-merge(si, barns_si[, c("id", "group")], "id")
row.names(si)<-si$SAMPLES
sample_data(data)<-si
data
diagdds = phyloseq_to_deseq2(data, ~ group)
names(sample_data(data))
head(si$group.x)
head(si$group.y)
diagdds = phyloseq_to_deseq2(data, ~ group.y)
str(si$group.y)
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
res = results(diagdds)
head(res)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab<-data.frame(sigtab)
dim(si)
sigtab<-merge(sigtab, tax, "row.name")
head(tax)
head(sigtab)
sigtab<-merge(sigtab, tax, "row.names")
head(sigtab)
plotMA(res)
cfsigtab = sigtab[sigtab[, "log2FoldChange"] < 0, ]
head(cfsigtabs)
head(cfsigtab)
length(unique(cfsigtab$genus))
(unique(cfsigtab$genus))
t100_2fold_sig_genera<-as.character(unique(cfsigtab$genus))
t100_2fold_sig_genera
saveRDS(t100_2fold_sig_genera, "sig_2fold_abunndant_t100vsAllOthers.RDS")
ls()
rm(data1.phy)
rm(data2.phy)
t100_2fold_sig.phy<-subset_taxa(data, genus %in% t100_2fold_sig_genera)
t100_2fold_sig.phy
min(sample_sums(t100_2fold_sig.phy))
plot_bar(t100_2fold_sig.phy, fill="phylum")
plot_bar(t100_2fold_sig.phy, x="category", fill="phylum")
plot_bar(t100_2fold_sig.phy, x="category", fill="phylum")+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")
t100_2fold_sig_genera<-as.character(unique(cfsigtab$genus))
t100_2fold_sig_otu<-as.character(unique(cfsigtab$Row.names))
length(t100_2fold_sig_otu)
t100_2fold_sig.phy<-subset_taxa(data, taxa_names %in% t100_2fold_sig_genera)
t100_2fold_sig.phy<-subset_taxa(data, taxa_names %in% t100_2fold_sig_otu)
head(taxa_names(data))
t100_2fold_sig.phy<-subset_taxa(data_1e5, taxa_names(data_1e5) %in% t100_2fold_sig_otu)
t100_2fold_sig.phy
min(sample_sums(t100_2fold_sig.phy))
plot_bar(t100_2fold_sig.phy, x="category", fill="phylum")+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")
dim(si)
data_1e5
sample_data(data_1e5)<-si
t100_2fold_sig.phy<-subset_taxa(data_1e5, taxa_names(data_1e5) %in% t100_2fold_sig_otu)
plot_bar(t100_2fold_sig.phy, x="category", fill="phylum")+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")
plot_bar
plot_bar(t100_2fold_sig.phy, x="category", fill="phylum", facet_grid=phylum~.)+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")
plot_bar(t100_2fold_sig.phy, x="category", fill="phylum", facet_grid=phylum~.)+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(scales="free")
plot_bar(t100_2fold_sig.phy, x="category", fill="phylum")+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~., scales="free")
str(sample_data(t100_2fold_sig.phy)$category)
sample_data(t100_2fold_sig.phy)$category<-factor(sample_data(t100_2fold_sig.phy)$category, levels=c("T0", "T0-10", "T10-50", "T50-75
sample_data(t100_2fold_sig.phy)$category<-factor(sample_data(t100_2fold_sig.phy)$category, levels=c("T0", "T0-10", "T10-50", "T50-75", "T75-T100", "T100"))
plot_bar(t100_2fold_sig.phy, x="category")+ geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+facet_grid(phylum~., scales="free")
plot_bar(t100_2fold_sig.phy, x="category", fill="genus")+ geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+facet_grid(phylum~., scales="free")
sample_data(data_1e5)$category<-factor(sample_data(data_1e5)$category, levels=c("T0", "T0-10", "T10-50", "T50-75", "T75-100", "T100"))
t100_2fold_sig.phy<-subset_taxa(data_1e5, taxa_names(data_1e5) %in% t100_2fold_sig_otu)
plot_bar(t100_2fold_sig.phy, x="category", fill="genus")+ geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+facet_grid(phylum~., scales="free")
plot_heatmap(t100_2fold_sig.phy, "category")
plot_heatmap
data("GlobalPatterns")
gpt <- subset_taxa(GlobalPatterns, Kingdom=="Bacteria")
gpt
names(sample_data(gpt)
)
head(sample_data(gpt)$SampleType)
rm(gpt)
rm(GlobalPatterns)
head(sample_data(t100_2fold_sig.phy)$category)
plot_heatmap(t100_2fold_sig.phy, sample.label="category")
totu<-t(data.frame(t100_2fold_sig.phy))
totu<-t(data.frame(otu_table(t100_2fold_sig.phy)))
si<-data.frame(sample_data(t100_2fold_sig.phy))
tax<-data.frame(tax_table(t100_2fold_sig.phy))
head(totu[, 1:4])
head(si[, 1:4])
plot(totu$OTU_1006 ~ si$Carbon)
class(totu)
plot(totu[, 1] ~ si$Carbon)
plot(totu[, 1] ~ si$TKN)
plot(totu[, 1] ~ si$Adj.MPR)
plot(si$Carbon ~ si$Adj.M)
plot_heatmap(t100_2fold_sig.phy, sample.label="category")
plot(totu[, 1] ~ si$Acetic.Acid)
head(tax)
plot(totu[, 2] ~ si$Acetic.Acid)
plot(totu[, 3] ~ si$Acetic.Acid)
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
ls()
library(phyloseq)
data_1e5
sample_names(data_1e5)
ls()
sample.variables(data_1e5)
dim(si)
data_1e5
library(ggplot2)
ggplot(si, aes(Observed, fill=Foaming.Status))+geom_boxplot
ggplot(si, aes(Observed, fill=Foaming.Status))+geom_boxplot()
ggplot(si, aes(x=myear, y=Observed, fill=Foaming.Status))+geom_boxplot()
ggplot(si, aes(x=myear, y=Observed, fill=as.factor(Foaming.Status)))+geom_boxplot()
ggplot(si, aes(x=id, y=Observed, fill=as.factor(Foaming.Status)))+geom_boxplot()
library(reshape2)
alpha
alpha<-melt(si[, c("SAMPLES", "Foaming.Status", "Observed", "Pielou")], id.vars="SAMPLES")
head(alpha)
alpha<-melt(si[, c("SAMPLES", "Foaming.Status", "Observed", "Pielou")], id.vars=c("SAMPLES", "Foaming.Status"))
head(alpha)
ggplot(alpha, aes(x=Foaming.Status, y=variable, fill=Foaming.Status)) + geom_boxplot()
ggplot(alpha, aes(x=Foaming.Status, y=variable, fill=as.factor(Foaming.Status))) + geom_boxplot()
ggplot(alpha, aes(x=Foaming.Status, y=variable, fill=as.factor(Foaming.Status))) + geom_boxplot()+facet_grid(variable~., scale="free")
head(alpha)
ggplot(alpha, aes(x=Foaming.Status, y=value, fill=as.factor(Foaming.Status))) + geom_boxplot()+facet_grid(variable~., scale="free")
ggplot(alpha, aes(x=Foaming.Status, y=value, fill=as.factor(Foaming.Status))) + geom_boxplot()+facet_grid(variable~., scale="free")+theme_bw()
ggplot(alpha, aes(value)) + geom_histogram(bins=10)+facet_grid(Foaming.Status ~ variable, scale="free")
ggplot(alpha, aes(value)) + geom_histogram(bins=30)+facet_grid(Foaming.Status ~ variable, scale="free")
wx <- function(d){
 w <- wilcox.test(
  # First vector (x)
    subset(d, Foaming.Status == "0", select = Pielou  )[,1], 
    subset(d, Foaming.Status == "1", select = Pielou  )[,1],
    na.rm=T  )
  # c(1,3) returns the Stat and the P-value (tweak that if you want something else)
  return(w[c(1,3)])
  }
wx(alpha)
wilcox.test(si$Observed[si$Foaming.Status=="0"], si$Observed[si$Foaming.Status=="1"]
)
head(si$Observed[si$Foaming.Status=="0"])
head(si$Observed)
length(si$Observed[si$Foaming.Status=="0"])
wilcox.test(si$Pielou[si$Foaming.Status=="0"], si$Pielou[si$Foaming.Status=="1"]
)
ls()
wilcox.test(si$Adj.MPR[si$Foaming.Status=="0"], si$Adj.MPR[si$Foaming.Status=="1"])
hist(si$Adj.MPR[si$Foaming.Status=="0"])
hist(si$Adj.MPR[si$Foaming.Status=="0"], "50")
hist(si$Adj.MPR[si$Foaming.Status=="0"], bins=10)
hist(si$[si$Foaming.Status=="1"], bins=10)
wilcox.test(si$Foam.Stability[si$Foaming.Status=="0"], si$Foam.Stability[si$Foaming.Status=="1"])
ls()
t100_2fold_sig_otu
length(t100_2fold_sig_otu)
write(t100_2fold_sig_otu, "t100_2fold_sig_otu_list.txt", quote=F)
write(t100_2fold_sig_otu, "t100_2fold_sig_otu_list.txt")
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
ls()
head(t100_genera_names)
rm(t100_genera_names)
rm(alpha)
rm(cfsigtab)
rm(D0.6)
rm(D0.7)
rm(dataset)
rm(t0_t100)
rm(t0_t100.dds)
library(phyloseq)
rm(t100_2fold_sig_genera)
head(totu.melt)
rm(totu.melt)
rm(test)
rm(test1)
head(wx)
rm(wx)
head(x)
rm(x)
ls()
ls(grepl("ds")
)
grepl("ds", ls)
grepl("ds", ls())
ls()[grepl("ds", ls())]
rm(ls()[grepl("ds", ls())])
rm(as.character(ls()[grepl("ds", ls())]))
rm(list=ls()[grepl("ds", ls())])
ls
ls()
rm(list=c('garciella', 'geoMeans', 'gm_mean', 'h', 'i', 'month', 'month_index', 'p', 'posigtab', 'res', 'results_hd', 'results_sp', 'rhos', 'rhos_melt', 'rumino'))
head(si_no_2b)
si(si_no_2b)
dim(si_no_2b)
dim(si
)
rm(si_no_2b)
ls()
rm(list=ls()[grepl("sp_", ls())])
ls()
t100_2fold_sig
t100_2fold_sig.phy
rm(test.psmelt)
rm(temp)
ls()
rm(bac.bac)
physeq
rm(physeq)
data
data_1e5
dim(si)
head(combined_cc)
rm(combined_cc)
ls
ls()
head(merged)
rm(merged)
head(env)
dim(env)
dim(si)
rm(env)
ls()
otu[1:5, 1:5])
otu[1:5, 1:5]
non_trt.phy
rm(non_trt.phy)
ls()
foam.height<-read.delim("foam.depth.txt", sep="\t", header=T)
head(foam.height)
test<-data.frame(do.call('rbind', strsplit(as.character(foam.height$Long_Sample_ID), "-", fixed=T)))
head(test)
test$X6<-data.frame(do.call('rbind', strsplit(as.character(test$X5), "", fixed=T)[, 1]
test$X6<-data.frame(do.call('rbind', strsplit(as.character(test$X5), "", fixed=T)))[, 1]
head(test)
dim(test)
test$new_long_id<-paste(test$X1, test$X2, test$X3, test$X4, test$X6, sep="-")
dim(foam.height)
foam.height<-cbind(foam.height, test)
head(foam.height)
head(si$Long_Sample_ID)
test<-data.frame(do.call('rbind', strsplit(as.character(si$Long_Sample_ID), "-", fixed=T)))
head(test)
names(si)
test<-data.frame(do.call('rbind', strsplit(as.character(si$layer), "-", fixed=T)))
head(test)
test<-data.frame(do.call('rbind', strsplit(as.character(si$layer), "", fixed=T)))
head(test)
si$sub_barn<-data.frame(do.call('rbind', strsplit(as.character(si$layer), "", fixed=T)))[, 1]
names(si)
si$new_long_id<-paste(si$state, si$some_number, si$full_date, si$id, si$sub_barn, sep="-")
head(si$new_long_id)
head(foam.height$new_long_id)
head(foam.height)
foam.height$X3<-as.numeric(as.character(foam.height$X3))
head(foam.height)
foam.height$new_long_id<-paste(foam.height$X1, foam.height$X2, foam.height$X3, foam.height$X4, foam.height$X6, sep="-")
head(foam.height)
head(foam.height$new_long_id)
head(si$new_long_id)
length(unique(si$new_long_id)))
length(unique(si$new_long_id))
length(unique(foam.height$new_long_id))
length(foam.height$new_long_id))
length(foam.height$new_long_id)
length((foam.height$new_long_id))
length((si$new_long_id))
length((si$Long_Sample_ID))
length(unique(si$Long_Sample_ID))
si$Long_Sample_ID[duplicated(si$Long_Sample_ID)]
si[si$Long_Sample_ID == "IL-3-031513-PTB1-1B", ]
data_1e5
data
totu<-t(data.frame(otu_table(data))
)
dim(totu)
totu[1:5, 1:5]
rowSums(totu["S_373", ])
rowSums(totu[c("S_373", "S_450"), ])
si_duprmed<-si[-"S_450", ]
si_duprmed<-si[-c("S_450"), ]
si_duprmed<-si[rownames(si) != "S_450", ]
dim(si_duprmed)
dim(si)
si_duprmed$new_long_id[!si_duprmed$new_long_id %in% foam.height$new_long_id]
si_duprmed$Long_Sample_ID[grepl("PTB5", si_duprmed$Long_Sample_ID)]
foam.height$Long_Sample_ID[grepl("PTB5", foam.height$Long_Sample_ID)]
head(foam.height)
names(si_duprmed)
dim(si_duprmed)
si_duprmed<-merge(si_duprmed, foam.height[, c(2:4, 11)])
dim(si_duprmed)
names(si_duprmed)
si_duprmed<-si[rownames(si) != "S_450", ]
names(si_duprmed)
si_duprmed<-merge(si_duprmed, foam.height[, c(2:4, 11)], "new_long_id", by.x=T)
si_duprmed<-merge(si_duprmed, foam.height[, c(2:4, 11)], "new_long_id", all.x=T)
dim(si_duprmed)
names(si)
names(si_duprmed)
head(si_duprmed[, c("Foaming.Status.x", "Foaming.Status.y")]
)
tail(si_duprmed[, c("Foaming.Status.x", "Foaming.Status.y")])
si_duprmed$Foaming.Status.y<-NULL
colnames(si_duprmed)[13]
colnames(si_duprmed)[13]<-"Foaming.Status"
dim(si_duprmed)
colnames(si_duprmed)[122]
colnames(si_duprmed)[122]<-"foam.type"
data
head(si_duprmed[, 1:5])
 write.table(si_duprmed, "fixed_inputs/sample_information_donotuseRemoved_rearranged_w_alpha_indices_time_index_foam_height_dupremoved.txt", sep="\t", quote=F, row.names=F)
data
row.names(si_duprmed)<-si_duprmed$SAMPLES
sample_data(data)<-si_duprmed
data
sample_data(data_1e5)<-si_duprmed
saveRDS(data, "RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/all_valid_samples_min_taxasums_5_min_seq_10k_physeq.RDS")
saveRDS(data_1e5, "RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/all_valid_samples_min_taxasums_5_min_seq_10k_1e5standardized_physeq.RDS")
si_duprmed<-si
dim(si_duprmed)
rm(si_duprmed)
si<-data.frame(sample_data(data))
dim(si)
names(si)
test<-si[, c("Foaming.Status", "surface.depth", "foam.type")]
head(test)
test$Foaming.Status<-factor(test$Foaming.Status, levles=c("0", "1"))
test$Foaming.Status<-factor(test$Foaming.Status, levels=c("0", "1"))
plot(test$surface.depth ~ test$Foaming.Status)
plot(test$foam.type ~ test$Foaming.Status)
plot(test$surface.depth ~ test$foam.type)
si$Foaming.Status<-factor(si$Foaming.Status, levels=c("0", "1"))
test<-si[, c("id", "Foaming.Status", "surface.depth", "foam.type")]
head(test)
unique(si$category)
si$category<-factor(si$category, levels=c("T0", "T0-10", "T10-50
si$category<-factor(si$category, levels=c("T0", "T0-10", "T10-50", "T50-75", "T75-100", "T100"))
head(si$category)
unique(si$category)
sample_data(data)<-si
sample_data(data_1e5)<-si
data
data_1e5
saveRDS(data, "RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/all_valid_samples_min_taxasums_5_min_seq_10k_physeq.RDS")
saveRDS(data_1e5, "RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/all_valid_samples_min_taxasums_5_min_seq_10k_1e5standardized_physeq.RDS")
test<-si[, c("id", "category", "Foaming.Status", "surface.depth", "foam.type")]
head(test)
plot(test$id ~ test$foam.type)
plot(test$surface.depth ~ test$id)
test<-si[, c("id", "category", "Foaming.Status", "surface.depth", "foam.type", "foaming.rate")]
head(test)
library(ggplot2)
ggplot(test, aes(x=category, y=surface.depth, fill=Foaming.Status)) + geom_point()
ggplot(test, aes(x=category, y=surface.depth, fill=as.factor(Foaming.Status))) + geom_point()
ggplot(test, aes(x=category, y=surface.depth, colors=Foaming.Status)) + geom_point()
str(test)
ggplot(test, aes(x=category, y=surface.depth, colors=Foaming.Status)) + geom_point(colors=Foaming.Status)
ggplot(test, aes(x=category, y=surface.depth, colors=Foaming.Status)) + geom_point(aes(colors=Foaming.Status))
ggplot(test, aes(x=category, y=surface.depth, colors=as.factor(Foaming.Status))) + geom_point(aes(colors=Foaming.Status))
ggplot(test, aes(x=category, y=surface.depth)) + geom_point(aes(colors=as.factor(Foaming.Status)))
ggplot(test, aes(surface.depth)) + geom_point(aes(colors=as.factor(Foaming.Status)))
ggplot(test, aes(id, surface.depth)) + geom_point(aes(colors=as.factor(Foaming.Status)))
ggplot(test, aes(id, surface.depth)) + geom_point(aes(fill=as.factor(Foaming.Status)))
ggplot(test, aes(id, surface.depth)) + geom_point(aes(fill=Foaming.Status))
p<-ggplot(test, aes(id, surface.depth))
p + geom_point(aes(colour=Foaming.Status))
p<-ggplot(test, aes(category, surface.depth))
p + geom_point(aes(colour=Foaming.Status))
p + geom_point(aes(colour=foam.type))
plot(si$Adj.MPR ~ si$surface.depth)
plot(si$Pielou ~ si$surface.depth)
plot(si$Observed ~ si$surface.depth)
plot(si$ ~ si$surface.depth)
ggplot(test, aes(surface.depth)) + geom_point(aes(color=foam.type))
ggplot(si, aes(surface.depth, Acetic.Acid)) + geom_point(aes(color=foam.type))
### adonis only ###
totu<-t(data.frame(otu_table(data_1e5)))
library(vegan)
### adonis only ###
totu<-t(data.frame(otu_table(data_1e5)))
si<-data.frame(sample_data(data_1e5))
cn<-names(si[, c(1:2, 7, 10:13, 15:20)]
df<-data.frame()
for (i in names(si[, cn])){
        tryCatch({
results.rela<-adonis(totu ~ as.factor(si[, i]), strata=si$id)
df1<-data.frame(i, results.rela$aov.tab[5][1, ], results.rela$aov.tab[6][1, ])
totu.trans<-decostand(totu, "pa")
results.pa<-adonis(totu.trans ~ as.factor(si[, i]), strata=si$id)
df2<-data.frame(i, results.pa$aov.tab[5][1, ], results.pa$aov.tab[6][1, ])
df.temp<-cbind(df1,df2)
df<-rbind(df, df.temp)
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
names(df)<-c("factor", "rela.adonis_R2", "rela.Pr(>F)", "factor", "pa.adonis_R2", "pa.Pr(>F)")
names(si[, c(1:2, 7, 10:13, 15:20)]
)
names(si)
cn<-names(si[, c(2, 115, 13:14, 18:21)])
df<-data.frame()
for (i in names(si[, cn])){
        tryCatch({
results.rela<-adonis(totu ~ as.factor(si[, i]), strata=si$id)
df1<-data.frame(i, results.rela$aov.tab[5][1, ], results.rela$aov.tab[6][1, ])
totu.trans<-decostand(totu, "pa")
results.pa<-adonis(totu.trans ~ as.factor(si[, i]), strata=si$id)
df2<-data.frame(i, results.pa$aov.tab[5][1, ], results.pa$aov.tab[6][1, ])
df.temp<-cbind(df1,df2)
df<-rbind(df, df.temp)
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
cn<-names(si[, c(2, 115, 13:14, 18:21, 118:122)])
df<-data.frame()
for (i in names(si[, cn])){
        tryCatch({
results.rela<-adonis(totu ~ as.factor(si[, i]), strata=si$id)
df1<-data.frame(i, results.rela$aov.tab[5][1, ], results.rela$aov.tab[6][1, ])
totu.trans<-decostand(totu, "pa")
results.pa<-adonis(totu.trans ~ as.factor(si[, i]), strata=si$id)
df2<-data.frame(i, results.pa$aov.tab[5][1, ], results.pa$aov.tab[6][1, ])
df.temp<-cbind(df1,df2)
df<-rbind(df, df.temp)
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
names(df)<-c("factor", "rela.adonis_R2", "rela.Pr(>F)", "factor", "pa.adonis_R2", "pa.Pr(>F)")
df
df[order(df[, c("relat.Pr(>F)", "rela.adonis_R2")]),]
df[order(df[, c(2:3)]),]
data_1e5
names(si)
totu<-t(data.frame(otu_table(data_1e5)))
si<-data.frame(sample_data(data_1e5))
cn<-names(si[, c(2, 115, 13:14, 118, 122, 18:21)])
df<-data.frame()
for (i in names(si[, cn])){
        tryCatch({
results.rela<-adonis(totu ~ as.factor(si[, i]), strata=si$id)
df1<-data.frame(i, results.rela$aov.tab[5][1, ], results.rela$aov.tab[6][1, ])
totu.trans<-decostand(totu, "pa")
results.pa<-adonis(totu.trans ~ as.factor(si[, i]), strata=si$id)
df2<-data.frame(i, results.pa$aov.tab[5][1, ], results.pa$aov.tab[6][1, ])
df.temp<-cbind(df1,df2)
df<-rbind(df, df.temp)
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
names(df)<-c("factor", "rela.adonis_R2", "rela.Pr(>F)", "factor", "pa.adonis_R2", "pa.Pr(>F)")
df
cn
write.table(df, "all_valid_donotuseremoved_min10k_mintax_5_dupremoved_adonis_strata_barns_summary.txt", sep="\t", quote=F, row.names=F)
df<-data.frame()
for (i in names(si[, cn])){
        tryCatch({
results.rela<-adonis(totu ~ as.factor(si[, i]))
df1<-data.frame(i, results.rela$aov.tab[5][1, ], results.rela$aov.tab[6][1, ])
totu.trans<-decostand(totu, "pa")
results.pa<-adonis(totu.trans ~ as.factor(si[, i]), strata=si$id)
df2<-data.frame(i, results.pa$aov.tab[5][1, ], results.pa$aov.tab[6][1, ])
df.temp<-cbind(df1,df2)
df<-rbind(df, df.temp)
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
names(df)<-c("factor", "rela.adonis_R2", "rela.Pr(>F)", "factor", "pa.adonis_R2", "pa.Pr(>F)")
df
write.table(df, "all_valid_donotuseremoved_min10k_mintax_5_dupremoved_adonis_summary.txt", sep="\t", quote=F, row.names=F)
totu<-t(data.frame(otu_table(data_1e5)))
si<-data.frame(sample_data(data_1e5))
cn<-names(si[, c(2, 115, 13:14, 118, 122, 18:21)])
df<-data.frame()
for (i in names(si[, cn])){
        tryCatch({
results.rela<-adonis(totu ~ as.factor(si[, i]))
df1<-data.frame(i, results.rela$aov.tab[5][1, ], results.rela$aov.tab[6][1, ])
totu.trans<-decostand(totu, "pa")
results.pa<-adonis(totu.trans ~ as.factor(si[, i]))
df2<-data.frame(i, results.pa$aov.tab[5][1, ], results.pa$aov.tab[6][1, ])
df.temp<-cbind(df1,df2)
df<-rbind(df, df.temp)
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
df
names(df)<-c("factor", "rela.adonis_R2", "rela.Pr(>F)", "factor", "pa.adonis_R2", "pa.Pr(>F)")
df
write.table(df, "all_valid_donotuseremoved_min10k_mintax_5_dupremoved_adonis_summary.txt", sep="\t", quote=F, row.names=F)
library(DESeq2)
data
names(si)
head(si$group.x)
head(si$group.y)
head(si[, c("group.y", "category")])
diagdds = phyloseq_to_deseq2(data, ~ group.y)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
head(res)
sigtab_t100<-subset(sigtab, log2FoldChange < 0)
dim(sigtab_t100)
t100_2fold_sig<-row.names(sigtab_t100)
write(t100_2fold_sig, "t100_2fold_sig_otu_list.txt")
head(sigtab_t100)
save.image("all_16s_physeq.RData")
savehistory(
"temp.R")
load("all_16s_physeq.RData")
ls()
library(phyloseq)
library(ggplot2)
data_1e5
library(vegan)
totu<-t(data.frame(otu_table(data_1e5)))
si<-data.frame(sample_data(data_1e5))
head(totu[, 1:5])
head(si[, 1:5])
head(si[, 1:20])
library(RColorBrewer)
colorCount = length(unique(si$Foaming.Status))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
test<-head(si[, c("month", "some_number)])
)
test<-head(si[, c("month", "some_number")])
test
test$x<-row.names(test)
aggregate(test[, 1:2], list(group=colnames(test)[3]), mean)
test
aggregate(test[, 1:2], list(group=x), mean)
aggregate(test, list(group=x), mean)
aggregate(test, list(x), mean)
aggregate(test[, 1:2], mean)
aggregate(test[, 1:2], fun=mean)
aggregate(test[, 1:2], FUN=mean)
aggregate(test[, 1:2], FUN=mean, by=x)
aggregate(test[, 1:2], FUN=mean, by="x")
str(test)
test<-data.frame(a=1:4, b=5:8, c=c(1,1,2,2))
test
test<-aggregate(test[, a:b], list(group=c), mean)
test<-aggregate(test[, 1:2], list(group=c), mean)
test<-aggregate(test[, 1:2], list(c), mean)
test<-aggregate(test[, 1:2], list(group=test$c), mean)
test
test<-data.frame(a=1:4, b=5:8, c=c(1,1,2,2))
\
test<-data.frame(a=1:4, b=5:8, c=c(1,1,2,2))
aggregate(test[, 1:2], list(group=test[, 3]), mean)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
data.mds<-metaMDS(data.trans, k=3, autotransform=FALSE)
for (i in names si[, c("Foaming.Status", "category")]){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, i, colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
data.mds<-metaMDS(totu, k=3, autotransform=FALSE)
for (i in names si[, c("Foaming.Status", "category")]){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, i, colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, i, colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
for (i in names(si[, c("Foaming.Status", "category")])){
print (i)
for (i in names(si[, c("Foaming.Status", "category")])){
print(i)}
ggplot.NMDS.ellipse
str(si[, c("Foaming.Status", "category")])
colors
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
data.mds<-metaMDS(totu, k=3, autotransform=FALSE)
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, i, colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, i, colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
plot(si$Pielou ~ si$surface.depth)
dev.off()
dev.off()
plot(si$Pielou ~ si$surface.depth)
quartz()
plot(si$Pielou ~ si$surface.depth)
ls()
t100_2fold_sig.phy
t100_2fold_sig.phy<-prune_taxa(taxa_names(data_1e5) %in% t100_2fold_sig_otu, data_1e5)
t100_2fold_sig.phy
length(t100_2fold_sig_otu)
dim(t100_2fold_sig)
length(t100_2fold_sig)
head(t100_2fold_sig_otu)
rm(t100_2fold_sig_otu)
dim(sigtab_t100)
t100_2fold_sig.phy<-prune_taxa(taxa_names(data_1e5) %in% t100_2fold_sig, data_1e5)
t100_2fold_sig.phy
t100_2fold_sig_1e5.phy<-t100_2fold_sig.phy
t100_2fold_sig.phy<-prune_taxa(taxa_names(data) %in% t100_2fold_sig, data)
t100_2fold_sig.phy
plot_bar(t100_2fold_sig.phy, x="category", fill="genus")+ geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+f
plot_bar(t100_2fold_sig.phy, x="category", fill="genus")+ geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+facet_grid(phylum~., scales="free")
plot_bar(t100_2fold_sig.phy, x="category", fill="phylum")+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack"
)+facet_grid(phylum~., scales="free")
plot_bar(t100_2fold_sig.phy, fill="phylum")
plot_bar(t100_2fold_sig.phy, fill="phylum")+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~., scales="free")
plot_bar(t100_2fold_sig.phy, fill="phylum")+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~category., scales="free")
plot_bar(t100_2fold_sig.phy, x="SAMPLES", group="category", fill="phylum")+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~category., scales="free")
plot_bar(t100_2fold_sig.phy, x="SAMPLES", fill="phylum")+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~category., scales="free")
names(si)
si$category
plot_bar(t100_2fold_sig.phy, x="SAMPLES", fill="phylum")+ geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(~category., scales="free")
t100_2fold_sig.phy
min(sample_sums(t100_2fold_sig.phy))
plot_bar(t100_2fold_sig.phy, fill="phylum", facet=.~category)
plot_bar(t100_2fold_sig.phy, fill="phylum") + facet_grid(~category., scales="free")
plot_bar(t100_2fold_sig.phy, fill="phylum") + facet_grid(.~category)
plot_bar(t100_2fold_sig.phy, fill="phylum") + facet_grid(phylum~category)
plot_bar(t100_2fold_sig.phy, fill="phylum") + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~category)
plot_bar(t100_2fold_sig.phy, fill="phylum") + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~category, scales="free")
plot_bar(t100_2fold_sig_1e5.phy, fill="phylum") + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~category, scales="free")
saveRDS(t100_2fold_sig.phy, "RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/t100_2fold_sig_physeq.RDS")
saveRDS(t100_2fold_sig_1e5.phy, "RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/t100_2fold_sig_1e5_physeq.RDS")
subset_taxa(data, taxa_names(data)=="OTU_885")
tax_table(subset_taxa(data, taxa_names(data)=="OTU_885"))
totu[1:5, 1:5]
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
ls()
t100_2fold_sig
t100_2fold_sig.phy
plot_bar(t100_2fold_sig.phy, fill="phylum") + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~category, scales="free")
plot_bar(t100_2fold_sig.phy, fill="surface.depth") + geom_bar(aes(color=surface.depth, fill=surface.depth), stat="identity", position="stack")+facet_grid(phylum~category, scales="free")
plot_bar(t100_2fold_sig.phy, fill="foam.type") + geom_bar(aes(color=foam.type, fill=foam.type), stat="identity", position="stack")+facet_grid(phylum~category, scales="free")
plot_bar(t100_2fold_sig.phy, fill="surface.depth") + geom_bar(aes(color=surface.depth, fill=surface.depth), stat="identity", position="stack")+facet_grid(phylum~category, scales="free")
names(si)
plot(si$time_index ~ surface.depth)
plot(si$time_index ~ si$surface.depth)
ggplot(si, aes(x=surface.depth, y=time_index)) + geom_point(aes(color=Foaming.Status))
ggplot(si, aes(x=surface.depth, y=Acetic.Acid)) + geom_point(aes(color=Foaming.Status))
plot(si$surface.depth ~ si$Pielou)
ggplot(si, aes(x=surface.depth, y=Pielou)) + geom_point(aes(color = Foaming.Status))
ggplot(si, aes(x=surface.depth, y=Pielou)) + geom_point(aes(color = foam.type))
ggplot(si, aes(x=surface.depth, y=time_index)) + geom_point(aes(color = foam.type))
ggplot(si, aes(x=surface.depth, y=Pielou)) + geom_point(aes(color = time_index))
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = time_index))
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = id))
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = category))
ls()
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = category))+facet_grid(category ~ .)
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = category))+facet_grid( ~ category)
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = category))+facet_wrap(category)
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = category))+facet_wrap(~category)
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = category))+facet_wrap(~category) + theme_bw()
test<-prune_taxa(taxa_names(data_1e5) == c("OTU_9287", "OTU_9294"), data_1e5)
test
tax_table(test)
plot(si$Foaming.Capacity ~ si$surface.depth)
plot(si$Foam.Stability ~ si$surface.depth)
ggplot(si, aes(x=surface.depth, y=Foam.Stability)) + geom_point(aes(color = foam.type))
ggplot(si, aes(x=Organic.N, y=Foam.Stability)) + geom_point(aes(color = foam.type))
ggplot(si, aes(x=surface.depth, y=Acetic.Acid)) + geom_point(aes(color = foam.type))
ggplot(si, aes(x=surface.depth, y=Acetic.Acid)) + geom_point(aes(color = category))
ggplot(si, aes(x=surface.depth, y=Acetic.Acid)) + geom_point(aes(color = category))+facet_wrap(~category) + theme_bw()
head(tax)
tax["OTU_13720", ]
tax[row.names(tax)=="OTU_13720", ]
tax[c("OTU_13720"), ]
tax[rownames(tax)=="OTU_13720", ]
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load(
"all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
ggplot(si, aes(x=surface.depth, y=Acetic.Acid)) + geom_point(aes(color = category))+facet_wrap(~category) + theme_bw()
dim(totu)
dim(si)
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
library(RColorBrewer)
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
library(vegan)
for (i in names(si[, c("Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=TKN)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~category) + theme_bw()
ggplot(si, aes(x=surface.depth, y=TKN)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=TKN)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=DDGS)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=Crude.Protein)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
names(si)
ggplot(si, aes(x=surface.depth, y=N.Take)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=(Crude.Protein+Lysine))) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=(Crude.Protein+Lysine)/Crude.Fiber)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=Lysine)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~category) + theme_bw()
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~Foaming.Status) + theme_bw()
ggplot(si, aes(x=surface.depth, y=Organic.N)) + geom_point(aes(color = Foaming.Status))+facet_wrap(~id) + theme_bw()
head(alpha)
ls()
head(rich)
ls()
head(totu[, 1:5])
dim(si)
names(si)
ls()
df<-envfit(si, data.mds, permu=999)
library(vegan)
df<-envfit(data.mds, si, permu=999)
df<-envfit(data.mds, si, permu=999, na.rm=T)
df<-envfit(data.mds, si[, 19:122], permu=999, na.rm=T)
data(varechem)
head(varechem)
data(varespec)
head(varespec)
si[1:5, 1:5])
si[1:5, 1:5]
rm(varechem)
rm(varespec)
ls()
df<-envfit(data.mds, si[, 19:20], permu=999, na.rm=T)
df
df<-envfit(data.mds, si[, 19:100], permu=999, na.rm=T)
df<-envfit(data.mds, si[, c("pH","time_index", "State", "surface.depth", "SCFA", "LCFA", "Adj.MPR", "Acetic.Acid")], permu=999, na.rm=T)
df<-envfit(data.mds, si[, c("pH","time_index", "state", "surface.depth", "SCFA", "LCFA", "Adj.MPR", "Acetic.Acid")], permu=999, na.rm=T)
df
df<-envfit(data.mds, si[, c("pH","time_index", "state", "surface.depth", "SCFA", "LCFA", "Adj.MPR", "Acetic.Acid", "TKN")], permu=999, na.rm=T)
df
df<-envfit(data.mds, si[, c("pH","time_index", "state", "surface.depth", "SCFA", "LCFA", "Adj.MPR", "Acetic.Acid", "TKN", "Crude.Protein", "Crude.Fiber", "Carbon")], permu=999, na.rm=T)
df
df<-envfit(data.mds, si[, c("pH","time_index", "state", "surface.depth", "month")], permu=999, na.rm=T)
df
df<-envfit(data.mds, si[, c("pH")], permu=999, na.rm=T)
df<-envfit(data.mds, si[, c("pH", "time_index")], permu=999, na.rm=T)
df
df<-envfit(data.mds, si[, c("month", "time_index")], permu=999, na.rm=T)
df
plot(data.mds, display="sites")
plot(df, p.max=0.1)
rquartz
quartz
quartz()
plot(data.mds, display="sites")
plot(df, p.max=0.1)
df<-envfit(data.mds, si[, c("month", "time_index", "surface.depth")], permu=999, na.rm=T)
df
plot(data.mds, display="sites")
plot(df, p.max=0.1)
df<-envfit(data.mds, si[, c("month", "time_index", "Carbon")], permu=999, na.rm=T)
df
plot(data.mds, display="sites")
plot(df, p.max=0.1)
df<-envfit(data.mds, si[, c("month", "time_index", "foaming.rate")], permu=999, na.rm=T)
df<-envfit(data.mds, si[, c("month", "time_index", "Acetic.Acid")], permu=999, na.rm=T)
df
df<-envfit(data.mds, si[, c("month", "time_index", "Adj.MPR")], permu=999, na.rm=T)
df
plot(data.mds, display="sites")
plot(df, p.max=0.1)
plot(si$surface.depth ~ si$id)
plot(si$surface.depth ~ si$category)
plot(si$surface.depth ~ si$Adj.MPR)
plot(si$surface.depth ~ si$MPR_slurry)
plot(si$surface.depth ~ si$Adj.MPR)
plot(si$surface.depth/si$Adj.MPR ~ si$Organic.N)
plot(si$surface.depth ~ si$Organic.N)
plot(si$surface.depth ~ si$Crude.Protein)
plot(si$surface.depth ~ si$Crude.Protein/si$Crude.Fiber)
plot(si$surface.depth ~ (si$Crude.Protein/si$Crude.Fiber))
plot((si$Crude.Protein/si$Crude.Fiber) ~ si$surface.depth)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
for (i in names(si[, c("id", "Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
for (i in names(si[, c("id", "Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
for (i in names(si[, c("id", "Foaming.Status", "category")])){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
colorCount = length(unique(si[, i]))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
pdf(paste("all_valid_1e5_by_",i,".pdf", sep=""))
p<-ggplot.NMDS.ellipse(data.mds, si[, i], colors)
print(p)
dev.off()
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library("phyloseq")
si$Long_Sample_ID[is.na(si$foam.type)]
test<-si$Long_Sample_ID[is.na(si$foam.type)]
test<-si[si$Long_Sample_ID %in% test, ]
library(plyr)
count(test, "Foaming.Status")
test<-subset_samples(data_1e5, !is.na(foam.type))
totu<-t(data.frame(otu_table(test)))
si<-data.frame(sample_data(test))
adonis(totu ~ si$Foaming.Status, strata=si$id)
library(vegan)
adonis(totu ~ si$Foaming.Status, strata=si$id)
cn<-names(si[, c(2, 115, 13:14, 118, 122, 18:21)])
cn
df<-data.frame()
for (i in names(si[, cn])){
        tryCatch({
results.rela<-adonis(totu ~ as.factor(si[, i]), strata=si$id)
df1<-data.frame(i, results.rela$aov.tab[5][1, ], results.rela$aov.tab[6][1, ])
totu.trans<-decostand(totu, "pa")
results.pa<-adonis(totu.trans ~ as.factor(si[, i]), strata=si$id)
df2<-data.frame(i, results.pa$aov.tab[5][1, ], results.pa$aov.tab[6][1, ])
df.temp<-cbind(df1,df2)
df<-rbind(df, df.temp)
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
names(df)<-c("factor", "rela.adonis_R2", "rela.Pr(>F)", "factor", "pa.adonis_R2", "pa.Pr(>F)")
head(df)
df
write.table(df, "all_valid_donotuseremoved_min10k_mintax_5_dupremoved_adonis_strata_barns_summary_by_foamtype.txt", sep="\t", quote=F, row.names=F)
unique(sample_variables(data))
unique(sample_variables(data)=="foam.type")
unique(sample_variables(data)$foam.type)
unique(sample_data(data)$foam.type)
df<-data.frame()
for (i in names(si[, cn])){
        tryCatch({
results.rela<-adonis(totu ~ as.factor(si[, i]))
df1<-data.frame(i, results.rela$aov.tab[5][1, ], results.rela$aov.tab[6][1, ])
totu.trans<-decostand(totu, "pa")
results.pa<-adonis(totu.trans ~ as.factor(si[, i]), strata=si$id)
df2<-data.frame(i, results.pa$aov.tab[5][1, ], results.pa$aov.tab[6][1, ])
df.temp<-cbind(df1,df2)
df<-rbind(df, df.temp)
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
names(df)<-c("factor", "rela.adonis_R2", "rela.Pr(>F)", "factor", "pa.adonis_R2", "pa.Pr(>F)")
df
write.table(df, "all_valid_donotuseremoved_min10k_mintax_5_dupremoved_adonis_summary_by_foamtype.txt", sep="\t", quote=F, row.names=F)
df<-data.frame()
for (i in names(si[, cn])){
        tryCatch({
results.rela<-adonis(totu ~ as.factor(si[, i]))
df1<-data.frame(i, results.rela$aov.tab[5][1, ], results.rela$aov.tab[6][1, ])
totu.trans<-decostand(totu, "pa")
results.pa<-adonis(totu.trans ~ as.factor(si[, i]))
df2<-data.frame(i, results.pa$aov.tab[5][1, ], results.pa$aov.tab[6][1, ])
df.temp<-cbind(df1,df2)
df<-rbind(df, df.temp)
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
names(df)<-c("factor", "rela.adonis_R2", "rela.Pr(>F)", "factor", "pa.adonis_R2", "pa.Pr(>F)")
df
write.table(df, "all_valid_donotuseremoved_min10k_mintax_5_dupremoved_adonis_summary_by_foamtype.txt", sep="\t", quote=F, row.names=F)
data.ftype<-subset_samples(data, !is.na(foam.type))
for (i in unique(sample_data(data.ftype)$foam.type){
     temp<-subset_samples(data.ftype, foam.type == i)
    n = nsamples(temp)
     temp.core<-filter_taxa(temp, function(x) sum(x >=1 ) = (1*length(x)), TRUE)
    print(i)
    print(temp.core)
data.ftype
for (i in unique(sample_data(data.ftype))$foam.type){
     temp<-subset_samples(data.ftype, foam.type == i)
    n = nsamples(temp)
     temp.core<-filter_taxa(temp, function(x) sum(x >=1 ) = (1*length(x)), TRUE)
    print(i)
    print(temp.core)
}
data.ftype<-subset_samples(data, !is.na(foam.type))
for (i in unique(sample_data(data.ftype))$foam.type){
     temp<-subset_samples(data.ftype, foam.type == i)
    n = nsamples(temp)
     temp.core<-filter_taxa(temp, function(x) sum(x >0 ) = (1*length(x)), TRUE)
    print(i)
    print(temp.core)
}
GP = filter_taxa(test, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
gp
GP
GP = filter_taxa(test, function(x) sum(x > 0) > (0.2*length(x)), TRUE)
GP
GP = filter_taxa(test, function(x) sum(x > 0) > (1*length(x)), TRUE)
GP = filter_taxa(test, function(x) sum(x > 0) = (1*length(x)), TRUE)
GP = filter_taxa(test, function(x) sum(x > 0) > (0.999*length(x)), TRUE)
GP
GP = filter_taxa(test, function(x) sum(x > 0) > (0.99999*length(x)), TRUE)
GP
GP = filter_taxa(test, function(x) sum(x > 0) > (0.99*length(x)), TRUE)
GP
GP = filter_taxa(test, function(x) sum(x > 0) >= (0.99*length(x)), TRUE)
GP = filter_taxa(test, function(x) sum(x > 0) >= (1*length(x)), TRUE)
GP
rm(GP)
for (i in unique(sample_data(data.ftype))$foam.type){
     temp<-subset_samples(data.ftype, foam.type == i)
    n = nsamples(temp)
     temp.core<-filter_taxa(temp, function(x) sum(x >0 ) >= (1*length(x)), TRUE)
    print(i)
    print(temp.core)
}
for (i in unique(sample_data(data.ftype)$foam.type)){
     temp<-subset_samples(data.ftype, foam.type == i)
    n = nsamples(temp)
     temp.core<-filter_taxa(temp, function(x) sum(x >0 ) >= (1*length(x)), TRUE)
    print(i)
    print(temp.core)
}
min(sample_sums(temp.core))
data.ftype<-subset_samples(data, !is.na(foam.type))
for (i in unique(sample_data(data.ftype)$foam.type)){
     temp<-subset_samples(data.ftype, foam.type == i)
    n = nsamples(temp)
     temp.core<-filter_taxa(temp, function(x) sum(x >0 ) >= (1*length(x)), TRUE)
    print(i)
    print(temp.core)
    saveRDS(temp, paste("data_ftype_", i, "_core.RDS", sep=""))
    print("DONE!")
}
colnames(si[, c(1:21, 108:120, 122)]])
colnames(si[, c(1:21, 108:120, 122)])
count(si, "foam.type")
data.ftype<-subset_samples(data, !is.na(foam.type))
for (i in unique(sample_data(data.ftype)$foam.type)){
     temp<-subset_samples(data.ftype, foam.type == i)
    n = nsamples(temp)
     temp.core<-filter_taxa(temp, function(x) sum(x >0 ) >= (1*length(x)), TRUE)
    print(i)
    print(temp.core)
    saveRDS(temp.core, paste("data_ftype_", i, "_core.RDS", sep=""))
    print("DONE!")
}
test<-readRDS("RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/data_ftype_C_core.RDS")
test
data.ftype
for (i in unique(sample_data(data.ftype)$foam.type)){
     temp<-subset_samples(data.ftype, foam.type == i)
    n = nsamples(temp)
     temp.core<-filter_taxa(temp, function(x) sum(x >0 ) >= (1*length(x)), TRUE)
    print(i)
    print(temp.core)
    saveRDS(temp.core, paste("data_ftype_", i, "_core.RDS", sep=""))
     print(temp.core)
    print("DONE!")
}
test<-readRDS("RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/data_ftype_C_core.RDS")
test
temp.core
saveRDS(temp.core, paste("data_ftype_test" , "_core.RDS", sep=""))
saveRDS(temp.core, "test.RDS")
test<-readRDS("test.RDS")
test
test<-readRDS("data_ftype_C_core.RDS")
test
for (i in unique(sample_data(data.ftype)$foam.type)){
     temp<-subset_samples(data.ftype, foam.type == i)
    n = nsamples(temp)
     temp.core<-filter_taxa(temp, function(x) sum(x >0 ) >= (1*length(x)), TRUE)
    print(i)
    saveRDS(temp.core, paste("data_ftype_", i, "_core.RDS", sep=""))
     print(temp.core)
    print("DONE!")
}
test<-readRDS("RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/data_ftype_C_core.RDS")
test
physeq<-subset_taxa(physeq, domain!="Archaea" & domain!="unclassified_Root")
physeq<-subset_taxa(test, domain!="Archaea" & domain!="unclassified_Root")
physeq
x<-rnorm(162)
y<-0+rnorm(162)
plot(y ~ x)
mean(taxa_sums(temp.core))
temp.core
mean(taxa_sums(temp.core))/162
min(taxa_sums(temp.core))/162
max(taxa_sums(temp.core))/162
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(ggplot2)
library(phyloseq)
ls()
data
data.ftype
data_1e5.ftype<-subset_samples(data_1e5, !is.na(foam.type))
data_1e5.ftype
min(taxa_sums(data_1e5.ftype)
)
ls
ls()
core_otu<-read.delim("foaming_status_cc/spearman_hoeffding_cc/raw_counts/core_by_foamtype/core_otus.txt", sep="\t", header=T)
head(core_otu)
data_1e5.ftype.core<-prune_taxa(taxa_names(data_1e5.ftype) %in% core_otu$core_otus, data_1e5.ftype
)
data_1e5.ftype.core
min(sample_sums(data_1e5.ftype.core)
)
plot_bar(data_1e5.ftype.core, fill="genus", facet="foam.type")
plot_bar(data_1e5.ftype.core, fill="genus") + geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack") + facet_grid(phylum ~ foam.type, scales="free")
plot_bar(data_1e5.ftype.core, fill="genus", group="time_index") + geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack") + facet_grid(phylum ~ foam.type, scales="free")
plot_bar(data_1e5.ftype.core, fill="genus") + geom_bar(aes(color=genus, fill=genus, group=time_index), stat="identity", position="stack") + facet_grid(phylum ~ foam.type, scales="free")
library(DESeq2)
data.ftype.core<-prune_taxa(taxa_names(data.ftype) %in% core_otu$core_otus, data.ftype)
data.ftype.core
sample_variables(data.ftype.core)
diagdds = phyloseq_to_deseq2(data, ~ group.y)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
head(sigtab)
summary(res)
summary(diagdds)
diagdds
res
dim(sigtab)
ls()
diagdds = phyloseq_to_deseq2(data.ftype.core, ~ group.y)
diagdds
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
res
head(sigtab)
dim(sigtab)
x<-tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
x
x=sort(x, T)
str(sigtab$phylum)
sigtab$phylum<-factor(as.character(sigtab$phylum, levels=names(x))
sigtab$phylum<-factor(as.character(sigtab$phylum), levels=names(x))
x<-tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x<-sort(x, T)
sigtab$genus<-factor(as.character(sigtab$genus), levels=names(x))
ggplot(sigtab, aes(y=Genus, x=log2FoldChange, color=Phylum)) + geom_vline(xintercept=0.0, color="black", size=1) + geom_point(size=6) + theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) 
ggplot(sigtab, aes(y=genus, x=log2FoldChange, color=phylum)) + geom_vline(xintercept=0.0, color="black", size=1) + geom_point(size=6) + theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) 
head(sigtab)
dim(si)
head(si[, c("SAMPLES", "id", "foam.type", "category", "group.y")]
)
diagdds = phyloseq_to_deseq2(data.ftype.core, ~ Foaming.Status)
diagdds = phyloseq_to_deseq2(data.ftype.core, ~ as.factor(Foaming.Status))
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
diagdds = phyloseq_to_deseq2(data.ftype.core, ~ foam.type)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
res
head(sigtab)
sigtab$test<-ifelse(sigtab$log2FoldChange > 0, "No_Foam", "Crust")
head(sigtab)
sigtab$test<-NULL
sigtab$foam.type<-ifelse(sigtab$log2FoldChange > 0, "No_Foam", "Crust")
sigtab$foam.type<-ifelse(sigtab$log2FoldChange > 0, "No_Foam", "Crust")
x = tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$phylum = factor(as.character(sigtab$phylum), levels=names(x))
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
head(sigtab)
ggplot(sigtab, aes(y=genus, x=log2FoldChange, color=phylum)) + 
  geom_point(size=6) + 
  facet_grid(~foam.type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
dim(sigtab)
ggplot(sigtab, aes(y=genus, x=log2FoldChange, color=phylum)) + 
  geom_point(size=3) + 
  facet_grid(~foam.type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
test<-subset_samples(data.ftype.core, foam.type != "C")
test
data.ftype.core
min(taxa_sums(test))
min(sample_sums(test))
test<-subset_samples(data.ftype.core, foam.type != "C")
diagdds = phyloseq_to_deseq2(test, ~ foam.type)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
sigtab$foam.type<-ifelse(sigtab$log2FoldChange > 0, "No_Foam", "Crust")
x = tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$phylum = factor(as.character(sigtab$phylum), levels=names(x))
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
ggplot(sigtab, aes(y=genus, x=log2FoldChange, color=phylum)) + 
  geom_point(size=3) + 
  facet_grid(~foam.type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
sigtab$foam.type<-NULL
res
sigtab$foam.type<-ifelse(sigtab$log2FoldChange > 0, "No_Foam", "Foam")
ggplot(sigtab, aes(y=genus, x=log2FoldChange, color=phylum)) + 
  geom_point(size=3) + 
  facet_grid(~foam.type) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
ggplot(sigtab, aes(y=genus, x=log2FoldChange, color=phylum)) + 
  geom_point(size=3) + 
  facet_grid(~foam.type, scale="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
test<-subset_samples(data.ftype.core, foam.type != "NF")
diagdds = phyloseq_to_deseq2(test, ~ foam.type)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
res
sigtab$foam.type<-ifelse(sigtab$log2FoldChange > 0, "Foam", "Crust")
x = tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$phylum = factor(as.character(sigtab$phylum), levels=names(x))
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
ggplot(sigtab, aes(y=genus, x=log2FoldChange, color=phylum)) + 
  geom_point(size=3) + 
  facet_grid(~foam.type, scale="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
test<-subset_samples(data.ftype.core, foam.type != "F")
diagdds = phyloseq_to_deseq2(test, ~ foam.type)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
res
sigtab$foam.type<-ifelse(sigtab$log2FoldChange > 0, "No_Foam", "Crust")
x = tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$phylum = factor(as.character(sigtab$phylum), levels=names(x))
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
ggplot(sigtab, aes(y=genus, x=log2FoldChange, color=phylum)) + 
  geom_point(size=3) + 
  facet_grid(~foam.type, scale="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
plot_heatmap(data_1e5.ftype.core)
plot_heatmap(data_1e5.ftype.core, "NMDS", "bray", "foam.type", "genus"))
plot_heatmap(data_1e5.ftype.core, "NMDS", "bray", "foam.type", "genus")
plot_heatmap(data_1e5.ftype.core, "NMDS", "bray", "foam.type", "genus", low="#000033", high="#FF3300"))
plot_heatmap(data_1e5.ftype.core, "NMDS", "bray", "foam.type", "genus", low="#000033", high="#FF3300")
plot_heatmap(data_1e5.ftype.core, "NMDS", "bray", "foam.type", low="#000033", high="#FF3300")
test<-tax_glom(data_1e5.ftype.core, "genus")
test.totu<-t(data.frame(otu_table(test)))
test.tax<-data.frame(tax_table(test))
test.si<-data.frame(sample_data(test))
head(test.totu[, 1:5])
head(test.tax)
colnames(test.totu)<-test.tax$genus
head(test.totu[, 1:5])
to_plot<-merge(test.totu, test.si[, c("SAMPLES", "time_index", "foam.type")], by.x="row.names", by.y="SAMPLES")
head(to_plot[, 1:6])
head(to_plot[, ])
library(reshape)
dim(to_plot)
to_plot.m<-melt(to_plot, id.vars=c("Row.names", "time_index", "foam.type"))
head(to_plot.m)
max(colSums(to_plot))
max(colSums(to_plot[, 2:26]))
test<-data.frame(colSums(to_plot[, 2:26]))
head(test)
colnames(test)[1]<-"total"
test$genus<-row.names(test)
head(test)
test$genus<-reorder(test$genus, -test$total)
str(test)
head(to_plot.m)
to_plot.m$variable<-factor(to_plot.m$variable, levels=levels(test$genus))
str(to_plot.m)
to_plot.m$foam.type<-factor(to_plot.m$foam.type, levels=c("NF", "C", "F"))
to_plot.m$Row.names<-reorder(to_plot.m$Row.names, to_plot.m$time_index)
str(to_plot.m)
ggplot(to_plot.m, aes(x=Row.names, y=variable, group=foam.type)) + geom_tile(aes(fill=value), color="white") + scale_fill_gradient(low="white", high="steelblue")
ggplot(to_plot.m, aes(x=Row.names, y=variable, group=foam.type)) + geom_tile(aes(fill=rescale(value)), color="white") + scale_fill_gradient(low="white", high="steelblue")
ggplot(to_plot.m, aes(x=Row.names, y=variable, group=foam.type)) + geom_tile(aes(fill=scale(value)), color="white") + scale_fill_gradient(low="white", high="steelblue")
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(value)), color="white") + scale_fill_gradient(low="white", high="steelblue") + facet_grid(~foam.type)
ggplot(to_plot.m, aes(x=Row.names, y=variable, group=foam.type)) + geom_tile(aes(fill=scale(value)), color="white") + scale_fill_gradient(low="#000033", high="#FF3300")
ggplot(to_plot.m, aes(x=Row.names, y=variable, group=foam.type)) + geom_tile(aes(fill=scale(value)), color="white") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type)
ggplot(to_plot.m, aes(x=Row.names, y=variable, group=foam.type)) + geom_tile(aes(fill=scale(value)), color="white") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free")
ggplot(to_plot.m, aes(x=Row.names, y=variable, group=foam.type)) + geom_tile(aes(fill=scale(value)), color="#000033") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free")
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=value), color="#000033") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free")
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=log(value)), color="#000033") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free")
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=t(scale(t(value)))), color="#000033") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free")
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=apply(scale(value), 2, mean)), color="#000033") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free")
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=apply(scale(value), 2, var)), color="#000033") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free")
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(value), color="#000033") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free")
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(value)), color="#000033") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free")
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(value)), color="#000033") + scale_fill_gradient(low="yellow", high="#FF3300") + facet_grid(~ foam.type, scales="free")
theme_change <- theme(
 plot.background = element_blank(),
 panel.grid.minor = element_blank(),
 panel.grid.major = element_blank(),
 panel.background = element_blank(),
 panel.border = element_blank(),
 axis.line = element_blank(),
 axis.ticks = element_blank(),
 axis.text.x = element_blank(),
axis.title.x = element_blank(),
)
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(value)), color="#000033") + scale_fill_gradient(low="yellow", high="#FF3300") + facet_grid(~ foam.type, scales="free") + theme(panel.border=element_black())
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(value)), color="#000033") + scale_fill_gradient(low="yellow", high="#FF3300") + facet_grid(~ foam.type, scales="free") + theme(panel.border=element_blank())
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(value)), color="#000033") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + theme(panel.border=element_blank())
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=value), color="black") + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + theme(panel.border=element_blank())
plot_heatmap
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=value), color=value) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + theme(panel.border=element_blank())
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=value, color=value)) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + theme(panel.border=element_blank())
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=value, color=value)) + scale_fill_gradient(low="#000033", high="#FF3300") + scale_color_gradient(low="#000033", high="#FF3300")+ facet_grid(~ foam.type, scales="free") + theme(panel.border=element_blank())
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=value, color=value)) + scale_fill_gradient(low="#000033", high="#FF3300") + scale_color_gradient(low="#000033", high="#FF3300")+ facet_grid(~ foam.type, scales="free") + theme(panel.border=element_blank()) + scale_x_continuous(expand=c(0,0)) + 
    scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=value, color=value)) + scale_fill_gradient(low="#000033", high="#FF3300") + scale_color_gradient(low="#000033", high="#FF3300")+ facet_grid(~ foam.type, scales="free") + theme(panel.border=element_blank()) + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=value)) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=value)) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(log(value+1))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(log(value+1)))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(value))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
head(to_plot.m)
test<-subset(to_plot.m, foam.type="F")
head(test)
dim(test)
test$group<-paste(test$Row.names, test$time_index, sep="_")
dim(test)
head(test)
test<-subset(to_plot.m, foam.type=="F")
test$group<-paste(test$Row.names, test$time_index, sep="_")
head(test)
length(unique(test$group)
)
dim(test)
test<-subset(test, variable == "unclassified_Porphyromonadaceae")
dim(test)
head(test)
test$group<-reorder(test$group, -test$value)
str(test)
head(to_plot.m)
to_plot.m$group<-paste(to_plot.m$Row.names, to_plot.m$time_index, sep="_")
head(to_plot.m)
to_plot.m$group<-factor(to_plot.m$group, levels=levels(test$group))
ggplot(to_plot.m, aes(x=group, y=variable)) + geom_tile(aes(fill=scale(value))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
head(to_plot.m)
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=scale(value))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=(value))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
library(vegan)
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, standardize))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "standardize"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "chi.square"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "log"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "hellinger"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "normalize"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "max"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "total"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "wisconsin"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "range"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "range"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0)) + tehme(axix.text.x=element_blank())
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "range"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0)) + theme(axis.text.x=element_blank())
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "range"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0)) + theme(axis.text.x=element_blank(), axis.x.ticks=element_blank())
ggplot(to_plot.m, aes(x=Row.names, y=variable)) + geom_tile(aes(fill=decostand(value, "range"))) + scale_fill_gradient(low="#000033", high="#FF3300") + facet_grid(~ foam.type, scales="free") + scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0)) + theme(axis.text.x=element_blank(), axis.ticks=element_blank())
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
ls()
dim(si)
data.mds
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_polygon_eclipse.R")
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_polygon_eclipse.R")
ggplot.NMDS.poly.ellipse(data.mds, si[, c("id", "foam.type")], colors)
library(ggplot2)
library(vegan)
ggplot.NMDS.poly.ellipse(data.mds, si[, c("id", "foam.type")], colors)
library(phyloseq)
ls()
data
data_1e5.ftype
totu<-t(data.frame(otu_table(data_1e5.ftype))
)
head(totu[, 1:5])
head(si[, 1:5])
data.mds<-metaMDS(totu, k=3, autotransform=FALSE)
library(RColorBrewer)
colorCount = length(unique(si$foam.type))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(colorCount)
ggplot.NMDS.poly.ellipse(data.mds, si[, c("id", "foam.type")], colors)
MDS1<-data.frame(scores(data.mds))$NMDS1
MDS2<-data.frame(scores(data.mds))$NMDS2
NMDS<-data.frame(MDS1,MDS2,si[, c("id", "foam.type")])
head(NMDS)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  df_ell <- data.frame()
  for(g in levels(NMDS[, 4])){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS[, 4]==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
  }
head(df_ell)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_polygon_eclipse.R")
 ggplot.NMDS.poly.ellipse(data.mds, si[, c("id", "foam.type")], colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_polygon_eclipse.R")
 ggplot.NMDS.poly.ellipse(data.mds, si[, c("id", "foam.type")], colors)
hull.data<-data.frame()
for (i in levels(NMDS[, 3])){
        temp<-subset(NMDS, NMDS[, 3] == i)
        temp.hull<-temp[chull(temp[, c(1:2)]),]
        hull.data<-rbind(hull.data, temp.hull)
}
head(hull.data)
head(df_ell)
ggplot() + geom_polygon(data=hull.data, aes_string(x="MDS1", y="MDS2", fill=colnames(hull.data[, 3]), color = "grey"), alpha=0.3)
dim(hull.data)
dim(NMDS)
ggplot(NMDS, aes(x=NMDS1, y=NMDS2)) + geom_polygon(data=hull.data, aes_string(x="MDS1", y="MDS2", fill=colnames(hull.data[, 3]), color = "grey"), alpha=0.3)
colors
ggplot(NMDS, aes(x=NMDS1, y=NMDS2)) + geom_polygon(data=hull.data, aes_string(x="MDS1", y="MDS2"), alpha=0.3)
ggplot(NMDS, aes(x=NMDS1, y=NMDS2)) + geom_polygon(data=hull.data, aes_string(x="MDS1", y="MDS2", color="black"), alpha=0.3)
ggplot(NMDS, aes(x=NMDS1, y=NMDS2)) + geom_polygon(data=hull.data, aes_string(x="MDS1", y="MDS2"), color="black", alpha=0.3)
ggplot(NMDS, aes(x=NMDS1, y=NMDS2))+geom_point(data = NMDS, aes(x=MDS1, y=MDS2, fill = Treatment, shape = Treatment) ,size=1,alpha=0.75) + geom_polygon(data=hull.data, aes_string(x="MDS1", y="MDS2"), color="black", alpha=0.3)
ggplot(NMDS, aes(x=NMDS1, y=NMDS2))+geom_point(data = NMDS, aes(x=MDS1, y=MDS2, shape = Treatment) ,size=1,alpha=0.75) + geom_polygon(data=hull.data, aes_string(x="MDS1", y="MDS2"), color="black", alpha=0.3)
ggplot(NMDS, aes(x=NMDS1, y=NMDS2))+geom_point(data = NMDS, aes(x=MDS1, y=MDS2) ,size=1,alpha=0.75) + geom_polygon(data=hull.data, aes_string(x="MDS1", y="MDS2"), color="black", alpha=0.3)
ggplot()+geom_point(data = NMDS, aes(x=MDS1, y=MDS2, color=id) ,size=1,alpha=0.75) + geom_polygon(data=hull.data, aes_string(x="MDS1", y="MDS2"), color="id", alpha=0.3)
ggplot()+geom_point(data = NMDS, aes(x=MDS1, y=MDS2, color=id) ,size=1,alpha=0.75) + geom_polygon(data=hull.data, aes(x=MDS1, y=MDS2), color=id, alpha=0.3)
ggplot()+geom_point(data = NMDS, aes(x=MDS1, y=MDS2, color=id) ,size=1,alpha=0.75) + geom_polygon(data=hull.data, aes(x=MDS1, y=MDS2, color=id), alpha=0.3)
ls()
si<-data.frame(sample_data(data.ftype))
dim(si)
str(si)
str(si[, "foam.type")
dim(si)
names(si)
str(si[, "foam.type"]
)
head(rich)
head(alpha)
library(reshape2)
library(plyr)
alpha<-data.frame(si[, c("Observed", "Pielou", "foam.type")])
head(alpha)
dim(alpha)
si$foam.type<-factor(si$foam.type, levels=c("NF", "C", "F"))
sample_data(data.ftype)<-si
alpha<-data.frame(si[, c("Observed", "Pielou", "foam.type")])
dim(alpha)
alpha.melt<-melt(alpha, id.vars=c("foam.type"))
head(alpha.melt)
tail(alpha.melt)
ggplot(alpha.melt, aes(x=foam.type, y=value, color=foam.type)) + geom_boxplot()
ggplot(alpha.melt, aes(x=foam.type, y=value, color=foam.type)) + geom_boxplot()+ facet_grid(variable~., scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value, color=foam.type)) + geom_points() + geom_boxplot()+ facet_grid(variable~., scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value, color=foam.type)) + geom_point(aes(color=foam.type)) + geom_boxplot(aes(color = "balck"))+ facet_grid(variable~., scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value, color=foam.type)) + geom_point(aes(color=foam.type)) + geom_boxplot(color = "balck")+ facet_grid(variable~., scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value, color=foam.type)) + geom_point(aes(color=foam.type)) + geom_boxplot(aes(color = balck))+ facet_grid(variable~., scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter()+ facet_grid(variable~., scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(notch=T)+ facet_grid(variable~., scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot(notch=T)+geom_jitter()+ facet_grid(variable~., scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot(notch=T)+geom_jitter()+ facet_grid(~variable, scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot(notch=T)+geom_jitter()+ facet_grid(~variable, scale.y="free")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot(notch=T)+geom_jitter()+ facet_grid(~variable, scale="free_y")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(width=0.5)+ facet_grid(variable~., scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(width=0.5)+ facet_wrap(variable, ncol=1,  scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(width=0.5)+ facet_wrap(~variable, ncol=1,  scale="free")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot(aes(color=foam.type))+geom_jitter(width=0.5)+ facet_wrap(~variable, ncol=1,  scale="free")+ theme_classic() + scale_color_brewer(palette="Dark2")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ facet_wrap(~variable, ncol=1,  scale="free")+ theme_classic() + scale_color_brewer(palette="Dark2")
alpha.melt$foam.type<-gsub("NF", "No Foam", alpha.melt$foam.type)
head(alpha.melt)
alpha.melt$foam.type<-gsub("C", "Crust", alpha.melt$foam.type)
alpha.melt$foam.type<-gsub("F", "Foam", alpha.melt$foam.type)
head(alpha.melt)
alpha.melt$foam.type<-gsub("No Foamoam", "No Foam", alpha.melt$foam.type)
head(alpha.melt)
str(alpha.melt)
alpha.melt$variable<-gsub("Observed", "Observed Richness", alpha.melt$variable)
alpha.melt$variable<-gsub("Pielou", "Pielou Evenness", alpha.melt$variable)
str(alpha.melt)
alpha.melt$foam.type<-factor(alpha.melt$foam.type, levels=c("No Foam", "Crust", "Foam"))
alpha.melt$variable<-factor(alpha.melt$variable, levels=c("Observed Richness", "Pielou Evenness"))
head(alpha.melt)
str(alpha.melt)
saveRDS(alpha.melt, "RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/alpha_indices_melt.RDS")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ facet_wrap(~variable, ncol=1,  scale="free")+ theme_classic() + scale_color_brewer(palette="Dark2")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic()+ facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic()+ facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2")+ theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_bw()+ facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_bw()+ facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_bw()+ facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12)
)
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_bw()+ facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(fill="white"), strip.text.x=element_text(size = 14, face="bold"))
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.title=element_blank(), legend.text=element_text(size = 12, face="bold"))
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none")
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1)
plot(lm(si$Observed ~ si$foam.type))
plot(lm(log(si$Observed) ~ si$foam.type))
anova(lm(log(si$Observed) ~ si$foam.type))
TukeyHSD(aov(log(si$Observed) ~ si$foam.type))
plot(lm(log(si$pielou) ~ si$foam.type))
plot(lm(log(si$Pielou) ~ si$foam.type))
plot(lm((si$Pielou) ~ si$foam.type))
anova(lm((si$Pielou) ~ si$foam.type))
TukeyHSD(aov((si$Pielou) ~ si$foam.type))
ggplot(alpha.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1)
ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment),size=3,alpha=0.75) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=2, linetype=5)+theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
X1
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
ggplot.NMDS.ellipse(data.mds, alpha$foam.type, colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, alpha$foam.type, colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, alpha$foam.type,)
ggplot.NMDS.ellipse(data.mds, alpha$foam.type,"Dark2")
colors = getPalette(8)
colors
colors<-colors[1:3]
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, alpha$foam.type,colors)
head(si$foam.type)
si$foam.type<-gsub("\\bF\\b", "Foam", si$foam.type)
head(si$foam.type)
si$foam.type<-gsub("\\bNF\\b", "No Foam", si$foam.type)
si$foam.type<-gsub("\\bC\\b", "Crust", si$foam.type)
si$foam.type<-factor(si$foam.type, levels=c("No Foam", "Crust", "Foam"))
head(si$foam.type)
ggplot.NMDS.ellipse(data.mds, alpha$foam.type,colors)
ggplot.NMDS.ellipse(data.mds, si$foam.type,colors)
plot(si$LCFA ~ si$foam.depth)
plot(si$LCFA ~ si$surface.depth)
foam<-data.frame(si[, c("LCFA", "SCFA", "surface.depth","foam.type")])
head(foam)
ggplot(foam, aes(x=surface.depth, y=LCFA)) + geom_point() + facet_wrap(foam.type, ncol=1)
ggplot(foam, aes(x=surface.depth, y=LCFA)) + geom_point() + facet_wrap(foam.type~ , ncol=1)
ggplot(foam, aes(x=surface.depth, y=LCFA)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(foam, aes(x=surface.depth, y=LCFA)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
foam<-merge(foam, si[, c("id", "pH")
foam<-merge(foam, si[, c("SAMPLES", "pH")], by.x="row.names", by.y="SAMPLES")
dim(foam)
head(foam)
hist(foam$pH)
ggplot(foam, aes(x=pH, y=LCFA)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
ggplot(si, aes(x=Calcium1, y=LCFA)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
ggplot(si, aes(x=calcium1, y=LCFA)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
names(si)
ggplot(si, aes(x=calcium.1, y=LCFA)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
ggplot(si, aes(x=Calcium.1, y=LCFA)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
plot(si$SCFA ~ si$foam.type)
plot(lm(si$SCFA ~ si$foam.type))
plot(lm(log(si$SCFA) ~ si$foam.type))
plot(lm((si$SCFA) ~ si$foam.type))
anova(lm((si$SCFA) ~ si$foam.type))
TukeyHSD(aov(((si$SCFA) ~ si$foam.type)))
ggplot(si, aes(x=Calcium, y=LCFA)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
ggplot(si, aes(x=Oleic.Acid, y=LCFA)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
plot(si$LCFA ~ si$Adj.MPR)
plot(si$Oleic.Acid ~ si$Adj.MPR)
ggplot(si, aes(x=Oleic.Acid, y=Adj.MPR)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
ggplot(si, aes(x=Oleic.Acid, y=Adj.MPR)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(x=Oleic.Acid, y=MPR_Slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(x=Oleic.Acid, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(x=Oleic.Acid, y=MPR_VS1)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(x=Oleic.Acid, y=Adj.MPR)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(x=LCFA, y=Adj.MPR)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(x=Oleic.Acid, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
plot(si$MRP_slurry ~ si$foam.type)
plot(si$MPR_slurry ~ si$foam.type)
plot(lm(si$MPR_slurry ~ si$foam.type))
plot(lm(log(si$MPR_slurry) ~ si$foam.type))
plot(lm(log(si$MPR_slurry+1) ~ si$foam.type))
plot(lm(log(si$MPR_slurry+0.1) ~ si$foam.type))
anova(lm(log(si$MPR_slurry+0.1) ~ si$foam.type))
anova(lm(log(si$MPR_slurry+1) ~ si$foam.type))
plot(lm(log(si$MPR_slurry+1) ~ si$foam.type))
plot(lm(log(si$MPR_slurry+0.001) ~ si$foam.type))
anova(lm(log(si$MPR_slurry+0.001) ~ si$foam.type))
TukeyHSD(log(si$MPR_slurry+0.001) ~ si$foam.type))
TukeyHSD(log(si$MPR_slurry+0.001) ~ si$foam.type)
TukeyHSD(aov(log(si$MPR_slurry+0.001) ~ si$foam.type))
plot(log(si$MPR_slurry+0.001) ~ si$foam.type))
plot(log(si$MPR_slurry+0.001) ~ si$foam.type)
plot(si$MPR_slurry ~ si$foam.type)
plot(si$SCFA ~ si$foam.type)
ggplot(si, aes(x=Oleic.Acid, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=Oleic.Acid, x=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(x=Oleic.Acid, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
length(unique(si$Oleic.Acid))
(si$Oleic.Acid)
plot(si$Oleic.Acid ~ si$foam.type)
plot(lm(si$Oleic.Acid ~ si$foam.type))
plot(lm(log(si$Oleic.Acid+0.001) ~ si$foam.type))
plot((log(si$Oleic.Acid+0.001) ~ si$foam.type))
anova(lm(log(si$Oleic.Acid+0.001) ~ si$foam.type))
plot(lm(si$pH ~ si$foam.type))
plot(lm(log(si$pH) ~ si$foam.type))
anova(lm(log(si$pH) ~ si$foam.type))
plot(log(si$pH) ~ si$foam.type))
plot(log(si$pH) ~ si$foam.type)
plot((si$pH) ~ si$Oleic.Acid)
ggplot(si, aes(y=Oleic.Acid, x=pH)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=Calcium, x=pH)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=Calcium.1, x=pH)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=Magnesium, x=pH)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=sum(Oleic.Acid, LCFA), x=MPR_Slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=sum(Oleic.Acid, LCFA), x=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
names(si)
si$test<-sum(si$LCFA, si$Oleic.Acid)
head(si$test)
ggplot(si, aes(y=test, x=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
head(si[, c("LCFA", "Oleic.Acid")])
head(si[, c("LCFA", "Oleic.Acid")], n=50)
head(si$test, n=50)
si$test<-(si$LCFA + si$Oleic.Acid)
head(si$test, n=50)
ggplot(si, aes(y=test, x=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=test, x=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale = "free")
ggplot(si, aes(x=test, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale = "free")
si$test[si$test > 4000]
si[si$test > 4000, ]
test<-subset(si, test >4000)
dim(test)
test
ggplot(si, aes(x=test, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale = "free")+ theme(axis.title.x=element_text("Free LCFA + Oleic.Acid")
)
ggplot(si, aes(x=test, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale = "free")+ labs(x="Free LCFA + Oleic Acid")
names(si)
si$test<-rowSums(si[, c(25:30, 32:34)])
head(si$test, n=50)
ggplot(si, aes(x=test, y=pH)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale = "free")
ggplot(si, aes(x=test, y=pH)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=test, x=pH)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=test, x=Oleic.Acid)) + geom_point() + facet_wrap(~foam.type , ncol=1)
plot(si$test ~ si$foam.type)
si$test<-rowSums(si[, c("LCFA", "Oleic.Acid")])
ggplot(si, aes(x=test, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale = "free")+ labs(x="Free LCFA + Oleic Acid")
ggplot(si, aes(x=test, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)+ labs(x="Free LCFA + Oleic Acid")
ggplot(si[, -"S_437"], aes(x=test, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)+ labs(x="Free LCFA + Oleic Acid")
ggplot(si[, -c("S_437")], aes(x=test, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)+ labs(x="Free LCFA + Oleic Acid")
ggplot(si[-c("S_437"), ], aes(x=test, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)+ labs(x="Free LCFA + Oleic Acid")
ggplot(si[-"S_437", ], aes(x=test, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)+ labs(x="Free LCFA + Oleic Acid")
test1<-si[-"S_437", ]
test1<-si[-c("S_437"), ]
dim(si["S_437", ])
plot(si$Ammonium.N ~ si$foam.type)
ggplot(si, aes(x=test, y=surface.depth)) + geom_point() + facet_wrap(~foam.type , ncol=1)+ labs(x="Free LCFA + Oleic Acid")
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
ggplot(si, aes(x=foam.type, y=MPR_slurry)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer(pallete="Dark2")
ggplot(si, aes(x=foam.type, y=MPR_slurry)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer()
ggplot(si, aes(x=foam.type, y=MPR_slurry)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer("Dark2)
ggplot(si, aes(x=foam.type, y=MPR_slurry)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer("Dark2")
ggplot(si, aes(x=foam.type, y=MPR_slurry)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer(palette="Dark2")
mpr_sfca<-data.frame(si[, c("SFCA", "MPR_slurry", "foam.type")])
mpr_sfca<-data.frame(si[, c("SCFA", "MPR_slurry", "foam.type")])
)
ggplot(si, aes(x=SCFA, y=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
ggplot(si, aes(y=SCFA, x=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1, scale="free")
ggplot(si, aes(y=SCFA, x=MPR_slurry)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(x=foam.type, y=MPR_slurry)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer("Dark2")
ggplot(si, aes(x=foam.type, y=MPR_slurry)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer(palette="Dark2")
ggplot(si, aes(x=foam.type, y=Surface.Tension)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer(palette="Dark2")
names(si)
ggplot(si, aes(x=foam.type, y=Viscosity)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer(palette="Dark2")
ggplot(si, aes(x=foam.type, y=Phophorus)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer(palette="Dark2")
plot(lm(si$Phophorus ~ si$foam.type)
)
plot(lm(log(si$Phophorus) ~ si$foam.type))
anova(lm(log(si$Phophorus) ~ si$foam.type))
TukeyHSD(aov(log(si$Phophorus) ~ si$foam.type))
ggplot(si, aes(y=SCFA, x=Phophorus)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=LCFA, x=Phophorus)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=MPR_slurry, x=Phophorus)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=LCFA, x=Phophorus)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=SCFA, x=Phophorus)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=Corn, x=Phophorus)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ggplot(si, aes(y=Crude.Fiber, x=Phophorus)) + geom_point() + facet_wrap(~foam.type , ncol=1)
ls()
head(mpr_sfca)
mpr_sfca.melt<-melt(mpr_sfca, id.vars="foam.type")
head(mpr_sfca.melt)
ggplot(mpr_sfca.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer(palette="Dark2") + facet_wrap(~variable, ncol=1, scale="free")
ggplot(mpr_sfca.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, line
type='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + facet
_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title
.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, fac
e="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rec
t(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(lege
theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1)
ggplot(mpr_sfca.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1)
pdf("../Manuscript/figures_and_tables/fig_SCFA_mprslurry.pdf")
theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1)
ggplot(mpr_sfca.melt, aes(x=foam.type, y=value)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + facet_wrap(~variable, ncol=1,  scale="free")+ scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1)
dev.off()
ggplot(si, aes(x=foam.type, y=Copper)) + geom_boxplot()+geom_jitter(aes(color=foam.type), width=0.5) + scale_color_brewer(palette="Dark2") + facet_wrap(~variable, ncol=1, scale="free")
plot(si$Copper ~ si$foam.type)
plot(si$Iron ~ si$foam.type)
plot(si$Foam.Stability ~ si$foam.type)
plot(lm(si$Foam.Stability ~ si$foam.type))
plot(lm(log(si$Foam.Stability) ~ si$foam.type))
plot(lm(log(si$Foam.Stability+1) ~ si$foam.type))
plot(lm(log(si$Foam.Stability+0.0001) ~ si$foam.type))
anova(lm(log(si$Foam.Stability+0.0001) ~ si$foam.type))
library(DESeq2)
data.ftype
data.ftype.core
diagdds = phyloseq_to_deseq2(data.ftype, ~ Foaming.Status)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
head(sigtab)
head(res)
foaming_deseq_res<-res
sigtab$foam.type<-ifelse(sigtab$log2FoldChange > 0, "Foam", "Crust+No Foam")
x = tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$phylum = factor(as.character(sigtab$phylum), levels=names(x))
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
ggplot(sigtab, aes(y=genus, x=log2FoldChange, color=phylum)) + 
  geom_point(size=3) + 
  facet_grid(~foam.type, scale="free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
head(sigtab)
foaming_sig<-subset(sigtab, log2FoldChange > 0)
data_1e5.foaming.sig<-prune_taxa(taxa_names(data_1e5) %in% row.names(foaming_sig), data_1e5)
data_1e5.foaming.sig
ls()
plot_bar(t100_2fold_sig.phy, fill="surface.depth") + geom_bar(aes(color=surface.depth, fill=sur
face.depth), stat="identity", position="stack")+facet_grid(phylum~category, scales="free")
plot_bar(data_1e5.foaming.sig, fill="phylum") + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~category, scales="free")
plot_bar(data_1e5.foaming.sig, fill="phylum") + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")+facet_grid(phylum~foam.type, scales="free")
names(sample_data(data_1e5)
)
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type, scales="free")
data_1e5.foaming.sig<-prune_taxa(taxa_names(data_1e5.ftype) %in% row.names(foaming_sig), data_1e5.ftype)
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type, scales="free")
save.image("all_16s_physeq.RData")
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type) + theme(strip.axis.x=element_text(angle=90)) 
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type) + theme(strip.text.x=element_text(angle=90)) 
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type) + theme(strip.text.y=element_text(angle=90)) 
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type) + theme(strip.text.y=element_text(angle=0)) 
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type, scale.x="free") + theme(strip.text.y=element_text(angle=0)) 
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type, scale="free_x") + theme(strip.text.y=element_text(angle=0)) 
ggplot(si, aes(x=pH, y=Phophorus)) + geom_point() + facet_wrap(foam.type ~ ., ncol=1)
ggplot(si, aes(x=pH, y=Phophorus)) + geom_point() + facet_wrap(foam.type ~ , ncol=1)
ggplot(si, aes(x=pH, y=Phophorus)) + geom_point() + facet_wrap(foam.type, ncol=1)
ggplot(si, aes(x=pH, y=Phophorus)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(si, aes(x=Zinc.Oxide, y=Foam.Stability)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(si, aes(x=Total.Salt, y=Foam.Stability)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(si, aes(x=Sodium, y=Foam.Stability)) + geom_point() + facet_wrap(~foam.type, ncol=1)
plot(si$Zinc ~ si$foam.type)
plot(si$Copper ~ si$foam.type)
plot(si$General.Copper ~ si$foam.type)
si$General.Copper
si$General.Zinc
plot(si$General.Zinc ~ si$foam.type)
plot(sum(si$General.Zinc, si$Zinc.Oxide) ~ si$foam.type)
plot(sum(si[, c("General.Zinc", "Zinc.Oxide")] ~ si$foam.type)
plot(sum(si[, c("General.Zinc", "Zinc.Oxide")]) ~ si$foam.type)
si$diet.zinc<-si$General.Zinc + si$Zinc.Oxide
plot(si$diet.zinc ~ si$foam.type)
plot(lm(si$diet.zinc ~ si$foam.type))
plot(lm(log(si$diet.zinc) ~ si$foam.type))
plot(lm(log(si$diet.zinc+0.001) ~ si$foam.type))
plot(lm(log(si$diet.zinc+1) ~ si$foam.type))
plot(lm(log(si$diet.zinc+0.1) ~ si$foam.type))
plot(lm(log(1/si$diet.zinc) ~ si$foam.type))
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type, scale.x="free") + theme(strip.text.y=element_text(angle=0)) 
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type, scale="free") + theme(strip.text.y=element_text(angle=0)) 
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type, scale.x="free") + theme(strip.text.y=element_text(angle=0)) 
plot_bar(data_1e5.foaming.sig) + geom_bar(stat="identity", position="stack")+facet_grid(phylum ~ foam.type, scale="free_x") + theme(strip.text.y=element_text(angle=0)) 
data_1e5.foaming.sig
ls()
head(foaming_sig)
head(foaming_sig[order(-foaming_sig$log2FoldChange), ])
test<-ddply(foaming_sig, .(phylum, foam.type), summarize, avg=mean(log2FoldChange))
dim(foaming_sig)
unique(foaming_sig$foam.type)
test<-ddply(foaming_sig, .(phylum), summarize, avg=mean(log2FoldChange))
count(foaming_sig)
count(foaming_sig, "phylum")
test<-ddply(foaming_sig, .(phylum), summarize, avg=mean(log2FoldChange))
head(foaming_sig
)
test<-ddply(foaming_sig, .(phylum), summarise, sum=mean(log2FoldChange))
test<-ddply(foaming_sig, .(phylum), summarise, sum=sum(log2FoldChange))
test<-ddply(foaming_sig, .(phylum), summarise, sum=sum(log2FoldChange))
dim(foaming_sig)
length(foaming_sig$log2FoldChange)
save.image("all_16s_physeq.RData")
saveRDS(data_1e5.foaming.sig, "RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/data_1e5_foam_sig_physeq.RDS")
saveRDS(data.foaming.sig, "RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/data_foam_sig_physeq.RDS")
ls()
saveRDS(foaming_sig, "RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/foaming_significant_abundant_otus.RDS")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
test<-si[, c("id", "foam.type", "time_index")]
head(test)
library(plyr)
count(test, c("id", "foam.type"))
count(test, c("time_index", "foam.type"))
test<-count(test, c("time_index", "foam.type"))
plot(test$freq ~ test$time_index)
library(ggplot)
library(ggplot2)
ggplot(test, aes(x=time_index, y=freq)) + geom_point(aes(color=foam.type), size=2)
head(test)
names(si)
str(si$myear)
si$myear<-reorder(si$myear, si$time_index)
str(si$myear)
test<-si[, c("id", "foam.type", "time_index", "myear")]
head(test)
test1<-count(test, c("time_index", "foam.type"))
head(test1)
test1<-count(test, c("time_index", "foam.type", "myear"))
head(test1)
test2<-count(test, "myear")
test2<-count(test, "myear")
head(test2)
test1<-merge(test1, test1, "myear")
head(test1)
dim(test1)
dim(test)
test1<-count(test, c("time_index", "foam.type", "myear"))
head(test1)
dim(test1)
dim(test2)
test1<-merge(test1, test2, "myear")
dim(test1)
head(test1)
str(test1)
test1$myear<-reorder(test1$myear, test1$time_index)
str(test1)
test1<-test1[order(test1$time_index), ]
str(test1)
head(test1)
dim(test1)
colnames(test1)[4:5]<-c("count", "total
colnames(test1)[4:5]<-c("count", "total")
head(test1)
test1$percent<-100*test1$count/test1$total
head(test1)
ggplot(test1, aes(x=myear, y=percent)) + geom_point(aes(color=foam.type), size=3)
library(RColorBrewer)
ggplot(test1, aes(x=myear, y=percent)) + geom_point(aes(color=foam.type), size=3)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))
ggplot(test1, aes(x=myear, y=percent)) + geom_point(aes(color=foam.type), size=3)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+geom_smooth(method="lm")
ggplot(test1, aes(x=myear, y=percent)) + geom_point(aes(color=foam.type), size=3)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+geom_smooth(method="lm", fill=NA))
ggplot(test1, aes(x=myear, y=percent)) + geom_point(aes(color=foam.type), size=3)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+geom_smooth(method="lm", fill=NA))=
ggplot(test1, aes(x=myear, y=percent)) + geom_point(aes(color=foam.type), size=3)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+geom_smooth(method="lm", fill=NA)
ggplot(test1, aes(x=myear, y=percent, color=foam.type)) + geom_point(size=3)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+geom_smooth(method="lm", fill=NA)
ggplot(test1, aes(x=myear, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))
ggplot(test1, aes(x=myear, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_date()
ggplot(test1, aes(x=myear, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_date(format="%b-%Y")
ggplot(test1, aes(x=myear, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + xlim(
head(test1)
head(si[, 3:4])
head(as.Date(test1$myear))
head(as.Date(si$full_date))
head(as.Date(si$full_date, "%m/%Y"))
head(as.Date(si$month, "%m))
head(as.Date(si$month, "%m"))
head(as.Date(as.numeric(si$month), "%m"))
str(si$month)
test1$date<-test1$myear
unique(test1$myear)
head(si$full_date)
head(si$Long_Sample_ID)
si$date<-data.frame(do.call('rbind', strsplit(as.character(si$Long_Sample_ID), "-", fixed=T)))[, 3]
head(si$date)
head(as.Date(si$date))
head(as.Date(as.character(si$date), "%Y%m%d")
)
str(si$date)
head(as.character(si$date))
head(as.numberic(as.character(si$date)))
head(as.numeric(as.character(si$date)))
head(as.Date(as.numeric(as.character(si$date))))
head(as.Date(as.numeric(as.character(si$date)), "%Y%m%d"))
si$year<-gsub("\\b12\\b", "2012", si$year)
si$year<-gsub("\\b13\\b", "2013", si$year)
head(formatC(si$month, width=2, format="d", flag="0"))
si$month2<-(formatC(si$month, width=2, format="d", flag="0"))
head(si[, c("month2", "year")]
)
si$date<-paste(si$year, si$month2, sep="-")
head(si$date)
head(as.Date(si$date), )
str(si$date)
head(si$full_date)
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + xlim(myear)
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_discrete(limits=seq(1, 13))
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_discrete(limits=seq(1, 14))
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_discrete(limits=seq(0, 14))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_discrete(limits=seq(0, 15))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(limits=seq(0, 14))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(limits=c(0, 14))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=1)+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 14))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 14), labels=test1$myear)+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=test1$myear)+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(test1$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(test1$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)")
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(test1$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(test1$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
head(test1)
nf
nf<-subset(test1, foam.type="No Foam")
plot(nf$time_index ~ nf$percent)
plot(nf$percent ~ nf$time_index)
nf<-subset(test1, foam.type=="No Foam")
plot(nf$percent ~ nf$time_index)
plot(lm(nf$percent ~ nf$time_index))
plot(lm(log(nf$percent) ~ nf$time_index))
predict.lm(lm(nf$percent ~ nf$time_index))
anova(lm(log(nf$percent) ~ nf$time_index))
predict(aov(log(nf$percent) ~ nf$time_index))
lm.ft=fucntion(y, x)
lm.ft=function(y, x)
return(lm(y~x)
lm.ft=function(y, x) lm(y~x)
summary(lm.ft(nf$percent, nf$time_index))
crust<-subset(test1, foam.type=="Crust")
head(crust)
summary(lm.ft(crust$percent, crust$time_index))
foam<-subset(test1, foam.type=="Foam")
summary(lm.ft(foam$percent, foam$time_index))
plot(log(foam$percent) ~ foam$time_index)
plot((foam$percent) ~ foam$time_index)
hist(foam$percent)
summary(nls(foam$percent~(a + foam$time_index^b), start=list(a=0.01, b=0), trace=T))
summary(nls(foam$percent~(a + foam$time_index^b), start=list(a=0, b=0), trace=T))
head(test1)
ft_lm_overtime<-test1
rm(test1)
summary(nls(foam$percent~(a + foam$time_index^b))
)
nls(foam$percent ~ foam$time_index^a+b)
nls(percent ~I(time_index^a) +b, data=foam, start=list(a=0, b=0), trace=T)
summary(nls(percent ~I(time_index^a) +b, data=foam, start=list(a=0, b=0), trace=T))
plot(nls(percent ~I(time_index^a) +b, data=foam, start=list(a=0, b=0), trace=T))
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
ls()
rm(ggplot.NMDS.poly.ellipse)
rm(ggplot.NMDS.ellipse)
data_1e5.foaming.sig
test<-psmelt(data_1e5.foaming.sig)
dim(test)
dim(foam)
head(foam)
uniq_relations.otus<-read.delim("foaming_status_cc/spearman_hoeffding_cc/raw_counts/core_by_foamtype/ftype_core_unique_relationship_otu_list.txt", sep="\t")
head(uniq_relations.otus)
data_1e5.uniq.relations<-prune_taxa(taxa_names(data_1e5.ftype) %in% uniq_relations.otus$x, data_1e5.ftype)
data_1e5.uniq.relations
min(sample_sums(data_1e5.uniq.relations))
head(test[, 1:9])
names(test)
library(plyr)
test<-ddply(test, .(phylum, SAMPLES, Phophorus), summarize, sum=sum(Abundance))
dim(test)
head(test)
length(unique(test$phylum))
foam_sig.psmelt<-psmelt(data_1e5.foaming.sig)
ls()
dim(foaming_deseq_res)
data_1e5.foaming.sig
test<-ddply(foam_sig.psmelt, .(phylum, SAMPLES, Phophorus), summarize, sum=sum(Abundance))
test<-ddply(foam_sig.psmelt, .(phylum, SAMPLES, Phophorus, foam.type), summarize, sum=sum(Abundance))
head(test1)
head(test)
test1<-subset(test, phylum=="Firmicutes")
head(test1)
plot(test1$sum ~ test1$Phophorus)
ggplot(test1, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type))
test1<-subset(test, phylum=="Proteobacteria")
ggplot(test1, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type))
test1<-subset(test, phylum=="Firmicutes")
ggplot(test1, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type))
ggplot(test, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type))+facet_grid(~phylum, scale="free")
ggplot(test, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type))+facet_wrap(~phylum,ncol=4, scale="free")
library(Hmisc)
rcorr(test1)
head(test1)
rcorr(as.matrix(test1[, c(3, 5)]), type="spearman")
hoeffd(as.matrix(test1[, c(3, 5)]), type="spearman")
hoeffd(as.matrix(test1[, c(3, 5)]))
ggplot(test1, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type))
ggplot(test1, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type), size=3)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + labs(x = "Manure Phosphorus Content", y=bquote("Relative Abundance of Firmicutes ("*1 %.% 10^-5*")")) + theme(aspect.ratio=1)+theme(legend.title=element_blank()) + scale_color_brewer(palette="Dark2")+ theme(legend.text=element_text(size=14))+theme(legend.justification=c(1,0), legend.position=c(1,0)
)
ggplot(test1, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type), size=3)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + labs(x = "Manure Phosphorus Content", y=bquote("Relative Abundance of Firmicutes ("*1 %.% 10^-5*")")) + theme(aspect.ratio=1)+theme(legend.title=element_blank()) + scale_color_brewer(palette="Dark2")+ theme(legend.text=element_text(size=14))+theme(legend.justification=c(1,0), legend.position=c(0.8,0))
ggplot(test1, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type), size=3)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + labs(x = "Manure Phosphorus Content", y=bquote("Relative Abundance of Firmicutes ("*1 %.% 10^-5*")")) + theme(aspect.ratio=1)+theme(legend.title=element_blank()) + scale_color_brewer(palette="Dark2")+ theme(legend.text=element_text(size=14))+theme(legend.justification=c(1,0), legend.position=c(1,0.8))
head(test1)
test1$foam.type<-gsub("\\bF\\b", "Foam", test1$foam.type)
test1$foam.type<-gsub("\\bC\\b", "Crust", test1$foam.type)
test1$foam.type<-gsub("\\bNF\\b", "No Foam", test1$foam.type)
head(test1)
test1$foam.type<-factor(test1$foam.type, levels=c("No Foam", "Crust", "Foam"))
ggplot(test1, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type), size=3)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + labs(x = "Manure Phosphorus Content", y=bquote("Relative Abundance of Firmicutes ("*1 %.% 10^-5*")")) + theme(aspect.ratio=1)+theme(legend.title=element_blank()) + scale_color_brewer(palette="Dark2")+ theme(legend.text=element_text(size=14))+theme(legend.justification=c(1,0), legend.position=c(1,0.8))
test<-ddply(foam_sig.psmelt, .(phylum, SAMPLES, Phophorus, foam.type, DDGS), summarize, sum=sum(Abundance))
ggplot(si, aes(x=DDGS, y=Phophorus)) + geom_point(aes(color=foam.type))
ggplot(test1, aes(x=Phophorus, y=sum)) + geom_point(aes(color=foam.type), size=3)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + labs(x = "Manure Phosphorus Content", y=bquote("Relative Abundance of Firmicutes ("*1 %.% 10^-5*")")) + theme(aspect.ratio=1)+theme(legend.title=element_blank()) + scale_color_brewer(palette="Dark2")+ theme(legend.text=element_text(size=14))+theme(legend.justification=c(1,0), legend.position=c(1,0.8))
ls()
data_1e5.uniq.relations
test.glom<-tax_glom(data_1e5.uniq.relations, "genus")
test.glom
otu<-data.frame(otu_table(test.glom))
tax.glom<-data.frame(tax_table(test.glom))
head(otu[, 1:5])
head(tax.glom[, 1:5])
row.names(otu)<-tax.glom$genus
totu<-data.frame(t(otu))
head(totu[, 1:5])
dim(totu)
totu<-merge(totu, si, by.x="row.name", by.y="SAMPLES")
totu<-merge(totu, si, by.x="row.names", by.y="SAMPLES")
dim(totu)
names(totu[, 1:27])
plot(totu$Tissierella ~ totu$SCFA)
plot(totu$Tissie ~ totu$Acetic.Acid)
ggplot(totu, aes(x=SCFA, y=Tissierella, color=foam.type)) + geom_point()
ggplot(totu, aes(x=SCFA, y=Clostridium.XI, color=foam.type)) + geom_point()
ggplot(totu, aes(x=SCFA, y=Clostridium.sensu.stricto, color=foam.type)) + geom_point()
ggplot(totu, aes(x=Clostridium.XI, y=Clostridium.sensu.stricto, color=foam.type)) + geom_point()
ggplot(totu, aes(x=Tissierella, y=Clostridium.sensu.stricto, color=foam.type)) + geom_point()
ggplot(totu, aes(x=SCFA, y=Lactobacillus, color=foam.type)) + geom_point()
ggplot(totu, aes(x=Tissierella, y=Lactobacillus, color=foam.type)) + geom_point()
ggplot(totu, aes(x=Tissierella, y=unclassified_Ruminococcaceae, color=foam.type)) + geom_point()
ggplot(totu, aes(x=Tissierella, y=Clostridium.sensu.stricto, color=foam.type)) + geom_point()
ggplot(totu, aes(x=Tissierella, y=unclassified_Ruminococcaceae, color=foam.type)) + geom_point()
ggplot(totu, aes(x=Tissierella, y=SCFA, color=foam.type)) + geom_point()
ggplot(totu, aes(y=Tissierella, x=SCFA, color=foam.type)) + geom_point()
save.image("all_16s_physeq.RData")
savehistory("temp.R)
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
ls()
head(foam_sig.psmelt)
head(foaming_sig)
ls()
data.foaming.sig<-prune_taxa(taxa_names(data.ftype) %in% row.names(foaming_sig), data.ftype)
data.foaming.sig
data_1e5.foaming.sig
ls()
foam_sig_1e5.psmelt<-foam_sig.psmelt
foam_sig.psmelt<-psmelt(data.foaming.sig)
dim(foam_sig.psmelt)
test<-subset(foam_sig.psmelt, phylum=="Firmicutes")
head(test[, 1:5])
head(test[, 1:10])
length(unique(tes$OTU))
length(unique(test$OTU))
library(plyr)
test<-ddply(test, .(Samples, Phophorus), summarize, total=sum(Abundance))
test<-ddply(test, .(Sample, Phophorus), summarize, total=sum(Abundance))
dim(test)
head9test)
head(test)
test<-ddply(test, .(Sample, Phophorus, foam.type), summarize, total=sum(Abundance))
head(test)
library(Hmisc)
rcorr(test$total, test$Phophorus, "spearman")
dim(foam_sig.psmelt)
test<-subset(foam_sig.psmelt, phylum=="Firmicutes")
test<-ddply(test, .(Sample, Phophorus, foam.type), summarize, total=sum(Abundance))
head(test)
dim(test)
head(test[, 1:10])
foam_sig.phyla<-ddply(foam_sig.psmelt, .(SAMPLES, Phophorus, foam.type), summarize, total=sum(Abundance))
foam_sig.phyla<-ddply(foam_sig.psmelt, .(SAMPLES, Phophorus), summarize, total=sum(Abundance))
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(plyr)
ls()
foam_sig.psmelt.phyla<-ddply(foam_sig.psmelt, .(SAMPLES, foam.type, Phophorus, phylum), summarise, total=sum(Abundance))
head(foam_sig.psmelt.phyla)
test<-subset(foam_sig.psmelt.phyla, phylum=="Firmicutes")
dim(test)
library(ggplot2)
test1<-subset(test, foam.type=="F")
library(Hmisc)
rcorr(test1$total, test1$Phophorus, "spearman")
rcorr(test$total, test$Phophorus, "spearman")
detach("package:Hmisc", unload=T)
foam_sig_1e5.psmelt.phyla<-ddply(foam_sig_1e5.psmelt, .(SAMPLES, foam.type, Phophorus, phylum), summarise, total=sum(Abundance))
head(foam_sig_1e5.psmelt.phyla)
test<-subset(foam_sig_1e5.psmelt.phyla, phylum=="Firmicutes")
head(test)
library(ggplot2)
ggplot(test, aes(x=Phophorus, y=total, color=foam.type)) + geom_point()
ggplot(test, aes(x=Phophorus, y=total, color=foam.type)) + geom_boxplot()
ggplot(test, aes(x=Phophorus, y=total, color=foam.type)) + geom_point()
test<-subset(foam_sig.psmelt.phyla, phylum=="Firmicutes")
ggplot(test, aes(x=Phophorus, y=total, color=foam.type)) + geom_point()
plot(rank(test$total) ~ rank(test$Phophorus))
test$y<-rank(test$total)
test$x<-rank(test$Phophorus)
ggplot(test, aes(x=x, y=y, color=foam.type)) + geom_point()
head(foam_sig.psmelt.phyla)
library(Hmisc)
test1<-ddply(foam_sig.psmelt.phyla, .(phylum), function(x) rcorr(as.matrix(x[, c(3,5)], "spearman"))
)
test1<-ddply(foam_sig.psmelt.phyla, .(phylum), function(x) rcorr(as.matrix(x[, c(3,5)]), "spearman")
)
test1<-ddply(foam_sig.psmelt.phyla, .(phylum), function(x) rcorr(as.matrix(x[, c("Phophorus", "total")]), "spearman")
)
test1<-ddply(foam_sig.psmelt.phyla, .(phylum), function(x) rcorr(as.matrix(x[, c("Phophorus", "total")]), "spearman", na.omit=T))
str(foam_sig.psmelt.phyla
)
test1<-dlply(foam_sig.psmelt.phyla, .(phylum), summarize, rho= rcorr(as.matrix(x[, c("Phophorus", "total")]), "spearman"))
test1<-dlply(foam_sig.psmelt.phyla, .(phylum), function(x) rcorr(as.matrix(x[, c("Phophorus", "total")]), "spearman", na.omit=T))
test1<-dlply(foam_sig.psmelt.phyla, .(phylum), function(x) rcorr(x[[3]], x[[5]], "spearman"))
head(test1)
test1<-ddply(foam_sig.psmelt.phyla, .(phylum), summarize, rho= rcorr(Phophorus, total, "spearman")$rho, p=rcorr(Phophorus, total, "spearman")$P)
test1<-plyr::ddply(foam_sig.psmelt.phyla, .(phylum), summarize, rho= rcorr(Phophorus, total, "spearman")$rho, p=rcorr(Phophorus, total, "spearman")$P)
test1<-plyr::ddply(foam_sig.psmelt.phyla, .(phylum), summarize, rho= rcorr(Phophorus, total, "spearman")$rho)
test1<-plyr::ddply(foam_sig.psmelt.phyla, .(phylum), summarize, rho=rcorr(Phophorus, total, "spearman")$rho)
test1<-plyr::ddply(foam_sig.psmelt.phyla, .(phylum), summarise, rho=rcorr(Phophorus, total, "spearman")$rho)
test1<-dlply(foam_sig.psmelt.phyla, .(phylum), function(x) rcorr(x[[3]], x[[5]], "spearman"))
test1
test1<-dlply(foam_sig.psmelt.phyla, .(phylum), function(x) rcorr(x[[3]], x[[5]], "spearman")$P)
head(test1)
test1<-ddply(foam_sig.psmelt.phyla, .(phylum), function(x) rcorr(x[[3]], x[[5]], "spearman")$P)
head(test1)
test1<-ddply(foam_sig.psmelt.phyla, .(phylum), function(x) rcorr(x[, c("Phophorus", "total")], "spearman")$P)
head(foam_sig.psmelt.phyla)
test1<-ddply(foam_sig.psmelt.phyla, .(phylum), function(x) !is.na(rcorr(x[, c("Phophorus", "total")], "spearman")$P$x))
head(foam_sig.psmelt[, 1:10])
test<-foam_sig.psmelt[, c("OTU", "SAMPLES", "foam.type", "Phophorus", "Abundance")]
head(test)
test1<-ddply(test, .(OTU), summarize, p=rcorr(as.matrix(cbind(Phophorus, Abundance)), "spearman")$P)
test1<-ddply(test, .(OTU), summarise, p=rcorr(as.matrix(cbind(Phophorus, Abundance)), "spearman")$P)
for (i in unique(test$OTU)){
tryCatch({
temp <- subset(test, OTU == i)
results_sp<-rcorr(as.matrix(temp[, 4:5]),type="spearman")
    results_hd<-hoeffd(as.matrix(temp[, 4:5]))
        
    #make two seperate objects for p-value and correlation coefficients
    rhos<-results_sp$r
            ds<-results_hd$D
        ds_ps<-results_hd$P
        
        # going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
        sp_melt<-na.omit(melt(sp_ps))
        ds_melt<-na.omit(melt(ds_ps))
        
        #creating a qvalue (adjusted pvalue) based on FDR
        sp_melt$spearman_qval<-p.adjust(sp_melt$value, "fdr")
        ds_melt$hoeffding_qval<-p.adjust(ds_melt$value, "fdr")
        #       sp_melt$spearman_qval<-fdrtool(sp_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        #       ds_melt$hoeffding_qval<-fdrtool(ds_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        
        #making column names more relevant
        names(sp_melt)[3]<-"spearman_pval"
        names(ds_melt)[3]<-"hoeffding_pval"
        
        # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
        sp_sub<-subset(sp_melt, spearman_qval < 0.05)
        ds_sub<-subset(ds_melt, hoeffding_qval < 0.05)
        
        # now melting the rhos, note the similarity between ps_melt and rhos_melt
        rhos_melt<-na.omit(melt(rhos))
        ds_melt<-na.omit(melt(ds))
        names(rhos_melt)[3]<-"rho"
        names(ds_melt)[3]<-"D"
        
        #merging together 
        sp_merged<-merge(sp_sub,rhos_melt,by=c("Var1","Var2"))
        ds_merged<-merge(ds_sub, ds_melt,by=c("Var1","Var2"))
        merged<-merge(sp_merged, ds_merged, by=c("Var1", "Var2"))
        
        merged$id<-i
        combined_barn_cc<-rbind(combined_barn_cc, merged)
        print(paste("finished ",i,sep=""))
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
library(reshape2)
for (i in unique(test$OTU)){
tryCatch({
temp <- subset(test, OTU == i)
results_sp<-rcorr(as.matrix(temp[, 4:5]),type="spearman")
    results_hd<-hoeffd(as.matrix(temp[, 4:5]))
        
    #make two seperate objects for p-value and correlation coefficients
    rhos<-results_sp$r
            ds<-results_hd$D
        ds_ps<-results_hd$P
        
        # going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
        sp_melt<-na.omit(melt(sp_ps))
        ds_melt<-na.omit(melt(ds_ps))
        
        #creating a qvalue (adjusted pvalue) based on FDR
        sp_melt$spearman_qval<-p.adjust(sp_melt$value, "fdr")
        ds_melt$hoeffding_qval<-p.adjust(ds_melt$value, "fdr")
        #       sp_melt$spearman_qval<-fdrtool(sp_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        #       ds_melt$hoeffding_qval<-fdrtool(ds_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        
        #making column names more relevant
        names(sp_melt)[3]<-"spearman_pval"
        names(ds_melt)[3]<-"hoeffding_pval"
        
        # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
        sp_sub<-subset(sp_melt, spearman_qval < 0.05)
        ds_sub<-subset(ds_melt, hoeffding_qval < 0.05)
        
        # now melting the rhos, note the similarity between ps_melt and rhos_melt
        rhos_melt<-na.omit(melt(rhos))
        ds_melt<-na.omit(melt(ds))
        names(rhos_melt)[3]<-"rho"
        names(ds_melt)[3]<-"D"
        
        #merging together 
        sp_merged<-merge(sp_sub,rhos_melt,by=c("Var1","Var2"))
        ds_merged<-merge(ds_sub, ds_melt,by=c("Var1","Var2"))
        merged<-merge(sp_merged, ds_merged, by=c("Var1", "Var2"))
        
        merged$id<-i
        combined_barn_cc<-rbind(combined_barn_cc, merged)
        print(paste("finished ",i,sep=""))
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
for (i in unique(test$OTU)){
tryCatch({
temp <- subset(test, OTU == i)
results_sp<-rcorr(as.matrix(temp[, 4:5]),type="spearman")
    results_hd<-hoeffd(as.matrix(temp[, 4:5]))
        
    #make two seperate objects for p-value and correlation coefficients
    rhos<-results_sp$r
    sp_ps<-results_sp$P
    ds<-results_hd$D
    ds_ps<-results_hd$P
        
        # going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
        sp_melt<-na.omit(melt(sp_ps))
        ds_melt<-na.omit(melt(ds_ps))
        
        #creating a qvalue (adjusted pvalue) based on FDR
        sp_melt$spearman_qval<-p.adjust(sp_melt$value, "fdr")
        ds_melt$hoeffding_qval<-p.adjust(ds_melt$value, "fdr")
        #       sp_melt$spearman_qval<-fdrtool(sp_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        #       ds_melt$hoeffding_qval<-fdrtool(ds_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        
        #making column names more relevant
        names(sp_melt)[3]<-"spearman_pval"
        names(ds_melt)[3]<-"hoeffding_pval"
        
        # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
        sp_sub<-subset(sp_melt, spearman_qval < 0.05)
        ds_sub<-subset(ds_melt, hoeffding_qval < 0.05)
        
        # now melting the rhos, note the similarity between ps_melt and rhos_melt
        rhos_melt<-na.omit(melt(rhos))
        ds_melt<-na.omit(melt(ds))
        names(rhos_melt)[3]<-"rho"
        names(ds_melt)[3]<-"D"
        
        #merging together 
        sp_merged<-merge(sp_sub,rhos_melt,by=c("Var1","Var2"))
        ds_merged<-merge(ds_sub, ds_melt,by=c("Var1","Var2"))
        merged<-merge(sp_merged, ds_merged, by=c("Var1", "Var2"))
        
        merged$id<-i
        combined_barn_cc<-rbind(combined_barn_cc, merged)
        print(paste("finished ",i,sep=""))
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
combined_barn_cc<-data.frame()
for (i in unique(test$OTU)){
tryCatch({
temp <- subset(test, OTU == i)
results_sp<-rcorr(as.matrix(temp[, 4:5]),type="spearman")
    results_hd<-hoeffd(as.matrix(temp[, 4:5]))
        
    #make two seperate objects for p-value and correlation coefficients
    rhos<-results_sp$r
    sp_ps<-results_sp$P
    ds<-results_hd$D
    ds_ps<-results_hd$P
        
        # going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
        sp_melt<-na.omit(melt(sp_ps))
        ds_melt<-na.omit(melt(ds_ps))
        
        #creating a qvalue (adjusted pvalue) based on FDR
        sp_melt$spearman_qval<-p.adjust(sp_melt$value, "fdr")
        ds_melt$hoeffding_qval<-p.adjust(ds_melt$value, "fdr")
        #       sp_melt$spearman_qval<-fdrtool(sp_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        #       ds_melt$hoeffding_qval<-fdrtool(ds_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
        
        #making column names more relevant
        names(sp_melt)[3]<-"spearman_pval"
        names(ds_melt)[3]<-"hoeffding_pval"
        
        # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
        sp_sub<-subset(sp_melt, spearman_qval < 0.05)
        ds_sub<-subset(ds_melt, hoeffding_qval < 0.05)
        
        # now melting the rhos, note the similarity between ps_melt and rhos_melt
        rhos_melt<-na.omit(melt(rhos))
        ds_melt<-na.omit(melt(ds))
        names(rhos_melt)[3]<-"rho"
        names(ds_melt)[3]<-"D"
        
        #merging together 
        sp_merged<-merge(sp_sub,rhos_melt,by=c("Var1","Var2"))
        ds_merged<-merge(ds_sub, ds_melt,by=c("Var1","Var2"))
        merged<-merge(sp_merged, ds_merged, by=c("Var1", "Var2"))
        
        merged$id<-i
        combined_barn_cc<-rbind(combined_barn_cc, merged)
        print(paste("finished ",i,sep=""))
        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
head(combined_barn_cc)
foam_sig_otu_rcorr<-combined_barn_cc
rm(combined_barn_cc)
dim(foam_sig_otu_rcorr)
foam_sig_otu_rcorr<-foam_sig_otu_rcorr[!duplicated(foam_sig_otu_rcorr[, 3:9]), ]
dim(foam_sig_otu_rcorr)
(foam_sig_otu_rcorr)
test<-subset(foam_sig_otu_rcorr, rho>=0.2)
test
foam_sig_rcorr_0.2<-foam_sig_1e5.psmelt[foam_sig_1e5.psmelt$OTU %in% test$id, ]
dim(foam_sig_rcorr_0.2)
head(foam_sig_rcorr_0.2[, 1:6])
unique(foam_sig_rcorr_0.2$foam.type)
unique(foam_sig_rcorr_0.2$phylum)
unique(foam_sig_rcorr_0.2$genus)
ggplot(foam_sig_rcorr_0.2, aes(x=Phophorus, y=Abundance)) + geom_point(aes(color=foam.type))
ggplot(foam_sig_rcorr_0.2, aes(x=Phophorus, y=Abundance)) + geom_point(aes(color=foam.type))+facet_wrap(~OTU, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=time_index, y=Abundance)) + geom_point(aes(color=foam.type))+facet_wrap(~OTU, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=time_index, y=Phophorus)) + geom_point(aes(color=foam.type))+facet_wrap(~OTU, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=time_index, y=Phophorus)) + geom_boxplot(aes(color=foam.type))+facet_wrap(~OTU, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=myear, y=Phophorus)) + geom_boxplot(aes(color=foam.type))+facet_wrap(~OTU, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=myear, y=Abundance)) + geom_boxplot(aes(color=foam.type))+facet_wrap(~OTU, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=myear, y=Abundance)) + geom_boxplot(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=myear, y=Abundance)) + geom_boxplot(aes(color=OTU))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=time_index, y=Abundance)) + geom_point(aes(color=OTU))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=Phophorus, y=Abundance)) + geom_point(aes(color=OTU))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=Phophorus, y=Abundance)) + geom_line(aes(color=OTU))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=Phophorus, y=Abundance)) + geom_point(aes(color=OTU))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=Phophorus, y=Abundance)) + geom_point(aes(color=OTU))+facet_grid(myear~foam.type, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=Phophorus, y=Abundance)) + geom_point(aes(color=OTU))+facet_grid(myear~foam.type, scale="free")
ggplot(foam_sig_rcorr_0.2, aes(x=Phophorus, y=Abundance)) + geom_point(aes(color=OTU))+facet_grid(foam.type~myear, scale="free")
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
ls()
data_1e5.ftype
totu<-data.frame(t(data.frame(otu_table(data_1e5.ftype)))
)
si<-data.frame(sample_data(data_1e5.ftype))
dim(si)
head(totu[, 1:5])
head(si[, 1:6])
plot(totu$OTU_274 ~ si$Surface.Tension)
test<-cbind(totu, si)
library(ggplot2)
ggplot(test, aes(y=OTU_274, x=Surface.Tension)) + geom_point(aes(color=foam.type))
ggplot(test, aes(y=OTU_1665, x=Surface.Tension)) + geom_point(aes(color=foam.type))
ggplot(test, aes(y=OTU_1610, x=LCFA)) + geom_point(aes(color=foam.type))
ggplot(test, aes(y=OTU_1610, x=LCFA)) + geom_point(aes(color=foam.type))+facet_wrap(foam.type~., ncol=1, scale="free")
ggplot(test, aes(y=OTU_1610, x=LCFA)) + geom_point(aes(color=foam.type))+facet_wrap(foam.type~, ncol=1, scale="free")
ggplot(test, aes(y=OTU_1610, x=LCFA)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1, scale="free")
head(tax)
tax["OTU_1610", ]
tax[OTU_1610, ]
ggplot(test, aes(y=rank(OTU_1610), x=rank(LCFA))) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(test, aes(y=OTU_1610, x=LCFA)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1, scale="free")
tax$otu<-row.names(tax)
head(tax)
tax[otu=="OTU_1610", ]
dim(tax)
tax<-data.frame(tax_table(data_1e5.ftype))
tax["OTU_1610", ]
ggplot(test, aes(y=OTU_1425, x=SCFA)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(test, aes(y=OTU_6930, x=SCFA)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(test, aes(y=OTU_6930, x=SCFA)) + geom_point(aes(x=SCFA, y=OTU_1425))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(test, aes(y=OTU_1425, x=OTU_6930)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(test, aes(y=OTU_379, x=OTU_507)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1, scale="free")
test1<-test[, c("OTU_1425", "OTU_6930", "SCFA")]
head(test1)
library(reshape3)
library(reshape2)
?melt
test1.melt<-melt(test1, id.vars=c("SCFA"))
head(test1.melt)
ggplot(test1.melt, aes(x=SCFA, y=value)) + geom_point(aes(color=variable)) + facet_wrap(~foam.type, ncol=1, scale="free")
test1<-test[, c("OTU_1425", "OTU_6930", "SCFA", "foam.type")]
test1.melt<-melt(test1, id.vars=c("SCFA", "foam.type"))
ggplot(test1.melt, aes(x=SCFA, y=value)) + geom_point(aes(color=variable)) + facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(test, aes(y=OTU_379, x=Choice.White.Grease)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1, scale="free")
names(si)
ggplot(test, aes(y=OTU_379, x=General.Fat+Choice.White.Grease)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(test, aes(y=OTU_379, x=General.Fat+Choice.White.Grease)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(y=OTU_595, x=General.Fat+Choice.White.Grease)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(y=OTU_595, x=Choice.White.Grease)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(y=rank(OTU_595), x=General.Fat+Choice.White.Grease)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(y=OTU_379, x=General.Fat+Choice.White.Grease)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1)
si[, c("General.Fat", "Choice.White.Grease")]
plot(si$Choice.White.Grease ~ si$foam.type)
ggplot(test, aes(y=OTU_432, x=Acetic.Acid)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(y=OTU_4886, x=Acetic.Acid)) + geom_point(aes(color=foam.type))+facet_wrap(~foam.type, ncol=1)
test<-readRDS("foaming_status_cc/spearman_hoeffding_cc/raw_counts/core_by_foamtype_new/nf_strong_acetic_acid_relationships.txt")
head(test)
ls()
dim(test)
test<-cbind(totu, si)
ace_rela<-readRDS("foaming_status_cc/spearman_hoeffding_cc/raw_counts/core_by_foamtype_new/nf_strong_acetic_acid_relationships.txt")
head(test[, 1:5])
test1<-test[, colnames(test) %in% ace_rela$Var2]
dim(test1)
head(test1)
test1<-data.frame(test1, test[, c("foam.type", "Acetic.Acid")])
head(test1)
test1.melt<-melt(test1, id.vars=c("foam.type", "Acetic.Acid"))
head(test1.melt)
dim(test1.melt)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value)) + geom_point(aes(color=variable)) + facet_wrap(~foam.type, ncol=1, scale="free")
test1.melt<-subset(test1.melt, grepl("OTU_431|OTU_260", variable))
ggplot(test1.melt, aes(x=Acetic.Acid, y=value)) + geom_point(aes(color=variable)) + facet_wrap(~foam.type, ncol=1, scale="free")
test1.melt<-melt(test1, id.vars=c("foam.type", "Acetic.Acid"))
test1.melt<-subset(test1.melt, grepl("OTU_432|OTU_260", variable))
ggplot(test1.melt, aes(x=Acetic.Acid, y=value)) + geom_point(aes(color=variable)) + facet_wrap(~foam.type, ncol=1, scale="free")
ggplot(test1.melt, aes(x=Acetic.Acid, y=value)) + geom_point(aes(color=variable)) + facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value)) + geom_point(aes(color=variable)) + geom_line()+facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value)) + geom_point(aes(color=variable)) + geom_smooth(method=nlm)+facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value)) + geom_point(aes(color=variable)) + geom_line(aes(color=variable))+facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value)) + geom_point(aes(color=variable)) + geom_smooth(method=lm)+facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value)) + geom_point(aes(color=variable)) + stat_smooth(method="loess", formula= y ~ x) +facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value)) + geom_point(aes(color=variable)) + stat_smooth(aes(method="loess", formula= y ~ x)) +facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ x, size=1) +facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1)
head(test1.melt)
test1.melt$foam.type<-gsub("\\bF\\b", "Foam", test1.melt$foam.type)
test1.melt$foam.type<-gsub("\\bC\\b", "Crust", test1.melt$foam.type)
test1.melt$foam.type<-gsub("\\bNF\\b", "No Foam", test1.melt$foam.type)
head(test1.melt)
test1.melt$foam.type<-factor(test1.melt$foam.type, levels=c("No Foam", "Crust", "Foam"))
head(test1.melt)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1)
tax["OTU_432", ]
tax["OTU_260", ]
test1.melt$variable<-gsub("\\bOTU_260\\b", "OTU_260: Clostridium sensu stricto", test1.melt$variable)
test1.melt$variable<-gsub("\\bOTU_432\\b", "OTU_432: Lactobacillus", test1.melt$variable)
head(test1.melt)
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw()
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Acetic Acid", y = "Relative Abundance (1E-5)")
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Acetic Acid", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
head(test1.melt)
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
colors = getPalette(4)
colors
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Acetic Acid", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Acetic Acid", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top") + scale_color_manual(values=colors[1:2])
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Acetic Acid", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top") + scale_color_manual(values=])
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Acetic Acid", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top") + scale_color_manual(values=colors[c(3, 4)])
ggplot(test1.melt, aes(x=Acetic.Acid, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Acetic Acid", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top") + scale_color_manual(values=colors[c(4, 2)])
test1<-test[, c("OTU_379", "Choice.White.Grease")]
test1<-test[, c("OTU_379", "Choice.White.Grease", "OTU_663", "foam.type")]
head(test1)
test1.melt<-melt(test1, id.vars=c("foam.type", "Choice.White.Grease"))
head(test1.melt)
ggplot(test1.melt, aes(x=Choice.White.Grease, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
ggplot(test, aes(x=OTU_379, y=OTU_663)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=Choice.White.Grease, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
ggplot(test, aes(x=OTU_379, y=General.Fat)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(y=OTU_379, x=General.Fat)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(y=OTU_379, x=Choice.White.Grease)) + geom_point() + facet_wrap(~foam.type, ncol=1)
test1<-test[, c("OTU_379", "Choice.White.Grease", "OTU_274", "foam.type")]
test1.melt<-melt(test1, id.vars=c("foam.type", "Choice.White.Grease"))
ggplot(test1.melt, aes(x=Choice.White.Grease, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
test1<-test[, c("OTU_379", "Choice.White.Grease", "OTU_274", "foam.type", "General.Fat")]
test1.melt<-melt(test1, id.vars=c("foam.type", "Choice.White.Grease", "General.Fat"))
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
ggplot(test, aes(y=OTU_663, x=LCFA)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(y=OTU_663, x=Choice.White.Grease)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(x=OTU_379, y=OTU_663)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
ggplot(test1.melt, aes(x=Choice.White.Grease, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
ggplot(test, aes(x=OTU_379, y=OTU_274)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(x=OTU_379, y=OTU_261)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(x=OTU_379, y=OTU_663)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(x=OTU_379, y=OTU_1006)) + geom_point() + facet_wrap(~foam.type, ncol=1)
test1<-test[, c("OTU_379", "Choice.White.Grease", "OTU_274", "foam.type", "OTU_595")]
test1.melt<-melt(test1, id.vars=c("foam.type", "Choice.White.Grease", "General.Fat"))
test1<-test[, c("OTU_379", "Choice.White.Grease", "OTU_274", "foam.type", "OTU_595", "General.Fat")]
test1.melt<-melt(test1, id.vars=c("foam.type", "Choice.White.Grease", "General.Fat"))
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
ggplot(test, aes(x=OTU_379, y=OTU_595)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(x=OTU_379, y=OTU_1006)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
plot(si$General.Fat ~ si$foam.type)
plot(lm(si$General.Fat ~ si$foam.type))
anova(lm(si$General.Fat ~ si$foam.type))
TukeyHSD(aov(si$General.Fat ~ si$foam.type))
TukeyHSD(aov(si$General.Fat ~ si$DDGS))
plot(si$General.Fat ~ si$DDGS)
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
ggplot(test, aes(x=OTU_379, y=corn)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(x=OTU_379, y=Corn)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(x=OTU_379, y=Salt)) + geom_point() + facet_wrap(~foam.type, ncol=1)
ggplot(test, aes(y=OTU_379, x=Salt)) + geom_point() + facet_wrap(~foam.type, ncol=1)
names(si)
si[,c("Salt", "Sodium")]
test1<-test[, c("OTU_379", "Choice.White.Grease", "OTU_274", "foam.type", "General.Fat")]
test1.melt<-melt(test1, id.vars=c("foam.type", "Choice.White.Grease", "General.Fat"))
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top")
test1.melt$foam.type<-gsub("\\bF\\b", "Foam", test1.melt$foam.type)
test1.melt$foam.type<-gsub("\\bC\\b", "Crust", test1.melt$foam.type)
test1.melt$foam.type<-gsub("\\bNF\\b", "No Foam", test1.melt$foam.type)
head(test1.melt)
test1.melt$foam.type<-factor(test1.melt$foam.type, levels=c("No Foam", "Crust", "Foam"))
head(test1.melt)
tax[, c("OTU_379", "OTU_274")]
tax[c("OTU_379", "OTU_274"),]
test1.melt$variable<-gsub("\\bOTU_379\\b", "OTU_379: unclassified_Ruminococcaceae", test1.melt$variable)
test1.melt$variable<-gsub("\\bOTU_274\\b", "OTU_274: unclassified_Ruminococcaceae", test1.melt$variable)
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top") + scale_color_manual(values=color[1:2])
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top") + scale_color_manual(values=colors[1:2])
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top") + scale_color_manual(values=colors[1:3])
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top") + scale_color_manual(values=colors[c(2, 3)])
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Choice.White.Grease", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top") + scale_color_manual(values=colors[c(1, 4)])
ggplot(test1.melt, aes(x=General.Fat, y=value, color=variable)) + geom_point() + stat_smooth(method="gam", formula= y ~ s(x, k=3), size=1) +facet_wrap(~foam.type, ncol=1) + theme_bw() + labs(x="Fat", y = "Relative Abundance (1E-5)") + theme(legend.title=element_blank(), legend.position = "top") + scale_color_manual(values=colors[c(1, 4)])
ggplot(test, aes(y=OTU_379, x=OTU_663)) + geom_point() + facet_wrap(~foam.type, ncol=1)
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(BayesFactor)
ls()
dim(si)
names(si)
plot(lm(DDGS ~ foam.type * id, data=si))
anova(lm(DDGS ~ foam.type * id, data=si))
plot(lm(DDGS ~ foam.type * date, data=si))
plot(lm(DDGS ~ foam.type * time_index, data=si))
anova(lm(DDGS ~ foam.type * time_index, data=si))
summary(aov(DDGS ~ foam.type*id + Error(time_index * id)), data=si)
head(si$DDGS)
summary(aov( DDGS ~ foam.type * id + Error(time_index * id)), data=si)
summary(aov( DDGS ~ foam.type * id + Error(id/(foam.type * time_index)), data=si)
)
summary(aov( DDGS ~ foam.type * time_index + Error(id/(foam.type * time_index)), data=si)
)
summary(aov( DDGS ~ time_index * id + (1|fo), data=si))
summary(aov( DDGS ~ time_index * id + Error(foam.type), data=si))
summary(aov( DDGS ~ foam.type * id + Error(time.index), data=si))
summary(aov( DDGS ~ foam.type * id + Error(time_index), data=si))
summary(aov( DDGS ~ foam.type * time_index + Error(id/(foam.type*time_index)), data=si))
bf<-anovaBF(DDGS ~ foam.type * id + time_index, data=si, whichRandom="id")
summary(aov( DDGS ~ foam.type * id + Error(time_index), data=si))
summary(aov( DDGS ~ foam.type * id + Error(id), data=si))
summary(aov( DDGS ~ foam.type * id + Error(time_index), data=si))
?Error
(aov( DDGS ~ foam.type * id + Error(time_index), data=si))
bf=anovaBF( DDGS ~ foam.type * id + time_index, data=si))
bf=anovaBF( DDGS ~ foam.type * id + time_index, data=si)
test<-si[, c("DDGS", "time_index", "id", "foam.type")]
head(test)
dim(test)
test<-subset(test, !is.na(DDGS))
dim(test)
bf=anovaBF( DDGS ~ foam.type * id + time_index, data=si)
bf=anovaBF( DDGS ~ foam.type * id + time_index, data=test)
str(test)
test<-si[, c("DDGS", "myear", "id", "foam.type")]
head(test)
dim(test)
test<-subset(test, !is.na(DDGS))
bf=anovaBF( DDGS ~ foam.type * id + myear, data=test)
str(test)
test$myear<-as.factor(test$myear)
str(test)
bf=anovaBF( DDGS ~ foam.type * id + myear, data=test)
bf
plot(bf)
bf1=anovaBF( DDGS ~ foam.type * id + myear, data=test, whichRandom="myear")
bf1
plot(bf1)
plot(bf1)
bf1=anovaBF( DDGS ~ foam.type * id + id*myear, data=test, whichRandom="myear")
bf1
bf1=anovaBF( DDGS ~ foam.type * id + id, data=test, whichRandom="myear")
bf1
bf1=anovaBF( DDGS ~ foam.type + id*foam.type, data=test, whichRandom="myear")
bf1
bf1=anovaBF( DDGS ~ foam.type + id*myear, data=test, whichRandom="myear")
bf1
plot(bf)
bfWithoutmyear = lmBF( DDGS ~ foam.type, data =test)
bfWithoutmyear
bfWithoutmyear = lmBF( DDGS ~ id*myear, data =test)
bf1=anovaBF( DDGS ~ foam.type + id*myear, data=test, whichRandom="id*myear")
bf1
plot(bf1)
plot(bf)
plot(bf1)
plot(lm(General.Fat ~ foam.type, data=si))
plot(lm(DDGS ~ foam.type, data=si))
summary(test)
head(xtabs(~foam.type + myear+id, test))
head(xtabs(~foam.type + myear, test))
head(xtabs(~foam.type + id, test))
library(lme4)
lmer(DDGS ~ foam.type + (1|myear) + (1|id), test)
summary(lmer(DDGS ~ foam.type + (1|myear) + (1|id), test))
summary(lmer(DDGS ~ foam.type + (1|myear) + (1|id), test, REML=0))
data(Penicillin)
head(Penicillin)
lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin))
lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin)
summary(lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin))
anova(lmer(DDGS ~ foam.type + (1|myear) + (1|id), test, REML=0))
dim(test)
head(test)
str(test)
summary(anova(lmer(DDGS ~ foam.type + (1|myear) + (1|id), test, REML=0)))
anova(lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin))
aov(lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin))
aov(lmer(DDGS ~ foam.type + (1|myear) + (1|id), test, REML=0))
lmer(DDGS ~ foam.type + (1|myear) + (1|id)+ (1|Treated.Status), si, REML=0))
lmer(DDGS ~ foam.type + (1|myear) + (1|id)+ (1|Treated.Status), si, REML=0)
coef(lmer(DDGS ~ foam.type + (1|myear) + (1|id)+ (1|Treated.Status), si, REML=0))
coef(summary(lmer(DDGS ~ foam.type + (1|myear) + (1|id)+ (1|Treated.Status), si, REML=0)))
coefs<-data.frame(coef(summary(lmer(DDGS ~ foam.type + (1|myear) + (1|id)+ (1|Treated.Status), si, REML=0))))
coefs$p.z<- 2* (1-pnorm(abs(coefs$t.value)))
coefs
df
bf
bf1=anovaBF( DDGS ~ foam.type + id*myear, data=test, whichRandom="id+myear")
bf1=anovaBF( DDGS ~ foam.type + myear, data=test, whichRandom="myear")
bf1
bf1=anovaBF( DDGS ~ foam.type + id*myear, data=test, whichRandom=c("id", "myear")
)
bf1
bf1=anovaBF( DDGS ~ foam.type + id + myear, data=test, whichRandom=c("id", "myear"))
bf1
head(xtab(~Treated.Status + myear, data=si))
head(xtabs(~Treated.Status + myear, data=si))
plotlmer(DDGS ~ foam.type + (1|myear) + (1|id)+ (1|Treated.Status), si, REML=0))
plotlmer(DDGS ~ foam.type + (1|myear) + (1|id)+ (1|Treated.Status), si, REML=0)
plot(lmer(DDGS ~ foam.type + (1|myear) + (1|id)+ (1|Treated.Status), si, REML=0))
plot(lmer(DDGS ~ foam.type + (1|myear) + (1|id)+ (1|Treated.Status) + (1|myear:id:Treated.Status), si, REML=0))
test<-data.frame(si$Crude.Fiber/si$Crude.Protein)
head(test)
names(test)<-"F_P"
head(test)
test[, c("id", "myear", "Treated.Status")]<-si[,c("id", "myear", "Treated.Status")]
head(test)
test[, c("id", "myear", "Treated.Status", "foam.type")]<-si[,c("id", "myear", "Treated.Status", "foam.type")]
head(test)
test<-test[! is.na(test$F_P), ]
dim(test)
plot(F_P ~ foam.type, data=test)
bf=anovaBF( F_P ~ foam.type + id + myear + Treated.Status, data=test, whichRandom=c("id", "myear", "Treated.Sattus"))
bf=anovaBF( F_P ~ foam.type + id + myear + Treated.Status, data=test, whichRandom=c("id", "myear", "Treated.Status"))
str(test)
test[, c("myear", "Treated.Status")] <- factor(test[, c("myear", "Treated.Status")])
test[, c("myear", "Treated.Status")] <- as.factor(test[, c("myear", "Treated.Status")])
test[, c("myear", "Treated.Status")] <- lapply(test[, c("myear", "Treated.Status")], factor)
str(test)
bf=anovaBF( F_P ~ foam.type + id + myear + Treated.Status, data=test, whichRandom=c("id", "myear", "Treated.Status"))
bf
bfmain<-lmBF(F_P ~ foam.type, data=test)
bfmain
bfrandom=lmBF( F_P ~ id + myear + Treated.Status, data=test, whichRandom=c("id", "myear", "Treated.Status"))
bfrandom
bf[4] / bfrandom
bf2 = bfmain/bfrandom
bf2
bf[4]/bf2
bf
bfall=c(bf, bf2)
bfall
bfall[4]
bfall[1]
bfall[1]/bfall[2]
plot(NDF ~ TDF, data=si)
plot(NDF ~ ADF, data=si)
plot(TDF ~ ADF, data=si)
plot(TDF ~ ADF+NDF, data=si)
plot(TDF ~ ADF*NDF, data=si)
plot(TDF ~ Corn, data=si)
plot(Corn ~ foam.type, data=si)
si$foam.type<-factor(si$foam.type, levels=c("NF", "C", "F"))
head(si$foam.type)
plot(Corn ~ foam.type, data=si)
par(mfrow=c(2, 2))
plot(Corn ~ foam.type, data=si)
plot(DDGS ~ foam.type, data=si)
plot(General.Fat ~ foam.type, data=si)
plot(Soybean.Meal ~ foam.type, data=si)
si$A_NDF<-si$ADF+si$NDF
plot(A_NDF ~ TDF, data=si)
head(si[, c("ADF", "NDF", "TDF"))
head(si[, c("ADF", "NDF", "TDF")])
plot(pH ~ foam.type, data=si)
si$A_NDF<-NULL
bf
si$Treated.Status
ls()
savehistory("temp.R")
load("all_16s_physeq.RData")
library(BayesFactor)
si$foam.type<-factor(si$foam.type, levels=c("NF", "C", "F"))
head(si$foam.type)
si$myear<-reorder(si$myear, si$time_index)
str(si$myear)
temp<-si[, -c(1, 3:21, 110, 114, 116:120)]
for (i in colnames(temp)){
sink("sample_information_bayes_factors.txt", append = TRUE)
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
test<-temp[, c(i, "myear", "foam.type", "id")]
colnames(test)[1]<-"to_test"
test<-test[! is.na(test$to_test), ]
test1<-anovaBF(to_test ~ foam.type + myear + id, data=test, whichRandom=c("myear", "id"))
print(test1[1])
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
for (i in colnames(temp)){
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
test<-temp[, c(i, "myear", "foam.type", "id")]
colnames(test)[1]<-"to_test"
test<-test[! is.na(test$to_test), ]
test1<-anovaBF(to_test ~ foam.type + myear + id, data=test, whichRandom=c("myear", "id"))
print(test1[1])
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
test1
i
head(test)
names(temp)
test<-temp[, c(pH, "myear", "foam.type", "id")]
test<-temp[, c("pH", "myear", "foam.type", "id")]
test<-test[! is.na(test$to_test), ]
test<-test[! is.na(test$pH), ]
bf<-anovaBF(pH ~ foam.type + myear + id, data=test, whichRandom=c("myear", "id"))
head(test)
test<-temp[, c("pH", "myear", "foam.type", "id")]
head(test)
test<-test[! is.na(test$pH), ]
head(test)
bf<-anovaBF(pH ~ foam.type + myear + id, data=test, whichRandom=c("myear", "id"))
bf
bf1<-lmBF(pH ~ foam.type, data=test)
bf2<-lmBF(pH ~ myear +id, data=test, whichRandom=c("myear", "id"))
bf_random<-bf1/bf2
bf_random
bfall<-c(bf, bf_random)
bfall
bf
bf1
bf2
bfall
bf[4]/bf2
bf[1]/bf2
bf/bf2
bf
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
bfs<-read.delim("sample_information_bayes_factors_order_restricted_H.table", sep="\t", header=F)
head(bfs)
dim(bfs)
names(bfs)<-c("value", "hypo", "bf")
dim(bfs)
head(bfs)
plot(value ~ bf, data=bfs)
barplot(value ~ bf, data=bfs)
str(bfs)
plot( bf ~ value, data=bfs)
plot( bf ~ log(value), data=bfs)
plot( log(bf) ~ value, data=bfs)
names(si)
head(meta)
meta[meta$genus == "Lysine", ]
meta[meta$genus == "General.Lysine", ]
dim(meta)
meta<-read.delim("foaming_status_cc/meta_w_genus_information_all_otus.txt", sep="\t")
dim(meta)
meta[meta$genus == "General.Lysine", ]
head(meta)
meta[meta$otu == "General.Lysine", ]
meta[meta$otu == "Lysine", ]
major_diet<-c("Corn", "DDGS", "Soybean.Meal", "Wheat.Midds", "Lysine", "General.Fat")
head(dfs)
head(bfs)
bfs_toplot<-bfs[ bfs$value %in% major_diet, ]
dim(bfs_toplot)
bfs_toplot
bfs_toplot$type<-"Diets"
bfs_toplot
bfs_toplot$type<-"Major diet ingredients"
major_diet<-c("Crude.Fiber", "Crude.Protein", "Starch", "NDF", "ADF", "TDF")
test<-bfs[ bfs$value %in% major_diet, ]
dim(test)
test
test$type<-"As formulated components"
test
bfs_toplot<-rbind(bfs_toplot, test)
dim(bfs_toplot)
bfs_toplot$genres<-"Diet"
si$Foaming.Capacity
isu_measurements<-c("TS", "VS", "MPR_slurry", "Acetic.Acid", "SCFA", "LCFA", "Manure.Temperature", "Manure.Depth", "Viscosity", "Surface.Tension", "Foaming.Capacity", "surface.depth")
test<-bfs[ bfs$value %in% isu_measurements, ]
dim(test)
length(isu_measurements)
head(test)
test$type<-"ISU manure characteristics"
test$genres<-"Manure Characteristics"
bfs_toplot<-rbind(bfs_toplot, test)
dim(bfs_toplot)
tail(bfs_toplot)
names(si)
midwest_measurements<-c("pH", "Mositure", "TKN", "Organic.N", "Carbon", "Calsium", "Copper", "Potassium", "Manganese", "Magnesium", "Sodium", "Sulfur", "Zinc", "Iron", "Phophorus")
length(midwest_measurements)
test<-bfs[ bfs$value %in% midwest_measurements, ]
dim(test)
test$value
midwest_measurements[6]
midwest_measurements[6]<-"Calcium"
test<-bfs[ bfs$value %in% midwest_measurements, ]
dim(test)
test$type<-"Midwest Laboratories manure characteristics"
test$genres<-"Manure Characteristics"
bfs_toplot<-rbind(bfs_toplot, test)
dim(bfs_toplot
)
head(bfs_toplot
)
unique(bfs_toplot$type)
unique(bfs_toplot$hypo)
library(ggplot2)
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo~.) + coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo~.) + scale_y_continuous(trans=log2_trans())+ coord_flip()
library(scales)
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo~.) + scale_y_continuous(trans=log2_trans())+ coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo~.) + scale_y_log10()+ coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo~., scale="free") + scale_y_continuous(trans=log2_trans())+ coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo~., scale_x="free") + scale_y_continuous(trans=log2_trans())+ coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo~., scales="free_x") + scale_y_continuous(trans=log2_trans())+ coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_wrap(~hypo, ncol=1, scales="free_x") + scale_y_continuous(trans=log2_trans())+ coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_wrap(~hypo, ncol=1) + scale_y_continuous(trans=log2_trans())+ coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_wrap(hypo~, ncol=1) + scale_y_continuous(trans=log2_trans())+ coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_wrap(hypo~., ncol=1) + scale_y_continuous(trans=log2_trans())+ coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_wrap(hypo~.) + scale_y_log10()+ coord_flip()
ggplot(bfs_toplot, aes(x=value, y=bf)) + geom_point(aes(color=genres)) + facet_wrap(~hypo, ncol=2) + scale_y_log10()+ coord_flip()
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_wrap(~hypo, ncol=2) + scale_x_log10()
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres) + scale_x_log10()
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres, scales= "free_y") + scale_x_log10()
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10()
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_wrap(hypo ~ genres, scales= "free") + scale_x_log10()
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_wrap(hypo ~ genres, scales= "free_y") + scale_x_log10()
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_wrap(hypo ~ genres, scales= "free_y", ncol=5) + scale_x_log10()
unique(bfs_toplot$type)
unique(bfs_toplot$genres)
str(bfs_toplot$genres)
bfs_toplot$genres<-factor(bfs_toplot$genres, levels=c("D
bfs_toplot$genres<-factor(bfs_toplot$genres, levels=c("Diet", "Manure Characteristics"))
str(bfs_toplot$genres)
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_wrap(hypo ~ genres, scales= "free_y", ncol=5) + scale_x_log10()
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10()
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10(limits=c(0.01, 1e24)
)
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10(breaks=c(0.01, 1, 10, 100))
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10()+geom_vline(xintercept=1)
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10()+geom_vline(xintercept=10)
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10()+geom_vline(xintercept=5)
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10()+geom_vline(xintercept=1)
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres)) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10()+geom_vline(xintercept=c(1, 5, 10))
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres), size=3) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10()+geom_vline(xintercept=c(1, 5, 10))
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres), size=3) + facet_grid(hypo ~ genres, scales= "free") + scale_x_log10()+geom_vline(xintercept=c(3.2, 10, 100))
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres), size=3) + facet_grid(hypo ~ type, scales= "free") + scale_x_log10()+geom_vline(xintercept=c(3.2, 10, 100))
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres), size=3) + facet_grid(hypo ~ type, scales= "free") + scale_x_log10()+geom_vline(xintercept=c(1, 3.2, 10, 100))
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres), size=3) + facet_grid(hypo ~ type, scales= "free") + scale_x_log10()+geom_vline(xintercept=c(1, 3, 20, 150))
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres), size=3) + facet_grid(hypo ~ type, scales= "free") + scale_x_log10()+geom_vline(xintercept=c(1, 3, 20, 150)) + theme_bw() + theme(strip.text.y = element_text(size=12, angle=0))
unique(bfs_toplot$type)
bfs_toplot$type<-factor(bfs_toplot$type, levels=c("Major diet ingredients", "As formulated components", "Midwest Laboratories manure characteristics", "ISU manure characteristics"))
ggplot(bfs_toplot, aes(y=value, x=bf)) + geom_point(aes(color=genres), size=3) + facet_grid(hypo ~ type, scales= "free") + scale_x_log10()+geom_vline(xintercept=c(1, 3, 20, 150)) + theme_bw() + theme(strip.text.y = element_text(size=12, angle=0))
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
ls(sorted=F)
library(phyloseq)
library(bayesfactors)
library(bayesfactor)
library(BayesFactor)
names(si)
temp<-si[, c("myear", "foam.type", "id", "Observed", "Chao1", "Shannon", "InvSimpson", "Pielou")]
for (i in colnames(temp)){
sink("diversity_bayes_factors_order_restricted_H.txt", append = TRUE)
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
test<-temp[, c(i, "myear", "foam.type", "id")]
colnames(test)[1]<-"to_test"
test<-test[! is.na(test$to_test), ]
bf1<-anovaBF(to_test ~ foam.type + myear + id, data=test, whichRandom=c("myear", "id"))
avg<-ddply(test, .(foam.type), summarize, avg=mean(to_test))
avg<-avg[order(avg$avg), ]
samples = posterior(bf1, iterations = 10000)
print(paste(avg[3, 1], avg[2,1], avg[1,1], sep=">"))
consistent = (samples[, paste("foam.type", avg[3,1], sep="-")] > samples[, paste("foam.type", avg[2,1], sep="-")]) & (samples[, paste("foam.type", avg[2,1], sep="-")] > samples[, paste("foam.type", avg[1,1], sep="-")])
N_consistent = sum(consistent)
bf_restriction_against_full = (N_consistent / 10000) / (1 / 6)
bf_restriction_against_full
bf_full_against_null = as.vector(bf1)
bf_restriction_against_null = bf_restriction_against_full * bf_full_against_null
print(bf_restriction_against_null)
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
library(plyr)
temp<-si[, c("myear", "foam.type", "id", "Observed", "Chao1", "Shannon", "InvSimpson", "Pielou")]
for (i in colnames(temp)){
sink("diversity_bayes_factors_order_restricted_H.txt", append = TRUE)
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
test<-temp[, c(i, "myear", "foam.type", "id")]
colnames(test)[1]<-"to_test"
test<-test[! is.na(test$to_test), ]
bf1<-anovaBF(to_test ~ foam.type + myear + id, data=test, whichRandom=c("myear", "id"))
avg<-ddply(test, .(foam.type), summarize, avg=mean(to_test))
avg<-avg[order(avg$avg), ]
samples = posterior(bf1, iterations = 10000)
print(paste(avg[3, 1], avg[2,1], avg[1,1], sep=">"))
consistent = (samples[, paste("foam.type", avg[3,1], sep="-")] > samples[, paste("foam.type", avg[2,1], sep="-")]) & (samples[, paste("foam.type", avg[2,1], sep="-")] > samples[, paste("foam.type", avg[1,1], sep="-")])
N_consistent = sum(consistent)
bf_restriction_against_full = (N_consistent / 10000) / (1 / 6)
bf_restriction_against_full
bf_full_against_null = as.vector(bf1)
bf_restriction_against_null = bf_restriction_against_full * bf_full_against_null
print(bf_restriction_against_null)
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
head(si[, 1:5])
ls()
data(data_1e5.ftype)
(data_1e5.ftype)
dim(totu)
head(totu[, 1:5])
head(otu_table(data_1e5.ftype)[, 1:5])
library(vegan)
totu_pcoa<-cmdscale(totu, dist="bray")
?cmdscale
totu.bray<-vegdist(totu, method="bray")
totu_pcoa<-cmdscale(totu.bray, k = 5)
head(totu_pcoa)
plot(totu_pcoa[, 1], totu_pcoa[, 2], type="n")
plot(totu_pcoa[, 1], totu_pcoa[, 2], type="p")
summary(totu_pcoa)
?cmdscale
?cmdscale
packageVersion("vegan")
totu_pcoa<-cmdscale(totu.bray, k = 5, eig =T)
head(totu_pcoa)
totu_pcoa$GOF
totu_pcoa$GOF[1]
totu.bray.gof <- data.frame()
for (i in c(1:10)){
totu_pcoa<-cmdscale(totu.bray, k = 1, eig =T)
temp<-cbind(i, totu_pcoa$GOF[1], totu_pcoa$GOF[2])
totu.bray.gof<-rbind(totu.bray.gof, temp)
}
totu.bray.gof
totu.bray.gof <- data.frame()
for (i in c(1:10)){
totu_pcoa<-cmdscale(totu.bray, k = i, eig =T)
temp<-cbind(i, totu_pcoa$GOF[1], totu_pcoa$GOF[2])
totu.bray.gof<-rbind(totu.bray.gof, temp)
}
totu.bray.gof
head(eigenvals(totu.bray))
head(eigenvals(totu_pcoa))
(eigenvals(totu_pcoa))
(eigenvals(totu_pcoa))
(eigenvals(totu_pcoa))/sum(eigenvals(totu_pcoa))
totu_pcoa_var_explained<-(eigenvals(totu_pcoa))/sum(eigenvals(totu_pcoa))
head(totu_pcoa_var_explained)
head(totu_pcoa$points)
head(totu_pcoa$eig)
totu_pcoa.pc1<-data.frame(totu_pcoa$points[, 1])
head(totu_pcoa.pc1)
colnames(totu_pcoa.pc1)[1]<-"PC1"
head(totu_pcoa.pc1)
merge(totu_pcoa.pc1)
dim(totu_pcoa.pc1)
head(si[, c("foam.type", "SAMPLES")])
totu_pcoa.pc1<-merge(totu_pcoa.pc1, si[, c("foam.type", "SAMPLES")], by="row.names")
dim(totu_pcoa.pc1)
head(totu_pcoa.pc1)
boxplot(foam.type ~ PC1, data=totu_pcoa.pc1)
plot(foam.type ~ PC1, data=totu_pcoa.pc1)
plot(PC1 ~ foam.type, data=totu_pcoa.pc1)
plot(totu_pcoa, type="p")
plot(totu_pcoa)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2)
p2<-ggplot(si, aes(x=foam.type, y=Pielou)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold")) +labs(y="Pielou Evenness") + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2)
plot_grid(p1, p2, labels=c("A", "B"), ncol=2, nrow=1)
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2)
p2<-ggplot(si, aes(x=foam.type, y=Pielou)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) +labs(y="Pielou Evenness") + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2)
plot_grid(p1, p2, labels=c("A", "B"), ncol=2, nrow=1)
plot_grid(p1, p2, labels=c("A", "B"), ncol=2, nrow=1)
si$foam.type<-gsub("\\bNF\\b", "No Foam", si$foam.type)
head(si$foam.type)
si$foam.type<-gsub("\\bF\\b", "Foam", si$foam.type)
si$foam.type<-gsub("\\bC\\b", "Crust", si$foam.type)
head(si$foam.type)
unique(si$foam.type)
str(si$foam.type)
si$foam.type<-factor(si$foam.type, levels=c("No Foam", "Crust", "Foam"))
str(si$foam.type)
head(si$foam.type)
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(x=2.5, y=2500, label="Richness hypothesis: No Foam > Crust > Foam")
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(x=2.5, y=2500, label=as.character(Richness hypothesis: No Foam > Crust > Foam))
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(x=2.5, y=2500, label=as.character(as.expression(Richness hypothesis: No Foam > Crust > Foam)))
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(x=2.5, y=2500, label=as.character(as.expression(Richness hypothesis: No Foam > Crust > Foam)))
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(x=2.5, y=2500, label=as.character("Richness hypothesis: No Foam > Crust > Foam"))
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(x=2.5, y=2500, label=as.character(as.expression("Richness hypothesis: No Foam > Crust > Foam")))
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(x=2.5, y=2500, label=as.character(as.expression("Richness hypothesis:" "No Foam" > "Crust" > "Foam")))
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(x=2.5, y=2500, label=as.character(as.expression("Richness hypothesis: No Foam" > "Crust" > "Foam")))
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(x=2.5, y=2500, label=as.character(as.expression("Richness hypothesis: No Foam > Crust > Foam")), parse=T, inherit.aes=F)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(x=2.5, y=2500, label=as.character("Richness hypothesis: No Foam > Crust > Foam"), parse=T, inherit.aes=F)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(aes(x=2.5, y=2500, label=as.character("Richness hypothesis: No Foam > Crust > Foam")), parse=T, inherit.aes=F)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(aes(x=2.5, y=2500, label="Richness hypothesis: No Foam > Crust > Foam"))
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(aes(x=2.5, y=2500, label="Richness hypothesis: No Foam > Crust > Foam"))
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(data=df1,aes(x = x, y = y,label=label), parse = TRUE, inherit.aes=FALSE)
p1
df1<-data.frame(x="2.5", y="2500", label=as.character("Richness hypothesis: No Foam > Crust > Foam"))
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(data=df1,aes(x = x, y = y,label=label), parse = TRUE, inherit.aes=FALSE)
p1
df1<-data.frame(x="Crust", y="2500", label=as.character("Richness hypothesis: No Foam > Crust > Foam"))
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(data=df1,aes(x = x, y = y,label=label), parse = TRUE, inherit.aes=FALSE)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + geom_text(aes(x=2.5, y=2500, label="Richness hypothesis: No Foam > Crust > Foam"))
p1
ggsave(p1, file="test.pdf")
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate(x=2.5, y=2500, label="Richness hypothesis: No Foam > Crust > Foam")
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Richness hypothesis: No Foam > Crust > Foam")
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam") + annotate("text", x=2.5, y=2300, label="Bayes Factor: 7070.722")
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam") + annotate("text", x=2.5, y=2400, label="Bayes Factor: 7071")
p2<-ggplot(si, aes(x=foam.type, y=Pielou)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) +labs(y="Pielou Evenness") + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2, y=0.8, label="Hypothesis: No Foam < Crust < Foam") + annotate("text", x=2, y=0.75, label="Bayes Factor: 153")
plot_grid(p1, p2, labels=c("A", "B"), ncol=2, nrow=1)
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam") + annotate("text", x=2.5, y=2400, label="Bayes Factor: 7071")
p2<-ggplot(si, aes(x=foam.type, y=Pielou)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) +labs(y="Pielou Evenness") + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2, y=0.77, label="Hypothesis: No Foam < Crust < Foam") + annotate("text", x=2, y=0.76, label="Bayes Factor: 153")
plot_grid(p1, p2, labels=c("A", "B"), ncol=2, nrow=1)
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam") + annotate("text", x=2.5, y=2400, label="Bayes Factor: 7071", size = 12)
p2<-ggplot(si, aes(x=foam.type, y=Pielou)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) +labs(y="Pielou Evenness") + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2, y=0.77, label="Hypothesis: No Foam < Crust < Foam") + annotate("text", x=2, y=0.76, label="Bayes Factor: 153", size = 12)
plot_grid(p1, p2, labels=c("A", "B"), ncol=2, nrow=1)
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", xmin=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam", size = 8) + annotate("text", xmin=2.5, y=2400, label="Bayes Factor: 7071", size = 8)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam", size = 8) + annotate("text", x=2.5, y=2400, label="Bayes Factor: 7071", size = 8)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam", size = 6) + annotate("text", x=2.5, y=2400, label="Bayes Factor: 7071", size = 6)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam", size = 6) + annotate("text", xend=2.5, y=2400, label="Bayes Factor: 7071", size = 6)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam", size = 6) + annotate("text", x=2.5, y=2400, label="Bayes Factor: 7071", size = 6, hjust=0)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam", size = 6) + annotate("text", x=2.5, y=2400, label="Bayes Factor: 7071", size = 6, hjust=1)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=2.5, y=2500, label="Hypothesis: No Foam > Crust > Foam", size = 5, hjust=1) + annotate("text", x=2.5, y=2400, label="Bayes Factor: 7071", size = 5, hjust=1)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=1.5, y=2500, label="Hypothesis: No Foam > Crust > Foam", size = 5, hjust=1) + annotate("text", x=1.5, y=2400, label="Bayes Factor: 7071", size = 5, hjust=0)
p1
p1<-ggplot(si, aes(x=foam.type, y=Observed)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) + labs(y="Observed Richness")+ theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=1.5, y=2500, label="Hypothesis: No Foam > Crust > Foam", size = 5, hjust=0) + annotate("text", x=1.5, y=2400, label="Bayes Factor: 7071", size = 5, hjust=0)
p1
p2<-ggplot(si, aes(x=foam.type, y=Pielou)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold")) + theme(axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face = "bold")) +labs(y="Pielou Evenness") + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=1, y=0.77, label="Hypothesis: No Foam < Crust < Foam", size=5, hjust=0) + annotate("text", x=1, y=0.76, label="Bayes Factor: 153", size = 5, hjust=0)
p2
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
p1<-ggplot(si, aes(x=foam.type, y=DDGS)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2)
p1
head(bfs_toplot)
p1<-ggplot(si, aes(x=foam.type, y=DDGS)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold") + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=1.5, y=65, label="Hypothesis: No Foam < Crust < Foam", size = 5, hjust=0) + annotate("text", x=1.5, y=63, label="Bayes Factor: 40", size = 5, hjust=0)
p1<-ggplot(si, aes(x=foam.type, y=DDGS)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=1.5, y=65, label="Hypothesis: No Foam < Crust < Foam", size = 5, hjust=0) + annotate("text", x=1.5, y=63, label="Bayes Factor: 40", size = 5, hjust=0)
p1
p1<-ggplot(si, aes(x=foam.type, y=DDGS)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=1, y=65, label="Hypothesis: No Foam < Crust < Foam", size = 5, hjust=0) + annotate("text", x=1, y=63, label="Bayes Factor: 40", size = 5, hjust=0)
p1
p1<-ggplot(si, aes(x=foam.type, y=DDGS)) + geom_boxplot()+geom_jitter(aes(color=foam.type),width=0.5)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + labs(y="DDGS (%)") + theme(strip.background=element_rect(color = "white", fill="white"), strip.text.x=element_text(size = 14, face="bold"))+theme(legend.position="none") +theme(aspect.ratio=1.2) + annotate("text", x=1, y=65, label="Hypothesis: No Foam < Crust < Foam", size = 5, hjust=0) + annotate("text", x=1, y=63, label="Bayes Factor: 40", size = 5, hjust=0)
p1
names(si)
head(si[, c("General.Fat", "Choice.White.Grease")])
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(ggthemes)
library(ggplot2)
library(RColorBrewer)
bank_slopes(si$Observed, si$foam.type)
bank_slopes(foam.type, Observed, data = si)
bank_slopes(si$foam.type, si$Observed)
head(test1)
test<-si[, c("id", "foam.type", "time_index", "myear")]
head(test)
 test1<-count(test, c("time_index", "foam.type"))
library(plyr)
 test1<-count(test, c("time_index", "foam.type"))
test2<-count(test, "myear")
head(test1)
head(test2)
test1<-count(test, c("time_index", "foam.type", "myear"))
head(test1)
test1<-merge(test1, test2, "myear")
head(test1)
test1$percent<-100*test1$freq.x/test1$freq.y
head(test1)
ggplot(test1, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(test1$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ft_collections<-test1
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_text(size=14, face="bold", angle = 45, vjust =1, hjust =1), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+stat_smooth(method="lm", formula="y~poly(x,2)", se=FALSE)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
lm.ft=function(y, x) lm(y~poly(x, 2))
head(nf)
summary(lm.ft(nf$percent, nf$time_index))
crust<-subset(test1, foam.type=="Crust")
head(crust)
summary(lm.ft(crust$percent, crust$time_index))
no_f<-subset(ft_collections, foam.type != "Foam")
head(no_f)
ggplot() + geom_point(data = no_f, aes(x=time_index, y=percent, color=foam.type), size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot() + geom_point(data = no_f, aes(x=time_index, y=percent, color=foam.type), size=3)+geom_smooth(data = no_f, method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot() + geom_point(data = no_f, aes(x=time_index, y=percent, color=foam.type), size=3)+geom_smooth(data = no_f, aes(color = foam.type), method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(group = Crust, method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+geom_smooth(aes(group = Crust), method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+stat_smooth(aes(group = Crust), method="lm", fill=NA)+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+stat_smooth(aes(group = Crust), method="lm", formula = "y ~ x", se = FALSE)
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+stat_smooth(aes(group = "Crust"), method="lm", formula = "y ~ x", se = FALSE)
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+stat_smooth(aes_string(group = "Crust"), method="lm", formula = "y ~ x", se = FALSE)
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+stat_smooth(aes(group = ft_collections$Crust), method="lm", formula = "y ~ x", se = FALSE)
ggplot(ft_collections, aes(x=time_index, y=percent, color=foam.type)) + geom_point(size=3)+stat_smooth(, method="lm", formula = "y ~ x", se = FALSE)
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_point(aes(color=foam.type), size=3)+stat_smooth(aes(group = Crust), method="lm", formula = "y ~ x", se = FALSE)
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE)
 display.brewer.all()
scale_color_brewer(palette="Dark2")
colors
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
head(f)
head(foam)
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
install.packages("ggthemes")
library(ggplot2)
library(RColorBrewer)
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe() + theme_tufte() + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
library(ggtheme)
library(ggthemes)
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe() + theme_tufte() + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe() + theme_tufte() + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3")
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe() + theme_tufte() + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe(lwd=2) + theme_tufte() + theme(axis.ticks.x = element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5), axis.ticks.length = unit(0.3, "cm"))+ scale_y_continuous(breaks = extended_range_breaks()(ft_collections$percent) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe(lwd=2) + theme_tufte() + theme(axis.ticks.x = element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5), axis.ticks.length = unit(0.3, "cm"))+ scale_y_continuous(breaks = extended_range_breaks()(ft_collections$percent)) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe(lwd=2) + theme_tufte() + theme(axis.ticks.x = element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5), axis.ticks.length = unit(0.3, "cm"))+ scale_y_continuous(breaks = extended_range_breaks()(round(ft_collections$percent, 2))) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe(lwd=2) + theme_tufte() + theme(axis.ticks.x = element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5), axis.ticks.length = unit(0.3, "cm"))+ scale_y_continuous(breaks = extended_range_breaks()(round(ft_collections$percent, 2))) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=30, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=30, face="bold"), axis.text.y=element_text(size=26, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=26, face="bold")) + scale_color_brewer(palette="Dark2")
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe(lwd=2) + theme_tufte() + theme(axis.ticks.x = element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5), axis.ticks.length = unit(0.3, "cm"))+ scale_y_continuous(breaks = extended_range_breaks()(round(ft_collections$percent, 2))) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2")
bank_slopes(crust$time_index, crust$percent)
bank_slopes(ft_collections$time_index, ft_collections$percent)
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe(lwd=2) + theme_tufte() + theme(axis.ticks.x = element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5), axis.ticks.length = unit(0.3, "cm"))+ scale_y_continuous(breaks = extended_range_breaks()(round(ft_collections$percent, 2))) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2") + theme(aspect.ratio = 0.85
)
ggplot(ft_collections, aes(x=time_index, y=percent)) + geom_rangeframe(lwd=2) + theme_tufte() + theme(axis.ticks.x = element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5), axis.ticks.length = unit(0.3, "cm"))+ scale_y_continuous(breaks = extended_range_breaks()(round(ft_collections$percent, 2))) + geom_point(aes(color=foam.type), size=3)+stat_smooth(data=crust, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#D95F02")+stat_smooth(data=nf, aes(x=time_index, y=percent), method="lm", formula = "y ~ x", se = FALSE, color = "#1B9E77") + stat_smooth(data=foam, aes(x=time_index, y=percent), method="lm", formula="y~poly(x,2)", se=FALSE, color = "#7570B3") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + scale_x_continuous(breaks=seq(1, 13), labels=levels(ft_collections$myear))+theme(legend.position=c(1,0), legend.justification=c(1,0))+ theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), axis.title.x=element_blank()) + labs(y="Percent of Samples Collected (%)") + theme(legend.text=element_text(size=12, face="bold")) + scale_color_brewer(palette="Dark2") + theme(aspect.ratio = 1)
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
lm.ft
summary(lm.ft(foam$percent, foam$time_index))
ls()
ls(sorted=F)
library(phyloseq)
data.ftype.core
head(sigtab)
length(unique(sigtab$genus))
dim(sigtab)
head(core_otu)
head(diagdds)
head(sample_data(data.ftype.core)$foam.type)
test<-subset_samples(data.ftype.core, foam.type != "F")
library(DESeq2)
diagdds = phyloseq_to_deseq2(test, ~ foam.type)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
head(res)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data.ftype.core)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
names(sample_data(data.ftype.core))
head(sample_data(data.ftype.core)$group.y)
names(si)
head(si$group.x)
head(si$group.y)
unique(si$group.y)
data_1e5.ftype.core
plot_bar(data_1e5.ftype.core, "genus", facet_grid=foam.type ~)
plot_bar(data_1e5.ftype.core, "genus", facet_grid=foam.type ~.)
plot_bar(data_1e5.ftype.core, "genus", fill = "foam.type")
plot_bar(data_1e5.ftype.core, "otu", facet_grid=foam.type ~.)
plot_bar(data_1e5.ftype.core, "tax", facet_grid=foam.type ~.)
head(tax)
plot_bar(data_1e5.ftype.core, "row.names", facet_grid=foam.type ~.)
plot_bar(data_1e5.ftype.core, "genus", fill = "foam.type") + geom_bar(stat = "identity", position = "dodge")
library(ggplot2)
plot_bar(data_1e5.ftype.core, "genus", fill = "foam.type") + geom_bar(stat = "identity", position = "dodge")
plot_bar(data_1e5.ftype.core, "genus") + geom_bar(aes(fill= foam.type) +stat = "identity", position = "dodge")
plot_bar(data_1e5.ftype.core, "genus") + geom_bar(aes(fill= foam.type), stat = "identity", position = "dodge")
ls()
plot_bar()
plot_bar
data_1e5.ftype.core.psmelt<-psmelt(data_1e5.ftype.core)
head(data_1e5.ftype.core.psmelt)
ggplot(data_1e5.ftype.core.psmelt, aes(x=genus, y=Abundance)) + geom_bar(aes(fill=foam.type), stat="identity", position="dodge")
data_1e5.ftype.core.psmelt$foam.type<-gsub("\\bF\\b", "Foam", data_1e5.ftype.core.psmelt$foam.type)
data_1e5.ftype.core.psmelt$foam.type<-gsub("\\bC\\b", "Crust", data_1e5.ftype.core.psmelt$foam.type)
data_1e5.ftype.core.psmelt$foam.type<-gsub("\\bNF\\b", "No Foam", data_1e5.ftype.core.psmelt$foam.type)
head(data_1e5.ftype.core.psmelt$foam.type)
unique(data_1e5.ftype.core.psmelt$foam.type)
data_1e5.ftype.core.psmelt$foam.type<-as.factor(data_1e5.ftype.core.psmelt$foam.type, levels=c("No Foam", "Crust", "Foam"))
data_1e5.ftype.core.psmelt$foam.type<-factor(data_1e5.ftype.core.psmelt$foam.type, levels=c("No Foam", "Crust", "Foam"))
head(data_1e5.ftype.core.psmelt$foam.type)
ggplot(data_1e5.ftype.core.psmelt, aes(x=genus, y=Abundance)) + geom_bar(aes(fill=foam.type), stat="identity", position="dodge")
ggplot(data_1e5.ftype.core.psmelt, aes(x=genus, y=Abundance)) + geom_point(aes(fill=foam.type))
ggplot(data_1e5.ftype.core.psmelt, aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))
ggplot(data_1e5.ftype.core.psmelt, aes(x=OTU, y=Abundance)) + geom_bar(aes(fill=foam.type), stat="identity", position="dodge") + facet_grid(foam.type ~.)
length(unique(data_1e5.ftype.core.psmelt$OTU))
length(unique(data_1e5.ftype.core.psmelt$genus))
ggplot(data_1e5.ftype.core.psmelt, aes(x=OTU, y=Abundance)) + geom_bar(aes(fill=foam.type), stat="identity", position="dodge") + facet_grid(foam.type ~ genus)
ggplot(data_1e5.ftype.core.psmelt, aes(x=OTU, y=Abundance)) + geom_bar(aes(fill=genus), stat="identity", position="dodge") + facet_grid(foam.type ~ .)
ggplot(data_1e5.ftype.core.psmelt, aes(x=genus, y=Abundance)) + geom_point(aes(fill=OTU)) + facet_grid(foam.type ~.)
ggplot(data_1e5.ftype.core.psmelt, aes(x=genus, y=Abundance)) + geom_point(aes(color=OTU)) + facet_grid(foam.type ~.)
ggplot(data_1e5.ftype.core.psmelt, aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type)) + scale_y_log10()
ggplot(data_1e5.ftype.core.psmelt, aes(x=OTU, y=Abundance)) + geom_boxplot(aes(fill=foam.type)) + scale_y_log10()
ggplot(data_1e5.ftype.core.psmelt, aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))
head(data_1e5.ftype.core.psmelt)
library(plyr)
data_1e5.ftype.core.otu.order<-ddply(data_1e5.ftype.core.psmelt, .(OTU), summarise, total = sum(Abundance))
head(data_1e5.ftype.core.otu.order)
data_1e5.ftype.core.otu.order$OTU <- reorder(data_1e5.ftype.core.otu.order$OTU, -data_1e5.ftype.core.otu.order$total)
str(data_1e5.ftype.core.otu.order)
head(data_1e5.ftype.core.otu.order[order(-data_1e5.ftype.core.otu.order$total),]
)
data_1e5.ftype.core.psmelt$OTU<-factor(data_1e5.ftype.core.psmelt$OTU, levels=data_1e5.ftype.core.otu.order$OTU)
head(data_1e5.ftype.core.psmelt$OTU)
str(data_1e5.ftype.core.psmelt$OTU)
ggplot(data_1e5.ftype.core.psmelt, aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))
ggplot(data_1e5.ftype.core.psmelt, aes(x=OTU, y=Abundance)) + geom_boxplot(aes(fill=foam.type))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(OTU)) %>% ggplot(aes(x=OTU, y=Abundance)) + geom_boxplot(aes(fill=foam.type))
library("magrittr")
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(OTU)) %>% ggplot(aes(x=OTU, y=Abundance)) + geom_boxplot(aes(fill=foam.type))
library(dplyr)
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(OTU)) %>% ggplot(aes(x=OTU, y=Abundance)) + geom_boxplot(aes(fill=foam.type))
head(data_1e5.ftype.core.otu.order[order(-data_1e5.ftype.core.otu.order$total),], n =40)
data_1e5.ftype.core.psmelt$OTU<-factor(data_1e5.ftype.core.psmelt$OTU, levels=levels(data_1e5.ftype.core.otu.order$OTU))
str(data_1e5.ftype.core.psmelt$OTU)
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(OTU)) %>% ggplot(aes(x=OTU, y=Abundance)) + geom_boxplot(aes(fill=foam.type))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(OTU)) %>% ggplot(aes(x=OTU, y=Abundance)) + geom_bar(aes(fill=foam.type), stat="identity", position="dodge")
data_1e5.ftype.core.genus.order<-ddply(data_1e5.ftype.core.psmelt, .(genus), summarise, total = sum(Abundance))
dim(data_1e5.ftype.core.genus.order)
str(data_1e5.ftype.core.genus.order)
data_1e5.ftype.core.genus.order$genus<-reorder(data_1e5.ftype.core.genus.order$genus, -data_1e5.ftype.core.genus.order$total)
str(data_1e5.ftype.core.genus.order)
data_1e5.ftype.core.psmelt$genus<-factor(data_1e5.ftype.core.psmelt$genus, levels=levels(data_1e5.ftype.core.genus.order$genus))
str(data_1e5.ftype.core.psmelt$genus)
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))+ facet_wrap(~ genus, nrow =1)
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))+ theme(axis.text.x=element_text(angle=45, hjust=1, vjust =1))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust =1))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))+ theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust =1))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) + scale_y_log10()
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
head( data_1e5.ftype.core.genus.order)
( data_1e5.ftype.core.genus.order)
834237.99/sum(data_1e5.ftype.core.genus.order$total)
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_bar(aes(fill=foam.type), stat="identity", position="dodge")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_violin(aes(fill=foam.type))+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_violin(aes(fill=foam.type))+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) + scale_y_log10()
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_boxplot(aes(fill=foam.type))+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5)) + scale_y_log10()
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_bar(aes(fill=foam.type), position="dodge")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, fill = foam.type)) + geom_bar(position="dodge")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance)) + geom_bar(aes(fill=foam.type), stat="identity", position="dodge")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
data_1e5.ftype.core.genus.ft<-ddply(data_1e5.ftype.core.psmelt, .(genus, foam.type), sumamrise, avg=mean(Abundance), se=sd(Abundance)/sqrt(length(Abundance)))
data_1e5.ftype.core.genus.ft<-ddply(data_1e5.ftype.core.psmelt, .(genus, foam.type), summarise, avg=mean(Abundance), se=sd(Abundance)/sqrt(length(Abundance)))
head(data_1e5.ftype.core.genus.ft)
data_1e5.ftype.core.genus.ft$genus<-factor(data_1e5.ftype.core.genus.ft$genus, levels=levels(data_1e5.ftype.core.genus.order$genus))
head(data_1e5.ftype.core.genus.ft)
data_1e5.ftype.core.genus.order$genus %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=avg)) + geom_bar(aes(fill = foam.type), stat="identity", position="dodge")
data_1e5.ftype.core.genus.order %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=avg)) + geom_bar(aes(fill = foam.type), stat="identity", position="dodge")
data_1e5.ftype.core.genus.ft %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=avg)) + geom_bar(aes(fill = foam.type), stat="identity", position="dodge")
p1<-data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, fill = foam.type)) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge")
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)
data_1e5.ftype.core.genus.order[order(-data_1e5.ftype.core.genus.order$total), ]
library(ggthemes)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangefrom() + tehme_tufte() + scale_y_continuous(breaks=extend_range_breaks()data_1e5.ftype.core.psmelt$Abundance)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangefrom() + tehme_tufte() + scale_y_continuous(breaks=extend_range_breaks()(data_1e5.ftype.core.psmelt$Abundance))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + tehme_tufte() + scale_y_continuous(breaks=extend_range_breaks()(data_1e5.ftype.core.psmelt$Abundance))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extend_range_breaks()(data_1e5.ftype.core.psmelt$Abundance))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(data_1e5.ftype.core.psmelt$Abundance))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte()
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(10000))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(c(0, 10000)))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) +  theme_tufte()
geom_rangeframe()
geom_rangeframe
head(data_1e5.ftype.core.genus.ft)
data_1e5.ftype.core.genus.ft.avg<-ddply(data_1e5.ftype.core.genus.ft, .(genus), summarise, ft.avg=mean(avg))
head(data_1e5.ftype.core.genus.ft.avg)
data_1e5.ftype.core.genus.ft.avg$genus<-reorder(data_1e5.ftype.core.genus.ft.avg$genus, -data_1e5.ftype.core.genus.ft.avg$ft.avg)
data_1e5.ftype.core.psmelt$genus<-factor(data_1e5.ftype.core.psmelt$genus, levels=levels(data_1e5.ftype.core.genus.ft.avg$genus))
str(data_1e5.ftype.core.psmelt$genus)
p1<-data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, fill = foam.type)) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte()
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte() + scale_y_continuous(break=0:10000)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=0:10000)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(data_1e5.ftype.core.genus.ft$avg)
)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(data_1e5.ftype.core.genus.ft$avg) + ylim(0, 10000)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(data_1e5.ftype.core.genus.ft$avg)) + ylim(0, 10000)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2) + geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(data_1e5.ftype.core.genus.ft$avg), limits = c(0, 10000))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)
data_1e5.ftype.core.genus.ft.avg[order(-data_1e5.ftype.core.genus.ft.avg$ft.avg), ]
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic()
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"))
library(RColorBrewer)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
par(oma=c(10,4,4,2))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
par(oma=c(10,10,10,10))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
p1 + stat_summary(fun.y = mean, geom="point", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
p1<-data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, color = foam.type))
p1 + stat_summary(fun.y = mean, geom="point", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + stat_smooth(method="lm")
p1<-data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, color = foam.type))
p1 + stat_summary(fun.y = mean, geom="point", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
p1 + stat_summary(fun.y = mean, geom="point") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
p1 + stat_summary(fun.y = mean, geom="point") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + stat_smooth(method="lm")
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
1<-data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, color = foam.type))
p1<-data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, color = foam.type))
p1 + stat_summary(fun.y = mean, geom="point") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + stat_smooth(method="lm")
p1<-data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
par(mar=c(4, 4.5, 2, 1))
par(oma=c(0, 0, 0, 0))
p1<-data_1e5.ftype.core.psmelt %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
par(mar=c(8, 4.5, 2, 1))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
par(mar=c(20, 4.5, 2, 1))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
par(mar=c(8, 4.5, 2, 1))
par(oma=c(5, 0, 0, 0))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)")
par(oma=c(10,10,10, 10))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=14, face="bold")) + ylab("Relative Abundance (1e-5)")
par(oma=c(20,20,20, 20))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=14, face="bold")) + ylab("Relative Abundance (1e-5)")
pdf("~/Box Sync/Manure Foaming/Manuscript/Figs_tables/core_otu_genus.pdf", width=20, height=10)
pdf("~/Box Sync/Manure Foaming/Manuscript/figures_and_tables/core_otu_genus.pdf", width=20, height=10)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=12, face="bold", angle = 45, hjust = 1, vjust = 1), axis.text.y=element_text(size=14, face="bold")) + ylab("Relative Abundance (1e-5)")
dev.off()
p1 + stat_summary(fun.y = mean, geom="point") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 90, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + stat_smooth(method="lm")
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=12, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=14, face="bold")) + ylab("Relative Abundance (1e-5)")
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=12, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=14, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=12, face="bold"), legend.position=c(1,0), legend.justification=c(1,0))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=12, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=14, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=12, face="bold"), legend.position=c(0.5,0), legend.justification=c(1,0))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=12, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=14, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=12, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
pdf("~/Box Sync/Manure Foaming/Manuscript/figures_and_tables/core_otu_genus.pdf", width=20, height=10)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
dev.off()
pdf("~/Box Sync/Manure Foaming/Manuscript/figures_and_tables/core_otu_genus.pdf", width=10, height=5)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
dev.off()
pdf("~/Box Sync/Manure Foaming/Manuscript/figures_and_tables/core_otu_genus.pdf", width=10, height=8)
dev.off()
pdf("~/Box Sync/Manure Foaming/Manuscript/figures_and_tables/core_otu_genus.pdf", width=10, height=8)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
dev.off()
pdf("~/Box Sync/Manure Foaming/Manuscript/figures_and_tables/core_otu_genus.pdf", width=10, height=6.6)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
dev.off()
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
head(bfs)
dim(bfs)
dim(bfall)
head(bfs)
fligner.test(si$SCFA, si$foam.type)
head(bfs_toplot)
bfs_toplot.bf20<-subset(bfs_toplot, bf > 20)
dim(bfs_toplot.bf20)
(bfs_toplot.bf20)
ggplot(bfs_toplot.bfs20, aes(x=value, y = bf)) + geom_point() + facet_grid(hypo ~ type)
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point() + facet_grid(hypo ~ type)
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point() + facet_grid(hypo ~ type, scales = "free")
ggplot(bfs_toplot.bf20, aes(y=value, a = bf)) + geom_point() + facet_grid(hypo ~ type, scales = "free")
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point() + facet_grid(type ~hypo, scales = "free")
plot(lm(si$DDGS ~ si$foam.type))
plot(lm(log(si$DDGS) ~ si$foam.type))
plot(lm(si$DDGS ~ si$foam.type))
plot(lm(log(si$DDGS) ~ si$foam.type))
plot(lm(log(si$DDGS+1) ~ si$foam.type))
anova(lm(si$DDGS ~ si$foam.type))
summary(aov(si$DDGS ~ si$foam.type))
plot(aov(si$DDGS ~ si$foam.type))
anova(lm(si$DDGS ~ si$foam.type | si$time_index))
head(fs)
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point() + facet_grid(type ~hypo, scales = "free")
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + facet_grid( ~hypo, scales = "free")
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_bar(aes(fill=type), stat = "identity") + facet_grid( ~hypo, scales = "free")
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_bar(aes(fill=type), stat = "identity") + facet_grid( ~hypo, scales = "free")
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_bar(aes(fill=type), stat = "identity") + facet_grid( ~hypo, scales = "free") +scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_bar(aes(fill=type), stat = "identity") + facet_wrap( ~hypo, nrow =1, scales = "free") +scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_bar(aes(fill=type), stat = "identity") + facet_wrap( ~hypo, nrow =1) +scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_bar(aes(fill=type), stat = "identity") + facet_wrap( ~hypo, nrow =1, scales="free_x") +scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_bar(aes(fill=type), stat = "identity") + facet_wrap( ~hypo, nrow =1, scales="free_x") +scale_y_log10() + theme_bw() + theme(axis.text.x=element_text(angle =90, hjsut =1, vjust =1, face="bold")
)
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_bar(aes(fill=type), stat = "identity") + facet_wrap( ~hypo, nrow =1, scales="free_x") +scale_y_log10() + theme_bw() + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold"))
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, nrow =1, scales="free_x") +scale_y_log10() + theme_bw() + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + geom_segment(aes(yend=bf), xend=0, colour="grey50")
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, nrow =1, scales="free_x") +scale_y_log10() + theme_bw() + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + geom_segment(aes(xend=bf), xend=0, colour="grey50")
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, nrow =1, scales="free_x") +scale_y_log10() + theme_bw() + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + geom_segment(aes(xend=value), xend=0, colour="grey50")
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, nrow =1, scales="free_x") +scale_y_log10() + theme_bw() + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + geom_segment(aes(yend=bf), xend=0, colour="grey50")
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, nrow =1, scales="free_x") +scale_x_log10() + theme_bw() + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + geom_segment(aes(yend=bf), xend=0, colour="grey50")
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, nrow =1, scales="free_x") + theme_bw() + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + geom_segment(aes(yend=bf), xend=0, colour="grey50")
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, nrow =1, scales="free_x") + theme_bw() + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + geom_segment(aes(yend=bf), xend=0, colour="grey50")
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, nrow =1, scales="free_x") + theme_bw() + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + geom_segment(aes(yend=bf), xend=0, colour="grey50") + scale_y_log10()
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, nrow =1, scales="free_x") + theme_bw() + geom_segment(aes(yend=value), xend=0, colour="grey50")
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, ncol =1, scales="free_y") + theme_bw() + geom_segment(aes(yend=value), xend=0, colour="grey50") + scale_x_log10()
bfs_toplot.bf20
bfs_toplot.bf20<-bfs_toplot.bf20[ bfs_toplot.bf20$value != "surface.depth", ]
bfs_toplot.bf20
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + facet_wrap( hypo ~., ncol =1, scales="free_y") + theme_bw() + geom_segment(aes(yend=value), xend=0, colour="grey50") + scale_x_log10()
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, ncol =1, scales="free_y") + theme_bw() + geom_segment(aes(yend=value), xend=0, colour="grey50") + scale_x_log10() + coord_flip()
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + theme_bw() + geom_segment(aes(yend=value), xend=0, colour="grey50") + scale_x_log10() + coord_flip()
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + theme_bw() + geom_segment(aes(yend=value), xend=0, colour="grey50") + scale_x_log10() + coord_flip() + facet_wrap( ~hypo, nrow =1)
ggplot(bfs_toplot.bf20, aes(y=value, x = bf)) + geom_point(aes(color=type)) + theme_bw() + geom_segment(aes(yend=value), xend=0, colour="grey50") + scale_x_log10() + coord_flip() + facet_wrap( ~hypo, nrow =1, scales="free_x")
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type)) + facet_wrap( ~hypo, nrow =1, scales="free_x") + theme_classic() + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type), size = 4) + facet_wrap( ~hypo, nrow =1, scales="free_x") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
head(bfs_toplot.bf20)
str(bfs_toplot.bf20$hypo)
bfs_toplot.bf20$hypo <- factor(bfs_toplot.bf20$hypo, levels=c("No Foam>Crust>Foam", "Crust>Foam>No Foam", "Foam>Crust>No Foam", "Foam>No Foam>Crust"))
str(bfs_toplot.bf20$hypo)
uniq(bfs_toplot.bf20$hypo)
unique(bfs_toplot.bf20$hypo)
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type), size = 4) + facet_wrap(generes ~hypo, nrow =1, scales="free_x") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
head(bfs_toplot.bf20)
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type), size = 4) + facet_wrap( ~hypo, nrow =1, scales="free_x") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type), size = 4) + facet_wrap(genres ~hypo, nrow =1, scales="free_x") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type), size = 4) + facet_grid(genres ~hypo, scales="free_x") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type), size = 4) + facet_grid(genres ~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_bar(aes(fill=type), stat="identity") + facet_grid(genres ~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_bar(aes(fill=type), stat="identity") + facet_grid(genres ~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()+theme_wsj() + scale_colour_wsj("colors6", "")
library(ggthemes)
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_bar(aes(fill=type), stat="identity") + facet_grid(genres ~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()+theme_wsj() + scale_colour_wsj("colors6", "")
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type), size = 4) + facet_grid(genres ~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10() + theme_wsj() + scale_colour_wsj("colors6", "")
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=type), size = 4) + facet_grid(genres ~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10() + theme_hc() + scale_colour_hc()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=genres), size = 4) + facet_grid(~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()+theme_hc() + scale_colour_hc()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=genres), size = 4) + facet_grid(~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()+theme_hc() + scale_colour_hc() + theme_tuftes()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=genres), size = 4) + facet_grid(~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()+theme_hc() + scale_colour_hc() + theme_tufte()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=genres), size = 4) + facet_grid(~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()+theme_hc() + scale_colour_hc() + theme_tufte() + geom_rangeframe()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=genres), size = 4) + facet_grid(~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()+theme_hc() + scale_colour_hc()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=genres), size = 4) + facet_grid(~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
bfs_toplot.bf20$value<-reorder(bfs_toplot.bf20$value, -bfs_toplot.bf20$bf)
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=genres), size = 4) + facet_grid(~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
bfs_toplot.bf20$value<-reorder(bfs_toplot.bf20$value, bfs_toplot.bf20$bf)
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=genres), size = 4) + facet_grid(~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=genres), size = 4) + facet_grid(~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =90, hjust =1, vjust =1, face="bold")) + scale_y_log10()+theme_hc() + scale_colour_hc()
colors
ggplot(bfs_toplot.bf20, aes(x=value, y = bf)) + geom_point(aes(color=genres), size = 4) + facet_grid(~hypo, scales="free") + theme_classic() + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + theme(axis.text.x=element_text(angle =45, hjust =1, vjust =1, face="bold")) + scale_y_log10() + theme_hc() + scale_colour_hc() + theme(strip.text.x = element_text(size=12, face="bold"), strip.background = element_rect(fill=c("#1B9E77", "#D95F02", "#7570B3", "#7570B3")) 
)
head(bfs_toplot.bf20)
write.table(bfs_toplot.bf20, file = "~/Box Sync/Manure Foaming/Manuscript/figures_and_tables/si_w_bf01_min20.txt", sep="\t", row.names=F)
p1 + stat_summary(fun.y = mean, geom="point") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 90, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + stat_smooth(method="lm")
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
head(data_1e5.ftype.core.genus.ft.avg)
hist(data_1e5.ftype.core.genus.ft.avg$ft.avg)
hist(data_1e5.ftype.core.genus.ft.avg$ft.avg, breaks=30)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
data_1e5.ftype.core.genus.ft.avg.top10<-head(data_1e5.ftype.core.genus.ft.avg, n=10)
data_1e5.ftype.core.genus.ft.avg.top10
data_1e5.ftype.core.genus.ft.avg.top10<-head(data_1e5.ftype.core.genus.ft.avg[order(-data_1e5.ftype.core.genus.ft.avg$ft.avg), ], n=10)
data_1e5.ftype.core.genus.ft.avg.top10
data_1e5.ftype.core.psmelt.top10<-data_1e5.ftype.core.psmelt[ data_1e5.ftype.core.psmelt$genus %in%  data_1e5.ftype.core.genus.ft.avg.top10$genus, ]
head(data_1e5.ftype.core.psmelt.top10)
p1<-data_1e5.ftype.core.psmelt.top10 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
library(dplyr)
p1<-data_1e5.ftype.core.psmelt.top10 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 45, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
par(oma=c(50, 0, 0, 0))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 45, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
par(oma=c(50, 50, 50, 50))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 45, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
par(mar=c(50, 50, 50, 50))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 45, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 45, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0)) +scale_x_continuous(expand = c(.1, .1))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 45, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0)) +scale_x_discrete(expand = c(1, 1))
pdf("test.pdf", width=100, height=20)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 45, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
dev.off()
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
dim(totu)
data.mds<-metaMDS(totu, k=8, autotransform=FALSE)
library(vegan)
data.mds<-metaMDS(totu, k=8, autotransform=FALSE)
head(totu[, 1:5])
head(rowSums(totu))
data.mds<-metaMDS(totu, k=3, autotransform=FALSE)
data.mds
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
data.mds<-metaMDS(totu, k=3, autotransform=FALSE)
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_eclipse.R")
ggplot.NMDS.ellipse(data.mds, si$foam.type, colors)
ggplot(si, aes(y=MPR_slurry, x=SCFA, color=foam.type)) + geom_point(size=2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold")) + theme(legend.title=element_blank()) + theme(legend.justification=c(1,0), legend.position=c(1,0.7)) + theme(legend.text = element_text(size = 12)) + xlab("MPR (Slurry) L/L-day") + ylab("SCFA (mg/L)")
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
ls(pattern="data_1e5.ftype.core"))
ls(pattern="data_1e5.ftype.core")
ls(pattern="data_1e5.ftype.core", sorted=F)
library(phyloseq)
data_1e5.ftype
data_1e5.ftype.core
dim(data_1e5.ftype.core.psmelt)
head(data_1e5.ftype.core.genus.order)
names(data_1e5.ftype.core.genus.ft)
names(data_1e5.ftype.core.psmelt)
head((data_1e5.ftype.core.psmelt[, c(1, 131)])
)
library(dplyr)
head((data_1e5.ftype.core.psmelt[, c(1:3,125, 131)])
)
 data_1e5.ftype.core.genus.ft<-data_1e5.ftype.core.psmelt %>% select("genus", "sample", "foam.type", "Abundance") %>% head
 data_1e5.ftype.core.genus.ft<-data_1e5.ftype.core.psmelt %>% select(genus, sample, foam.type, Abundance) %>% head
 data_1e5.ftype.core.genus.ft<-data_1e5.ftype.core.psmelt %>% select(genus, Sample, foam.type, Abundance) %>% head
data_1e5.ftype.core.psmelt %>% select(genus, Sample, foam.type, Abundance) %>% head
 data_1e5.ftype.core.genus.ft<-data_1e5.ftype.core.psmelt %>% select(genus, Sample, foam.type, Abundance) %>% summarise(genus, total = sum(Abundance))
 data_1e5.ftype.core.genus.ft<-data_1e5.ftype.core.psmelt %>% select(genus, Sample, foam.type, Abundance) %>% summarise(c(genus, Sample, foam.type), total = sum(Abundance))%
 data_1e5.ftype.core.genus.ft<-data_1e5.ftype.core.psmelt %>% select(genus, Sample, foam.type, Abundance) %>% group_by(genus, Sample, foam.type) %>%  summarise(total = sum(Abundance)) %>% group_by(genus, foam.type
head(data_1e5.ftype.core.genus.ft)
str(data_1e5.ftype.core.genus.ft)
class(data_1e5.ftype.core.genus.ft)
data_1e5.ftype.core.genus.ft.avg<-data_1e5.ftype.core.genus.ft %>% group_by(genus, foam.type) %>% summarise(avg=mean(total))
head(data_1e5.ftype.core.genus.ft.avg)
data_1e5.ftype.core.genus.order<-data_1e5.ftype.core.genus.ft.avg %>% group_by(genus) %>% summarise(ft_total=sum(avg))
head(data_1e5.ftype.core.genus.order)
data_1e5.ftype.core.genus.order$genus<-reorder(data_1e5.ftype.core.genus.order$genus, -data_1e5.ftype.core.genus.order$ft_total)
head(data_1e5.ftype.core.genus.order)
str(data_1e5.ftype.core.genus.order)
head(data_1e5.ftype.core.genus.order)
head(data_1e5.ftype.core.genus.ft)
data_1e5.ftype.core.genus.ft
data_1e5.ftype.core.genus.ft$genus<-factor(data_1e5.ftype.core.genus.ft$genus, levels=levels(data_1e5.ftype.core.genus.order$genus)
)
str(data_1e5.ftype.core.genus.ft)
library(ggplot2)
p1<-data_1e5.ftype.core.genus.ft %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
library(RColorBrewer)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
library(baysfactor)
library(bayesfactor)
library(bayesfactors)
library(BayesFactor)
unique(data_1e5.ftype.core.genus.ft$genus)
for (i in unique(data_1e5.ftype.core.genus.ft$genus)){
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
test<-subset(data_1e5.ftype.core.genus.ft, genus == i)
dim(test)
}
head(test)
for (i in unique(data_1e5.ftype.core.genus.ft$genus)){
test<-subset(data_1e5.ftype.core.genus.ft, genus == i)
print(dim(test)
)
}
for (i in unique(data_1e5.ftype.core.genus.ft$genus)){
test<-subset(data_1e5.ftype.core.genus.ft, genus == i)
print(head(test))
}
avg<-test %>% group_by(foam.type) %>% summarise(ft_avg=mean(total))
head(avg)
avg<-avg[order(avg$ft_avg), ]
head(avg)
for (i in unique(data_1e5.ftype.core.genus.ft$genus)){
sink("bacteria_core_genus_bayes_factors_order_restricted_H.txt", append = TRUE)
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
test<-subset(data_1e5.ftype.core.genus.ft, genus == i)
bf1<-anovaBF(total ~ foam.type, data=test)
avg<-test %>% group_by(foam.type) %>% summarise(ft_avg=mean(total))
avg<-avg[order(avg$ft_avg), ]
samples = posterior(bf1, iterations = 10000)
print(paste(avg[3, 1], avg[2,1], avg[1,1], sep=">"))
consistent = (samples[, paste("foam.type", avg[3,1], sep="-")] > samples[, paste("foam.type", avg[2,1], sep="-")]) & (samples[, paste("foam.type", avg[2,1], sep="-")] > samples[, paste("foam.type", avg[1,1], sep="-")])
N_consistent = sum(consistent)
bf_restriction_against_full = (N_consistent / 10000) / (1 / 6)
bf_restriction_against_full
bf_full_against_null = as.vector(bf1)
bf_restriction_against_null = bf_restriction_against_full * bf_full_against_null
print(bf_restriction_against_null)
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
head(test)
head(data.frame(test))
str(data.frame(test)$total)
str(test$total)
bf1<-anovaBF(total ~ foam.type, data=test)
bf1<-anovaBF(total ~ foam.type, data=data.frame(test))
for (i in unique(data_1e5.ftype.core.genus.ft$genus)){
sink("bacteria_core_genus_bayes_factors_order_restricted_H.txt", append = TRUE)
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
test<-subset(data_1e5.ftype.core.genus.ft, genus == i)
bf1<-anovaBF(total ~ foam.type, data=data.frame(test))
avg<-test %>% group_by(foam.type) %>% summarise(ft_avg=mean(total))
avg<-avg[order(avg$ft_avg), ]
samples = posterior(bf1, iterations = 10000)
print(paste(avg[3, 1], avg[2,1], avg[1,1], sep=">"))
consistent = (samples[, paste("foam.type", avg[3,1], sep="-")] > samples[, paste("foam.type", avg[2,1], sep="-")]) & (samples[, paste("foam.type", avg[2,1], sep="-")] > samples[, paste("foam.type", avg[1,1], sep="-")])
N_consistent = sum(consistent)
bf_restriction_against_full = (N_consistent / 10000) / (1 / 6)
bf_restriction_against_full
bf_full_against_null = as.vector(bf1)
bf_restriction_against_null = bf_restriction_against_full * bf_full_against_null
print(bf_restriction_against_null)
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
avg
avg<-data.frame(avg[order(avg$ft_avg), ])
avg
print(paste(avg[3, 1], avg[2,1], avg[1,1], sep=">"))
for (i in unique(data_1e5.ftype.core.genus.ft$genus)){
sink("bacteria_core_genus_bayes_factors_order_restricted_H.txt", append = TRUE)
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
test<-subset(data_1e5.ftype.core.genus.ft, genus == i)
bf1<-anovaBF(total ~ foam.type, data=data.frame(test))
avg<-test %>% group_by(foam.type) %>% summarise(ft_avg=mean(total))
avg<-data.frame(avg[order(avg$ft_avg), ])
samples = posterior(bf1, iterations = 10000)
print(paste(avg[3, 1], avg[2,1], avg[1,1], sep=">"))
consistent = (samples[, paste("foam.type", avg[3,1], sep="-")] > samples[, paste("foam.type", avg[2,1], sep="-")]) & (samples[, paste("foam.type", avg[2,1], sep="-")] > samples[, paste("foam.type", avg[1,1], sep="-")])
N_consistent = sum(consistent)
bf_restriction_against_full = (N_consistent / 10000) / (1 / 6)
bf_restriction_against_full
bf_full_against_null = as.vector(bf1)
bf_restriction_against_null = bf_restriction_against_full * bf_full_against_null
print(bf_restriction_against_null)
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
bac_genus_bf<-read.delim("bacteria_core_genus_bayes_factors_order_restricted_H.table", sep="\t", header=F)
bac_genus_bf
names(bac_genus_bf)<-c("genus", "H1", "B01")
bac_genus_bf$genus<-factor(bac_genus_bf$genus, levels=levels(data_1e5.ftype.core.genus.order$genus))
str(bac_genus_bf)
unique(bac_genus_bf$H1)
bac_genus_bf<-bac_genus_bf[order(bac_genus_bf$H1), ]
(bac_genus_bf)
bac_genus_bf_min20<-subset(bac_genus_bf, B01 >= 20)
dim(bac_genus_bf)
dim(bac_genus_bf_min20)
bac_genus_bf_min100<-subset(bac_genus_bf, B01 >= 100)
dim(bac_genus_bf_min100)
(bac_genus_bf_min100)
data_1e5.ftype.core.genus.order
5304.543/sum(data_1e5.ftype.core.genus.order$ft_total)
bfs
(bac_genus_bf_min100)
head(data_1e5.ftype.core.genus.ft)
data_1e5.ftype.core.genus.ft.bf100<-data_1e5.ftype.core.genus.ft %>% filter( genus %in% bac_genus_bf_min100$genus)
dim(data_1e5.ftype.core.genus.ft.bf100)
unique(data_1e5.ftype.core.genus.ft.bf100$genus)
head(data_1e5.ftype.core.genus.ft.bf100)
test<-merge(data_1e5.ftype.core.genus.ft.bf100, ac_genus_bf_min100, "genus")
test<-merge(data_1e5.ftype.core.genus.ft.bf100, bac_genus_bf_min100, "genus")
head(test)
dim(data_1e5.ftype.core.genus.ft.bf100)
dim(test)
data_1e5.ftype.core.genus.ft.bf100<-merge(data_1e5.ftype.core.genus.ft.bf100, bac_genus_bf_min100, "genus")
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=Abundance, color = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, color = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0)) + facet_grid(~H1)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0)) + facet_grid(~H1, scales="free")
p1 + stat_summary(fun.y = mean, geom="point") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 90, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + stat_smooth(method="lm")
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=foam.type, size=total))
p1 + stat_summary(fun.y = mean, geom="point") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_color_brewer(palette="Dark2") + theme(axis.title.x=element_blank()) + theme(axis.text.x=element_text(size=12, face="bold", angle = 90, hjust = 1, vjust = 1), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + stat_smooth(method="lm")
p1 + geom_point()
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free")
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free", space="free")
head(data_1e5.ftype.core.genus.ft.bf100)
unique(data_1e5.ftype.core.genus.ft.bf100$H1)
data_1e5.ftype.core.genus.ft.bf100$H1<-factor(data_1e5.ftype.core.genus.ft.bf100$H1, levels=c("No Foam>Crust>Foam", "Crust>No Foam>Foam", "Foam>Crust>No Foam")
)
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free", space="free")
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free", space="free") + theme(strip.text.x=element_text(angle=45))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free", space="free") + theme(strip.text.x=element_text(angle=90))
p1 + stat_summary(fun.y = mean, geom="boxplot", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free", space="free")
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type)) + geom_boxplot()
p1 + stat_summary(fun.y = mean, position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free", space="free") + theme(strip.text.x=element_text(angle=90))
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type)) + geom_violin()
p1 + stat_summary(fun.y = mean, position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free", space="free") + theme(strip.text.x=element_text(angle=90))
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type)) + geom_boxplot()
p1
p1 + stat_summary(fun.y = mean, position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free", space="free") + theme(strip.text.x=element_text(angle=90))
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free", space="free") + theme(strip.text.x=element_text(angle=90))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))  + facet_grid(~H1, scales="free", space="free") + theme(strip.text.x=element_text(angle=90, face="bold", size=14))
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
ls()
data_1e5.ftype
data_1e5.ftype.core
data_1e5.ftype.core.samplesum<-data.frame(sample_sums(data_1e5.ftype.core))
head(data_1e5.ftype.core.samplesum
)
data_ftype_F_core.RDS
data_ftype_F_core
data.ftype.f.core <- readRDS("RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/data_ftype_F_core.RDS")
data.ftype.f.core
data.ftype.c.core <- readRDS("RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/data_ftype_C_core.RDS")
data.ftype.nf.core <- readRDS("RDS_inputs_MINtaxasums5_donotuseRemoved_MINseq10k/data_ftype_NF_core.RDS")
f.core.otu<-data.frame(taxa_names(data.ftype.f.core))
head(f.core.otu)
names(f.core.otu)<-"core"
c.core.otu<-data.frame(taxa_names(data.ftype.c.core))
names(c.core.otu)<-"core"
nf.core.otu<-data.frame(taxa_names(data.ftype.nf.core))
names(nf.core.otu)<-"core"
head(nf.core.otu)
test<-rbind(f.core.otu, c.core.otu, nf.core.otu)
head(test)
length(unique(test$core))
core_otu$core_otus
test$core
unique(test$core)
head(test)
dim(test)
dim(f.core.otu)
data.ftype.f.core
head(core_otu)
dim(core_otu)
core_otu<-data.frame(unique(test$core))
dim(core_otu)
head(core_otu)
names(core_otu)<-"core_otus"
test<-subset_samples(data.ftype, foam.type=="F")
test
test<-prune_taxa(taxa_sums(test)>0, test)
test
test<-filter_taxa(test, function(x) sum(x>0) > (1*length(x)), T)
test<-filter_taxa(test, function(x) sum(x=length(x)), T)
test<-filter_taxa(test, function(x) x=length(x), T)
?filter_taxa
test<-filter_taxa(test, function(x) sum(x>=1) > (1*length(x)), T)
test
test1<-filter_taxa(test, function(x) sum(x>0) > (0.5*length(x)), T)
test1
test1<-filter_taxa(test, function(x) sum(x>0) > (0.8*length(x)), T)
test1
test1<-filter_taxa(test, function(x) sum(x>0) > (1.0*length(x)), T)
test1<-filter_taxa(test, function(x) sum(x>0) == (1.0*length(x)), T)
test1
data.ftype.f.core
data_1e5.ftype.core<-prune_taxa(taxa_names(data_1e5.ftype) %in% core_otu$core_otus, data_1e5.ftype
)
data_1e5.ftype.core
data.ftype.core
data.ftype.core<-prune_taxa(taxa_names(data.ftype) %in% core_otu$core_otus, data.ftype)
data.ftype.core
data_1e5.ftype.core.psmelt<-psmelt(data_1e5.ftype.core)
head(data_1e5.ftype.core.psmelt)
data_1e5.ftype.core.psmelt$foam.type<-gsub("\\bF\\b", "Foam", data_1e5.ftype.core.psmelt$foam.type)
data_1e5.ftype.core.psmelt$foam.type<-gsub("\\bC\\b", "Crust", data_1e5.ftype.core.psmelt$foam.type)
data_1e5.ftype.core.psmelt$foam.type<-gsub("\\bNF\\b", "No Foam", data_1e5.ftype.core.psmelt$foam.type)
head(data_1e5.ftype.core.psmelt$foam.type)
unique(data_1e5.ftype.core.psmelt$foam.type)
data_1e5.ftype.core.psmelt$foam.type<-factor(data_1e5.ftype.core.psmelt$foam.type, levels=c("No Foam", "Crust", "Foam"))
head(data_1e5.ftype.core.psmelt$foam.type)
data_1e5.ftype.core.genus.ft <- data_1e5.ftype.core.psmelt %>% select(genus, Sample, foam.type, Abundance) %>% group_by(genus, Sample, foam.type) %>%  summarise(total = sum(Abundance))
head(data_1e5.ftype.core.genus.ft)
data_1e5.ftype.core.genus.ft.avg<-data_1e5.ftype.core.genus.ft %>% group_by(genus, foam.type) %>% summarise(avg=mean(total))
head(data_1e5.ftype.core.genus.ft.avg)
data_1e5.ftype.core.genus.order<-data_1e5.ftype.core.genus.ft.avg %>% group_by(genus) %>% summarise(ft_total=sum(avg))
head(data_1e5.ftype.core.genus.order)
data_1e5.ftype.core.genus.order$genus<-reorder(data_1e5.ftype.core.genus.order$genus, -data_1e5.ftype.core.genus.order$ft_total)
head(data_1e5.ftype.core.genus.order)
str(data_1e5.ftype.core.genus.order)
data_1e5.ftype.core.genus.ft$genus<-factor(data_1e5.ftype.core.genus.ft$genus, levels=levels(data_1e5.ftype.core.genus.order$genus))
head(data_1e5.ftype.core.genus.ft)
str(data_1e5.ftype.core.genus.ft)
p1<-data_1e5.ftype.core.genus.ft %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
library(RColorBrewer)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0))
library(BayesFactor)
for (i in unique(data_1e5.ftype.core.genus.ft$genus)){
sink("bacteria_core_genus_bayes_factors_order_restricted_H.txt", append = TRUE)
tryCatch({
print(c("Processing treatment: ", i, "!!!!!!!!!!!!"))
test<-subset(data_1e5.ftype.core.genus.ft, genus == i)
bf1<-anovaBF(total ~ foam.type, data=data.frame(test))
avg<-test %>% group_by(foam.type) %>% summarise(ft_avg=mean(total))
avg<-data.frame(avg[order(avg$ft_avg), ])
samples = posterior(bf1, iterations = 10000)
print(paste(avg[3, 1], avg[2,1], avg[1,1], sep=">"))
consistent = (samples[, paste("foam.type", avg[3,1], sep="-")] > samples[, paste("foam.type", avg[2,1], sep="-")]) & (samples[, paste("foam.type", avg[2,1], sep="-")] > samples[, paste("foam.type", avg[1,1], sep="-")])
N_consistent = sum(consistent)
bf_restriction_against_full = (N_consistent / 10000) / (1 / 6)
bf_restriction_against_full
bf_full_against_null = as.vector(bf1)
bf_restriction_against_null = bf_restriction_against_full * bf_full_against_null
print(bf_restriction_against_null)
}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
sink()
}
bac_genus_bf<-read.delim("bacteria_core_genus_bayes_factors_order_restricted_H.table", sep="\t", header=F)
bac_genus_bf
names(bac_genus_bf)<-c("genus", "H1", "B01")
bac_genus_bf$genus<-factor(bac_genus_bf$genus, levels=levels(data_1e5.ftype.core.genus.order$genus))
str(bac_genus_bf)
unique(bac_genus_bf$H1)
bac_genus_bf<-bac_genus_bf[order(bac_genus_bf$H1), ]
(bac_genus_bf)
bac_genus_bf_min20<-subset(bac_genus_bf, B01 >= 20)
dim(bac_genus_bf)
dim(bac_genus_bf_min20)
bac_genus_bf_min100<-subset(bac_genus_bf, B01 >= 100)
dim(bac_genus_bf_min100)
(bac_genus_bf_min100)
head(data_1e5.ftype.core.genus.ft)
length(unique(data_1e5.ftype.core.genus.ft))
length(unique(data_1e5.ftype.core.genus.ft$genus))
data_1e5.ftype.core.genus.ft.bf100<-data_1e5.ftype.core.genus.ft %>% filter( genus %in% bac_genus_bf_min100$genus)
dim(data_1e5.ftype.core.genus.ft.bf100)
unique(data_1e5.ftype.core.genus.ft.bf100$genus)
head(data_1e5.ftype.core.genus.ft.bf100)
data_1e5.ftype.core.genus.ft.bf100<-merge(data_1e5.ftype.core.genus.ft.bf100, bac_genus_bf_min100, "genus")
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0)) + facet_grid(~H1, scales="free", space="free") + theme(strip.text.x=element_text(angle=90, face="bold", size=14))
unique(data_1e5.ftype.core.genus.ft.bf100$H1)
data_1e5.ftype.core.genus.ft.bf100$H1<-factor(data_1e5.ftype.core.genus.ft.bf100$H1, levels=c("No Foam>Crust>Foam", "Crust>No Foam>Foam", "Crust>Foam>No Foam", "Foam>Crust>No Foam", "Foam>No Foam>Crust")
)
p1<-data_1e5.ftype.core.genus.ft.bf100 %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=12, face="bold")) + ylab("Relative Abundance (1e-5)") + theme(legend.title = element_blank(), legend.text=element_text(size=10, face="bold"), legend.position=c(1,0.5), legend.justification=c(1,0)) + facet_grid(~H1, scales="free", space="free") + theme(strip.text.x=element_text(angle=90, face="bold", size=14))
save.image("all_16s_physeq.RData")
data_1e5.ftype.core.samplesum<-data.frame(sample_sums(data_1e5.ftype.core))
head(data_1e5.ftype.core.samplesum
)
data_1e5.ftype.samplesum<-data.frame(sample_sums(data_1e5.ftype))
head(data_1e5.ftype.samplesum)
tail(data_1e5.ftype.samplesum)
tail(data_1e5.ftype.core.samplesum)
data_1e5.ftype.core.samplesum$total<-data_1e5.ftype.samplesum$sample_sums.data_1e5.ftype.
tail(data_1e5.ftype.core.samplesum)
data_1e5.ftype.core.samplesum$core.percent<-data_1e5.ftype.core.samplesum$sample_sums.data_1e5.ftype.core. * 100/data_1e5.ftype.core.samplesum$total
tail(data_1e5.ftype.core.samplesum)
head(data_1e5.ftype.core.samplesum)
colnames(data_1e5.ftype.core.samplesum)[1]<-"core.sum"
head(data_1e5.ftype.core.samplesum)
names(si)
dim(si)
dim(data_1e5.ftype.core.samplesum)
data_1e5.ftype.core.samplesum<-merge(data_1e5.ftype.core.samplesum, si[, c("SAMPLES", "time_index", "foam.type", "myear")], by.x="row.names", by.y="SAMPLES")
dim(data_1e5.ftype.core.samplesum)
data_1e5.ftype.core.samplesum$foam.type <- factor(data_1e5.ftype.core.samplesum<-merge(data_1e5.ftype.core.samplesum, si[, c("SAMPLES", "time_index", "foam.type", "myear")], by.x="row.names", by.y="SAMPLES")
data_1e5.ftype.core.samplesum$foam.type <- factor(data_1e5.ftype.core.samplesum$foam.type, levels=c("No Foam", "Crust", "Foam"))
dim(data_1e5.ftype.core.samplesum)
head(data_1e5.ftype.core.samplesum)
str(data_1e5.ftype.core.samplesum)
ggplot(data_1e5.ftype.core.samplesum, aes(x = Row.names, y = core.percent, fill = myear)) + geom_bar(stat="identity") + facet_grid(~foam.type) 
ggplot(data_1e5.ftype.core.samplesum, aes(x = Row.names, y = core.percent)) + geom_bar(stat="identity") + facet_grid(~foam.type) 
ggplot(data_1e5.ftype.core.samplesum, aes(y = core.percent)) + geom_bar(stat="identity") + facet_grid(~foam.type) 
ggplot(data_1e5.ftype.core.samplesum, aes(x = myeasr, y = core.percent)) + geom_bar(stat="identity") + facet_grid(~foam.type) 
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent)) + geom_bar(stat="identity") + facet_grid(~foam.type) 
max(data_1e5.ftype.core.samplesum$core.percent)
data_1e5.ftype.core.samplesum$Row.names <- as.factor(data_1e5.ftype.core.samplesum$Row.names)
str(data_1e5.ftype.core.samplesum)
data_1e5.ftype.core.samplesum$Row.names <- reorder(data_1e5.ftype.core.samplesum$Row.names, data_1e5.ftype.core.samplesum$time_index)
str(data_1e5.ftype.core.samplesum)
ggplot(data_1e5.ftype.core.samplesum, aes(x = Row.names, y = core.percent)) + geom_bar(stat="identity") + facet_grid(~foam.type) 
ggplot(data_1e5.ftype.core.samplesum, aes(x = Row.names, y = core.percent)) + geom_bar(stat="identity") + facet_grid(~foam.type, scale="free") 
data_1e5.ftype.core.samplesum.timeavg<-data_1e5.ftype.core.samplesum %>% group_by(myear) %>% summarise(avg_month=mean(core.percent))
head(data_1e5.ftype.core.samplesum.timeavg)
data_1e5.ftype.core.samplesum.timeavg<-data_1e5.ftype.core.samplesum %>% group_by(myear, foam.type) %>% summarise(avg_month=mean(core.percent))
head(data_1e5.ftype.core.samplesum.timeavg)
ggplot(data_1e5.ftype.core.samplesum.timeavg, aes(x = myear, y = avg_month)) + geom_bar(stat="identity") + facet_grid(~foam.type, scale="free") 
ggplot(data_1e5.ftype.core.samplesum, aes(x = Row.names, y = core.percent)) + geom_boxplot()+ facet_grid(~foam.type, scale="free") 
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent)) + geom_boxplot()+ facet_grid(~foam.type, scale="free") 
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent)) + geom_violin()+ facet_grid(~foam.type, scale="free") 
library(ggtheme)
library(ggthemes)
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent)) + geom_violin()+ facet_grid(~foam.type, scale="free") + theme_tufte()
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent)) + geom_tufteboxplot()+ facet_grid(~foam.type, scale="free") + theme_tufte()
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ theme_tufte()
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ theme_bw() + theme(axis.text.x = element_text(angle=45, face="bold") + scale_fill_brewer(pallete="Dark2")
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ theme_bw() + theme(axis.text.x = element_text(angle=45, face="bold")) + scale_fill_brewer(pallete="Dark2")
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ theme_bw() + theme(axis.text.x = element_text(angle=45, face="bold")) + scale_fill_brewer(palette="Dark2")
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ theme_bw() + theme(axis.text.x = element_text(angle=45, face="bold")) + scale_fill_brewer(palette="Dark2") + geom_rangeframe() + theme_tufte()
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(data_1e5.ftype.core.samplesum$core.percent)) + theme(axis.text.x = element_text(face="bold")) + scale_fill_brewer(palette="Dark2")
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent), 2)) + theme(axis.text.x = element_text(face="bold")) + scale_fill_brewer(palette="Dark2")
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent, 2)) + theme(axis.text.x = element_text(face="bold")) + scale_fill_brewer(palette="Dark2")
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent, 2))) + theme(axis.text.x = element_text(face="bold")) + scale_fill_brewer(palette="Dark2")
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent, 2))) + theme(axis.text.x = element_text(size =12, face="bold"), title.text.x=element_blank(), axis.text.y=element_text(size=12, face="bold"), title.text.y=element_text(size=14, face="bold"))+ ylab("Relative Abundance of Core Communities (%)") + scale_fill_brewer(palette="Dark2")
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent, 2))) + theme(axis.text.x = element_text(size =12, face="bold"), axis.title.x=element_blank(), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+ ylab("Relative Abundance of Core Communities (%)") + scale_fill_brewer(palette="Dark2")
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent, 2))) + theme(axis.text.x = element_text(size =12, face="bold", angle=45, hjust=1, vjust=1), axis.title.x=element_blank(), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+ ylab("Relative Abundance of Core Communities (%)") + scale_fill_brewer(palette="Dark2")
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent, 2))) + theme(axis.text.x = element_text(size =12, face="bold", angle=45, hjust=1, vjust=1), axis.title.x=element_blank(), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+ ylab("Proportion of Core Communities (%)") + scale_fill_brewer(palette="Dark2") + theme(legend.title = element_blank(), legend.text=element_text(size=12, face="bold"), legend.position=c(0.2,0.2), legend.justification=c(1,0))
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent, 2))) + theme(axis.text.x = element_text(size =12, face="bold", angle=45, hjust=1, vjust=1), axis.title.x=element_blank(), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+ ylab("Proportion of Core Communities (%)") + scale_fill_brewer(palette="Dark2") + theme(legend.title = element_blank(), legend.text=element_text(size=12, face="bold"), legend.position=c(0.2,0.1), legend.justification=c(1,0))
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent, 2))) + theme(axis.text.x = element_text(size =12, face="bold", angle=45, hjust=1, vjust=1), axis.title.x=element_blank(), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+ ylab("Proportion of Core Communities (%)") + scale_fill_brewer(palette="Dark2") + theme(legend.title = element_blank(), legend.text=element_text(size=12, face="bold"), legend.position=c(0.2,0), legend.justification=c(1,0))
save.images("all_16s_physeq.RData")
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(dplyr)
library(phyloseq)
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent, 2))) + theme(axis.text.x = element_text(size =12, face="bold", angle=45, hjust=1, vjust=1), axis.title.x=element_blank(), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+ ylab("Proportion of Core Communities (%)") + scale_fill_brewer(palette="Dark2") + theme(legend.title = element_blank(), legend.text=element_text(size=12, face="bold"), legend.position=c(0.2,0), legend.justification=c(1,0))
ggplot(data_1e5.ftype.core.samplesum, aes(x = myear, y = core.percent, fill = foam.type)) + geom_boxplot()+ geom_rangeframe() + theme_tufte() + scale_y_continuous(breaks=extended_range_breaks()(round(data_1e5.ftype.core.samplesum$core.percent, 2))) + theme(axis.text.x = element_text(size =12, face="bold", angle=45, hjust=1, vjust=1), axis.title.x=element_blank(), axis.text.y=element_text(size=12, face="bold"), axis.title.y=element_text(size=14, face="bold"))+ ylab("Proportion of Bacteria Core Communities (%)") + scale_fill_brewer(palette="Dark2") + theme(legend.title = element_blank(), legend.text=element_text(size=12, face="bold"), legend.position=c(0.2,0), legend.justification=c(1,0))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Methanogen Relative Abundance (1e-5)") + facet_wrap(~H1, scales="free", ncol=2) + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Methanogen Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
bac_41samples<-read.delim("../mcra/bac_41samples_share_mcra_info_ids.txt", sep="\t")
head(bac_41samples)
dim(bac_41samples)
length(bac_41samples$Long_Sample_ID[bac_41samples$Long_Sample_ID %in% si$Long_Sample_ID])
ls()
dim(core_otu)
head(core_otu)
length(unique(core_otu$core_otus)
)
data_1e5.ftype.core
length(unique(data_1e5.ftype.core.psmelt$OTU))
head(bac_41samples)
head(data_1e5.ftype.core.psmelt$Long_Sample_ID)
bac_41samp_1e5.ftype.core.psmelt<-data_1e5.ftype.core.psmelt[ data_1e5.ftype.core.psmelt$Long_Sample_ID %in% bac_41samples$Long_Sample_ID, ]
dim(bac_41samp_1e5.ftype.core.psmelt)
length(unique(bac_41samp_1e5.ftype.core.psmelt$Long_Sample_ID))
length(unique(bac_41samp_1e5.ftype.core.psmelt$OTU))
saveRDS(bac_41samp_1e5.ftype.core.psmelt, "../mcra/bac_41samples_w_mcra_info_1e5.ftype.core.psmelt.RDS")
data.ftype.core
data_1e5.ftype.core
bac_41samp.ftype.core<-subset_samples(data.ftype.core, Long_Sample_ID %in% bac_41samples$Long_Sample_ID)
bac_41samp.ftype.core
bac_41samp.ftype.core<-prune_taxa(taxa_sums(bac_41samp.ftype.core)>0, bac_41samp.ftype.core)
bac_41samp.ftype.core
saveRDS(bac_41samp.ftype.core, "../mcra/bac_41samples_w_mcra_info.ftype.core.phy.RDS")
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
source("~/Documents/repos/R_code/R_functions/ggplot_nmds_ellipse_and_arrow.R")
ggplot.NMDS.ellipse.arrow
colors
dim(data.mds)
df
ls()
names(si)
dim(si)
si.envfit<-si[, c("MPR_slurry", "Lysine")]
head(si.envfit)
si.envfit<-si[, c("MPR_slurry", "Lysine", "SAMPLES")]
head(si.envfit)
source("~/Documents/repos/R_code/R_functions/nmds_envfit_arrow_extraction.R")
envfit.df<-mds.envfit.arrows(data.mds, si.envfit, "SAMPLES")
dim(envfit.df)
head(envfit.df)
ggplot.NMDS.ellipse.arrow(data.mds, envfit.df, si$foam.type, colors)
ggplot.NMDS.ellipse.arrow(data.mds, envfit.df[1, ], si$foam.type, colors)
pdf("../Manuscript/figures_and_tables/bacteria_ellipse_mpr_arrow.pdf", height=10, width=10)
ggplot.NMDS.ellipse.arrow(data.mds, envfit.df[1, ], si$foam.type, colors)
dev.off()
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
dim(data_1e5.ftype.core.genus.ft)
head(data_1e5.ftype.core.genus.ft)
p1
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
library(ggplot2(
library(ggplot2)
library(ggthemes)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_sdl, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_sdl, geom="errorbar", mult=1, position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + geom_boxplot()+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_sdl, geom="errorbar", mult=1, position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_sdl, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_sdl, geom="errorbar", mult=1, position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + geom_boxplot()+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
head(data_1e5.ftype.core.genus.ft.bf100)
dim(data_1e5.ftype.core.genus.ft.bf100)
unique(data_1e5.ftype.core.genus.ft.bf100$H1
)
data_1e5.ftype.core.genus.ft.bf100<-data_1e5.ftype.core.genus.ft.bf100[order(-data_1e5.ftype.core.genus.ft.bf100$H1), ]
dim(data_1e5.ftype.core.genus.ft.bf100)
head(data_1e5.ftype.core.genus.ft.bf100)
data_1e5.ftype.core.genus.ft.bf100<-data_1e5.ftype.core.genus.ft.bf100[order(-data_1e5.ftype.core.genus.ft.bf100$foam.type), ]
data_1e5.ftype.core.genus.ft.bf100<-data_1e5.ftype.core.genus.ft.bf100[order(-data_1e5.ftype.core.genus.ft.bf100$total), ]
dim(data_1e5.ftype.core.genus.ft.bf100)
head(data_1e5.ftype.core.genus.ft.bf100, n= 40)
data_1e5.ftype.core.genus.ft.bf100$most_abundant <- data_1e5.ftype.core.genus.ft.bf100$H1
head(data_1e5.ftype.core.genus.ft.bf100)
data_1e5.ftype.core.genus.ft.bf100$most_abundant<-gsub("\\bCrust>No Foam>Foam\\b", "Crust", data_1e5.ftype.core.genus.ft.bf100$most_abundant)
data_1e5.ftype.core.genus.ft.bf100$most_abundant<-gsub("\\bCrust>Foam>No Foam\\b", "Crust", data_1e5.ftype.core.genus.ft.bf100$most_abundant)
data_1e5.ftype.core.genus.ft.bf100$most_abundant<-gsub("\\bFoam>Crust>No Foam\\b", "Foam", data_1e5.ftype.core.genus.ft.bf100$most_abundant)
data_1e5.ftype.core.genus.ft.bf100$most_abundant<-gsub("\\bFoam>No Foam>Crust\\b", "Foam", data_1e5.ftype.core.genus.ft.bf100$most_abundant)
data_1e5.ftype.core.genus.ft.bf100$most_abundant<-gsub("\\bNo Foam>Crust>Foam\\b", "No-foam", data_1e5.ftype.core.genus.ft.bf100$most_abundant)
data_1e5.ftype.core.genus.ft.bf100$most_abundant<-gsub("\\bNo Foam>Foam>Crust\\b", "No-foam", data_1e5.ftype.core.genus.ft.bf100$most_abundant)
head(data_1e5.ftype.core.genus.ft.bf100)
unique( data_1e5.ftype.core.genus.ft.bf100$most_abundant)
dim( data_1e5.ftype.core.genus.ft.bf100)
head(data_1e5.ftype.core.genus.ft)
head(data_1e5.ftype.core.psmelt)
dim(data_1e5.ftype.core.psmelt)
library(plyr)
count(data_1e5.ftype.core.genus.ft, "most_abundant")
names(data_1e5.ftype.core.genus.ft)
names(data_1e5.ftype.core.genus.ft.bf100)
count(data_1e5.ftype.core.genus.ft.bf100, "most_abundant")
count(data_1e5.ftype.core.genus.ft.bf100, c("Sample", "most_abundant"))
head(data_1e5.ftype.core.genus.ft.bf100)
head(data_1e5.ftype.core.genus.ft)
head(bac_genus_bf_min20)
library(dplyr)
data_1e5.ftype.core.genus.ft.bf20<-data_1e5.ftype.core.genus.ft %>% filter( genus %in% bac_genus_bf_min20$genus)
head(data_1e5.ftype.core.genus.ft.bf20)
data_1e5.ftype.core.genus.ft.bf20<-data_1e5.ftype.core.genus.ft %>% filter( genus %in% bac_genus_bf_min20$genus) %>% data.frame()
head(data_1e5.ftype.core.genus.ft.bf20)
data_1e5.ftype.core.genus.ft.bf20<-merge(data_1e5.ftype.core.genus.ft.bf20, bac_genus_bf_min20, "genus")
head(data_1e5.ftype.core.genus.ft.bf20)
data_1e5.ftype.core.genus.ft.bf20$most_abundant<-data_1e5.ftype.core.genus.ft.bf20$H1
head(data_1e5.ftype.core.genus.ft.bf20)
dim(data_1e5.ftype.core.genus.ft.bf20)
dim(data_1e5.ftype.core.genus.ft.bf100)
dim(bac_genus_bf_min20)
dim(bac_genus_bf_min100)
dim(bac_genus_bf)
(bac_genus_bf)
rm(data_1e5.ftype.core.genus.ft.bf20)
head(data_1e5.ftype.core.genus.ft.bf100)
plyr::count(data_1e5.ftype.core.genus.ft.bf100, c("foam.type", "most_abundant"))
p1<-data_1e5.ftype.core.genus.ft[1:15, ] %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~H1, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~most_abundant, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
str(data_1e5.ftype.core.genus.ft.bf100)
data_1e5.ftype.core.genus.ft.bf100$most_abundant<-factor(data_1e5.ftype.core.genus.ft.bf100$most_abundant, levels=c(
"No-foam", "Crust", "Foam"))
head(data_1e5.ftype.core.genus.ft.bf100)
p1<-data_1e5.ftype.core.genus.ft.bf100[1:15, ] %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~most_abundant, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
length(unique(data_1e5.ftype.core.genus.ft.bf100$genus))
head(unique(data_1e5.ftype.core.genus.ft.bf100$genus))
bf100.top10.genus<-unique(data_1e5.ftype.core.genus.ft.bf100$genus)
bf100.top10.genus
bf100.top10.genus<-unique(data_1e5.ftype.core.genus.ft.bf100$genus)[1:10]
bf100.top10.genus
p1<-data_1e5.ftype.core.genus.ft.bf100[ data_1e5.ftype.core.genus.ft.bf100$genus %in% bf100.top10.genus, ] %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~most_abundant, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
bac_genus_bf
head(data_1e5.ftype.core.genus.ft.bf100)
bf100.top10.genus<-ddply(data_1e5.ftype.core.genus.ft.bf100, .(genus), summarise, all_total = sum(total))
head(bf100.top10.genus)
(bf100.top10.genus)
p1<-data_1e5.ftype.core.genus.ft.bf100[ data_1e5.ftype.core.genus.ft.bf100$genus %in% bf100.top10.genus$genus, ] %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~most_abundant, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
p1<-data_1e5.ftype.core.genus.ft.bf100[ data_1e5.ftype.core.genus.ft.bf100$genus %in% bf100.top10.genus$genus[1:10], ] %>% ungroup() %>% arrange(as.integer(genus)) %>% ggplot(aes(x=genus, y=total, fill = foam.type))
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~most_abundant, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
pdf("../Manuscript/figures_and_tables/bacteria_top10_abundant_most_significant_genera.pdf", height=10, width=15)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~most_abundant, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
dev.off()
pdf("../Manuscript/figures_and_tables/bacteria_top10_abundant_most_significant_genera.pdf", height=8, width=10)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_normal, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~most_abundant, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
dev.off()
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(ggplot2)
pdf("../Manuscript/figures_and_tables/bacteria_top10_abundant_most_significant_genera.pdf", height=8, width=10)
p1 + stat_summary(fun.y = mean, geom="bar", position="dodge") + stat_summary(fun.data = mean_cl_boot, geom="errorbar", position = position_dodge(width = 0.90), width = 0.2)+ theme_classic() +theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + scale_fill_brewer(palette="Dark2") + theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14, face="bold")) + theme(axis.text.x=element_text(size=14, face="bold", angle = 90, hjust = 1, vjust = 0.5), axis.text.y=element_text(angle = 90, hjust =0.5, vjust =0.5, size=12, face="bold")) + ylab("Bacteria Relative Abundance (1e-5)") + facet_grid(~most_abundant, scales="free", space="free") + theme(legend.position = "none")+ theme(strip.text.x=element_text(angle=90, face="bold", size=14))
dev.off()
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
source("/Users/fanyang/Documents/repos/R_code/R_functions/ordisurf_extraction.R")
mpr.env<-data.frame(si[, c("SAMPLES", "MPR_slurry")])
mpr.ordisf<-ordi.sf(data.mds, mpr.env, "SAMPLES")
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
ls()
colors
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, colors)
ft.colors<-colors
pdf("../Manuscript/figures_and_tables/bac_ellipse_mpr_ordisurf.pdf", height=8, width=8)
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
dev.off()
save.image("all_16s_physeq.RData")
savehistory("temp.R")

load("all_16s_physeq.RData")
library(phyloseq)
ls(pattern = "core")
dim(c.core.otu)
dim(f.core.otu)
dim(nf.core.otu)
m(core_otu)
head(nf.core.otu)
install.packages("VennDiagram")
library(Vennerable)
source("http://bioconductor.org/biocLite.R")
biocLite("graph")
install.packages("Vennerable", repos = "http://R-Forge.R-project.org")
biocLite("RBGL")
install.packages("Vennerable", repos = "http://R-Forge.R-project.org")
library(Vennerable)
ls(pattern = "core")
core.otu.venn <- list(nf.core.otu$core, c.core.otu$core, f.core.otu$core)
class(core.otu.venn)
Vcore.otu.venn <- Venn(core.otu.venn)
Vcore.otu.venn
head(core.otu.venn)
core.otu.venn <- list(No-foam = nf.core.otu$core, Crust = c.core.otu$core, Foam = f.core.otu$core)
core.otu.venn <- list(No-foam=c(nf.core.otu$core), Crust = c(c.core.otu$core), Foam = c(f.core.otu$core))
head(core.otu.venn)
names(core.otu.venn)
head(nf.core.otu)
names(core.otu.venn) <- c("No-foam", "Crust", "Foam")
head(core.otu.venn)
core.otu.venn$No-foam
core.otu.venn$`No-foam`
Vcore.otu.venn <- Venn(core.otu.venn)
Vcore.otu.venn
plot(Vcore.otu.venn)
plot(Vcore.otu.venn, doWeight =F)
colors
plot(Vcore.otu.venn, doWeight =F, col = colors)
?plot
plot(Vcore.otu.venn, doWeight =F, Faces=F)
plot(Vcore.otu.venn, doWeight =F)
plot(Vcore.otu.venn, doWeight =F, faces = F)
?plotVenn
library(VennDiagram)
venn.diagram(core.otu.venn)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.pdf")
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, cat.position = 0.5, margin=0.2)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, cat.position = 1, margin=0.2)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, cat.pos = 0, margin=0.2)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, cat.pos = 1, margin=0.2)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2)
?venn.diagram
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2, cat.dist = 0.1)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2, cat.dist = 0.07)
colors
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.1, offset = 0.05)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.1, offset = 0.1)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2, cat.dist = 0.02)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2, cat.dist = 0.05)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2, cat.pos = c(-0.1, 0, 0.1))
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2, cat.dist = 0.03)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2, cat.dist = 0.05)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=2, margin=0.2, cat.dist = 0.07)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=1.5, margin=0.2)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=1.5, margin=0.2, ext.text =F)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, fill = colors, alpha=0.5, cex=1.5, cat.cex=1.5, margin=0.2, ext.text =T)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.pos = 0)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1, lwd = c(2, 2, 2))
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1, lwd = 3)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1, lwd = 5)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1, lwd = 5, cat.pos = c(0.1, 0.05, 0.05)
)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1, lwd = 5, cat.dist = c(0.1, 0.05, 0.05))
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1, lwd = 5, cat.dist = c(0.01, 0.01, 0.01))
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1, lwd = 5, cat.dist = c(0.03, 0.02, 0.02))
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1, lwd = 5, cat.dist = c(0.05, 0.04, 0.04))
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1, lwd = 5, cat.dist = c(0.06, 0.04, 0.03))
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=1.5, margin=0.1, lwd = 5, cat.dist = c(0.06, 0.04, 0.03), cat.cex=2)
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=2, margin=0.1, lwd = 5, cat.dist = c(0.06, 0.04, 0.03))
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=2, margin=0.1, lwd = 5, cat.dist = c(0.07, 0.05, 0.04))
venn.diagram(core.otu.venn, filename= "~/Downloads/test.tiff", col=colors, cex=1.5, cat.cex=2, margin=0.1, lwd = 5, cat.dist = c(0.08, 0.06, 0.04))
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
source("/Users/fanyang/Documents/repos/R_code/R_functions/ordisurf_extraction.R")
head(mpr.ordisf)
mpr.ordisf<-ordi.sf(data.mds, mpr.env, "SAMPLES")
head(mpr.ordisf)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
pdf("../Manuscript/figures_and_tables/bac_ellipse_mpr_ordisurf.pdf", height=8, width=8)
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
dev.off()
min(mpr.ordisf$z)
max(mpr.ordisf$z)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
pdf("../Manuscript/figures_and_tables/bac_ellipse_mpr_ordisurf.pdf", height=8, width=8)
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
dev.off()
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
source("/Users/fanyang/Documents/repos/R_code/R_functions/ggplot_nmds_ordisurf.R")
pdf("../Manuscript/figures_and_tables/bac_ellipse_mpr_ordisurf.pdf", height=8, width=8)
ggplot.NMDS.ordisurf(data.mds, mpr.ordisf, si$foam.type, ft.colors)
dev.off()
save.image("all_16s_physeq.RData")
savehistory("temp.R")
load("all_16s_physeq.RData")
library(phyloseq)
library(plyr)
dim(si)
count(si, "id")
count(si, "foam.type")
C <- read.delim("~/Documents/repos/pit_foaming_manuscript/Text_outputs/data.taxmin5.sample_goods_coverage.txt")
head(C)
dim(si)
dim(C)
si<-merge(si, C, "SAMPLES")
dim(si)
ddply(si, .(foam.type), summarise, avg=mean(sample_sums))
ddply(si, .(foam.type), summarise, avg=median(sample_sums))
ddply(si, .(foam.type), summarise, avg=mean(C))
savehistory("temp.R")
load("all_16s_physeq.RData")
library(ggplot2)
ggplot(otu_432_1665, aes(y=Abundance, x=SCFA, shape = genus, linetype=genus)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(2, 15)) +
 facet_grid(~foam.type) + 
 geom_smooth(method="lm", formula=y ~ poly(x,2), se=F, color = "black") + 
 theme_bw() + 
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
library(plyr)
names(otu_432_1665)
ddply(otu_432_1665, .(SCFA, foam.type, OTU, genus), summarise, function(x) boxplot.stats(x$Abundance)$out)
ddply(otu_432_1665, .(SCFA, foam.type, OTU, genus), function(x) boxplot.stats(x$Abundance)$out)
ddply(otu_432_1665, .(foam.type, OTU, genus), function(x) boxplot.stats(x$Abundance)$out)
boxplot.stats(unique(otu_432_1665$SCFA))$out
boxplot.stats(unique(otu_432_1665[otu_432_1665$foam.type == "No-foam", ]$SCFA))$out
boxplot.stats(unique(otu_432_1665[otu_432_1665$foam.type == "Crust", ]$SCFA))$out
boxplot.stats(unique(otu_432_1665[otu_432_1665$foam.type == "Foam", ]$SCFA))$out
ddply(otu_432_1665, .(foam.type, OTU, genus), function(x) boxplot.stats(x$Abundance))
head(otu_432_1665$foam.type
)
head(otu_432_1665$OTU)
unique(otu_432_1665$OTU)
ddply(otu_432_1665, .(foam.type, OTU, genus), function(x) boxplot.stats(Abundance, x)$out)
ddply(otu_432_1665, .(foam.type, OTU, genus), function(x) boxplot.stats(x)$out)
ddply(otu_432_1665, .(foam.type, OTU, genus), function(Abundance) boxplot.stats(Abundance)$out)
head(otu_432_1665$Abundance)
ddply(otu_432_1665, .(foam.type, OTU), summarise, outliers = boxplot.stats(Abundance)$out)
ddply(otu_432_1665, .(foam.type, genus), summarize, outliers = boxplot.stats(Abundance)$out)
dim(otu_432_1665)
otu_432_1665.outliers <- ddply(otu_432_1665, .(foam.type, genus), summarize, outliers = boxplot.stats(Abundance)$out)
dim(otu_432_1665.outliers)
scfa.outliers <- ddply(otu_432_1665, .(foam.type, genus), summarize, outliers = boxplot.stats(SCFA)$out)
scfa.outliers
dim(otu_432_1665)
names(otu_432_1665)
otu_432_1665_toplot<-otu_432_1665[, c("OTU", "SAMPLES", "foam.type", "genus", "SCFA","Abundance")]
head( otu_432_1665_toplot)
otu_432_1665_toplot$outlier<-ifelse(otu_432_1665_toplot$SCFA %in% scfa.outliers$outliers, "Outliers", ifelse(otu_432_1665_toplot$Abundance %in% otu_432_1665.outliers$outliers, "Outliers", "Non-outliers"))
head( otu_432_1665_toplot)
str( otu_432_1665_toplot)
unique( otu_432_1665_toplot$outlier)
myformula <- formula(Abundance~poly(SCFA, 2))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers", ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(Abundance~poly(SCFA, 2))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA), ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(Abundance~poly(SCFA))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA), ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(Abundance~poly(SCFA, 2))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA), ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(Abundance~poly(SCFA))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA) & otu_432_1665_toplot$OTU == "OTU_432", ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(Abundance~poly(SCFA, 2))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA) & otu_432_1665_toplot$OTU == "OTU_1665", ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
head(otu_432_1665_toplot)
str(otu_432_1665_toplot)
unique(otu_432_1665_toplot$outlier)
otu_432_1665_toplot$outlier<-factor(otu_432_1665_toplot$outlier, levels=c("Outliers", "Non-outliers")
)
head(otu_432_1665_toplot)
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_grid(~genus) + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665", ], method="lm", formula=y ~ poly(x), se=F) + 
 theme_bw() + 
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
dim(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA) & otu_432_1665_toplot$OTU == "OTU_1665", ])
otu_432_1665_toplot$outlier<-factor(otu_432_1665_toplot$outlier, levels=c("Non-outliers", "Outliers"))
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_grid(~genus) + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_bw() + 
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
otu_432_1665_toplot$facet<-ifelse(otu_432_1665_toplot$Abundance > 20000, "two", ifelse(otu_432_1665_toplot$OTU=="OTU_432" & otu_432_1665_toplot$Abundance >5000, "two", "one"))
head(otu_432_1665_toplot$facet<-ifelse(otu_432_1665_toplot$Abundance > 20000, "two", ifelse(otu_432_1665_toplot$OTU=="OTU_432" & otu_432_1665_toplot$Abundance >5000, "two", "one"))
head(otu_432_1665_toplot)
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_grid(facet~genus) + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_bw() + 
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
otu_432_1665_toplot$facet<-factor(otu_432_1665_toplot$facet, levels=c("two", "one"))
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_grid(facet~genus) + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_bw() + 
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_grid(facet~genus, scale = "free_y") + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_bw() + 
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_grid(genus~., scale = "free_y") + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_bw() + 
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_wrap(~genus, ncol = 2, scale = "free_y") + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_wrap(~genus, ncol = 2, scale = "free_y") + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.8))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_wrap(~genus, ncol = 2, scale = "free_y") + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
pdf("../Manuscript/figures_and_tables/SCFA_vs_lacto_turicibacter.pdf", width=10,height=6)
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_wrap(~genus, ncol = 2, scale = "free_y") + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
dev.off()
pdf("../Manuscript/figures_and_tables/SCFA_vs_lacto_turicibacter.pdf", width=10,height=6)
ggplot(otu_432_1665_toplot, aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_wrap(~genus, ncol = 2, scale = "free_y") + 
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
dev.off()
names(otu_15145)
otu_15145_toplot<-otu_15145[, c("OTU", "SAMPLES", "LCFA" "genus")]
otu_15145_toplot<-otu_15145[, c("OTU", "SAMPLES", "LCFA", "genus")]
lcfa.outliers <- ddply(otu_15145_toplot, .(foam.type, genus), summarize, outliers = boxplot.stats(LCFA)$out)
otu_15145.outliers <- ddply(otu_15145_toplot, .(foam.type, genus), summarize, outliers = boxplot.stats(Abundance)$out)
otu_15145_toplot<-otu_15145[, c("OTU", "SAMPLES", "LCFA", "genus", "foam.type")]
lcfa.outliers <- ddply(otu_15145_toplot, .(foam.type, genus), summarize, outliers = boxplot.stats(LCFA)$out)
otu_15145.outliers <- ddply(otu_15145_toplot, .(foam.type, genus), summarize, outliers = boxplot.stats(Abundance)$out)
otu_15145_toplot$outlier<-ifelse(otu_15145_toplot$LCFA %in% lcfa.outliers$outliers, "Outliers", ifelse(otu_15145_toplot$Abundance %in% otu_15145.outliers$outliers, "Outliers", "Non-outliers"))
otu_15145_toplot<-otu_15145[, c("OTU", "SAMPLES", "LCFA", "genus", "foam.type", "Abundance")]
lcfa.outliers <- ddply(otu_15145_toplot, .(foam.type, genus), summarize, outliers = boxplot.stats(LCFA)$out)
otu_15145.outliers <- ddply(otu_15145_toplot, .(foam.type, genus), summarize, outliers = boxplot.stats(Abundance)$out)
otu_15145_toplot$outlier<-ifelse(otu_15145_toplot$LCFA %in% lcfa.outliers$outliers, "Outliers", ifelse(otu_15145_toplot$Abundance %in% otu_15145.outliers$outliers, "Outliers", "Non-outliers"))
head(otu_15145_toplot)
otu_15145_toplot$outlier<-factor(otu_15145_toplot$outlier, levels=c("Non-outliers", "Outliers"))
head(otu_15145_toplot)
myformula <- formula(Abundance~poly(LCFA))
models<-dlply(otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers" & !is.na(otu_15145_toplot$LCFA), ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(Abundance~poly(LCFA, 2))
models<-dlply(otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers" & !is.na(otu_15145_toplot$LCFA), ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, shape = outlier, color = foam.ty[e)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, shape = outlier, color = foam.ty[e)) + 
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
myformula <- formula(LCFA~poly(Abundance, 2))
models<-dlply(otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers" & !is.na(otu_15145_toplot$LCFA), ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(LCFA~poly(Abundance))
models<-dlply(otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers" & !is.na(otu_15145_toplot$LCFA), ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
dim(otu_15145_toplot)
dim(otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers" & !is.na(otu_15145_toplot$LCFA), ])
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
# geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
otu_15145_toplot$facet<-ifelse(otu_15145_toplot$LCFA > 3000, "two", "one")
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
# geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(~facet, scale="free")+
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
otu_15145_toplot<-otu_15145_toplot[! is.na(otu_15145_toplot$LCFA, ]
otu_15145_toplot<-otu_15145_toplot[! is.na(otu_15145_toplot$LCFA), ]
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
# geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(~facet, scale="free")+
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(~facet, scale="free")+
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
#geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(~foam.type, scale="free")+
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
#geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(~foam.type, scale="free")+
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
myformula <- formula(LCFA~poly(Abundance, 2))
models<-dlply(otu_15145_toplot, .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(Abundance~poly(LCFA, 2))
models<-dlply(otu_15145_toplot, .(foam.type, genus), function(x) lm(myformula, x))
l_ply(models, summary, .print=T)
boxplot.stats(otu_15145_toplot$Abundance)$out
boxplot.stats(otu_15145_toplot$LCFA)$out
lcfa.outliers
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(~foam.type, scale="free")+
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
save.image("all_16s_physeq.RData")
savehistory("temp.R")
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU=="OTU_432", aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ x, se=F) +  
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU=="OTU_432", ], aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ x, se=F) +  
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU=="OTU_432", ], aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) +  
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
myformula <- formula(Abundance~poly(SCFA, 2))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA) & otu_432_1665_toplot$OTU == "OTU_432", ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(SCFA~poly(Abundance, 2))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA) & otu_432_1665_toplot$OTU == "OTU_432", ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU=="OTU_432", ], aes(x=Abundance, y=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) +  
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
dev.off()
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU=="OTU_432", ], aes(x=Abundance, y=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) +  
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
myformula <- formula(Abundance~poly(SCFA, 2))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA) & otu_432_1665_toplot$OTU == "OTU_432", ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU=="OTU_432", ], aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) +  
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
pdf("../Manuscript/figures_and_tables/SCFA_vs_lacto.pdf", width=6,height=6)
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU=="OTU_432", ], aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) +  
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
dev.off()
pdf("../Manuscript/figures_and_tables/SCFA_vs_lacto.pdf", width=8,height=8)
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU=="OTU_432", ], aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) +  
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
dev.off()
pdf("../Manuscript/figures_and_tables/SCFA_vs_lacto.pdf", width=6,height=6)
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU=="OTU_432", ], aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) +  
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.2, 0.7))+
ylab("Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
dev.off()
ggplot(otu_432_1665_toplot, aes(x=Abundance, y=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 facet_wrap(~genus, ncol = 2, scale = "free_y") + 
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
xlab("Relative Abundance (1e-5)") + ylab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665", ], aes(x=Abundance, y=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.1, 0.7))+
xlab("Relative Abundance (1e-5)") + ylab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665", ], aes(x=Abundance, y=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.7))+
xlab("Relative Abundance (1e-5)") + ylab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
pdf("../Manuscript/figures_and_tables/turicibacter_vs_SCFA.pdf", width=6,height=6)
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665", ], aes(x=Abundance, y=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Relative Abundance (1e-5)") + ylab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
dev.off()
pdf("../Manuscript/figures_and_tables/turicibacter_vs_SCFA.pdf", width=6,height=6)
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665", ], aes(x=Abundance, y=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
  geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_1665" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Turicibacter Relative Abundance (1e-5)") + ylab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
dev.off()
pdf("../Manuscript/figures_and_tables/SCFA_vs_lacto.pdf", width=6,height=6)
ggplot(otu_432_1665_toplot[otu_432_1665_toplot$OTU=="OTU_432", ], aes(y=Abundance, x=SCFA, shape = outlier, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 geom_smooth(data = otu_432_1665_toplot[otu_432_1665_toplot$OTU == "OTU_432" & otu_432_1665_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) +  
 theme_classic() + 
 theme(aspect.ratio =1) +
 scale_color_manual(values = colors) +
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.2, 0.7))+
ylab("Lactobacillus Relative Abundance (1e-5)") + xlab("SCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
dev.off()
head(otu_15145_toplot)
outliers<-dlply(otu_15145_toplot, .(foam.type, genus), function(x) compute.bagplot(cbind(x$LCFA, x$Abundance), na.rm=T))
outliers<-data.frame(unlist(outliers))
outliers<-subset(outliers, grepl("pxy.outlier", row.names(outliers)))
dim(outliers)
outliers
ggplot(otu_15145_toplot, aes(y=Abundance, x=LCFA, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(~foam.type, scale="free")+
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Relative Abundance (1e-5)") + xlab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
otu_15145_toplot$outlier<-ifelse(grepl("\\b76.|\\b3577.|\\b630.|\\b141.|\\b12680.|\\b1820.", otu_15145_toplot$LCFA) & grepl("\\b116.|\\b1453.|\\b337.|\\b114.|\\b108.|\\b426.", otu_15145_toplot$Abundance, "Outliers", "Non-outliers")
)
otu_15145_toplot$outlier<-ifelse(grepl("\\b76.|\\b3577.|\\b630.|\\b141.|\\b12680.|\\b1820.", otu_15145_toplot$LCFA) & grepl("\\b116.|\\b1453.|\\b337.|\\b114.|\\b108.|\\b426.", otu_15145_toplot$Abundance), "Outliers", "Non-outliers")
head(otu_15145_toplot)
otu_15145_toplot[otu_15145_toplot$outlier == "Outliers", ]
otu_15145_toplot$outlier<-ifelse(grepl("\\b76.|\\b3577.|\\b630.|\\b141.|\\b12680.|\\b1820.", otu_15145_toplot$LCFA) | grepl("\\b116.|\\b1453.|\\b337.|\\b114.|\\b108.|\\b426.", otu_15145_toplot$Abundance), "Outliers", "Non-outliers")
otu_15145_toplot[otu_15145_toplot$outlier == "Outliers", ]
outliers
outliers<-dlply(otu_15145_toplot, .(foam.type, genus), function(x) compute.bagplot(cbind(x$LCFA, x$Abundance), na.rm=T))
outliers$No-foam.unclassified_Firmicutes
outliers$`No-foam.unclassified_Firmicutes`
outliers$`No-foam.unclassified_Firmicutes`$pxy.outliers
otu_15145_toplot$outlier<-factor(otu_15145_toplot$outlier, levels=c("Non-outliers", "Outliers"))
str(otu_15145_toplot)
myformula <- formula(Abundance~poly(LCFA, 2))
models<-dlply(otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers" & !is.na(otu_15145_toplot$LCFA), ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(LCFA~poly(Abundance, 2))
models<-dlply(otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers" & !is.na(otu_15145_toplot$LCFA), ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.3, 0.9))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
head(otu_15145_toplot)
otu_15145_toplot[otu_15145_toplot$facet, ]
otu_15145_toplot[otu_15145_toplot$facet=="two", ]
otu_15145_toplot$facet<-ifelse(otu_15145_toplot$LCFA > 10000, "two", "one")
otu_15145_toplot$facet<-factor(otu_15145_toplot$facet, levels = c("two", "one"))
str(otu_15145_toplot)
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(facet~., scale="free") +
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
otu_15145_toplot$facet<-ifelse(otu_15145_toplot$LCFA > 3000, "two", "one")
otu_15145_toplot$facet<-factor(otu_15145_toplot$facet, levels = c("two", "one"))
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(facet~., scale="free") +
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(facet~., scale="free") +
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.y = element_blank())
pdf("../Manuscript/figures_and_tables/uc_firmicutes_vs_lcfa.pdf", width=6,height=12)
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(facet~., scale="free") +
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.y = element_blank())
dev.off()
pdf("../Manuscript/figures_and_tables/uc_firmicutes_vs_lcfa.pdf", width=8,height=10)
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(facet~., scale="free") +
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.y = element_blank())
dev.off()
pdf("../Manuscript/figures_and_tables/uc_firmicutes_vs_lcfa.pdf", width=8,height=10)
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(facet~., scale="free", space="free_y") +
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.y = element_blank())
dev.off()
pdf("../Manuscript/figures_and_tables/uc_firmicutes_vs_lcfa.pdf", width=8,height=10)
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
 facet_grid(facet~., scale="free") +
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.y = element_blank())
dev.off()
pdf("../Manuscript/figures_and_tables/uc_firmicutes_vs_lcfa.pdf", width=8,height=10)
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
#facet_grid(facet~., scale="free") +
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.y = element_blank())
dev.off()
pdf("../Manuscript/figures_and_tables/uc_firmicutes_vs_lcfa.pdf", width=8,height=10)
ggplot(otu_15145_toplot, aes(x=Abundance, y=LCFA, color = foam.type, shape = outlier)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
#geom_smooth(method="lm", formula=y ~ poly(x, 2), se=F)+
geom_smooth(data=otu_15145_toplot[otu_15145_toplot$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_classic() + 
facet_grid(facet~., scale="free") +
 scale_color_manual(values=colors)+
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.7, 0.8))+
xlab("Unclassified Firmicutes Relative Abundance (1e-5)") + ylab("LCFA (µg/g)")+
theme(strip.text.y = element_blank())
dev.off()
dim(otu_15145_1665)
head(otu_15145_1665)
outliers<-dlply(otu_15145_1665, .(foam.type, genus), function(x) compute.bagplot(cbind(x$unclassified_Firmicutes, x$Turicibacter), na.rm=T))
outliers<-data.frame(unlist(outliers))
outliers<-dlply(otu_15145_1665, .(foam.type), function(x) compute.bagplot(cbind(x$unclassified_Firmicutes, x$Turicibacter), na.rm=T))
outliers<-data.frame(unlist(outliers))
outliers<-subset(outliers, grepl("pxy.outlier", row.names(outliers)))
dim(outliers)
head(outliers)
outliers$ftype<-data.frame(do.call('rbind', strsplit(as.character(row.names(outliers), ".", fixed=T)))[, 1]
outliers$ftype<-data.frame(do.call('rbind', strsplit(as.character(row.names(outliers), ".", fixed=T))))[, 1]
outliers$ftype<-data.frame(do.call('rbind', strsplit(as.character(row.names(outliers), split =".", fixed=T))))[, 1]
outliers$ftype<-data.frame(do.call('rbind', strsplit(as.character(row.names(outliers)), ".", fixed=T)))[, 1]
head(outliers)
unique(outliers$ftype)
outliers<-dlply(otu_15145_1665, .(foam.type), function(x) compute.bagplot(cbind(x$unclassified_Firmicutes, x$Turicibacter), na.rm=T))
outliers$`No-foam`$pxy.outlier
outliers.df<-rbind(outliers$`No-foam`$pxy.outlier, outliers$Crust$pxy.outlier, outliers$Foam$pxy.outlier)
dim(outliers.df)
head(outliers.df)
otu_15145_1665$outlier<-ifelse(otu_15145_1665$unclassified_Firmicutes %in% outliers.df[, 1] & otu_15145_1665$Turicibacter %in% outliers.df[, 2], "Outliers", "Non-outliers")
head(otu_15145_1665)
otu_15145_1665[otu_15145_1665$outlier == "Outliers", ]
otu_15145_1665$outlier<-factor(otu_15145_1665$outlier, levels=c("Non-outliers", "Outliers"))
myformula <- formula(unclassified_Firmicutes ~ poly(Turicibacter, 2))
models<-dlply(otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], .(foam.type), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(Turicibacter ~ poly(unclassified_Firmicutes, 2))
models<-dlply(otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], .(foam.type), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(Turicibacter ~ poly(unclassified_Firmicutes))
models<-dlply(otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], .(foam.type), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(unclassified_Firmicutes ~ poly(Turicibacter))
models<-dlply(otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], .(foam.type), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
ggplot(otu_15145_1665, aes(y=unclassified_Firmicutes, x=Turicibacter, fill=foam.type, shape = outlier, color=foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
 geom_smooth(data = otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x, 2), se=F) + 
 theme_bw() + 
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
xlab("Turicibacter Relative Abundance (1e-5)") + ylab("Unclassified_Firmicutes Relative abundance (1e-5)")+
theme(aspect.ratio=1)+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_1665, aes(y=unclassified_Firmicutes, x=Turicibacter, fill=foam.type, shape = outlier, color=foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
 geom_smooth(data = otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
 theme_classic() + 
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
xlab("Turicibacter Relative Abundance (1e-5)") + ylab("Unclassified_Firmicutes Relative abundance (1e-5)")+
theme(aspect.ratio=1)+
theme(strip.text.x = element_text(size=15, face="bold"))
ggplot(otu_15145_1665, aes(x=unclassified_Firmicutes, y=Turicibacter, fill=foam.type, shape = outlier, color=foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
 geom_smooth(data = otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
 theme_classic() + 
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
xlab("Turicibacter Relative Abundance (1e-5)") + ylab("Unclassified_Firmicutes Relative abundance (1e-5)")+
theme(aspect.ratio=1)+
theme(strip.text.x = element_text(size=15, face="bold"))
myformula <- formula(Turicibacter ~ poly(unclassified_Firmicutes))
models<-dlply(otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], .(foam.type), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(Turicibacter ~ poly(unclassified_Firmicutes, 2))
models<-dlply(otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], .(foam.type), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
ggplot(otu_15145_1665, aes(x=unclassified_Firmicutes, y=Turicibacter, fill=foam.type, shape = outlier, color=foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
 geom_smooth(data = otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
 theme_classic() + 
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Turicibacter Relative Abundance (1e-5)") + xlab("Unclassified Firmicutes Relative abundance (1e-5)")+
theme(aspect.ratio=1)+
theme(strip.text.x = element_text(size=15, face="bold"))
pdf("../Manuscript/figures_and_tables/turicibacter_vs_uc_firmicutes.pdf", width=8,height=8)
ggplot(otu_15145_1665, aes(x=unclassified_Firmicutes, y=Turicibacter, fill=foam.type, shape = outlier, color=foam.type)) + 
 geom_point(size = 3) +
 scale_shape_manual(values=c(16:18)) +
 scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
 geom_smooth(data = otu_15145_1665[otu_15145_1665$outlier == "Non-outliers", ], method="lm", formula=y ~ poly(x), se=F) + 
 theme_classic() + 
theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
theme(legend.title=element_blank(),legend.text=element_text(size=15), legend.position=c(0.13, 0.9))+
ylab("Turicibacter Relative Abundance (1e-5)") + xlab("Unclassified Firmicutes Relative abundance (1e-5)")+
theme(aspect.ratio=1)+
theme(strip.text.x = element_text(size=15, face="bold"))
dev.off()
myformula <- formula(SCFA~poly(Abundance))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA) & otu_432_1665_toplot$OTU == "OTU_1665", ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
myformula <- formula(SCFA~poly(Abundance, 2))
models<-dlply(otu_432_1665_toplot[otu_432_1665_toplot$outlier == "Non-outliers" & !is.na(otu_432_1665_toplot$SCFA) & otu_432_1665_toplot$OTU == "OTU_1665", ], .(foam.type, genus), function(x) lm(myformula, x))
ldply(models, coef)
l_ply(models, summary, .print=T)
save.image("all_16s_physeq.RData")
savehistory("temp.R")

plot_bar(data_1e5.ftype, "SAMPLES", fill="phylum") + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")
