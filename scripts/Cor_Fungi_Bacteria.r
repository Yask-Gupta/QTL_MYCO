library(pheatmap)
library(reshape2)


################# correlate mycobiome and standing community of bacteria (microbiome DNA) ################## 
rdp.microDNA <- read.table(file="../raw_data/AIL_micro_DNA.txt",header=T,row.names=1,sep="\t")[,-c(1)]
rdp.microDNA$name <- gsub("unidentified","Un",rdp.microDNA$name)
rdp.microDNA.data <- rdp.microDNA[rdp.microDNA$rank %in% c("phylum","class","order","family","genus"),]
rdp.microDNA.data$rank <- as.factor(as.character(rdp.microDNA.data$rank))
colnames(rdp.microDNA.data) <- gsub("X|.fasta|.nonchimeras","",colnames(rdp.microDNA.data))
ncmb.phen <- read.table(file="/home/yask/artem_data/MicrMycoGmice/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.microDNA.count <- rdp.microDNA.data[,colnames(rdp.microDNA.data) %in% rownames(ncmb.phen)]
ncmb.cage <- read.table(file="/home/yask/artem_data/MicrMycoGmice/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.microDNA.count <- rdp.microDNA.count[,colnames(rdp.microDNA.count) %in% rownames(ncmb.cage)]
microDNA.taxa.count <- list()
for(i in levels(rdp.microDNA.data$rank)){
microDNA.taxa.count[[i]] <- rdp.microDNA.count[which(rdp.microDNA.data$rank == i),]
}

rdp.myco <- read.table(file="../raw_data/AIL_myco_hier.tab",header=T,row.names=1,sep="\t")[,-c(1)]
rdp.myco$name <- gsub("unidentified","Un",rdp.myco$name)
rdp.myco.data <- rdp.myco[rdp.myco$rank %in% c("phylum","class","order","family","genus"),]
rdp.myco.data$rank <- as.factor(as.character(rdp.myco.data$rank))
colnames(rdp.myco.data) <- gsub("X|.fasta","",colnames(rdp.myco.data))
ncmb.phen <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.myco.count <- rdp.myco.data[,colnames(rdp.myco.data) %in% rownames(ncmb.phen)]
ncmb.cage <- read.table(file="/home/yask/artem_data/MicrMycoGmice/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.myco.count <- rdp.myco.count[,colnames(rdp.myco.count) %in% rownames(ncmb.cage)]
myco.taxa.count <- list()
for(i in levels(rdp.myco.data$rank)){
myco.taxa.count[[i]] <- rdp.myco.count[which(rdp.myco.data$rank == i),]
}

myco.g <- myco.taxa.count$genus
rownames(myco.g) <- rdp.myco$name[match(rownames(myco.taxa.count$genus),rownames(rdp.myco))]

microDNA.g <- microDNA.taxa.count$genus
rownames(microDNA.g) <- rdp.microDNA$name[match(rownames(microDNA.taxa.count$genus),rownames(rdp.microDNA))]

myco.g <- myco.g[,colnames(myco.g) %in% colnames(microDNA.g)]
microDNA.g <- microDNA.g[,match(colnames(myco.g),colnames(microDNA.g))]
myco.microDNA.g <- rbind(myco.g,microDNA.g)
myco.microDNA.g <- myco.microDNA.g[which(apply(myco.microDNA.g,1,function(x){sum(x > 0)}) >= 0.1*dim(myco.microDNA.g)[2]),]
rownames(myco.microDNA.g) <- gsub(" ","_",rownames(myco.microDNA.g))
write.table(myco.microDNA.g,file="Myco_MicroDNA.tsv",sep="\t",quote=F)

####### run fastspar on output file #######
######## this is bash commands for fastspar #######
#fastspar --otu_table Myco_MicroDNA.tsv --correlation Myco_MicroDNA_median_correlation.tsv --covariance Myco_MicroDNA_median_covariance.tsv
#mkdir Myco_MicroDNA_bootstrap_counts
#fastspar_bootstrap --otu_table Myco_MicroDNA.tsv --number 1000 --prefix Myco_MicroDNA_bootstrap_counts/Myco_MicroDNA
#mkdir Myco_MicroDNA_bootstrap_correlation

#for i in {0..999}
#do
#	fastspar --otu_table Myco_MicroDNA_bootstrap_counts/Myco_MicroDNA_${i}.tsv --correlation Myco_MicroDNA_bootstrap_correlation/Myco_MicroDNA_cor_${i}.tsv --covariance Myco_MicroDNA_bootstrap_correlation/Myco_MicroDNA_cov_${i}.tsv -i 5 -t 5
#done
#fastspar_pvalues --otu_table Myco_MicroDNA.tsv --correlation Myco_MicroDNA_median_correlation.tsv --prefix Myco_MicroDNA_bootstrap_correlation/Myco_MicroDNA_cor_ --permutations 1000 --outfile Myco_MicroDNA_pvalues.tsv


cor.myco.microDNA <- read.table(file="Myco_MicroDNA_median_correlation.tsv",comment.char="",header=T,sep="\t",row.names=1)
cor.myco.microDNA.p <- read.table(file="Myco_MicroDNA_pvalues.tsv",comment.char="",header=T,sep="\t",row.names=1)
all.p.microDNA <- as.numeric(as.matrix(cor.myco.microDNA.p))
p.thresh.microDNA <- max(all.p.microDNA[which(p.adjust(all.p.microDNA,method="BH") < 0.05)])
rdp.myco <- read.table(file="../raw_data/AIL_myco_hier.tab",header=T,row.names=1,sep="\t")[,-c(1)]
rdp.myco$name <- gsub("unidentified","Un",rdp.myco$name)
rdp.myco.data <- rdp.myco[rdp.myco$rank %in% c("phylum","class","order","family","genus"),]
rdp.myco.data$rank <- as.factor(as.character(rdp.myco.data$rank))
colnames(rdp.myco.data) <- gsub("X|.fasta","",colnames(rdp.myco.data))
ncmb.phen <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.myco.count <- rdp.myco.data[,colnames(rdp.myco.data) %in% rownames(ncmb.phen)]
ncmb.cage <- read.table(file="../raw_data/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.myco.count <- rdp.myco.count[,colnames(rdp.myco.count) %in% rownames(ncmb.cage)]
myco.taxa.count <- list()
for(i in levels(rdp.myco.data$rank)){
myco.taxa.count[[i]] <- rdp.myco.count[which(rdp.myco.data$rank == i),]
}

myco.g <- myco.taxa.count$genus
rownames(myco.g) <- rdp.myco$name[match(rownames(myco.taxa.count$genus),rownames(rdp.myco))]
idx.myco <- rownames(cor.myco.microDNA) %in% rownames(myco.g)
cor.myco.microDNA <- cor.myco.microDNA[idx.myco,!idx.myco]
cor.myco.microDNA.p <- cor.myco.microDNA.p[idx.myco,!idx.myco]
idx.myco.sig <- which(apply(cor.myco.microDNA.p,1,function(x){sum(x < p.thresh.microDNA)}) > 0)
idx.microDNA.sig <- which(apply(cor.myco.microDNA.p,2,function(x){sum(x < p.thresh.microDNA)}) > 0)
cor.myco.microDNA <- cor.myco.microDNA[idx.myco.sig,idx.microDNA.sig]
cor.myco.microDNA.p <- cor.myco.microDNA.p[idx.myco.sig,idx.microDNA.sig]
cor.myco.microDNA.1 <- cor.myco.microDNA
cor.myco.microDNA.1[cor.myco.microDNA.p > p.thresh.microDNA] <- 0
pdf("cor_myco_microDNA.pdf")
pheatmap(t(cor.myco.microDNA.1),scale="none", color=colorRampPalette(c('dark blue','white','dark red'))(7),cex=0.8)
dev.off()
save(cor.myco.microDNA,cor.myco.microDNA.p,cor.myco.microDNA.1,file="cor.myco.microDNA.RData")

cor.myco.microDNA.melt <- melt(as.matrix(cor.myco.microDNA))
cor.myco.microDNA.p.melt <- melt(as.matrix(cor.myco.microDNA.p))
cor.myco.microDNA.df <- cor.myco.microDNA.melt[which(cor.myco.microDNA.p.melt$value <= p.thresh.microDNA),]
colnames(cor.myco.microDNA.df) <- c("Fungi","Bacteria","rho")
write.table(cor.myco.microDNA.df,file="cor.myco.microDNA.cyto.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(cor.myco.microDNA.df[paste0(cor.myco.microDNA.df$Fungi,"-",cor.myco.microDNA.df$Bacteria) %in% paste0(cor.myco.microRNA.df$Fungi,"-",cor.myco.microRNA.df$Bacteria),],file="cor.myco.microDNA.microRNA.cyto.txt",sep="\t",quote=F)


################# correlate mycobiome and standing community of bacteria (microbiome DNA) ################## 

rdp.microRNA <- read.table(file="../raw_data/AIL_micro_RNA.txt",header=T,row.names=1,sep="\t")[,-c(1)]
rdp.microRNA$name <- gsub("unidentified","Un",rdp.microRNA$name)
rdp.microRNA.data <- rdp.microRNA[rdp.microRNA$rank %in% c("phylum","class","order","family","genus"),]
rdp.microRNA.data$rank <- as.factor(as.character(rdp.microRNA.data$rank))
colnames(rdp.microRNA.data) <- gsub("X|.fasta|.nonchimeras","",colnames(rdp.microRNA.data))
ncmb.phen <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.microRNA.count <- rdp.microRNA.data[,colnames(rdp.microRNA.data) %in% rownames(ncmb.phen)]
ncmb.cage <- read.table(file="/home/yask/artem_data/MicrMycoGmice/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.microRNA.count <- rdp.microRNA.count[,colnames(rdp.microRNA.count) %in% rownames(ncmb.cage)]
microRNA.taxa.count <- list()
for(i in levels(rdp.microRNA.data$rank)){
microRNA.taxa.count[[i]] <- rdp.microRNA.count[which(rdp.microRNA.data$rank == i),]
}

rdp.myco <- read.table(file="../raw_data/AIL_myco_hier.tab",header=T,row.names=1,sep="\t")[,-c(1)]
rdp.myco$name <- gsub("unidentified","Un",rdp.myco$name)
rdp.myco.data <- rdp.myco[rdp.myco$rank %in% c("phylum","class","order","family","genus"),]
rdp.myco.data$rank <- as.factor(as.character(rdp.myco.data$rank))
colnames(rdp.myco.data) <- gsub("X|.fasta","",colnames(rdp.myco.data))
ncmb.phen <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.myco.count <- rdp.myco.data[,colnames(rdp.myco.data) %in% rownames(ncmb.phen)]
ncmb.cage <- read.table(file="../raw_data/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.myco.count <- rdp.myco.count[,colnames(rdp.myco.count) %in% rownames(ncmb.cage)]
myco.taxa.count <- list()
for(i in levels(rdp.myco.data$rank)){
myco.taxa.count[[i]] <- rdp.myco.count[which(rdp.myco.data$rank == i),]
}

myco.g <- myco.taxa.count$genus
rownames(myco.g) <- rdp.myco$name[match(rownames(myco.taxa.count$genus),rownames(rdp.myco))]

microRNA.g <- microRNA.taxa.count$genus
rownames(microRNA.g) <- rdp.microRNA$name[match(rownames(microRNA.taxa.count$genus),rownames(rdp.microRNA))]

myco.g <- myco.g[,colnames(myco.g) %in% colnames(microRNA.g)]
microRNA.g <- microRNA.g[,match(colnames(myco.g),colnames(microRNA.g))]
myco.microRNA.g <- rbind(myco.g,microRNA.g)
myco.microRNA.g <- myco.microRNA.g[which(apply(myco.microRNA.g,1,function(x){sum(x > 0)}) >= 0.1*dim(myco.microRNA.g)[2]),]
rownames(myco.microRNA.g) <- gsub(" ","_",rownames(myco.microRNA.g))
write.table(myco.microRNA.g,file="Myco_MicroRNA.tsv",sep="\t",quote=F)

###### bash command similar to DNA #####
#fastspar --otu_table Myco_MicroRNA.tsv --correlation Myco_MicroRNA_median_correlation.tsv --covariance Myco_MicroRNA_median_covariance.tsv
#mkdir Myco_MicroRNA_bootstrap_counts
#fastspar_bootstrap --otu_table Myco_MicroRNA.tsv --number 1000 --prefix Myco_MicroRNA_bootstrap_counts/Myco_MicroRNA
#mkdir Myco_MicroRNA_bootstrap_correlation

#for i in {0..999}
#do
#	fastspar --otu_table Myco_MicroRNA_bootstrap_counts/Myco_MicroRNA_${i}.tsv --correlation Myco_MicroRNA_bootstrap_correlation/Myco_MicroRNA_cor_${i}.tsv --covariance Myco_MicroRNA_bootstrap_correlation/Myco_MicroRNA_cov_${i}.tsv -i 5 -t 5
#done
#fastspar_pvalues --otu_table Myco_MicroRNA.tsv --correlation Myco_MicroRNA_median_correlation.tsv --prefix Myco_MicroRNA_bootstrap_correlation/Myco_MicroRNA_cor_ --permutations 1000 --outfile Myco_MicroRNA_pvalues.tsv
###############################################

library(pheatmap)
cor.myco.microRNA <- read.table(file="Myco_MicroRNA_median_correlation.tsv",comment.char="",header=T,sep="\t",row.names=1)
cor.myco.microRNA.p <- read.table(file="Myco_MicroRNA_pvalues.tsv",comment.char="",header=T,sep="\t",row.names=1)
all.p.microRNA <- as.numeric(as.matrix(cor.myco.microRNA.p))
p.thresh.microRNA <- max(all.p.microRNA[which(p.adjust(all.p.microRNA,method="BH") < 0.05)])
rdp.myco <- read.table(file="../raw_data/AIL_myco_hier.tab",header=T,row.names=1,sep="\t")[,-c(1)]
rdp.myco$name <- gsub("unidentified","Un",rdp.myco$name)
rdp.myco.data <- rdp.myco[rdp.myco$rank %in% c("phylum","class","order","family","genus"),]
rdp.myco.data$rank <- as.factor(as.character(rdp.myco.data$rank))
colnames(rdp.myco.data) <- gsub("X|.fasta","",colnames(rdp.myco.data))
ncmb.phen <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.myco.count <- rdp.myco.data[,colnames(rdp.myco.data) %in% rownames(ncmb.phen)]
ncmb.cage <- read.table(file="../raw_data/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.myco.count <- rdp.myco.count[,colnames(rdp.myco.count) %in% rownames(ncmb.cage)]
myco.taxa.count <- list()
for(i in levels(rdp.myco.data$rank)){
myco.taxa.count[[i]] <- rdp.myco.count[which(rdp.myco.data$rank == i),]
}

myco.g <- myco.taxa.count$genus
rownames(myco.g) <- rdp.myco$name[match(rownames(myco.taxa.count$genus),rownames(rdp.myco))]
idx.myco <- rownames(cor.myco.microRNA) %in% rownames(myco.g)
cor.myco.microRNA <- cor.myco.microRNA[idx.myco,!idx.myco]
cor.myco.microRNA.p <- cor.myco.microRNA.p[idx.myco,!idx.myco]
idx.myco.sig <- which(apply(cor.myco.microRNA.p,1,function(x){sum(x < p.thresh.microRNA)}) > 0)
idx.microRNA.sig <- which(apply(cor.myco.microRNA.p,2,function(x){sum(x < p.thresh.microRNA)}) > 0)
cor.myco.microRNA <- cor.myco.microRNA[idx.myco.sig,idx.microRNA.sig]
cor.myco.microRNA.p <- cor.myco.microRNA.p[idx.myco.sig,idx.microRNA.sig]
cor.myco.microRNA.1 <- cor.myco.microRNA
cor.myco.microRNA.1[cor.myco.microRNA.p > p.thresh.microRNA] <- 0
pdf("cor_myco_microRNA.pdf")
pheatmap(t(cor.myco.microRNA.1),scale="none", color=colorRampPalette(c('dark blue','white','dark red'))(7),cex=0.8)
dev.off()
save(cor.myco.microRNA,cor.myco.microRNA.p,cor.myco.microRNA.1,file="cor.myco.microRNA.RData")

library(reshape2)
cor.myco.microRNA.melt <- melt(as.matrix(cor.myco.microRNA))
cor.myco.microRNA.p.melt <- melt(as.matrix(cor.myco.microRNA.p))
cor.myco.microRNA.df <- cor.myco.microRNA.melt[which(cor.myco.microRNA.p.melt$value <= p.thresh.microRNA),]
colnames(cor.myco.microRNA.df) <- c("Fungi","Bacteria","rho")
write.table(cor.myco.microRNA.df,file="cor.myco.microRNA.cyto.txt",sep="\t",col.names=T,row.names=F,quote=F)

####### these interaction were visualized in cytoscape ######### 