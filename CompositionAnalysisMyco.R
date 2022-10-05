library(indicspecies)
library(pheatmap)
library(phyloseq)
library(reshape2)
library(ggplot2)
library(vegan)
library(ggpubr)
library("MANOVA.RM")

## code for formatting RDP heirarchial file

rdp.myco <- read.table(file="../raw_data/AIL_myco_hier.tab",header=T,row.names=1,sep="\t")[,-c(1)]
rdp.myco$name <- gsub("unidentified","Un",rdp.myco$name)
rdp.myco.data <- rdp.myco[rdp.myco$rank %in% c("phylum","class","order","family","genus"),]
rdp.myco.data$rank <- as.factor(as.character(rdp.myco.data$rank))
colnames(rdp.myco.data) <- gsub("X|.fasta","",colnames(rdp.myco.data))
ncmb.phen <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.myco.count <- rdp.myco.data[,colnames(rdp.myco.data) %in% rownames(ncmb.phen)]
ncmb.cage <- read.table(file="../raw_data/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.myco.count <- rdp.myco.count[,colnames(rdp.myco.count) %in% rownames(ncmb.cage)]
a <- NULL;
for(i in levels(rdp.myco.data$rank)){
  a <- c(a,names(which(apply(rdp.myco.count[which(rdp.myco.data$rank == i),],2,sum) < 5000)))
}
a <- unique(a)
rdp.myco.count <- rdp.myco.count[,!colnames(rdp.myco.count) %in% a]
for(i in levels(rdp.myco.data$rank)){
  a <- c(a,names(which(apply(rdp.myco.count[which(rdp.myco.data$rank == i),],2,sum) < 5000)))
}
a <- unique(a)
rdp.myco.count <- rdp.myco.count[,!colnames(rdp.myco.count) %in% a]
ncmb.cage.myco <- ncmb.cage[match(colnames(rdp.myco.count),rownames(ncmb.cage)),]
ncmb.cage.myco$cage <- as.factor(as.character(ncmb.cage.myco$cage))
ncmb.phen.myco <- ncmb.phen[match(colnames(rdp.myco.count),rownames(ncmb.phen)),]
myco.taxa.count <- list()
for(i in levels(rdp.myco.data$rank)){
  myco.taxa.count[[i]] <- rdp.myco.count[which(rdp.myco.data$rank == i),]
  myco.taxa.count[[i]] <- myco.taxa.count[[i]][which(apply(myco.taxa.count[[i]],1,function(x){length(which(x >=5))} >= 20)),]
}

## estimate Phylum abundance from rdp output

myco.phylum <- myco.taxa.count$phylum
rownames(myco.phylum) <- as.character(rdp.myco$name[match(rownames(myco.phylum),rownames(rdp.myco))])
phylum.melt <- melt(apply(myco.phylum,2,function(x){x/sum(x)}))
colnames(phylum.melt) = c("taxa","Diet","rl")
phylum.melt$Diet <- ncmb.phen.myco$Diet[match(phylum.melt$Diet,rownames(ncmb.phen.myco))]
phylum.melt$Diet <- factor(ifelse(phylum.melt$Diet == 0,"Cal",ifelse(phylum.melt$Diet == 1,"Con","Wes")),levels=c("Cal","Con","Wes"))

## estimate Genus abundance from rdp output

myco.genus <- myco.taxa.count$genus
rownames(myco.genus) <- as.character(rdp.myco$name[match(rownames(myco.genus),rownames(rdp.myco))])
myco.genus.rl <- apply(myco.genus,2,function(x){x/sum(x)})
myco.genus.10 <- names(sort(apply(apply(myco.genus,2,function(x){x/sum(x)}),1,mean),decreasing=T)[1:10])
myco.genus.others <- apply(myco.genus.rl[!(rownames(myco.genus.rl) %in% myco.genus.10),],2,sum)
genus.melt <- melt(rbind(myco.genus.rl[rownames(myco.genus.rl) %in% myco.genus.10,],Others=myco.genus.others))
colnames(genus.melt) <- c("taxa","Diet","rl")
genus.melt$Diet <- ncmb.phen.myco$Diet[match(genus.melt$Diet,rownames(ncmb.phen.myco))]
genus.melt$Diet <- factor(ifelse(genus.melt$Diet == 0,"Cal",ifelse(genus.melt$Diet == 1,"Con","Wes")),levels=c("Cal","Con","Wes"))
genus.melt$taxa <- factor(genus.melt$taxa,levels=c(myco.genus.10,"Others"))
genus.melt$taxa <- gsub("_1_1","",genus.melt$taxa)

## combine data and plot rdp abundances

taxa.melt <- rbind(phylum.melt,genus.melt)
pdf("../output/taxa_plot_myco.pdf",height=4,width=5)
ggplot(taxa.melt, aes(x=taxa, y=rl*100,fill=Diet)) + scale_fill_manual(values=c("chartreuse4", "cornflowerblue", "firebrick1")) + geom_boxplot(lwd=0.2,outlier.size=0.2) +labs(x="Phylum   Genus", y = "Relative Abundance (%)") + theme_classic() + theme(axis.text.x = element_text(size=10,angle = 60, hjust = 1,color="black"),axis.text.y = element_text(size=10,color="black"),legend.title = element_text(size=10),legend.text = element_text(size=10),axis.title.x = element_text(size=12,face="bold"),axis.title.y = element_text(size=10),legend.position="top") + ylim(0, 100) + geom_vline(xintercept = nrow(myco.phylum)+0.5,linetype = "dotted",color="black",size=0.8)
dev.off()

# OTU level analysis
##load data
myco.sp <- read.table(file=paste0("../raw_data/otu_table_even5k.txt"),skip=1,header=T,row.names=1,comment.char="",sep="\t")
colnames(myco.sp) <- gsub("X","",colnames(myco.sp))
t.myco.sp <- t(myco.sp)
ncmb.phen <- read.table(paste0(inpath,"../raw_data/NCMB_Sept2018.tab"),sep="\t",header=T,row.names=1)
myco.sp <- myco.sp[,colnames(myco.sp) %in% rownames(ncmb.phen)]
myco.sp <- myco.sp[which(apply(myco.sp,1,sum) > 0),]
ncmb.phen.myco <- ncmb.phen[match(colnames(myco.sp),rownames(ncmb.phen)),]
ncmb.myco <- data.frame(Generation=as.factor(paste("G",ncmb.phen.myco$Generation,sep="")),Sex = as.factor(ifelse(ncmb.phen.myco$Sex == 0,"M","F")),Diet = as.factor(ifelse(ncmb.phen.myco$Diet == 0,"Cal",ifelse(ncmb.phen.myco$Diet == 1,"Con","Wes"))))
rownames(ncmb.myco) <- rownames(ncmb.phen.myco)

## Indicator species analysis

t.myco.sp <- t(myco.sp)
t.myco.sp <- t.myco.sp[,which(apply(t.myco.sp,2,function(x){length(which(x > 5))}) > 1)]
indval.diet <- multipatt(as.data.frame(t.myco.sp),cluster=ncmb.myco$Diet,control = how(nperm=999),func="IndVal.g",duleg=TRUE)
indval.gen <- multipatt(as.data.frame(t.myco.sp),cluster=ncmb.myco$Generation,control = how(nperm=999),func="IndVal.g",duleg=TRUE)
write.table(file="../output/diet_indi_myco.txt",subset(indval.diet$sign,p.value < 0.01),sep="\t",row.names=T,col.names=T,quote=F)
write.table(file="../output/gen_indi_myco.txt",subset(indval.gen$sign,p.value < 0.01),sep="\t",row.names=T,col.names=T,quote=F)
indval.diet.p <- subset(indval.diet$sign,p.value < 0.01)
indval.gen.p <- subset(indval.gen$sign,p.value < 0.01)
indval.diet.only <- indval.diet.p[!rownames(indval.diet.p) %in% rownames(indval.gen.p),]
#indval.diet.only <- indval.diet.only[which(apply(indval.diet.only,1,function(x){sum(x[1:3])}) == 1),]
t.myco.diet <- t.myco.sp[,match(rownames(indval.diet.only),colnames(t.myco.sp))]

# use ncbi blast align representative sequence (Myco_repseqs.fasta) against database of ITS from fungi
## output file : NCBI_blast_ITS_OTU.txt
# Make uniq hit based on ncbi output from each OTU make_uniqhits.pl
## output file : HitTable_modified.csv
# extract all the taxanomy for all the NCBI BLAST hits : hits_annot.txt
# Use assign_taxa.pl to annotate OTUs

t.myco.diet.annot <- read.table(file="../output/Indi_otu_sp_annot2.tab",header=T,sep="\t",stringsAsFactors=F)[,c("OTUID","species")]
t.myco.diet.annot <- t.myco.diet.annot[t.myco.diet.annot$OTUID %in% colnames(t.myco.diet),]
t.myco.diet <- t.myco.diet[,match(t.myco.diet.annot$OTUID,colnames(t.myco.diet))]
t.myco.diet.annot$species <- as.factor(t.myco.diet.annot$species)
pheatmap.sp.in <- t(aggregate(t(aggregate(t(t.myco.diet),list(t.myco.diet.annot$species),FUN=sum)[,-1]),list(ncmb.myco$Diet),FUN=mean)[-1])
rownames(pheatmap.sp.in) <- levels(t.myco.diet.annot$species)
colnames(pheatmap.sp.in) <- levels(ncmb.myco$Diet)
pheatmap(pheatmap.sp.in,scale="row",cluster_cols = FALSE,color=colorRampPalette(c("white","grey","orange"))(10),cellwidth=12,cellheight=12,file="Indispecies_Myco.pdf")


###calculate alpha diversity


myco.phylo <- phyloseq(sample_data(ncmb.myco),otu_table(myco.sp,taxa_are_rows = TRUE))
myco.phylo.alpha <- estimate_richness(myco.phylo)
myco.phylo.alpha.norm <- as.data.frame(apply(myco.phylo.alpha,2,function(x){(x-(min(x)))/((max(x))-(min(x)))}))
myco.phylo.alpha.norm$Diet <- ncmb.myco$Diet
myco.phylo.alpha.melt <- melt(myco.phylo.alpha.norm,variable.name = "Diet")
colnames(myco.phylo.alpha.melt) <- c("Diet","Index","value")
pdf("../output/alpha_diversity_myco.pdf")
ggplot(myco.phylo.alpha.melt, aes(x=Index, y=value,fill=Diet)) + scale_fill_manual(values=c("chartreuse4", "cornflowerblue", "firebrick1")) + geom_boxplot(outlier.size=0.5) +labs(x="alpha-diversity metrics", y = "Scaled Indices") + theme_classic() + theme(axis.text.x = element_text(size=12,angle = 60, hjust = 1, face="bold"),axis.text.y = element_text(size=12,face="bold"),legend.title = element_text(size=12,face="bold"),legend.text = element_text(size=12,face="bold"),axis.title.x = element_text(size=14, face="bold"),axis.title.y = element_text(size=14, face="bold"))
dev.off()

###calculate beta diversity
## Capscale Analysis
myco.cap <- capscale(t(myco.sp) ~ ncmb.myco$Diet + Condition(ncmb.myco$Generation),distance = "bray")
myco.cap.scores <- as.data.frame(scores(myco.cap,scaling = "sites", correlation = TRUE,hill=TRUE)$sites)
myco.cap.scores$Diet <- ncmb.myco$Diet
myco.cap1 <- round(as.numeric(eigenvals(myco.cap)['CAP1']),digit=2)
myco.cap2 <- round(as.numeric(eigenvals(myco.cap)['CAP2']),digit=2)
pdf("../output/beta_diversity_myco.pdf")
ggscatter(myco.cap.scores, x = "CAP1", y = "CAP2", color = "Diet",palette=c("chartreuse4", "cornflowerblue", "firebrick1"),ellipse=TRUE,conf.int = TRUE, ellipse.level = 0.5,conf.int.level=0.95,ellipse.type="norm",star.plot=T,star.plot.lwd=0.1,star.plot.lty=3,size=0.1,ellipse.alpha=0.3,ellipse.border.remove=TRUE,xlab=paste("CAP1 [",myco.cap1,"%]",sep=""),ylab=paste("CAP2 [",myco.cap2,"%]",sep=""),font.label=c(12,"bold","black"))
dev.off()

#### Test for significance ###
anova(myco.cap,by="margin",model="reduced")
myco.manova <- MANOVA.wide(cbind(CAP1,CAP2) ~ Diet,data=myco.cap.scores, iter = 1000, CPU = 1)
simCI(myco.manova, contrast = "pairwise", type = "Tukey",base=2)





