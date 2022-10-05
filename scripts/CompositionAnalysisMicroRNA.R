library(reshape2)
library(ggplot2)
library(phyloseq)
library(vegan)
library(ggpubr)
library(indicspecies)
library(pheatmap)
library(ape)
library(ggtree)
library(vegan)
library("MANOVA.RM")

## code for formatting RDP heirarchial file
inDIR="../raw_data/"
rdp.microRNA <- read.table(file=paste0(inDIR,"/AIL_micro_RNA.txt"),header=T,row.names=1,sep="\t")[,-c(1)]
rdp.microRNA$name <- gsub("unidentified","Un",rdp.microRNA$name)
rdp.microRNA.data <- rdp.microRNA[rdp.microRNA$rank %in% c("phylum","class","order","family","genus"),]
rdp.microRNA.data$rank <- as.factor(as.character(rdp.microRNA.data$rank))
colnames(rdp.microRNA.data) <- gsub("X|.fasta|.nonchimeras","",colnames(rdp.microRNA.data))
ncmb.phen <- read.table(file=file=paste0(inDIR,"/NCMB_Sept2018.tab"),sep="\t",header=T,row.names=1)
rdp.microRNA.count <- rdp.microRNA.data[,colnames(rdp.microRNA.data) %in% rownames(ncmb.phen)]
ncmb.cage <- read.table(file=file=paste0(inDIR,"/ncmb_cage_reformat.tsv"),sep="\t",header=T,row.names=1)
rdp.microRNA.count <- rdp.microRNA.count[,colnames(rdp.microRNA.count) %in% rownames(ncmb.cage)]
a <- NULL;
for(i in levels(rdp.microRNA.data$rank)){
a <- c(a,names(which(apply(rdp.microRNA.count[which(rdp.microRNA.data$rank == i),],2,sum) > 5000)))
}
a <- unique(a)
rdp.microRNA.count <- rdp.microRNA.count[,(colnames(rdp.microRNA.count) %in% a)]
ncmb.cage.microRNA <- ncmb.cage[match(colnames(rdp.microRNA.count),rownames(ncmb.cage)),]
ncmb.cage.microRNA$cage <- as.factor(as.character(ncmb.cage.microRNA$cage))
ncmb.phen.microRNA <- ncmb.phen[match(colnames(rdp.microRNA.count),rownames(ncmb.phen)),]
microRNA.taxa.count <- list()
for(i in levels(rdp.microRNA.data$rank)){
microRNA.taxa.count[[i]] <- rdp.microRNA.count[which(rdp.microRNA.data$rank == i),]
microRNA.taxa.count[[i]] <- microRNA.taxa.count[[i]][which(apply(microRNA.taxa.count[[i]],1,function(x){length(which(x >=5))} >= 20)),]
}

## estimate Phylum abundance

microRNA.phylum <- microRNA.taxa.count$phylum
rownames(microRNA.phylum) <- as.character(rdp.microRNA$name[match(rownames(microRNA.phylum),rownames(rdp.microRNA))])
phylum.melt <- melt(apply(microRNA.phylum,2,function(x){x/sum(x)}))
colnames(phylum.melt) = c("taxa","Diet","rl")
phylum.melt$Diet <- ncmb.phen.microRNA$Diet[match(phylum.melt$Diet,rownames(ncmb.phen.microRNA))]
phylum.melt$Diet <- factor(ifelse(phylum.melt$Diet == 0,"Cal",ifelse(phylum.melt$Diet == 1,"Con","Wes")),levels=c("Cal","Con","Wes"))

## estimate Genus abundance

microRNA.genus <- microRNA.taxa.count$genus
rownames(microRNA.genus) <- as.character(rdp.microRNA$name[match(rownames(microRNA.genus),rownames(rdp.microRNA))])
microRNA.genus.rl <- apply(microRNA.genus,2,function(x){x/sum(x)})
microRNA.genus.10 <- names(sort(apply(apply(microRNA.genus,2,function(x){x/sum(x)}),1,mean),decreasing=T)[1:10])
microRNA.genus.others <- apply(microRNA.genus.rl[!(rownames(microRNA.genus.rl) %in% microRNA.genus.10),],2,sum)
genus.melt <- melt(rbind(microRNA.genus.rl[rownames(microRNA.genus.rl) %in% microRNA.genus.10,],Others=microRNA.genus.others))
colnames(genus.melt) <- c("taxa","Diet","rl")
genus.melt$Diet <- ncmb.phen.microRNA$Diet[match(genus.melt$Diet,rownames(ncmb.phen.microRNA))]
genus.melt$Diet <- factor(ifelse(genus.melt$Diet == 0,"Cal",ifelse(genus.melt$Diet == 1,"Con","Wes")),levels=c("Cal","Con","Wes"))
genus.melt$taxa <- factor(genus.melt$taxa,levels=c(microRNA.genus.10,"Others"))
genus.melt$taxa <- gsub("_1_1","",genus.melt$taxa)

## combine data and plot

taxa.melt <- rbind(phylum.melt,genus.melt)
pdf("taxa_plot.pdf",height=4,width=5)
ggplot(taxa.melt, aes(x=taxa, y=rl*100,fill=Diet)) + scale_fill_manual(values=c("green", "cornflowerblue", "red")) + geom_boxplot(lwd=0.2,outlier.size=0.2) +labs(x="Phylum   Genus", y = "Relative Abundance (%)") + theme_classic() + theme(axis.text.x = element_text(size=10,angle = 60, hjust = 1,color="black"),axis.text.y = element_text(size=10,color="black"),legend.title = element_text(size=10),legend.text = element_text(size=10),axis.title.x = element_text(size=12,face="bold"),axis.title.y = element_text(size=10),legend.position="top") + ylim(0, 100) + geom_vline(xintercept = nrow(microRNA.phylum)+0.5,linetype = "dotted",color="black",size=0.8)
dev.off()


##### create data for lefse #####
options(scipen=999)
microRNA.rl <- list()
for(i in levels(rdp.microRNA.data$rank)){
	microRNA.rl[[i]] <- apply(microRNA.taxa.count[[i]],2,function(x){round(x/sum(x),digit=8)})
}
microRNA.rl.data <- rbind(microRNA.rl$phylum,microRNA.rl$class,microRNA.rl$order,microRNA.rl$family,microRNA.rl$genus)
rdp.microRNA.lin <- read.table(file="../raw_data/AIL_micro_RNA.txt",header=T,row.names=1,sep="\t",stringsAsFactors=F)[,c(1:2)]
rownames(microRNA.rl.data) <- gsub("Root;rootrank;|\\|$","",gsub(";domain;|;phylum;|;class;|;order;|;family;|;genus;","|",rdp.microRNA.lin$lineage[match(rownames(microRNA.rl.data),rownames(rdp.microRNA.lin))]))

ncmb.microRNA <- data.frame(Generation=as.factor(paste("G",ncmb.phen.microRNA$Generation,sep="")),Sex = as.factor(ifelse(ncmb.phen.microRNA$Sex == 0,"M","F")),Diet = as.factor(ifelse(ncmb.phen.microRNA$Diet == 0,"Cal",ifelse(ncmb.phen.microRNA$Diet == 1,"Con","Wes"))))
write.table(rbind(t(ncmb.microRNA),microRNA.rl.data),file="MicroRNA_lefse.txt",sep="\t",col.names=T,row.names=T,quote=F)

##command for lefse
#/home/yask/software/nsegata-lefse-82605a2ae7b7/format_input.py MicroRNA_lefse_sex.txt MicroRNA_lefse_sex.in -c 1 -o 1000000
#/home/yask/software/nsegata-lefse-82605a2ae7b7/run_lefse.py MicroRNA_lefse_gen.in MicroRNA_lefse_gen.res

#~/software/graphlan_commit_6ca8735/export2graphlan/export2graphlan.py -i MicroRNA_lefse_diet.in -o MicroRNA_lefse_diet_only.res -t diet_tree.txt -a diet_annot.txt --title "AIL MICROBIOTA(RNA) DIET" --annotations 2,3,4,5 --external_annotations 6 --fname_row 0 --skip_rows 1,2 --ftop 200
#~/software/graphlan_commit_6ca8735/graphlan_annotate.py --annot diet_annot.txt diet_tree.txt diet_outtree.txt
#~/software/graphlan_commit_6ca8735/graphlan.py --dpi 300 --size 7.0 --format pdf diet_outtree.txt diet_outtree.pdf

######################################################################################################################
norm_data <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}


for(i in names(microRNA.taxa.count)){
  rownames(microRNA.taxa.count[[i]]) <- rdp.microRNA.data$name[match(rownames(microRNA.taxa.count[[i]]),rownames(rdp.microRNA.data))]
}

microRNA.taxa.rl = microRNA.taxa.count
for(i in names(microRNA.taxa.count)){
  microRNA.taxa.rl[[i]] = apply(microRNA.taxa.count[[i]],2,function(x){x/sum(x)})
}

microRNA.taxa.count.df <- NULL
microRNA.taxa.rl.df <- NULL
taxlevels = c( "phylum","class","order","family","genus")
for(i in taxlevels){
  if(i == "phylum"){
    microRNA.taxa.count.df <- microRNA.taxa.count[[i]]
    microRNA.taxa.count.df <- cbind(Taxname = rownames(microRNA.taxa.count[[i]]),
                                    class = i, microRNA.taxa.count.df)
    microRNA.taxa.rl.df <- as.data.frame(microRNA.taxa.rl[[i]])
    microRNA.taxa.rl.df <- cbind(Taxname = rownames(microRNA.taxa.rl[[i]]),
                                 class = i, microRNA.taxa.rl.df)
  }else{
    df_tmp_ct <- microRNA.taxa.count[[i]]
    df_tmp_ct <- cbind(Taxname = rownames(microRNA.taxa.count[[i]]),
                       class = i, df_tmp_ct)
    df_tmp_rl <- as.data.frame(microRNA.taxa.rl[[i]])
    df_tmp_rl <- cbind(Taxname = rownames(microRNA.taxa.rl[[i]]),
                       class = i, df_tmp_rl)
    microRNA.taxa.count.df <- rbind(microRNA.taxa.count.df,df_tmp_ct)
    microRNA.taxa.rl.df <- rbind(microRNA.taxa.rl.df,df_tmp_rl)
    
  }
}
write.table(microRNA.taxa.count.df, file = "microRNA_RDP_count.txt",col.names = T,row.names = F,sep = "\t",quote = F)
write.table(microRNA.taxa.rl.df, file = "microRNA_RDP_relAbund.txt",col.names = T,row.names = F,sep = "\t",quote = F)

ncmb.phen.microRNA.conf <- ncmb.phen.microRNA[,c("Generation","Sex","Diet")]
ncmb.phen.microRNA.conf$Diet2 = factor(ifelse(ncmb.phen.microRNA.conf$Diet == 0,"cal",
                                              ifelse(ncmb.phen.microRNA.conf$Diet == 1, "con","wes")),
                                       levels=c("cal","con","wes"))

p_out_all <- NULL
for(cl in taxlevels){
  test_1 <- NULL; test_out <- NULL
  test_1 <- apply(microRNA.taxa.rl[[cl]],1,function(x){
    x = norm_data(x)
    pval = coef(summary(glm(x ~ Generation + Sex + Diet,data = ncmb.phen.microRNA.conf)))[-1,4]
    pval = data.frame(Phen = names(pval),P = as.numeric(pval))
    pval_pair = melt(pairwise.t.test(x, ncmb.phen.microRNA.conf$Diet2)$p.value)
    pval_pair = subset(pval_pair, !is.na(value))
    pval_all = data.frame(Phen=paste0(pval_pair$Var2,"-",pval_pair$Var1),P = pval_pair$value)
    pval_all <-rbind(pval,pval_all)
  })
  
  for(i in 1:length(test_1)){
    if(i == 1){
      test_out <- data.frame(matrix(NA,ncol = nrow(test_1[[i]]),nrow = length(test_1)))
      rownames(test_out) <- names(test_1)
      colnames(test_out) <- c(test_1[[i]]$Phen)
      test_out[i,] <- test_1[[i]]$P
    }else{
      test_out[i,] <- test_1[[i]]$P
    }
  }
  test_out <- cbind(test_out,as.data.frame(apply(test_out,2,function(x){p.adjust(x,method = "BH")})))
  test_out <- cbind(TaxName= rownames(test_out),class=cl,test_out)
  if(cl == 1){
    p_out_all <- test_out
  }else{
    p_out_all <- rbind(p_out_all,test_out)
  }
}
write.table(p_out_all, file = "microRNA_RDP_Pvals.txt",col.names = T,row.names = F,sep = "\t",quote = F)

##load data
microRNA.sp <- read.table(file="../raw_data/RNAotus10K.txt",skip=1,header=T,row.names=1,comment.char="",sep="\t")
colnames(microRNA.sp) <- gsub("X","",colnames(microRNA.sp))
t.microRNA.sp <- t(microRNA.sp)
ncmb.phen <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
microRNA.sp <- microRNA.sp[,colnames(microRNA.sp) %in% rownames(ncmb.phen)]
microRNA.sp <- microRNA.sp[which(apply(microRNA.sp,1,sum) > 0),]
ncmb.phen.microRNA <- ncmb.phen[match(colnames(microRNA.sp),rownames(ncmb.phen)),]
ncmb.microRNA <- data.frame(Generation=as.factor(paste("G",ncmb.phen.microRNA$Generation,sep="")),Sex = as.factor(ifelse(ncmb.phen.microRNA$Sex == 0,"M","F")),Diet = as.factor(ifelse(ncmb.phen.microRNA$Diet == 0,"Cal",ifelse(ncmb.phen.microRNA$Diet == 1,"Con","Wes"))))
rownames(ncmb.microRNA) <- rownames(ncmb.phen.microRNA)

## Indicator species analysis

t.microRNA.sp <- t(microRNA.sp)
t.microRNA.sp <- t.microRNA.sp[,which(apply(t.microRNA.sp,2,function(x){length(which(x > 10))}) > 1)]
indval.diet <- multipatt(as.data.frame(t(microRNA.sp)),cluster=ncmb.microRNA$Diet,control = how(nperm=999),func="IndVal.g",duleg=TRUE)
indval.gen <- multipatt(as.data.frame(t(microRNA.sp)),cluster=ncmb.microRNA$Generation,control = how(nperm=999),func="IndVal.g",duleg=TRUE)
indval.sex <- multipatt(as.data.frame(t(microRNA.sp)),cluster=ncmb.microRNA$Sex,control = how(nperm=999),func="IndVal.g",duleg=TRUE)
write.table(file="diet_indi_microRNA.txt",subset(indval.diet$sign,p.value < 0.01),sep="\t",row.names=T,col.names=T,quote=F)
write.table(file="gen_indi_microRNA.txt",subset(indval.gen$sign,p.value < 0.01),sep="\t",row.names=T,col.names=T,quote=F)
write.table(file="sex_indi_microRNA.txt",subset(indval.gen$sign,p.value < 0.01),sep="\t",row.names=T,col.names=T,quote=F)
indval.diet.p <- subset(indval.diet$sign,p.value < 0.01)
indval.gen.p <- subset(indval.gen$sign,p.value < 0.01)
indval.sex.p <- subset(indval.gen$sign,p.value < 0.01)
exclude.rownames <- unique(c(rownames(indval.sex.p),rownames(indval.gen.p)))
t.microRNA.diet <- t.microRNA.sp[,match(rownames(indval.diet.p)[!(rownames(indval.diet.p) %in% exclude.rownames)],colnames(t.microRNA.sp))]
indval.diet.only <- indval.diet.p[!(rownames(indval.diet.p) %in% exclude.rownames),]
t.microRNA.diet <- t.microRNA.sp[,match(rownames(indval.diet.only),colnames(t.microRNA.sp))]
micro.sp.annot <- read.table(file="all.otusRNA.mseq",header=T,skip=1,comment.char="",sep="\t")[,c("X.query","Species")]
micro.sp.annot$X.query <- unlist(lapply(strsplit(as.character(micro.sp.annot$X.query),split=";"),function(x){x <- x[1]}))
t.microRNA.diet.annot <- micro.sp.annot[match(colnames(t.microRNA.diet),micro.sp.annot$X.query),]
t.microRNA.diet.annot$Species <- as.factor(as.character(t.microRNA.diet.annot$Species))
pheatmap.sp.in <- t(aggregate(t(aggregate(t(t.microRNA.diet),list(t.microRNA.diet.annot$Species),FUN=sum)[,-1]),list(ncmb.microRNA$Diet),FUN=mean)[-1])
rownames(pheatmap.sp.in) <- levels(t.microRNA.diet.annot$Species)
colnames(pheatmap.sp.in) <- levels(ncmb.microRNA$Diet)
pheatmap(pheatmap.sp.in,scale="row",cluster_cols = FALSE,color=colorRampPalette(c("white","grey","orange"))(10),cellwidth=12,cellheight=12,file="Indispecies_MicroRNA.pdf")


micro.sp.annot2 <- read.table(file="../raw_data/all.otusRNA.mseq",header=T,skip=1,comment.char="",sep="\t",stringsAsFactors=F)[,c("NCBIv2.2b.Kingdom","Phylum","Class","Order","Family","Genus","Species")]
micro.sp.annot2 <- micro.sp.annot2[match(rownames(pheatmap.sp.in),micro.sp.annot2$Species),]
rownames(micro.sp.annot2) <- rownames(pheatmap.sp.in);
sp.names <- unlist(lapply(strsplit(micro.sp.annot2$Species," "),function(x){
	genus <- x[1];
	species <- paste(x[2:length(x)],collapse =" ")
	genus <- gsub("\\[|\\]","",genus);
	g.abbr <- paste(unlist(strsplit(genus,"")[[1]][1]),".",sep="");
	paste(g.abbr,species,sep=" ")
}))
ggtree.data.in = as.phylo.hclust(hclust(taxa2dist(x=micro.sp.annot2)))
gr.indisp <- apply(pheatmap.sp.in,1,which.max)
groups.diet <- list(Cal=names(gr.indisp)[which(gr.indisp == 1)],Con=names(gr.indisp)[which(gr.indisp == 2)],Wes=names(gr.indisp)[which(gr.indisp == 3)])
ggtree.data.in <- groupOTU(ggtree.data.in,groups.diet)
ggtree.data.in$tip.label <- sp.names
ggtree.data.in$edge.length <- log1p(ggtree.data.in$edge.length)
pdf("Indispecies_MicroRNA_tree.pdf",width = 15, height = 15)
ggtree(ggtree.data.in, aes(color=group),layout="circular",ladderize=F) + geom_tiplab2(aes(angle=angle)) + scale_color_manual(values=c("green", "cornflowerblue","red")) + theme(legend.position="right")
dev.off()

###calculate alpha diversity 

microRNA.phylo <- phyloseq(sample_data(ncmb.microRNA),otu_table(microRNA.sp,taxa_are_rows = TRUE))
microRNA.phylo.alpha <- estimate_richness(microRNA.phylo)
microRNA.phylo.alpha.norm <- as.data.frame(apply(microRNA.phylo.alpha,2,function(x){(x-(min(x)))/((max(x))-(min(x)))}))
microRNA.phylo.alpha.norm$Diet <- ncmb.microRNA$Diet
microRNA.phylo.alpha.melt <- melt(microRNA.phylo.alpha.norm,variable.name = "Diet")
colnames(microRNA.phylo.alpha.melt) <- c("Diet","Index","value")
pdf("alpha_diversity.pdf")
ggplot(microRNA.phylo.alpha.melt, aes(x=Index, y=value,fill=Diet)) + scale_fill_manual(values=c("green", "cornflowerblue", "red")) + geom_boxplot(outlier.size=0.5) +labs(x="alpha-diversity metrics", y = "Scaled Indices") + theme_classic() + theme(axis.text.x = element_text(size=12,angle = 60, hjust = 1, face="bold"),axis.text.y = element_text(size=12,face="bold"),legend.title = element_text(size=12,face="bold"),legend.text = element_text(size=12,face="bold"),axis.title.x = element_text(size=14, face="bold"),axis.title.y = element_text(size=14, face="bold"))
dev.off()

###calculate beta diversity
## Capscale Analysis
microRNA.cap <- capscale(t(microRNA.sp) ~ ncmb.microRNA$Diet + Condition(ncmb.microRNA$Generation),distance = "bray")
microRNA.cap.scores <- as.data.frame(scores(microRNA.cap,scaling = "sites", correlation = TRUE,hill=TRUE)$sites)
microRNA.cap.scores$Diet <- ncmb.microRNA$Diet
microRNA.cap1 <- round(as.numeric(eigenvals(microRNA.cap)['CAP1']),digit=2)
microRNA.cap2 <- round(as.numeric(eigenvals(microRNA.cap)['CAP2']),digit=2)
pdf("beta_diversity.pdf")
ggscatter(microRNA.cap.scores, x = "CAP1", y = "CAP2", color = "Diet",palette=c("green","cornflowerblue","red"),ellipse=TRUE,conf.int = TRUE, ellipse.level = 0.5,conf.int.level=0.95,ellipse.type="norm",star.plot=T,star.plot.lwd=0.1,star.plot.lty=3,size=0.1,ellipse.alpha=0.3,ellipse.border.remove=TRUE,xlab=paste("CAP1 [",microRNA.cap1,"%]",sep=""),ylab=paste("CAP2 [",microRNA.cap2,"%]",sep=""),font.label=c(12,"bold","black"))
dev.off()

#### Test for significance ###
anova(microRNA.cap,by="margin",model="reduced")
microRNA.manova <- MANOVA.wide(cbind(CAP1,CAP2) ~ Diet,data=microRNA.cap.scores, iter = 1000, CPU = 1)
simCI(microRNA.manova, contrast = "pairwise", type = "Tukey",base=2)



