library(reshape2)
library(ggplot2)
library(phyloseq)
library(vegan)
library(ggpubr)
library("MANOVA.RM")
library(indicspecies)
library(pheatmap)
## code for formatting RDP heirarchial file
inDIR="../raw_data/"
rdp.microDNA <- read.table(file=paste0(inDIR,"/AIL_micro_DNA.txt"),header=T,row.names=1,sep="\t")[,-c(1)]
rdp.microDNA$name <- gsub("unidentified","Un",rdp.microDNA$name)
rdp.microDNA.data <- rdp.microDNA[rdp.microDNA$rank %in% c("phylum","class","order","family","genus"),]
rdp.microDNA.data$rank <- as.factor(as.character(rdp.microDNA.data$rank))
colnames(rdp.microDNA.data) <- gsub("X|.fasta|.nonchimeras","",colnames(rdp.microDNA.data))
ncmb.phen <- read.table(file=paste0(inDIR,"/NCMB_Sept2018.tab"),sep="\t",header=T,row.names=1)
rdp.microDNA.count <- rdp.microDNA.data[,colnames(rdp.microDNA.data) %in% rownames(ncmb.phen)]
ncmb.cage <- read.table(file=paste0(inDIR,"/ncmb_cage_reformat.tsv"),sep="\t",header=T,row.names=1)
rdp.microDNA.count <- rdp.microDNA.count[,colnames(rdp.microDNA.count) %in% rownames(ncmb.cage)]
a <- NULL;
for(i in levels(rdp.microDNA.data$rank)){
a <- c(a,names(which(apply(rdp.microDNA.count[which(rdp.microDNA.data$rank == i),],2,sum) < 5000)))
}
a <- unique(a)
rdp.microDNA.count <- rdp.microDNA.count[,!colnames(rdp.microDNA.count) %in% a]
for(i in levels(rdp.microDNA.data$rank)){
a <- c(a,names(which(apply(rdp.microDNA.count[which(rdp.microDNA.data$rank == i),],2,sum) < 5000)))
}
a <- unique(a)
rdp.microDNA.count <- rdp.microDNA.count[,!colnames(rdp.microDNA.count) %in% a]
ncmb.cage.microDNA <- ncmb.cage[match(colnames(rdp.microDNA.count),rownames(ncmb.cage)),]
ncmb.cage.microDNA$cage <- as.factor(as.character(ncmb.cage.microDNA$cage))
ncmb.phen.microDNA <- ncmb.phen[match(colnames(rdp.microDNA.count),rownames(ncmb.phen)),]
microDNA.taxa.count <- list()
for(i in levels(rdp.microDNA.data$rank)){
microDNA.taxa.count[[i]] <- rdp.microDNA.count[which(rdp.microDNA.data$rank == i),]
microDNA.taxa.count[[i]] <- microDNA.taxa.count[[i]][which(apply(microDNA.taxa.count[[i]],1,function(x){length(which(x >=5))} >= 20)),]
}

## estimate Phylum abundance

microDNA.phylum <- microDNA.taxa.count$phylum
rownames(microDNA.phylum) <- as.character(rdp.microDNA$name[match(rownames(microDNA.phylum),rownames(rdp.microDNA))])
phylum.melt <- melt(apply(microDNA.phylum,2,function(x){x/sum(x)}))
colnames(phylum.melt) = c("taxa","Diet","rl")
phylum.melt$Diet <- ncmb.phen.microDNA$Diet[match(phylum.melt$Diet,rownames(ncmb.phen.microDNA))]
phylum.melt$Diet <- factor(ifelse(phylum.melt$Diet == 0,"Cal",ifelse(phylum.melt$Diet == 1,"Con","Wes")),levels=c("Cal","Con","Wes"))

## estimate Genus abundance

microDNA.genus <- microDNA.taxa.count$genus
rownames(microDNA.genus) <- as.character(rdp.microDNA$name[match(rownames(microDNA.genus),rownames(rdp.microDNA))])
microDNA.genus.rl <- apply(microDNA.genus,2,function(x){x/sum(x)})
microDNA.genus.10 <- names(sort(apply(apply(microDNA.genus,2,function(x){x/sum(x)}),1,mean),decreasing=T)[1:10])
microDNA.genus.others <- apply(microDNA.genus.rl[!(rownames(microDNA.genus.rl) %in% microDNA.genus.10),],2,sum)
genus.melt <- melt(rbind(microDNA.genus.rl[rownames(microDNA.genus.rl) %in% microDNA.genus.10,],Others=microDNA.genus.others))
colnames(genus.melt) <- c("taxa","Diet","rl")
genus.melt$Diet <- ncmb.phen.microDNA$Diet[match(genus.melt$Diet,rownames(ncmb.phen.microDNA))]
genus.melt$Diet <- factor(ifelse(genus.melt$Diet == 0,"Cal",ifelse(genus.melt$Diet == 1,"Con","Wes")),levels=c("Cal","Con","Wes"))
genus.melt$taxa <- factor(genus.melt$taxa,levels=c(microDNA.genus.10,"Others"))
genus.melt$taxa <- gsub("_1_1","",genus.melt$taxa)

## combine data and plot

taxa.melt <- rbind(phylum.melt,genus.melt)
pdf("taxa_plot.pdf",height=4,width=5)
ggplot(taxa.melt, aes(x=taxa, y=rl*100,fill=Diet)) + scale_fill_manual(values=c("chartreuse4", "cornflowerblue", "firebrick1")) + geom_boxplot(lwd=0.2,outlier.size=0.2) +labs(x="Phylum   Genus", y = "Relative Abundance (%)") + theme_classic() + theme(axis.text.x = element_text(size=10,angle = 60, hjust = 1,color="black"),axis.text.y = element_text(size=10,color="black"),legend.title = element_text(size=10),legend.text = element_text(size=10),axis.title.x = element_text(size=12,face="bold"),axis.title.y = element_text(size=10),legend.position="top") + ylim(0, 100) + geom_vline(xintercept = nrow(microDNA.phylum)+0.5,linetype = "dotted",color="black",size=0.8)
dev.off()


options(scipen=999)

microDNA.rl <- list()
for(i in levels(rdp.microDNA.data$rank)){
	microDNA.rl[[i]] <- apply(microDNA.taxa.count[[i]],2,function(x){round(x/sum(x),digit=8)})
}
microDNA.rl.data <- rbind(microDNA.rl$phylum,microDNA.rl$class,microDNA.rl$order,microDNA.rl$family,microDNA.rl$genus)
rdp.microDNA.lin <- read.table(file="../raw_data/AIL_micro_RNA.txt",header=T,row.names=1,sep="\t",stringsAsFactors=F)[,c(1:2)]
rownames(microDNA.rl.data) <- gsub("Root;rootrank;|\\|$","",gsub(";domain;|;phylum;|;class;|;order;|;family;|;genus;","|",rdp.microDNA.lin$lineage[match(rownames(microDNA.rl.data),rownames(rdp.microDNA.lin))]))

ncmb.microDNA <- data.frame(Generation=as.factor(paste("G",ncmb.phen.microDNA$Generation,sep="")),Sex = as.factor(ifelse(ncmb.phen.microDNA$Sex == 0,"M","F")),Diet = as.factor(ifelse(ncmb.phen.microDNA$Diet == 0,"Cal",ifelse(ncmb.phen.microDNA$Diet == 1,"Con","Wes"))))
write.table(rbind(t(ncmb.microDNA),microDNA.rl.data),file="MicroDNA_lefse.txt",sep="\t",col.names=T,row.names=T,quote=F)



##command for lefse
#/home/yask/software/nsegata-lefse-82605a2ae7b7/format_input.py MicroRNA_lefse_sex.txt MicroRNA_lefse_sex.in -c 1 -o 1000000
#/home/yask/software/nsegata-lefse-82605a2ae7b7/run_lefse.py MicroRNA_lefse_gen.in MicroRNA_lefse_gen.res

#~/software/graphlan_commit_6ca8735/export2graphlan/export2graphlan.py -i MicroDNA_lefse_Diet.in -o MicroDNA_lefse_diet_only.res -t diet_tree.txt -a diet_annot.txt --title "AIL MICROBIOTA(RNA) DIET" --annotations 2,3,4,5 --external_annotations 6 --fname_row 0 --skip_rows 1,2 --ftop 200
#~/software/graphlan_commit_6ca8735/graphlan_annotate.py --annot diet_annot.txt diet_tree.txt diet_outtree.txt
#~/software/graphlan_commit_6ca8735/graphlan.py --dpi 300 --size 7.0 --format pdf diet_outtree.txt diet_outtree.pdf



#############################################################################################


norm_data <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}


for(i in names(microDNA.taxa.count)){
  rownames(microDNA.taxa.count[[i]]) <- rdp.microDNA.data$name[match(rownames(microDNA.taxa.count[[i]]),rownames(rdp.microDNA.data))]
}

microDNA.taxa.rl = microDNA.taxa.count
for(i in names(microDNA.taxa.count)){
  microDNA.taxa.rl[[i]] = apply(microDNA.taxa.count[[i]],2,function(x){x/sum(x)})
}

microDNA.taxa.count.df <- NULL
microDNA.taxa.rl.df <- NULL
taxlevels = c( "phylum","class","order","family","genus")
for(i in taxlevels){
  if(i == "phylum"){
    microDNA.taxa.count.df <- microDNA.taxa.count[[i]]
    microDNA.taxa.count.df <- cbind(Taxname = rownames(microDNA.taxa.count[[i]]),
                                    class = i, microDNA.taxa.count.df)
    microDNA.taxa.rl.df <- as.data.frame(microDNA.taxa.rl[[i]])
    microDNA.taxa.rl.df <- cbind(Taxname = rownames(microDNA.taxa.rl[[i]]),
                                 class = i, microDNA.taxa.rl.df)
  }else{
    df_tmp_ct <- microDNA.taxa.count[[i]]
    df_tmp_ct <- cbind(Taxname = rownames(microDNA.taxa.count[[i]]),
                       class = i, df_tmp_ct)
    df_tmp_rl <- as.data.frame(microDNA.taxa.rl[[i]])
    df_tmp_rl <- cbind(Taxname = rownames(microDNA.taxa.rl[[i]]),
                       class = i, df_tmp_rl)
    microDNA.taxa.count.df <- rbind(microDNA.taxa.count.df,df_tmp_ct)
    microDNA.taxa.rl.df <- rbind(microDNA.taxa.rl.df,df_tmp_rl)
    
  }
}
write.table(microDNA.taxa.count.df, file = "~/microDNA_RDP_count.txt",col.names = T,row.names = F,sep = "\t",quote = F)
write.table(microDNA.taxa.rl.df, file = "~/microDNA_RDP_relAbund.txt",col.names = T,row.names = F,sep = "\t",quote = F)

ncmb.phen.microDNA.conf <- ncmb.phen.microDNA[,c("Generation","Sex","Diet")]
ncmb.phen.microDNA.conf$Diet2 = factor(ifelse(ncmb.phen.microDNA.conf$Diet == 0,"cal",
                                              ifelse(ncmb.phen.microDNA.conf$Diet == 1, "con","wes")),
                                       levels=c("cal","con","wes"))

p_out_all <- NULL
for(cl in taxlevels){
  test_1 <- NULL; test_out <- NULL
  test_1 <- apply(microDNA.taxa.rl[[cl]],1,function(x){
    x = norm_data(x)
    pval = coef(summary(glm(x ~ Generation + Sex + Diet,data = ncmb.phen.microDNA.conf)))[-1,4]
    pval = data.frame(Phen = names(pval),P = as.numeric(pval))
    pval_pair = melt(pairwise.t.test(x, ncmb.phen.microDNA.conf$Diet2)$p.value)
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
write.table(p_out_all, file = "~/microDNA_RDP_Pvals.txt",col.names = T,row.names = F,sep = "\t",quote = F)

##load data
microDNA.sp <- read.table(file="../raw_data/OtusDNA_10k.txt",skip=1,header=T,row.names=1,comment.char="",sep="\t")
colnames(microDNA.sp) <- gsub("X","",colnames(microDNA.sp))
t.microDNA.sp <- t(microDNA.sp)
ncmb.phen <- read.table(file="../NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
microDNA.sp <- microDNA.sp[,colnames(microDNA.sp) %in% rownames(ncmb.phen)]
microDNA.sp <- microDNA.sp[which(apply(microDNA.sp,1,sum) > 0),]
ncmb.phen.microDNA <- ncmb.phen[match(colnames(microDNA.sp),rownames(ncmb.phen)),]
ncmb.microDNA <- data.frame(Generation=as.factor(paste("G",ncmb.phen.microDNA$Generation,sep="")),Sex = as.factor(ifelse(ncmb.phen.microDNA$Sex == 0,"M","F")),Diet = as.factor(ifelse(ncmb.phen.microDNA$Diet == 0,"Cal",ifelse(ncmb.phen.microDNA$Diet == 1,"Con","Wes"))))
rownames(ncmb.microDNA) <- rownames(ncmb.phen.microDNA)

## Indicator species analysis

t.microDNA.sp <- t(microDNA.sp)
t.microDNA.sp <- t.microDNA.sp[,which(apply(t.microDNA.sp,2,function(x){length(which(x > 10))}) > 1)]
indval.diet <- multipatt(as.data.frame(t(microDNA.sp)),cluster=ncmb.microDNA$Diet,control = how(nperm=999),func="IndVal.g",duleg=TRUE)
indval.gen <- multipatt(as.data.frame(t(microDNA.sp)),cluster=ncmb.microDNA$Generation,control = how(nperm=999),func="IndVal.g",duleg=TRUE)
indval.sex <- multipatt(as.data.frame(t(microDNA.sp)),cluster=ncmb.microDNA$Sex,control = how(nperm=999),func="IndVal.g",duleg=TRUE)
write.table(file="diet_indi_microDNA.txt",subset(indval.diet$sign,p.value < 0.01),sep="\t",row.names=T,col.names=T,quote=F)
write.table(file="gen_indi_microDNA.txt",subset(indval.gen$sign,p.value < 0.01),sep="\t",row.names=T,col.names=T,quote=F)
write.table(file="sex_indi_microDNA.txt",subset(indval.gen$sign,p.value < 0.01),sep="\t",row.names=T,col.names=T,quote=F)
indval.diet.p <- subset(indval.diet$sign,p.value < 0.01)
indval.gen.p <- subset(indval.gen$sign,p.value < 0.01)
indval.sex.p <- subset(indval.gen$sign,p.value < 0.01)
exclude.rownames <- unique(c(rownames(indval.sex.p),rownames(indval.gen.p)))
t.microDNA.diet <- t.microDNA.sp[,match(rownames(indval.diet.p)[!(rownames(indval.diet.p) %in% exclude.rownames)],colnames(t.microDNA.sp))]
indval.diet.only <- indval.diet.p[!(rownames(indval.diet.p) %in% exclude.rownames),]
t.microDNA.diet <- t.microDNA.sp[,match(rownames(indval.diet.only),colnames(t.microDNA.sp))]
micro.sp.annot <- read.table(file="all.otusDNA.mseq",header=T,skip=1,comment.char="",sep="\t")[,c("X.query","Species")]
micro.sp.annot$X.query <- unlist(lapply(strsplit(as.character(micro.sp.annot$X.query),split=";"),function(x){x <- x[1]}))
t.microDNA.diet.annot <- micro.sp.annot[match(colnames(t.microDNA.diet),micro.sp.annot$X.query),]
t.microDNA.diet.annot$Species <- as.factor(as.character(t.microDNA.diet.annot$Species))
pheatmap.sp.in <- t(aggregate(t(aggregate(t(t.microDNA.diet),list(t.microDNA.diet.annot$Species),FUN=sum)[,-1]),list(ncmb.microDNA$Diet),FUN=mean)[-1])
rownames(pheatmap.sp.in) <- levels(t.microDNA.diet.annot$Species)
colnames(pheatmap.sp.in) <- levels(ncmb.microDNA$Diet)
pheatmap(pheatmap.sp.in,scale="row",cluster_cols = FALSE,color=colorRampPalette(c("white","grey","orange"))(10),cellwidth=12,cellheight=12,file="Indispecies_MicroDNA.pdf")


micro.sp.annot2 <- read.table(file="../raw_data/all.otusDNA.mseq",header=T,skip=1,comment.char="",sep="\t",stringsAsFactors=F)[,c("NCBIv2.2b.Kingdom","Phylum","Class","Order","Family","Genus","Species")]
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
pdf("Indispecies_MicroDNA_tree.pdf",width = 15, height = 15)
ggtree(ggtree.data.in, aes(color=group),layout="circular",ladderize=F) + geom_tiplab2(aes(angle=angle)) + scale_color_manual(values=c("green", "cornflowerblue","red")) + theme(legend.position="right")
dev.off()


###calculate alpha diversity 

microDNA.phylo <- phyloseq(sample_data(ncmb.microDNA),otu_table(microDNA.sp,taxa_are_rows = TRUE))
microDNA.phylo.alpha <- estimate_richness(microDNA.phylo)
microDNA.phylo.alpha.norm <- as.data.frame(apply(microDNA.phylo.alpha,2,function(x){(x-(min(x)))/((max(x))-(min(x)))}))
microDNA.phylo.alpha.norm$Diet <- ncmb.microDNA$Diet
microDNA.phylo.alpha.melt <- melt(microDNA.phylo.alpha.norm,variable.name = "Diet")
colnames(microDNA.phylo.alpha.melt) <- c("Diet","Index","value")
pdf("alpha_diversity.pdf")
ggplot(microDNA.phylo.alpha.melt, aes(x=Index, y=value,fill=Diet)) + scale_fill_manual(values=c("green", "cornflowerblue", "red")) + geom_boxplot(outlier.size=0.5) +labs(x="alpha-diversity metrics", y = "Scaled Indices") + theme_classic() + theme(axis.text.x = element_text(size=12,angle = 60, hjust = 1, face="bold"),axis.text.y = element_text(size=12,face="bold"),legend.title = element_text(size=12,face="bold"),legend.text = element_text(size=12,face="bold"),axis.title.x = element_text(size=14, face="bold"),axis.title.y = element_text(size=14, face="bold"))
dev.off()

###calculate beta diversity
## Capscale Analysis
t.microDNA.sp <- t(microDNA.sp)
microDNA.cap <- capscale(t.microDNA.sp ~ ncmb.microDNA$Diet + Condition(as.factor(as.character(paste(ncmb.microDNA$Generation,ncmb.microDNA$Sex,sep="+")))),distance = "bray")
microDNA.cap.scores <- as.data.frame(scores(microDNA.cap,scaling = "sites", correlation = TRUE,hill=TRUE)$sites)
microDNA.cap.scores$Diet <- ncmb.microDNA$Diet
microDNA.cap1 <- round(as.numeric(eigenvals(microDNA.cap)['CAP1']),digit=2)
microDNA.cap2 <- round(as.numeric(eigenvals(microDNA.cap)['CAP2']),digit=2)
pdf("beta_diversity.pdf")
ggscatter(microDNA.cap.scores, x = "CAP1", y = "CAP2", color = "Diet",palette=c("green","cornflowerblue","red"),ellipse=TRUE,conf.int = TRUE, ellipse.level = 0.5,conf.int.level=0.95,ellipse.type="norm",star.plot=T,star.plot.lwd=0.1,star.plot.lty=3,size=0.1,ellipse.alpha=0.3,ellipse.border.remove=TRUE,xlab=paste("CAP1 [",microDNA.cap1,"%]",sep=""),ylab=paste("CAP2 [",microDNA.cap2,"%]",sep=""),font.label=c(12,"bold","black"))
dev.off()

#### Test for significance ###

anova(microDNA.cap,by="margin",model="reduced")
microDNA.manova <- MANOVA.wide(cbind(CAP1,CAP2) ~ Diet,data=microDNA.cap.scores, iter = 1000, CPU = 1)
simCI(microDNA.manova, contrast = "pairwise", type = "Tukey",base=2)



