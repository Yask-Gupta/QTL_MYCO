library("phyloseq")
library(DOQTL)

gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}


############ QTL phylum to Genus ##############
rdp.microRNA <- read.table(file="../raw_data/AIL_micro_RNA.txt",header=T,row.names=1,sep="\t")[,-c(1,2)]
rdp.microRNA.data <- rdp.microRNA[rdp.microRNA$rank %in% c("phylum","class","order","family","genus"),]
rdp.microRNA.data$rank <- as.factor(as.character(rdp.microRNA.data$rank))
colnames(rdp.microRNA.data) <- gsub("X|.nonchimeras.fasta","",colnames(rdp.microRNA.data))

ncmb.phen <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.microRNA.count <- rdp.microRNA.data[,colnames(rdp.microRNA.data) %in% rownames(ncmb.phen)]

ncmb.cage <- read.table(file="../raw_data/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.microRNA.count <- rdp.microRNA.count[,colnames(rdp.microRNA.count) %in% rownames(ncmb.cage)]

a <- NULL;
for(i in levels(rdp.microRNA.data$rank)){
	a <- c(a,names(which(apply(rdp.microRNA.count[which(rdp.microRNA.data$rank == i),],2,sum) > 5000)))
}
a <- unique(a)


rdp.microRNA.count <- rdp.microRNA.count[,colnames(rdp.microRNA.count) %in% a]

ncmb.cage.microRNA <- ncmb.cage[match(colnames(rdp.microRNA.count),rownames(ncmb.cage)),]
ncmb.cage.microRNA$cage <- as.factor(as.character(ncmb.cage.microRNA$cage))

ncmb.phen.microRNA <- ncmb.phen[match(colnames(rdp.microRNA.count),rownames(ncmb.phen)),]


library(DESeq2)
library(lme4)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

microRNA.taxa.count <- list()
dds.taxa.microRNA <- list()
geoMeans <- list()
vst.taxa.microRNA <- list()
lmer.vst.taxa.microRNA <- list()

for(i in levels(rdp.microRNA.data$rank)){
	microRNA.taxa.count[[i]] <- rdp.microRNA.count[which(rdp.microRNA.data$rank == i),]
	microRNA.taxa.count[[i]] <- microRNA.taxa.count[[i]][which(apply(microRNA.taxa.count[[i]],1,function(x){length(which(x >=5))} >= 20)),]
	dds.taxa.microRNA[[i]] <- DESeqDataSetFromMatrix(countData=microRNA.taxa.count[[i]], colData = ncmb.phen.microRNA, design =~1)
	geoMeans[[i]] <- apply(counts(dds.taxa.microRNA[[i]]), 1, gm_mean)
	dds.taxa.microRNA[[i]] <- estimateSizeFactors(dds.taxa.microRNA[[i]], geoMeans = geoMeans[[i]])
	dds.taxa.microRNA[[i]] = estimateDispersions(dds.taxa.microRNA[[i]])
	vst.taxa.microRNA[[i]] <- getVarianceStabilizedData(dds.taxa.microRNA[[i]])
	lmer.vst.taxa.microRNA[[i]] <- t(apply(vst.taxa.microRNA[[i]],1,function(x){residuals(lmer(x ~ ncmb.phen.microRNA$Generation + (1|ncmb.cage.microRNA$cage)))}))
}

library(DOQTL)
load("../robj/QTLRelInputHglm.RData")
array.G.mice.microRNA <- array.G.mice[match(colnames(vst.taxa.microRNA$phylum),dimnames(array.G.mice)[[1]]),,]
ncmb.addcovar.microRNA <- ncmb.phen.microRNA[,c("Sex","Diet")]
rm(array.G.mice)
gc()

kin.microRNA <- kinship.probs(array.G.mice.microRNA)

eqtl.taxa.microRNA = list()
for(i in levels(rdp.microRNA.data$rank)){
	eqtl.taxa.microRNA[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.microRNA[[i]]),probs=array.G.mice.microRNA,K=kin.microRNA,addcovar=ncmb.addcovar.microRNA,snps=snp.info),digit=4)
}
save(rdp.microRNA.data,eqtl.taxa.microRNA,kin.microRNA,lmer.vst.taxa.microRNA,array.G.mice.microRNA,ncmb.addcovar.microRNA,snp.info,file="Add_eQTL_rdp_microRNA.RData")

array.G.mice.microRNA.diet <- array(NA, c(dim(array.G.mice.microRNA)[1],dim(array.G.mice.microRNA)[2]*2, dim(array.G.mice.microRNA)[3]))
dimnames(array.G.mice.microRNA.diet)[[1]] <- dimnames(array.G.mice.microRNA)[[1]]
dimnames(array.G.mice.microRNA.diet)[[2]] <- c(dimnames(array.G.mice.microRNA)[[2]],paste(dimnames(array.G.mice.microRNA)[[2]],"D",sep=""))
dimnames(array.G.mice.microRNA.diet)[[3]] <- dimnames(array.G.mice.microRNA)[[3]]
for(i in 1:dim(array.G.mice.microRNA)[3]){
 array.G.mice.microRNA.diet[,,i] <- cbind(array.G.mice.microRNA[,,i],array.G.mice.microRNA[,,i]*ncmb.addcovar.microRNA$Diet)
}

eqtl.taxa.microRNA.full = list()
eqtl.taxa.microRNA = list()
eqtl.taxa.microRNA.diet = list()

for(i in levels(rdp.microRNA.data$rank)){
	eqtl.taxa.microRNA[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.microRNA[[i]]),probs=array.G.mice.microRNA,K=kin.microRNA,addcovar=as.matrix(ncmb.addcovar.microRNA),snps=snp.info),digit=4)
	eqtl.taxa.microRNA.full[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.microRNA[[i]]),probs=array.G.mice.microRNA.diet,K=kin.microRNA,addcovar=as.matrix(ncmb.addcovar.microRNA),snps=snp.info),digit=4)
	eqtl.taxa.microRNA.diet[[i]] <- eqtl.taxa.microRNA.full[[i]]-eqtl.taxa.microRNA[[i]]
}
save(rdp.microRNA.data,eqtl.taxa.microRNA.diet,kin.microRNA,lmer.vst.taxa.microRNA,array.G.mice.microRNA.diet,ncmb.addcovar.microRNA,snp.info,file="Int_eQTL_rdp_microRNA.diet.RData")

array.G.mice.microRNA.sex <- array(NA, c(dim(array.G.mice.microRNA)[1],dim(array.G.mice.microRNA)[2]*2, dim(array.G.mice.microRNA)[3]))
dimnames(array.G.mice.microRNA.sex)[[1]] <- dimnames(array.G.mice.microRNA)[[1]]
dimnames(array.G.mice.microRNA.sex)[[2]] <- c(dimnames(array.G.mice.microRNA)[[2]],paste(dimnames(array.G.mice.microRNA)[[2]],"S",sep=""))
dimnames(array.G.mice.microRNA.sex)[[3]] <- dimnames(array.G.mice.microRNA)[[3]]
for(i in 1:dim(array.G.mice.microRNA)[3]){
 array.G.mice.microRNA.sex[,,i] <- cbind(array.G.mice.microRNA[,,i],array.G.mice.microRNA[,,i]*ncmb.addcovar.microRNA$Sex)
}

eqtl.taxa.microRNA.full = list()
eqtl.taxa.microRNA = list()
eqtl.taxa.microRNA.sex = list()

for(i in levels(rdp.microRNA.data$rank)){
	eqtl.taxa.microRNA[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.microRNA[[i]]),probs=array.G.mice.microRNA,K=kin.microRNA,addcovar=as.matrix(ncmb.addcovar.microRNA),snps=snp.info),digit=4)
	eqtl.taxa.microRNA.full[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.microRNA[[i]]),probs=array.G.mice.microRNA.sex,K=kin.microRNA,addcovar=as.matrix(ncmb.addcovar.microRNA),snps=snp.info),digit=4)
	eqtl.taxa.microRNA.sex[[i]] <- eqtl.taxa.microRNA.full[[i]]-eqtl.taxa.microRNA[[i]]
}
save(rdp.microRNA.data,eqtl.taxa.microRNA.sex,kin.microRNA,lmer.vst.taxa.microRNA,array.G.mice.microRNA.sex,ncmb.addcovar.microRNA,snp.info,file="Int_eQTL_rdp_microRNA.sex.RData")

############# QTL OTU (species) ##################

micro.rna.biome <- import_biom(BIOMfilename="../raw_data/RNAotus.biom")
ncmb.data <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
ncmb.cage <- read.table(file="../raw_data/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
ncmb.cage <- ncmb.cage[rownames(ncmb.cage) %in% rownames(ncmb.data),]
ncmb.data <- ncmb.data[match(rownames(ncmb.cage),rownames(ncmb.data)),]
ncmb.info <- data.frame(Generation=ncmb.data$Generation,Sex=ncmb.cage$Gender,Diet=ncmb.data$Diet,cage=as.factor(as.character(ncmb.cage$cage)))
rownames(ncmb.info) <- rownames(ncmb.data)
ncmb.info <- ncmb.info[match(colnames(micro.rna.biome),rownames(ncmb.info)),]
ncmb.info$cage <- as.factor(as.character(ncmb.info$cage))
micro.rna.prune = filter_taxa(micro.rna.biome, function(x) sum(x >= 5) > (0.2*length(x)), TRUE)
library(DESeq2)
micro.rna.dds <- DESeqDataSetFromMatrix(countData=as.matrix(as.data.frame(micro.rna.prune)), colData = ncmb.info, design =~1)
geoMeans.rna <- apply(counts(micro.rna.dds), 1, gm_mean)
micro.rna.dds <- estimateSizeFactors(micro.rna.dds, geoMeans = geoMeans.rna)
micro.rna.dds <- estimateDispersions(micro.rna.dds,fitType='local')
micro.rna.vst <- getVarianceStabilizedData(micro.rna.dds)
library(lme4)
micro.rna.lmer <- t(apply(micro.rna.vst,1,function(x){residuals(lmer(x ~ ncmb.info$Generation + (1|ncmb.info$cage)))}))

load("../robj/QTLRelInputHglm.RData")
array.G.mice.add <- array.G.mice[match(colnames(micro.rna.lmer),dimnames(array.G.mice)[[1]]),,]
ncmb.addcovar.micro <- ncmb.info[,c("Sex","Diet")]
ncmb.addcovar.micro$Sex = ifelse(ncmb.addcovar.micro$Sex == "m",0,1)

kin.micro <- kinship.probs(array.G.mice.add)

array.G.mice.sex <- array(NA, c(dim(array.G.mice.add)[1],dim(array.G.mice.add)[2]*2, dim(array.G.mice.add)[3]))
dimnames(array.G.mice.sex)[[1]] <- dimnames(array.G.mice.add)[[1]]
dimnames(array.G.mice.sex)[[2]] <- c(dimnames(array.G.mice.add)[[2]],paste(dimnames(array.G.mice.add)[[2]],"S",sep=""))
dimnames(array.G.mice.sex)[[3]] <- dimnames(array.G.mice.add)[[3]]
for(i in 1:dim(array.G.mice.add)[3]){
 array.G.mice.sex[,,i] <- cbind(array.G.mice.add[,,i],array.G.mice.add[,,i]*ncmb.addcovar.micro$Sex)
}

array.G.mice.diet <- array(NA, c(dim(array.G.mice.add)[1],dim(array.G.mice.add)[2]*2, dim(array.G.mice.add)[3]))
dimnames(array.G.mice.diet)[[1]] <- dimnames(array.G.mice.add)[[1]]
dimnames(array.G.mice.diet)[[2]] <- c(dimnames(array.G.mice.add)[[2]],paste(dimnames(array.G.mice.add)[[2]],"D",sep=""))
dimnames(array.G.mice.diet)[[3]] <- dimnames(array.G.mice.add)[[3]]
for(i in 1:dim(array.G.mice.add)[3]){
 array.G.mice.diet[,,i] <- cbind(array.G.mice.add[,,i],array.G.mice.add[,,i]*ncmb.addcovar.micro$Diet)
}

save(micro.rna.lmer,kin.micro,array.G.mice.add,ncmb.addcovar.micro,snp.info,file="Add_eQTL_species_microRNA.RData")
save(micro.rna.lmer,kin.micro,array.G.mice.add,array.G.mice.sex,ncmb.addcovar.micro,snp.info,file="IntSex_eQTL_species_microRNA.RData")
save(micro.rna.lmer,kin.micro,array.G.mice.add,array.G.mice.diet,ncmb.addcovar.micro,snp.info,file="IntDiet_eQTL_species_microRNA.RData")

############## permutation computatitonaly expensive #########

##### Additive model ###########
dir.create("spQTL")
### copy all rdata file spQTL folder

sink('/scratch/analysis-output.txt')
inpath <- "/SpQTL/"
library(DOQTL)
load("SpQTL/Add_eQTL_species_microRNA.RData")

eqtl.micro.RNA.add <- round(scanone.eqtl(expr=t(micro.RNA.lmer),probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)

nperm <- 1000

eqtl.micro.RNA.add.perm <- matrix(0,ncol=nrow(micro.RNA.lmer),nrow=nperm)
colnames(eqtl.micro.RNA.add.perm) <- rownames(micro.RNA.lmer)
rownames(eqtl.micro.RNA.add.perm) <- paste("Perm",c(1:nperm),sep="")

for(j in 1:nperm){
	eqtl.micro.RNA.add.tmp <- NULL
	micro.RNA.lmer.sample <- NULL
	micro.RNA.lmer.sample <- apply(micro.RNA.lmer,1,sample)
	rownames(micro.RNA.lmer.sample) <- colnames(micro.RNA.lmer)
	eqtl.micro.RNA.add.tmp <- round(scanone.eqtl(expr=micro.RNA.lmer.sample,probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
	eqtl.micro.RNA.add.perm[j,] <- apply(eqtl.micro.RNA.add.tmp,2,max)
}
outfile=paste(inpath,"Add_eQTL_species_microRNA_perm.RData",sep="")
save(eqtl.micro.RNA.add,eqtl.micro.RNA.add.perm,file=outfile)

##### Int Sex model ###########

#!/usr/bin/env Rscript
library(DOQTL)
sink('/scratch/IntsexRNAanalysis-output.txt')
inpath <- "/SpQTL/"
library(DOQTL)
load("/SpQTL/IntSex_eQTL_species_microRNA.RData")
eqtl.micro.RNA.add <- round(scanone.eqtl(expr=t(micro.RNA.lmer),probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
eqtl.micro.RNA.full <- round(scanone.eqtl(expr=t(micro.RNA.lmer),probs=array.G.mice.sex,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
eqtl.micro.RNA.sex <- eqtl.micro.RNA.full-eqtl.micro.RNA.add


nperm=1000

eqtl.micro.RNA.sex.perm <- matrix(0,ncol=nrow(micro.RNA.lmer),nrow=nperm)
colnames(eqtl.micro.RNA.sex.perm) <- rownames(micro.RNA.lmer)
rownames(eqtl.micro.RNA.sex.perm) <- paste("Perm",c(1:nperm),sep="")

for(j in 1:nperm){
	eqtl.micro.RNA.add.tmp <- NULL
	eqtl.micro.RNA.full.tmp <- NULL
	eqtl.micro.RNA.sex.tmp <- NULL
	micro.RNA.lmer.sample <- NULL
	micro.RNA.lmer.sample <- apply(micro.RNA.lmer,1,sample)
	rownames(micro.RNA.lmer.sample) <- colnames(micro.RNA.lmer)
	eqtl.micro.RNA.add.tmp <- round(scanone.eqtl(expr=micro.RNA.lmer.sample,probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
	eqtl.micro.RNA.full.tmp <- round(scanone.eqtl(expr=micro.RNA.lmer.sample,probs=array.G.mice.sex,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
	eqtl.micro.RNA.sex.tmp <- eqtl.micro.RNA.full.tmp-eqtl.micro.RNA.add.tmp
	eqtl.micro.RNA.sex.perm[j,] <- apply(eqtl.micro.RNA.sex.tmp,2,max)
}
outfile=paste(inpath,"IntSex_eQTL_species_microRNA_perm.RData",sep="")
save(eqtl.micro.RNA.sex,eqtl.micro.RNA.sex.perm,file=outfile)
sink()

##### Diet Sex model ###########
library(DOQTL)
sink('/scratch/IntdietRNAanalysis-output.txt')
inpath <- "SpQTL"
load("/IntDiet_eQTL_species_microRNA.RData")
eqtl.micro.RNA.add <- round(scanone.eqtl(expr=t(micro.RNA.lmer),probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
eqtl.micro.RNA.full <- round(scanone.eqtl(expr=t(micro.RNA.lmer),probs=array.G.mice.diet,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
eqtl.micro.RNA.diet <- eqtl.micro.RNA.full-eqtl.micro.RNA.add


nperm=1000

eqtl.micro.RNA.diet.perm <- matrix(0,ncol=nrow(micro.RNA.lmer),nrow=nperm)
colnames(eqtl.micro.RNA.diet.perm) <- rownames(micro.RNA.lmer)
rownames(eqtl.micro.RNA.diet.perm) <- paste("Perm",c(1:nperm),sep="")

for(j in 1:nperm){
	eqtl.micro.RNA.add.tmp <- NULL
	eqtl.micro.RNA.full.tmp <- NULL
	eqtl.micro.RNA.diet.tmp <- NULL
	micro.RNA.lmer.sample <- NULL
	micro.RNA.lmer.sample <- apply(micro.RNA.lmer,1,sample)
	rownames(micro.RNA.lmer.sample) <- colnames(micro.RNA.lmer)
	eqtl.micro.RNA.add.tmp <- round(scanone.eqtl(expr=micro.RNA.lmer.sample,probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
	eqtl.micro.RNA.full.tmp <- round(scanone.eqtl(expr=micro.RNA.lmer.sample,probs=array.G.mice.diet,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
	eqtl.micro.RNA.diet.tmp <- eqtl.micro.RNA.full.tmp-eqtl.micro.RNA.add.tmp
	eqtl.micro.RNA.diet.perm[j,] <- apply(eqtl.micro.RNA.diet.tmp,2,max)
}
outfile=paste(inpath,"IntDiet_eQTL_species_microRNA_perm.RData",sep="")
save(eqtl.micro.RNA.diet,eqtl.micro.RNA.diet.perm,file=outfile)
sink()

###### function for calculating percentage of variance ######
pov <- function(datavst = NULL,ncmbInfo=NULL){
	library(lme4)
	pov.taxa <- apply(datavst,1,function(x){
		data.tax <- data.frame(taxa=x,sex=as.factor(ncmbInfo$Sex),Diet=as.factor(ncmbInfo$Diet),Gen=as.factor(ncmbInfo$Generation),cage=as.factor(ncmbInfo$cage))
		lm.tax <- lm(taxa ~ sex + Gen + Diet, data=data.tax)
		af.tax <- anova(lm.tax)
		afss.tax <- af.tax$"Sum Sq"
		PctExp.tax=afss.tax/sum(afss.tax)*100
		names(PctExp.tax) <- rownames(af.tax)
		data.tax$residTax <- as.numeric(residuals(lm.tax))
		lmer.cage.tax <- lmer(residTax ~ 1 + (1|cage),data=data.tax)
		pct.cage <- as.data.frame(VarCorr(lmer.cage.tax))[1,4]/(as.data.frame(VarCorr(lmer.cage.tax))[1,4]+as.data.frame(VarCorr(lmer.cage.tax))[2,4])*100
		names(pct.cage) <- "cage"
		j <-  c(PctExp.tax,pct.cage)[-4]
		return(j)
	})
	pov.df <- data.frame(taxa=rownames(datavst),t(round(pov.taxa,digit=4)))
	return(pov.df)
}



######## for calculating % of variance in species RNA ###########
library("phyloseq")
gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
micro.rna.biome <- import_biom(BIOMfilename="../raw_data/RNAotus.biom")
ncmb.data <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
ncmb.cage <- read.table(file="../raw_data/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
ncmb.cage <- ncmb.cage[rownames(ncmb.cage) %in% rownames(ncmb.data),]
ncmb.data <- ncmb.data[match(rownames(ncmb.cage),rownames(ncmb.data)),]
ncmb.info <- data.frame(Generation=ncmb.data$Generation,Sex=ncmb.cage$Gender,Diet=ncmb.data$Diet,cage=as.factor(as.character(ncmb.cage$cage)))
rownames(ncmb.info) <- rownames(ncmb.data)
ncmb.info <- ncmb.info[match(colnames(micro.rna.biome),rownames(ncmb.info)),]
ncmb.info$cage <- as.factor(as.character(ncmb.info$cage))
micro.rna.prune = filter_taxa(micro.rna.biome, function(x) sum(x >= 5) > (0.2*length(x)), TRUE)
library(DESeq2)
micro.rna.dds <- DESeqDataSetFromMatrix(countData=as.matrix(as.data.frame(micro.rna.prune)), colData = ncmb.info, design =~1)
geoMeans.rna <- apply(counts(micro.rna.dds), 1, gm_mean)
micro.rna.dds <- estimateSizeFactors(micro.rna.dds, geoMeans = geoMeans.rna)
micro.rna.dds <- estimateDispersions(micro.rna.dds,fitType='local')
micro.rna.vst <- getVarianceStabilizedData(micro.rna.dds)
micro.rna.vst.pov <- pov(datavst = micro.rna.vst,ncmbInfo=ncmb.info)
write.table(micro.rna.vst.pov,file="MicroRNA_species_pov.tab",sep="\t",row.names=F,col.names=T,quote=F)



######## for calculating % of variance in phylum to genus rdp RNA ###########
rdp.microRNA <- read.table(file="../raw_data/AIL_micro_RNA.txt",header=T,row.names=1,sep="\t")[,-c(1,2)]
rdp.microRNA.data <- rdp.microRNA[rdp.microRNA$rank %in% c("phylum","class","order","family","genus"),]
rdp.microRNA.data$rank <- as.factor(as.character(rdp.microRNA.data$rank))
colnames(rdp.microRNA.data) <- gsub("X|.nonchimeras.fasta","",colnames(rdp.microRNA.data))

ncmb.phen <- read.table(file="/home/yask/artem_data/MicrMycoGmice/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.microRNA.count <- rdp.microRNA.data[,colnames(rdp.microRNA.data) %in% rownames(ncmb.phen)]

ncmb.cage <- read.table(file="/home/yask/artem_data/MicrMycoGmice/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.microRNA.count <- rdp.microRNA.count[,colnames(rdp.microRNA.count) %in% rownames(ncmb.cage)]

a <- NULL;
for(i in levels(rdp.microRNA.data$rank)){
	a <- c(a,names(which(apply(rdp.microRNA.count[which(rdp.microRNA.data$rank == i),],2,sum) > 5000)))
}
a <- unique(a)


rdp.microRNA.count <- rdp.microRNA.count[,colnames(rdp.microRNA.count) %in% a]

ncmb.cage.microRNA <- ncmb.cage[match(colnames(rdp.microRNA.count),rownames(ncmb.cage)),]
ncmb.cage.microRNA$cage <- as.factor(as.character(ncmb.cage.microRNA$cage))

ncmb.phen.microRNA <- ncmb.phen[match(colnames(rdp.microRNA.count),rownames(ncmb.phen)),]

ncmb.cov.microRNA <- data.frame(Generation=as.factor(ncmb.phen.microRNA$Generation),Sex=as.factor(ncmb.phen.microRNA$Sex),Diet=as.factor(ncmb.phen.microRNA$Diet),cage=as.factor(ncmb.cage.microRNA$cage))

rownames(ncmb.cov.microRNA) <- rownames(ncmb.phen.microRNA)



library(DESeq2)
library(lme4)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

microRNA.taxa.count <- list()
dds.taxa.microRNA <- list()
geoMeans <- list()
vst.taxa.microRNA <- list()
pov.taxa.microRNA <- list()

for(i in levels(rdp.microRNA.data$rank)){
	microRNA.taxa.count[[i]] <- rdp.microRNA.count[which(rdp.microRNA.data$rank == i),]
	microRNA.taxa.count[[i]] <- microRNA.taxa.count[[i]][which(apply(microRNA.taxa.count[[i]],1,function(x){length(which(x >=5))} >= 20)),]
	dds.taxa.microRNA[[i]] <- DESeqDataSetFromMatrix(countData=microRNA.taxa.count[[i]], colData = ncmb.phen.microRNA, design =~1)
	geoMeans[[i]] <- apply(counts(dds.taxa.microRNA[[i]]), 1, gm_mean)
	dds.taxa.microRNA[[i]] <- estimateSizeFactors(dds.taxa.microRNA[[i]], geoMeans = geoMeans[[i]])
	dds.taxa.microRNA[[i]] <- estimateDispersions(dds.taxa.microRNA[[i]])
	vst.taxa.microRNA[[i]] <- getVarianceStabilizedData(dds.taxa.microRNA[[i]])
	pov.taxa.microRNA[[i]] <- pov(datavst = vst.taxa.microRNA[[i]],ncmbInfo=ncmb.cov.microRNA)
	write.table(pov.taxa.microRNA[[i]],file=paste("MicroRNA_",i,"_pov.tab",sep=""),sep="\t",row.names=F,col.names=T,quote=F)
	
}

