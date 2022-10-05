library(DESeq2)
library(lme4)
library(DOQTL)
library(phyloseq)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
############ phylum to genus QTL based on rdp ############

rdp.microDNA <- read.table(file="../raw_data/AIL_micro_DNA.txt",header=T,row.names=1,sep="\t")[,-c(1,2)]
rdp.microDNA.data <- rdp.microDNA[rdp.microDNA$rank %in% c("phylum","class","order","family","genus"),]
rdp.microDNA.data$rank <- as.factor(as.character(rdp.microDNA.data$rank))
colnames(rdp.microDNA.data) <- gsub("X|.nonchimeras.fasta","",colnames(rdp.microDNA.data))

ncmb.phen <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.microDNA.count <- rdp.microDNA.data[,colnames(rdp.microDNA.data) %in% rownames(ncmb.phen)]

ncmb.cage <- read.table(file="../raw_data/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.microDNA.count <- rdp.microDNA.count[,colnames(rdp.microDNA.count) %in% rownames(ncmb.cage)]

a <- NULL;
for(i in levels(rdp.microDNA.data$rank)){
	a <- c(a,names(which(apply(rdp.microDNA.count[which(rdp.microDNA.data$rank == i),],2,sum) > 5000)))
}
a <- unique(a)


rdp.microDNA.count <- rdp.microDNA.count[,colnames(rdp.microDNA.count) %in% a]

ncmb.cage.microDNA <- ncmb.cage[match(colnames(rdp.microDNA.count),rownames(ncmb.cage)),]
ncmb.cage.microDNA$cage <- as.factor(as.character(ncmb.cage.microDNA$cage))

ncmb.phen.microDNA <- ncmb.phen[match(colnames(rdp.microDNA.count),rownames(ncmb.phen)),]


microDNA.taxa.count <- list()
dds.taxa.microDNA <- list()
geoMeans <- list()
vst.taxa.microDNA <- list()
lmer.vst.taxa.microDNA <- list()

for(i in levels(rdp.microDNA.data$rank)){
	microDNA.taxa.count[[i]] <- rdp.microDNA.count[which(rdp.microDNA.data$rank == i),]
	microDNA.taxa.count[[i]] <- microDNA.taxa.count[[i]][which(apply(microDNA.taxa.count[[i]],1,function(x){length(which(x >=5))} >= 20)),]
	dds.taxa.microDNA[[i]] <- DESeqDataSetFromMatrix(countData=microDNA.taxa.count[[i]], colData = ncmb.phen.microDNA, design =~1)
	geoMeans[[i]] <- apply(counts(dds.taxa.microDNA[[i]]), 1, gm_mean)
	dds.taxa.microDNA[[i]] <- estimateSizeFactors(dds.taxa.microDNA[[i]], geoMeans = geoMeans[[i]])
	dds.taxa.microDNA[[i]] = estimateDispersions(dds.taxa.microDNA[[i]])
	vst.taxa.microDNA[[i]] <- getVarianceStabilizedData(dds.taxa.microDNA[[i]])
	lmer.vst.taxa.microDNA[[i]] <- t(apply(vst.taxa.microDNA[[i]],1,function(x){residuals(lmer(x ~ ncmb.phen.microDNA$Generation + (1|ncmb.cage.microDNA$cage)))}))
}

load("../robj/QTLRelInputHglm.RData")
array.G.mice.microDNA <- array.G.mice[match(colnames(vst.taxa.microDNA$phylum),dimnames(array.G.mice)[[1]]),,]
ncmb.addcovar.microDNA <- ncmb.phen.microDNA[,c("Sex","Diet")]
rm(array.G.mice)
gc()

kin.microDNA <- kinship.probs(array.G.mice.microDNA)

eqtl.taxa.microDNA = list()
for(i in levels(rdp.microDNA.data$rank)){
	eqtl.taxa.microDNA[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.microDNA[[i]]),probs=array.G.mice.microDNA,K=kin.microDNA,addcovar=ncmb.addcovar.microDNA,snps=snp.info),digit=4)
}
save(rdp.microDNA.data,eqtl.taxa.microDNA,kin.microDNA,lmer.vst.taxa.microDNA,array.G.mice.microDNA,ncmb.addcovar.microDNA,snp.info,file="Add_eQTL_rdp_microDNA.RData")

array.G.mice.microDNA.diet <- array(NA, c(dim(array.G.mice.microDNA)[1],dim(array.G.mice.microDNA)[2]*2, dim(array.G.mice.microDNA)[3]))
dimnames(array.G.mice.microDNA.diet)[[1]] <- dimnames(array.G.mice.microDNA)[[1]]
dimnames(array.G.mice.microDNA.diet)[[2]] <- c(dimnames(array.G.mice.microDNA)[[2]],paste(dimnames(array.G.mice.microDNA)[[2]],"D",sep=""))
dimnames(array.G.mice.microDNA.diet)[[3]] <- dimnames(array.G.mice.microDNA)[[3]]
for(i in 1:dim(array.G.mice.microDNA)[3]){
 array.G.mice.microDNA.diet[,,i] <- cbind(array.G.mice.microDNA[,,i],array.G.mice.microDNA[,,i]*ncmb.addcovar.microDNA$Diet)
}

eqtl.taxa.microDNA.full = list()
eqtl.taxa.microDNA = list()
eqtl.taxa.microDNA.diet = list()

for(i in levels(rdp.microDNA.data$rank)){
	eqtl.taxa.microDNA[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.microDNA[[i]]),probs=array.G.mice.microDNA,K=kin.microDNA,addcovar=as.matrix(ncmb.addcovar.microDNA),snps=snp.info),digit=4)
	eqtl.taxa.microDNA.full[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.microDNA[[i]]),probs=array.G.mice.microDNA.diet,K=kin.microDNA,addcovar=as.matrix(ncmb.addcovar.microDNA),snps=snp.info),digit=4)
	eqtl.taxa.microDNA.diet[[i]] <- eqtl.taxa.microDNA.full[[i]]-eqtl.taxa.microDNA[[i]]
}
save(rdp.microDNA.data,eqtl.taxa.microDNA.diet,kin.microDNA,lmer.vst.taxa.microDNA,array.G.mice.microDNA.diet,ncmb.addcovar.microDNA,snp.info,file="Int_eQTL_rdp_microDNA.diet.RData")

array.G.mice.microDNA.sex <- array(NA, c(dim(array.G.mice.microDNA)[1],dim(array.G.mice.microDNA)[2]*2, dim(array.G.mice.microDNA)[3]))
dimnames(array.G.mice.microDNA.sex)[[1]] <- dimnames(array.G.mice.microDNA)[[1]]
dimnames(array.G.mice.microDNA.sex)[[2]] <- c(dimnames(array.G.mice.microDNA)[[2]],paste(dimnames(array.G.mice.microDNA)[[2]],"S",sep=""))
dimnames(array.G.mice.microDNA.sex)[[3]] <- dimnames(array.G.mice.microDNA)[[3]]
for(i in 1:dim(array.G.mice.microDNA)[3]){
 array.G.mice.microDNA.sex[,,i] <- cbind(array.G.mice.microDNA[,,i],array.G.mice.microDNA[,,i]*ncmb.addcovar.microDNA$Sex)
}

eqtl.taxa.microDNA.full = list()
eqtl.taxa.microDNA = list()
eqtl.taxa.microDNA.sex = list()

for(i in levels(rdp.microDNA.data$rank)){
	eqtl.taxa.microDNA[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.microDNA[[i]]),probs=array.G.mice.microDNA,K=kin.microDNA,addcovar=as.matrix(ncmb.addcovar.microDNA),snps=snp.info),digit=4)
	eqtl.taxa.microDNA.full[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.microDNA[[i]]),probs=array.G.mice.microDNA.sex,K=kin.microDNA,addcovar=as.matrix(ncmb.addcovar.microDNA),snps=snp.info),digit=4)
	eqtl.taxa.microDNA.sex[[i]] <- eqtl.taxa.microDNA.full[[i]]-eqtl.taxa.microDNA[[i]]
}
save(rdp.microDNA.data,eqtl.taxa.microDNA.sex,kin.microDNA,lmer.vst.taxa.microDNA,array.G.mice.microDNA.sex,ncmb.addcovar.microDNA,snp.info,file="Int_eQTL_rdp_microDNA.sex.RData")

######################## QTL for species ###########

micro.dna.biome <- import_biom(BIOMfilename="../raw_data/OtusDNA.biom")
ncmb.data <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
ncmb.cage <- read.table(file="../raw_data/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
ncmb.cage <- ncmb.cage[rownames(ncmb.cage) %in% rownames(ncmb.data),]
ncmb.data <- ncmb.data[match(rownames(ncmb.cage),rownames(ncmb.data)),]
ncmb.info <- data.frame(Generation=ncmb.data$Generation,Sex=ncmb.cage$Gender,Diet=ncmb.data$Diet,cage=as.factor(as.character(ncmb.cage$cage)))
rownames(ncmb.info) <- rownames(ncmb.data)
ncmb.info <- ncmb.info[match(colnames(micro.dna.biome),rownames(ncmb.info)),]
ncmb.info$cage <- as.factor(as.character(ncmb.info$cage))
micro.dna.prune = filter_taxa(micro.dna.biome, function(x) sum(x >= 5) > (0.2*length(x)), TRUE)
library(DESeq2)
micro.dna.dds <- DESeqDataSetFromMatrix(countData=as.matrix(as.data.frame(micro.dna.prune)), colData = ncmb.info, design =~1)
geoMeans.dna <- apply(counts(micro.dna.dds), 1, gm_mean)
micro.dna.dds <- estimateSizeFactors(micro.dna.dds, geoMeans = geoMeans.dna)
micro.dna.dds <- estimateDispersions(micro.dna.dds,fitType='local')
micro.dna.vst <- getVarianceStabilizedData(micro.dna.dds)
library(lme4)
micro.dna.lmer <- t(apply(micro.dna.vst,1,function(x){residuals(lmer(x ~ ncmb.info$Generation + (1|ncmb.info$cage)))}))


load("../robj/QTLRelInputHglm.RData")
array.G.mice.add <- array.G.mice[match(colnames(micro.dna.lmer),dimnames(array.G.mice)[[1]]),,]
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

save(micro.dna.lmer,kin.micro,array.G.mice.add,ncmb.addcovar.micro,snp.info,file="Add_eQTL_species_microDNA.RData")
save(micro.dna.lmer,kin.micro,array.G.mice.sex,ncmb.addcovar.micro,snp.info,file="IntSex_eQTL_species_microDNA.RData")
save(micro.dna.lmer,kin.micro,array.G.mice.diet,ncmb.addcovar.micro,snp.info,file="IntDiet_eQTL_species_microDNA.RData")


############## permutation computatitonaly expensive #########

##### Additive model ###########
dir.create("spQTL")
### copy all rdata file spQTL folder

sink('/scratch/analysis-output.txt')
inpath <- "/SpQTL/"
library(DOQTL)
load("SpQTL/Add_eQTL_species_microDNA.RData")

eqtl.micro.dna.add <- round(scanone.eqtl(expr=t(micro.dna.lmer),probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)

nperm <- 1000

eqtl.micro.dna.add.perm <- matrix(0,ncol=nrow(micro.dna.lmer),nrow=nperm)
colnames(eqtl.micro.dna.add.perm) <- rownames(micro.dna.lmer)
rownames(eqtl.micro.dna.add.perm) <- paste("Perm",c(1:nperm),sep="")

for(j in 1:nperm){
	eqtl.micro.dna.add.tmp <- NULL
	micro.dna.lmer.sample <- NULL
	micro.dna.lmer.sample <- apply(micro.dna.lmer,1,sample)
	rownames(micro.dna.lmer.sample) <- colnames(micro.dna.lmer)
	eqtl.micro.dna.add.tmp <- round(scanone.eqtl(expr=micro.dna.lmer.sample,probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
	eqtl.micro.dna.add.perm[j,] <- apply(eqtl.micro.dna.add.tmp,2,max)
}
outfile=paste(inpath,"Add_eQTL_species_microDNA_perm.RData",sep="")
save(eqtl.micro.dna.add,eqtl.micro.dna.add.perm,file=outfile)

##### Int Sex model ###########

#!/usr/bin/env Rscript
library(DOQTL)
sink('/scratch/IntsexDNAanalysis-output.txt')
inpath <- "/SpQTL/"
library(DOQTL)
load("/SpQTL/IntSex_eQTL_species_microDNA.RData")
eqtl.micro.dna.add <- round(scanone.eqtl(expr=t(micro.dna.lmer),probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
eqtl.micro.dna.full <- round(scanone.eqtl(expr=t(micro.dna.lmer),probs=array.G.mice.sex,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
eqtl.micro.dna.sex <- eqtl.micro.dna.full-eqtl.micro.dna.add


nperm=1000

eqtl.micro.dna.sex.perm <- matrix(0,ncol=nrow(micro.dna.lmer),nrow=nperm)
colnames(eqtl.micro.dna.sex.perm) <- rownames(micro.dna.lmer)
rownames(eqtl.micro.dna.sex.perm) <- paste("Perm",c(1:nperm),sep="")

for(j in 1:nperm){
	eqtl.micro.dna.add.tmp <- NULL
	eqtl.micro.dna.full.tmp <- NULL
	eqtl.micro.dna.sex.tmp <- NULL
	micro.dna.lmer.sample <- NULL
	micro.dna.lmer.sample <- apply(micro.dna.lmer,1,sample)
	rownames(micro.dna.lmer.sample) <- colnames(micro.dna.lmer)
	eqtl.micro.dna.add.tmp <- round(scanone.eqtl(expr=micro.dna.lmer.sample,probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
	eqtl.micro.dna.full.tmp <- round(scanone.eqtl(expr=micro.dna.lmer.sample,probs=array.G.mice.sex,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
	eqtl.micro.dna.sex.tmp <- eqtl.micro.dna.full.tmp-eqtl.micro.dna.add.tmp
	eqtl.micro.dna.sex.perm[j,] <- apply(eqtl.micro.dna.sex.tmp,2,max)
}
outfile=paste(inpath,"IntSex_eQTL_species_microDNA_perm.RData",sep="")
save(eqtl.micro.dna.sex,eqtl.micro.dna.sex.perm,file=outfile)
sink()

##### Diet Sex model ###########
library(DOQTL)
sink('/scratch/IntdietDNAanalysis-output.txt')
inpath <- "SpQTL"
load("/IntDiet_eQTL_species_microDNA.RData")
eqtl.micro.dna.add <- round(scanone.eqtl(expr=t(micro.dna.lmer),probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
eqtl.micro.dna.full <- round(scanone.eqtl(expr=t(micro.dna.lmer),probs=array.G.mice.diet,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
eqtl.micro.dna.diet <- eqtl.micro.dna.full-eqtl.micro.dna.add


nperm=1000

eqtl.micro.dna.diet.perm <- matrix(0,ncol=nrow(micro.dna.lmer),nrow=nperm)
colnames(eqtl.micro.dna.diet.perm) <- rownames(micro.dna.lmer)
rownames(eqtl.micro.dna.diet.perm) <- paste("Perm",c(1:nperm),sep="")

for(j in 1:nperm){
	eqtl.micro.dna.add.tmp <- NULL
	eqtl.micro.dna.full.tmp <- NULL
	eqtl.micro.dna.diet.tmp <- NULL
	micro.dna.lmer.sample <- NULL
	micro.dna.lmer.sample <- apply(micro.dna.lmer,1,sample)
	rownames(micro.dna.lmer.sample) <- colnames(micro.dna.lmer)
	eqtl.micro.dna.add.tmp <- round(scanone.eqtl(expr=micro.dna.lmer.sample,probs=array.G.mice.add,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
	eqtl.micro.dna.full.tmp <- round(scanone.eqtl(expr=micro.dna.lmer.sample,probs=array.G.mice.diet,K=kin.micro,addcovar=as.matrix(ncmb.addcovar.micro),snps=snp.info),digit=4)
	eqtl.micro.dna.diet.tmp <- eqtl.micro.dna.full.tmp-eqtl.micro.dna.add.tmp
	eqtl.micro.dna.diet.perm[j,] <- apply(eqtl.micro.dna.diet.tmp,2,max)
}
outfile=paste(inpath,"IntDiet_eQTL_species_microDNA_perm.RData",sep="")
save(eqtl.micro.dna.diet,eqtl.micro.dna.diet.perm,file=outfile)
sink()

###### function for calculating percentage of variance ######
pov <- function(datavst = NULL,ncmbInfo=NULL){
	library(lme4)
	pov.taxa <- apply(datavst,1,function(x){
		data.tax <- data.frame(taxa=x,sex=as.factor(ncmbInfo$Sex),Diet=as.factor(ncmbInfo$Diet),Gen=as.factor(ncmbInfo$Generation),cage=as.factor(ncmbInfo$cage))
		lmer.cage.tax <- lmer(taxa ~ 1 + (1|cage),data=data.tax)
		var.cage <- as.data.frame(VarCorr(lmer.cage.tax))[1,4]
		var.resid <- as.data.frame(VarCorr(lmer.cage.tax))[2,4]
		pct.cage <- (var.cage/(var.cage+var.resid))*100
		data.tax$residTax <- as.numeric(residuals(lmer.cage.tax))
		names(pct.cage) <- "cage"
		lm.tax <- lm(residTax ~ sex + Gen + Diet, data=data.tax)
		af.tax <- anova(lm.tax)
		afss.tax <- af.tax$"Sum Sq"
		PctExp.tax=afss.tax/sum(afss.tax)*100
		names(PctExp.tax) <- rownames(af.tax)
		j <-  c(PctExp.tax,pct.cage)[-4]
		return(j)
	})
	pov.df <- data.frame(taxa=rownames(datavst),t(round(pov.taxa,digit=4)))
	return(pov.df)
}


######## for calculating % of variance in species DNA ###########
library("phyloseq")
gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
micro.dna.biome <- import_biom(BIOMfilename="/home/yask/artem_data/MicrMycoGmice/MicrobiomeAIL/OtusDNA.biom")
ncmb.data <- read.table(file="/home/yask/artem_data/MicrMycoGmice/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
ncmb.cage <- read.table(file="/home/yask/artem_data/MicrMycoGmice/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
ncmb.cage <- ncmb.cage[rownames(ncmb.cage) %in% rownames(ncmb.data),]
ncmb.data <- ncmb.data[match(rownames(ncmb.cage),rownames(ncmb.data)),]
ncmb.info <- data.frame(Generation=ncmb.data$Generation,Sex=ncmb.cage$Gender,Diet=ncmb.data$Diet,cage=as.factor(as.character(ncmb.cage$cage)))
rownames(ncmb.info) <- rownames(ncmb.data)
ncmb.info <- ncmb.info[match(colnames(micro.dna.biome),rownames(ncmb.info)),]
ncmb.info$cage <- as.factor(as.character(ncmb.info$cage))
micro.dna.prune = filter_taxa(micro.dna.biome, function(x) sum(x >= 5) > (0.2*length(x)), TRUE)
library(DESeq2)
micro.dna.dds <- DESeqDataSetFromMatrix(countData=as.matrix(as.data.frame(micro.dna.prune)), colData = ncmb.info, design =~1)
geoMeans.dna <- apply(counts(micro.dna.dds), 1, gm_mean)
micro.dna.dds <- estimateSizeFactors(micro.dna.dds, geoMeans = geoMeans.dna)
micro.dna.dds <- estimateDispersions(micro.dna.dds,fitType='local')
micro.dna.vst <- getVarianceStabilizedData(micro.dna.dds)

micro.dna.vst.pov <- pov(datavst = micro.dna.vst,ncmbInfo=ncmb.info)
write.table(micro.dna.vst.pov,file="MicroDNA_species_pov.tab",sep="\t",row.names=F,col.names=T,quote=F)



######## for calculating % of variance in phylum to genus rdp DNA ###########

rdp.microDNA <- read.table(file="../raw_data/AIL_micro_DNA.txt",header=T,row.names=1,sep="\t")[,-c(1,2)]
rdp.microDNA.data <- rdp.microDNA[rdp.microDNA$rank %in% c("phylum","class","order","family","genus"),]
rdp.microDNA.data$rank <- as.factor(as.character(rdp.microDNA.data$rank))
colnames(rdp.microDNA.data) <- gsub("X|.nonchimeras.fasta","",colnames(rdp.microDNA.data))

ncmb.phen <- read.table(file="../raw_data/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
rdp.microDNA.count <- rdp.microDNA.data[,colnames(rdp.microDNA.data) %in% rownames(ncmb.phen)]

ncmb.cage <- read.table(file="../raw_data/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
rdp.microDNA.count <- rdp.microDNA.count[,colnames(rdp.microDNA.count) %in% rownames(ncmb.cage)]

a <- NULL;
for(i in levels(rdp.microDNA.data$rank)){
	a <- c(a,names(which(apply(rdp.microDNA.count[which(rdp.microDNA.data$rank == i),],2,sum) > 5000)))
}
a <- unique(a)


rdp.microDNA.count <- rdp.microDNA.count[,colnames(rdp.microDNA.count) %in% a]

ncmb.cage.microDNA <- ncmb.cage[match(colnames(rdp.microDNA.count),rownames(ncmb.cage)),]
ncmb.cage.microDNA$cage <- as.factor(as.character(ncmb.cage.microDNA$cage))

ncmb.phen.microDNA <- ncmb.phen[match(colnames(rdp.microDNA.count),rownames(ncmb.phen)),]

ncmb.cov.microDNA <- data.frame(Generation=as.factor(ncmb.phen.microDNA$Generation),Sex=as.factor(ncmb.phen.microDNA$Sex),Diet=as.factor(ncmb.phen.microDNA$Diet),cage=as.factor(ncmb.cage.microDNA$cage))

rownames(ncmb.cov.microDNA) <- rownames(ncmb.phen.microDNA)

library(DESeq2)
library(lme4)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

microDNA.taxa.count <- list()
dds.taxa.microDNA <- list()
geoMeans <- list()
vst.taxa.microDNA <- list()
pov.taxa.microDNA <- list()

for(i in levels(rdp.microDNA.data$rank)){
	microDNA.taxa.count[[i]] <- rdp.microDNA.count[which(rdp.microDNA.data$rank == i),]
	microDNA.taxa.count[[i]] <- microDNA.taxa.count[[i]][which(apply(microDNA.taxa.count[[i]],1,function(x){length(which(x >=5))} >= 20)),]
	dds.taxa.microDNA[[i]] <- DESeqDataSetFromMatrix(countData=microDNA.taxa.count[[i]], colData = ncmb.phen.microDNA, design =~1)
	geoMeans[[i]] <- apply(counts(dds.taxa.microDNA[[i]]), 1, gm_mean)
	dds.taxa.microDNA[[i]] <- estimateSizeFactors(dds.taxa.microDNA[[i]], geoMeans = geoMeans[[i]])
	dds.taxa.microDNA[[i]] <- estimateDispersions(dds.taxa.microDNA[[i]])
	vst.taxa.microDNA[[i]] <- getVarianceStabilizedData(dds.taxa.microDNA[[i]])
	pov.taxa.microDNA[[i]] <- pov(datavst = vst.taxa.microDNA[[i]],ncmbInfo=ncmb.cov.microDNA)
	write.table(pov.taxa.microDNA[[i]],file=paste("MicroDNA_",i,"_pov.tab",sep=""),sep="\t",row.names=F,col.names=T,quote=F)
	
}




