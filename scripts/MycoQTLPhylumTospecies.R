library("DESeq2")
library("phyloseq")
library(DOQTL)

###### Mycobiome OTU QTL ##########

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
myco.biome <- import_biom(BIOMfilename="../raw_data/otu_table.biom")
myco.biome.dds <- phyloseq_to_deseq2(myco.biome)
ncmb.data = read.table(file="../raw_data/NCMB_Sept2018.tab",header=T,row.names=1)
myco.biome.5k = prune_samples(sample_names(myco.biome.5k) != "NC", myco.biome.5k)
ncmb.data.myco <- ncmb.data[match(sample_names(myco.biome.5k),rownames(ncmb.data)),]
ncmb.data.myco$Diet <- factor(ifelse(ncmb.data.myco$Diet == 0,"Cal",ifelse(ncmb.data.myco$Diet == 1,"Con","Wes")),levels=c("Cal","Con","Wes"))
ncmb.data.myco$Sex <- factor(ifelse(ncmb.data.myco$Sex == 0,"M","F"),levels=c("M","F"))
sample_data(myco.biome.5k) <- ncmb.data.myco
myco.prune = filter_taxa(myco.biome.5k, function(x) sum(x >= 5) > (0.2*length(x)), TRUE)
diagdds = phyloseq_to_deseq2(myco.prune, ~ Diet)
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
myco.vst <- varianceStabilizingTransformation(diagdds,fitType='local')
myco.vst.data <- assay(myco.vst)

ncmb.info <- read.table(file="../raw_data/NcmbInfo.txt",header=T,sep="\t")
ncmb.info.myco <- ncmb.info[match(colnames(myco.vst),ncmb.info$Auto.ID),]
ncmb.info.myco <- as.factor(as.character(ncmb.info.myco$Käfignummer))

library(lmer)
library(lme4)

myco.vst.res <- myco.vst.data;
for(i in 1:nrow(myco.vst.data)){
myco.vst.res[i,] <- as.numeric(residuals(lmer(as.numeric(myco.vst.data[i,]) ~ myco.vst$Generation + (1|ncmb.info.myco$cage))))
}

load("../robj/QTLRelInputHglm.RData")
idx.myco <- match(colnames(myco.vst.res),rownames(array.G.mice[,,1]))
array.G.mice.myco <- array.G.mice[idx.myco,,]
kin.Gmice <- kinship.probs(array.G.mice.myco)
myco.eqtl <- scanone.eqtl(t(myco.vst.res),probs=array.G.mice.myco,K=kin.Gmice,addcovar=as.matrix(ncmb.addcovar))
save(myco.eqtl,file="../robj/speciesMycoeQTL.RData")

################# qtl phylum to genus rdp #############

rdp.myco <- read.table(file="../raw_data/AIL_myco_hier.tab",header=T,row.names=1,sep="\t")[,-c(1,2)]
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


library(DESeq2)
library(lme4)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

myco.taxa.count <- list()
dds.taxa.myco <- list()
geoMeans <- list()
vst.taxa.myco <- list()
lmer.vst.taxa.myco <- list()

for(i in levels(rdp.myco.data$rank)){
	myco.taxa.count[[i]] <- rdp.myco.count[which(rdp.myco.data$rank == i),]
	myco.taxa.count[[i]] <- myco.taxa.count[[i]][which(apply(myco.taxa.count[[i]],1,function(x){length(which(x >=5))} >= 20)),]
	dds.taxa.myco[[i]] <- DESeqDataSetFromMatrix(countData=myco.taxa.count[[i]], colData = ncmb.phen.myco, design =~1)
	geoMeans[[i]] <- apply(counts(dds.taxa.myco[[i]]), 1, gm_mean)
	dds.taxa.myco[[i]] <- estimateSizeFactors(dds.taxa.myco[[i]], geoMeans = geoMeans[[i]])
	dds.taxa.myco[[i]] = estimateDispersions(dds.taxa.myco[[i]])
	vst.taxa.myco[[i]] <- getVarianceStabilizedData(dds.taxa.myco[[i]])
	lmer.vst.taxa.myco[[i]] <- t(apply(vst.taxa.myco[[i]],1,function(x){residuals(lmer(x ~ ncmb.phen.myco$Generation + (1|ncmb.cage.myco$cage)))}))
}

library(DOQTL)
load("../robj/QTLRelInputHglm.RData")
array.G.mice.myco <- array.G.mice[match(colnames(vst.taxa.myco$phylum),dimnames(array.G.mice)[[1]]),,]
ncmb.addcovar.myco <- ncmb.phen.myco[,c("Sex","Diet")]
rm(array.G.mice)
gc()

kin.myco <- kinship.probs(array.G.mice.myco)

eqtl.taxa.myco = list()
for(i in levels(rdp.myco.data$rank)){
	eqtl.taxa.myco[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.myco[[i]]),probs=array.G.mice.myco,K=kin.myco,addcovar=ncmb.addcovar.myco,snps=snp.info),digit=4)
}
save(rdp.myco.data,eqtl.taxa.myco,kin.myco,lmer.vst.taxa.myco,array.G.mice.myco,ncmb.addcovar.myco,snp.info,file="../robj/Add_eQTL_rdp_myco.RData")

q()

library(DOQTL)
load("/home/yask/artem_data/MicrMycoGmice/Add_eQTL_rdp_myco.RData")
eqtl.taxa.myco.perm <- list()
for(i in levels(rdp.myco.data$rank)){
	eqtl.taxa.myco.perm[[i]] <- matrix(0,ncol=nrow(lmer.vst.taxa.myco[[i]]),nrow=1000)
	colnames(eqtl.taxa.myco.perm[[i]]) <- rownames(lmer.vst.taxa.myco[[i]])
	rownames(eqtl.taxa.myco.perm[[i]]) <- paste("Perm",1:1000,sep="")
	for(j in 1:1000){
 		lmer.tmp.taxa.myco <- apply(lmer.vst.taxa.myco[[i]],1,sample)
 		rownames(lmer.tmp.taxa.myco) <- colnames(lmer.vst.taxa.myco[[i]])
		eqtl.taxa.myco.perm[[i]][j,] <- apply(round(scanone.eqtl(expr=lmer.tmp.taxa.myco,probs=array.G.mice.myco,K=kin.myco,addcovar=ncmb.addcovar.myco,snps=snp.info),digit=4),2,max)
		print(paste(j,"Permuations Done for",i,sep=" "))
	}
}
save(lmer.vst.taxa.myco,eqtl.taxa.myco.perm,eqtl.taxa.myco,file="../robj/Add_eQTL_rdp_myco_perm.RData")

library(DOQTL)

load("../robj/Add_eQTL_rdp_myco.RData")
array.G.mice.myco.diet <- array(NA, c(dim(array.G.mice.myco)[1],dim(array.G.mice.myco)[2]*2, dim(array.G.mice.myco)[3]))
dimnames(array.G.mice.myco.diet)[[1]] <- dimnames(array.G.mice.myco)[[1]]
dimnames(array.G.mice.myco.diet)[[2]] <- c(dimnames(array.G.mice.myco)[[2]],paste(dimnames(array.G.mice.myco)[[2]],"D",sep=""))
dimnames(array.G.mice.myco.diet)[[3]] <- dimnames(array.G.mice.myco)[[3]]
for(i in 1:dim(array.G.mice.myco)[3]){
 array.G.mice.myco.diet[,,i] <- cbind(array.G.mice.myco[,,i],array.G.mice.myco[,,i]*ncmb.addcovar.myco$Diet)
}

eqtl.taxa.myco.diet = list()
for(i in levels(rdp.myco.data$rank)){
	eqtl.taxa.myco.diet[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.myco[[i]]),probs=array.G.mice.myco.diet,K=kin.myco,addcovar=as.matrix(ncmb.addcovar.myco),snps=snp.info),digit=4)
}
save(rdp.myco.data,eqtl.taxa.myco.diet,kin.myco,lmer.vst.taxa.myco,array.G.mice.myco.diet,ncmb.addcovar.myco,snp.info,file="../robj/Int_eQTL_rdp_myco.diet.RData")

library(DOQTL)
load("../robj/Add_eQTL_rdp_myco.RData")
load("../robj/Int_eQTL_rdp_myco.diet.RData")

eqtl.taxa.myco.full = list()
eqtl.taxa.myco = list()
eqtl.taxa.myco.diet = list()

for(i in levels(rdp.myco.data$rank)){
	eqtl.taxa.myco[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.myco[[i]]),probs=array.G.mice.myco,K=kin.myco,addcovar=as.matrix(ncmb.addcovar.myco),snps=snp.info),digit=4)
	eqtl.taxa.myco.full[[i]] <- round(scanone.eqtl(expr=t(lmer.vst.taxa.myco[[i]]),probs=array.G.mice.myco.diet,K=kin.myco,addcovar=as.matrix(ncmb.addcovar.myco),snps=snp.info),digit=4)
	eqtl.taxa.myco.diet[[i]] <- eqtl.taxa.myco.full[[i]]-eqtl.taxa.myco[[i]]
}

eqtl.taxa.myco.diet.perm <- list()

for(i in levels(rdp.myco.data$rank)){
	eqtl.taxa.myco.diet.perm[[i]] <- matrix(0,ncol=nrow(lmer.vst.taxa.myco[[i]]),nrow=1000)
	colnames(eqtl.taxa.myco.diet.perm[[i]]) <- rownames(lmer.vst.taxa.myco[[i]])
	rownames(eqtl.taxa.myco.diet.perm[[i]]) <- paste("Perm",1:1000,sep="")
	for(j in 1:1000){
		eqtl.taxa.myco.perm  <- NULL
		eqtl.taxa.myco.full.perm <- NULL
		lmer.tmp.taxa.myco <- NULL
 		lmer.tmp.taxa.myco <- apply(lmer.vst.taxa.myco[[i]],1,sample)
		rownames(lmer.tmp.taxa.myco) <- colnames(lmer.vst.taxa.myco[[i]])
		eqtl.taxa.myco.perm <- round(scanone.eqtl(expr=lmer.tmp.taxa.myco,probs=array.G.mice.myco,K=kin.myco,addcovar=as.matrix(ncmb.addcovar.myco),snps=snp.info),digit=4)
		eqtl.taxa.myco.full.perm <- round(scanone.eqtl(expr=lmer.tmp.taxa.myco,probs=array.G.mice.myco.diet,K=kin.myco,addcovar=as.matrix(ncmb.addcovar.myco),snps=snp.info),digit=4)
		eqtl.taxa.myco.diet.perm[[i]][j,] <- apply((eqtl.taxa.myco.full.perm-eqtl.taxa.myco.perm),2,max)
		print(paste(j,"Permuations Done for",i,sep=" "))
	}
}

save(eqtl.taxa.myco.diet,eqtl.taxa.myco.diet.perm,file="../robj/Int_diet_eQTL_rdp_myco_perm.RData")
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



######## for calculating % of variance in species RNA ###########

lnames <- load("../robj/speciesMycoeQTL.RData")

ncmb.data <- read.table(file="/home/yask/artem_data/MicrMycoGmice/NCMB_Sept2018.tab",sep="\t",header=T,row.names=1)
ncmb.cage <- read.table(file="/home/yask/artem_data/MicrMycoGmice/ncmb_cage_reformat.tsv",sep="\t",header=T,row.names=1)
ncmb.cage <- ncmb.cage[rownames(ncmb.cage) %in% rownames(ncmb.data),]
ncmb.data <- ncmb.data[match(rownames(ncmb.cage),rownames(ncmb.data)),]
ncmb.info <- data.frame(Generation=ncmb.data$Generation,Sex=ncmb.cage$Gender,Diet=ncmb.data$Diet,cage=as.factor(as.character(ncmb.cage$cage)))
rownames(ncmb.info) <- rownames(ncmb.data)
ncmb.info <- ncmb.info[match(colnames(myco.vst.data),rownames(ncmb.info)),]
ncmb.info$cage <- as.factor(as.character(ncmb.info$cage))

myco.dna.vst.pov <- pov(datavst = myco.vst.data,ncmbInfo=ncmb.info)
write.table(myco.dna.vst.pov,file="MycoDNA_species_pov.tab",sep="\t",row.names=F,col.names=T,quote=F)



######## for calculating % of variance in phylum to genus rdp RNA ###########

rdp.myco <- read.table(file="../raw_data/AIL_myco_hier.tab",header=T,row.names=1,sep="\t")[,-c(1,2)]
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

ncmb.cov.myco <- data.frame(Generation=as.factor(ncmb.phen.myco$Generation),Sex=as.factor(ncmb.phen.myco$Sex),Diet=as.factor(ncmb.phen.myco$Diet),cage=as.factor(ncmb.cage.myco$cage))

rownames(ncmb.cov.myco) <- rownames(ncmb.phen.myco)

myco.taxa.count <- list()
dds.taxa.myco <- list()
geoMeans <- list()
vst.taxa.myco <- list()
pov.taxa.myco <- list()

for(i in levels(rdp.myco.data$rank)){
	myco.taxa.count[[i]] <- rdp.myco.count[which(rdp.myco.data$rank == i),]
	myco.taxa.count[[i]] <- myco.taxa.count[[i]][which(apply(myco.taxa.count[[i]],1,function(x){length(which(x >=5))} >= 20)),]
	dds.taxa.myco[[i]] <- DESeqDataSetFromMatrix(countData=myco.taxa.count[[i]], colData = ncmb.phen.myco, design =~1)
	geoMeans[[i]] <- apply(counts(dds.taxa.myco[[i]]), 1, gm_mean)
	dds.taxa.myco[[i]] <- estimateSizeFactors(dds.taxa.myco[[i]], geoMeans = geoMeans[[i]])
	dds.taxa.myco[[i]] <- estimateDispersions(dds.taxa.myco[[i]])
	vst.taxa.myco[[i]] <- getVarianceStabilizedData(dds.taxa.myco[[i]])
	pov.taxa.myco[[i]] <- pov(datavst = vst.taxa.myco[[i]],ncmbInfo=ncmb.cov.myco)
	write.table(pov.taxa.myco[[i]],file=paste("MycoDNA_",i,"_pov.tab",sep=""),sep="\t",row.names=F,col.names=T,quote=F)
	
}
