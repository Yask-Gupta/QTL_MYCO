# QTL_MYCO 
The repository contains all the source code and script required for reproducing the compositional data and QTL data for Mycobiome (Fungi) and Microbiome (Bacteria, DNA and RNA) from AIL mice.
# Overview
**Impact of host genetics and diet on the murine intestinal mycobiome.**

Yask Gupta, Anna Lara Ernst, Artem Vorobyev, Foteini Beltsiou, Detlef Zillikens, Katja Bieber, Simone Sanna-Cherchi, Angela M. Christiano, Christian D. Sadik, Ralf J. Ludwig and Tanya Sezin

**SUMMARY**

The mammalian gut is a home to a diverse microbial ecosystem, whose composition affects various physiological and pathophysiological traits of the host organism. Recently, next generation sequencing-based meta-genomic approaches demonstrated how the interplay of host genetics, microbiome, and environmental factors shape complex traits and clinical outcomes. Much of the data on microbial communities within mammalian hosts has been derived from the analysis of bacterial communities, whereas other communities, such as fungi, have been understudied. More recently, sequencing of the ribosomal DNA internal transcribed spacer (ITS) of nuclear DNA has helped to reveal the biological role of fungi and their regulation by host genetics. Here, we provide in vivo evidence in mice that fungal communities are regulated by host genetics. Using whole genome sequencing and genotyping, we mapped QTL associated with various fungal species to single genes in the host. Moreover, we show that diet and its’ interaction with host genetics alter the composition of fungi in outbred mice, and we identified fungal indicator species associated with different dietary regimes. Collectively, we uncovered an association of the intestinal fungal community with host genetics, as well as a regulatory role of diet on this ecological niche.  

# Documentation
All the script for the analysis are in script folder. 

If starting from Fastq files, the users needs to run PIPITS pipeline (for Fungi) and QIIME/vsearch (for bacteria) to identify OTUs. 

For phylum to genus level the users need to use RDP classifier. 

The output of these pipelines have been incorporated within the repository in raw_data folder as BIOM files (OTUs) and as heir file (RDP, Phylum to Genus).

For QTL mapping the HAPPY R packages prosterior probabilities (PProbs) array (NSNPs x Nsamples x PProbs[4 strains]) are required. This has been formatted as R object and deposited elsewhere (https://drive.google.com/file/d/1o4bPHlPkxsFReoqgnc65oBAFb3wbEN8m/view?usp=sharing) due to limited strorage in GitHub.

For correlation analysis Fastspar need to be install within the system https://github.com/scwatts/fastspar

## System requirements
### OS Requirements
The following code has been test on Ubuntu 16.04.
For performing permuation in QTL analysis server (CentOS) from University of Lubeck has been used.
### Dependencies

To run the R scripts the several packages need to be installed within R. These include:
+ indicspecies (v1.7.12)
+ pheatmap (v1.0.12)
+ phyloseq (v1.40.0)
+ reshape2 (v1.4.4)
+ ggplot2 (v3.3.6)
+ vegan (v2.6.2)
+ ggpubr (v0.4.0)
+ MANOVA.RM (v0.5.3)
+ DOQTL (v1.19.0) : https://rdrr.io/bioc/DOQTL/
+ DESeq2 (v1.36.0)
+ lme4 (v1.1.30)

To run perl scripts Perl should be installed in the system.
## Demo and Usage
The repository can be cloned in local machine using :

git clone https://github.com/Yask-Gupta/QTL_MYCO.git

Please copy the RData file (QTLRelInputHglm.RData) from the link provided above to Robj folder for execution of the code.

### Compositional analysis 
Compositional analysis of Fungi can be performed using script **"CompositionAnalysisMyco.R"** in script folder.

Compositional analysis of Microbiome (DNA) can be performed using script **"CompositionAnalysisMicroDNA.R"** in script folder.

Compositional analysis of Microbiome (RNA) can be performed using script **"CompositionAnalysisMicroRNA.R"** in script folder.

All the output file and figures will be generated in output folder.

#### QTL analysis
This analysis is computationally expensive due to permuation. The analysis requires minimum of 128 GB ram and it was executed on the server.

QTL analysis of Fungi can be performed using script **"MycoQTLPhylumTospecies.R"** in script folder.

QTL analysis of Microbiome (DNA) can be performed using script **"MicroDNAPhylumTospecies.R"** in script folder.

QTL analysis of Microbiome (RNA) can be performed using script **"MicroRNAPhylumTospecies.R"** in script folder.

The analysis will produce several RData file which will be stored in cuurent directory. Permutations and calculation of phenotypic variation has been coded within the scripts. 
