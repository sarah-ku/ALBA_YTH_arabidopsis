library(tidyverse)
library(GenomicRanges)
library(gridExtra)

#########SEQUENCE########
library(BSgenome.Athaliana.TAIR.TAIR9)
seqlevels(BSgenome.Athaliana.TAIR.TAIR9) <- c("1","2","3","4","5","Pt","Mt")

########ARAPORT GTF#######
gtf <- rtracklayer::import("./Araport11_GFF3_genes_transposons.201606.gtf")
seqlevels(gtf) <- c("1","2","3","4","5","Pt","Mt")
gtfGR <- gtf

########Gene descriptions#######
gene.desc <- read.csv("./gene_descriptions_arabidopsis.tsv",sep="\t",header=F)
row.names(gene.desc) <- gene.desc$V1
load(paste0("./gene_ids.Rdat"))

#########ALBA GRANGE OBJECTS########
load(file="./ALBAGR_with_meta.Rdat")
load(file="./ALBAGR_R_with_meta.Rdat")
ALBAGR_R_STR <- ALBAGR_R[ALBAGR_R$score>20]


########ECT2 iCLIP#######
load(file="./iCLIP_with_meta.Rdat")


########nanopore m6A#######
load(file="./parker19_nanoGR.Rdat")
nanopore_genes <- unique(nanoGR$gene)
load(file="./parker19_micGR.Rdat")

ALBA_genes <- unique(ALBAGR$gene)
ECT_genes <- unique(iCLIP[[3]]$gene)
ECT_genes <- unique(ECT_genes[!(ECT_genes=="")])

#HYPERTRIBE RESULTs FOR ALBA2 AND ALBA4
load(file="./ALBA2_results_design.Rdat")
design_ALBA2 <- design
load("./ALBA4_TRIBE_results_design.Rdat")

####motifs GRANGE######
load(file="./motGR.Rdat")


#####GENE SET LISTS########
load(file="./ALBA_project_SET_LIST.Rdat")