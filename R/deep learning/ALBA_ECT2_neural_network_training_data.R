library(jsonlite)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(BSgenome)
library(GenomicRanges)
library(rsample)

##########################################################################
################GENERATE TIBBLE FOR SITES/LABELS##########################
##########################################################################

load(file="m6AGR_combined.Rdat")
m6AGR <- allGRE
NT <- 100
m6AGR$near_iclip <- countOverlaps(m6AGR+NT,ECT2GR)>0
m6AGR$near_alba <- countOverlaps(m6AGR+NT,ALBAGR)>0
m6AGR$label <- as.numeric(as.factor(m6AGR$label))

CENTER_GR <- m6AGR 
CENTER_GR <- CENTER_GR + 300
CENTER_GR <- CENTER_GR[as.vector(strand(CENTER_GR)) %in% c("+","-")]

seqs <- getSeq(BSgenome.Athaliana.TAIR.TAIR9,CENTER_GR)

data_tibble <- tibble(index = 1:length(CENTER_GR))
data_tibble <- data_tibble %>% dplyr::mutate(loc=paste0(seqnames(CENTER_GR),"_",end(CENTER_GR)),
                                             strand = as.vector(strand(CENTER_GR)),
                                             label_ALBA = as.numeric(CENTER_GR$near_alba),
                                             label_ECT2 = as.numeric(CENTER_GR$near_iclip),
                                             label_both  = as.numeric(CENTER_GR$near_alba & CENTER_GR$near_iclip),
                                             gene = as.vector(CENTER_GR$gene))

gene_names <- unique(data_tibble$gene)
gene_group <- sample(1:5,length(gene_names),replace=T)
names(gene_group) <- gene_names
data_tibble$gene_group <- gene_group[as.vector(data_tibble$gene)]

##########################################################################
################GENERATE TIBBLE FOR SITES/LABELS##########################
##########################################################################

#load(file="data_tibble_MO_fold.Rdat")

#leave out one gene group at a time
set.seed(49089)
split_tibble <- data_tibble %>% group_vfold_cv(group = gene_group)

for(i in 1:nrow(split_tibble))
{
  label_split_TRAIN <- training(split_tibble$splits[[i]]) %>% pull(label_ALBA)
  label_split_TEST <- testing(split_tibble$splits[[i]]) %>% pull(label_ALBA)
  
  write_json(label_split_TRAIN, paste0("model_MO_train_label_SPLIT_ALBA_fold_",i,".json"))
  write_json(label_split_TEST, paste0("model_MO_test_label_SPLIT_ALBA_fold_",i,".json"))
  
  label_split_TRAIN <- training(split_tibble$splits[[i]]) %>% pull(label_ECT2)
  label_split_TEST <- testing(split_tibble$splits[[i]]) %>% pull(label_ECT2)
  
  write_json(label_split_TRAIN, paste0("model_MO_train_label_SPLIT_ECT2_fold_",i,".json"))
  write_json(label_split_TEST, paste0("model_MO_test_label_SPLIT_ECT2_fold_",i,".json"))
  
  label_split_TRAIN <- training(split_tibble$splits[[i]]) %>% pull(label_both)
  label_split_TEST <- testing(split_tibble$splits[[i]]) %>% pull(label_both)
  
  write_json(label_split_TRAIN, paste0("model_MO_train_label_SPLIT_both_fold_",i,".json"))
  write_json(label_split_TEST, paste0("model_MO_test_label_SPLIT_both_fold_",i,".json"))
  
  seqs_split_TRAIN <- seqs[training(split_tibble$splits[[i]]) %>% pull(index)]
  seqs_split_TEST <- seqs[testing(split_tibble$splits[[i]]) %>% pull(index)]
  
  writeXStringSet(seqs_split_TRAIN, paste0("model_MO_iclip_model_sequence_fasta_train_fold_SPLIT_",i,".fasta"), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
  
  writeXStringSet(seqs_split_TEST, paste0("model_MO_iclip_model_sequence_fasta_test_fold_SPLIT_",i,".fasta"), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
}