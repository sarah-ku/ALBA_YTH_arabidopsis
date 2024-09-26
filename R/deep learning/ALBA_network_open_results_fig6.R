library(rsample)
library(tidyverse)
library(zoo)
library(seqLogo)
library(Biostrings)
library(jsonlite)
library(ROCR)
library(universalmotif)


source("./ALBA_FUNCTIONS.R")

#########################
load(file="./preds.Rdat")
load(file="./data_tibble_splits_MO_fold.Rdat")

#data_tibble
testing_tibble <- do.call(rbind,(split_tibble %>% dplyr::mutate(test=purrr::map2(splits,id,~testing(.x) %>% mutate(fold=.y))))$test)
ftable(testing_tibble$fold,testing_tibble$gene_group)

names(preds) <- unique(testing_tibble$fold)
testing_tibble <- testing_tibble %>% mutate(num = 1:nrow(testing_tibble),
                                            OOF = map2(fold,num,~preds[[.x]][.y,]))

testing_tibble$OOF_ECT2 <- unlist(lapply(testing_tibble$OOF,function(x) x[1]))
testing_tibble$OOF_ALBA <- unlist(lapply(testing_tibble$OOF,function(x) x[2]))
testing_tibble$OOF_BOTH <- unlist(lapply(testing_tibble$OOF,function(x) x[3]))

ect2_fold <- rep(0,5)
alba_fold <- rep(0,5)
both_fold <- rep(0,5)

for(k in 1:5)
{
  myfd <- unique(testing_tibble$fold)[k]
  apred <- ROCR::prediction(preds[[k]][testing_tibble$fold==myfd,3],testing_tibble$label_both[testing_tibble$fold==myfd])
  pboth <- (performance(apred,measure="auc")@y.values[[1]])
  
  apred <- ROCR::prediction(preds[[k]][testing_tibble$fold==myfd,2],testing_tibble$label_ALBA[testing_tibble$fold==myfd])
  palba <- (performance(apred,measure="auc")@y.values[[1]])
  
  apred <- ROCR::prediction(preds[[k]][testing_tibble$fold==myfd,1],testing_tibble$label_ECT2[testing_tibble$fold==myfd])
  pect2 <- (performance(apred,measure="auc")@y.values[[1]])
  
  both_fold[k] <- pboth
  alba_fold[k] <- palba
  ect2_fold[k] <- pect2
  
  print(paste(c(pboth,palba,pect2)))
}


seq_list <- list()
for(i in 1:5)
{
  seq_list[[i]] <- Biostrings::readDNAStringSet(paste0("./model_MO_iclip_model_sequence_fasta_test_fold_SPLIT_",i,".fasta"))
}

seqs <- c(seq_list[[1]],seq_list[[2]],seq_list[[3]],seq_list[[4]],seq_list[[5]])

load(file="motif_list.Rdat")

df <- data.frame("ect2"=preds[[k]][,1],"alba"=preds[[k]][,2])

#############################################################
##DENSITY PLOT########################################
#############################################################

library(viridis)

TS <- ggplot(df, aes(x=ect2, y=alba) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + theme_bw() + 
  scale_fill_distiller(palette=2, direction=1)
  
ggsave(plot=TS,filename="density_ML_predictions.pdf",width=8,height=7)

#############################################################
##SUPPLEMENTARY TABLE########################################
#############################################################
load(file="./m6A_site_network_set.Rdat")

testing_tibble %>% arrange(index)

SUPP_TABLE <- data.frame(chr = seqnames(m6AGR),
pos = start(m6AGR),
gene = testing_tibble$gene,
name = m6AGR$name,
experiment = m6AGR$experiment,
probability = m6AGR$prediction,
prediction_ECT2 = testing_tibble$OOF_ECT2,
prediction_ALBA4 = testing_tibble$OOF_ALBA) %>% arrange(chr,desc(experiment))

write_csv(SUPP_TABLE,"./supplementary_m6A_sites.csv")
