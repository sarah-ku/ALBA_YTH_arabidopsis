library(abind)
library(reticulate)
np <- import("numpy")

features_1 <- list()
for(fold in 1:5)
{
  npz1 <- np$load(paste0("FEATURES_1_10CV_fold_",fold,".npz"))
  features_1[[fold]] <- npz1$f[["arr_0"]]
}
features_1 <- abind(features_1)
dim(features_1)


preds <- list()
for(fold in 1:5)
{
  p1 <- np$load(paste0("PREDICTIONS_MO_10CV_fold_",fold,".npz"))
  p1 <- p1$f[["arr_0"]][,,1]
  p1 <- t(p1)
  preds[[fold]] <- p1
}


grad_1 <- list()
grad_2 <- list()
grad_3 <- list()
for(fold in 1:5)
{
  npz1 <- np$load(paste0("GRADS_1_MO_10CV_fold_",fold,".npz"))
  grad_1[[fold]] <- npz1$f[["arr_0"]]
  
  npz1 <- np$load(paste0("GRADS_2_MO_10CV_fold_",fold,".npz"))
  grad_2[[fold]] <- npz1$f[["arr_0"]]
  
  npz1 <- np$load(paste0("GRADS_3_MO_10CV_fold_",fold,".npz"))
  grad_3[[fold]] <- npz1$f[["arr_0"]]
}

grad_1 <- abind(grad_1,along=1)
grad_2 <- abind(grad_2,along=1)
grad_3 <- abind(grad_3,along=1)
grad_list <- list(grad_1,grad_2,grad_3)

filters_1 <- list()
for(fold in 1:5)
{
  npz1 <- np$load(paste0("WEIGHTS_CONV_1_10CV_fold_",fold,".npz"))
  filters_1[[fold]] <- npz1$f[["arr_0"]]
}
filters_1 <- abind(filters_1)
dim(filters_1)

library(Biostrings)
library(universalmotif)

seq_list <- list()
for(i in 1:5)
{
  seq_list[[i]] <- Biostrings::readDNAStringSet(paste0("model_MO_iclip_model_sequence_fasta_test_fold_SPLIT_",i,".fasta"))
}

seqs <- c(seq_list[[1]],seq_list[[2]],seq_list[[3]],seq_list[[4]],seq_list[[5]])

seq_mat <- as.matrix(seqs)

nt_freqs <- as.vector(table(as.vector(seq_mat)))/sum(as.vector(table(as.vector(seq_mat))))
names(nt_freqs) <- c("A","T","C","G")

motList <- list()
for(i in 1:(64*5))
{
  print(i)
  pool_size <- 1
  r_mat <- features_1[,,i]
  filter_size <- 8
  cut_off <- .5
  
  motList[[i]] <- returnMotif(r_mat,pool_size = pool_size,MOT_ID = i,cut_off = 0,max_instances=2500,filter_size = filter_size)
}


#save everything
save(motList,file="motif_list.Rdat")
save(filters_1,file="filters.Rdat")
save(grad_list,file="grad_list.Rdat")
save(preds,file="preds.Rdat")
save(features_1,file="features.Rdat")
