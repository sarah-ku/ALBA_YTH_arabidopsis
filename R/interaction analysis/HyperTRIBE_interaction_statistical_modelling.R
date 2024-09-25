###########################################################################################
###########################################################################################

###REQUIRED OBJECTS FROM ZENODO
load(file = "design_with_DEA_interaction_fig7.Rdat")
load(file = "design_with_results_interaction_fig7.Rdat")
load(file="counts_HyperTRIBE_interaction_fig7.Rdat") #counts.mat
load(file="abundance_HyperTRIBE_interaction_fig7.Rdat") #tpm.mat

###########################################################################################
###########################################################################################
########################### FUNCTIONS TO RUN ANALYSES #####################################
###########################################################################################
###########################################################################################

make_data_int <- function(x){
  data_list <- design$data_list[[x]]
  adar_control <-  design$adar_control[[x]]
  adar_treatment <- design$adar_treatment[[x]]
  data_int <- list("data_list"=data_list,"adar_control"=adar_control,"adar_treatment"=adar_treatment)
  return(data_int)
}

makeCountData <- function(data_int)
{
  #data_int
  adar <- c(data_int$adar_control,data_int$adar_treatment)
  positions_1 <- row.names(data_int$data_list[[1]])
  locsP <- locsGR[positions_1]
  
  #get A+G and G counts
  atg <- locsP[locsP$ref=="A"]$names
  data_atg <- lapply(data_int$data_list,function(x) x[atg,])
  mexp <- do.call(cbind,(lapply(data_atg[1:10],function(x) x[,"A"]+x[,"G"])))
  meG <- do.call(cbind,(lapply(data_atg[1:10],function(x) x[,"G"])))
  
  #get T+C and C counts
  ttc <- locsP[locsP$ref=="T"]$names
  data_atg <- lapply(data_int$data_list,function(x) x[ttc,])
  mexpTC <- do.call(cbind,(lapply(data_atg[1:10],function(x) x[,"T"]+x[,"C"])))
  meC <- do.call(cbind,(lapply(data_atg[1:10],function(x) x[,"C"])))
  
  #join together A>G and T>C candidates
  mexp <- as_tibble(mexp,rownames="pos") %>% full_join(as_tibble(mexpTC,rownames="pos"))
  meG <- as_tibble(meG,rownames="pos") %>% full_join(as_tibble(meC,rownames="pos"))
  
  condition <- c(rep("contrl",5),rep("treat",5))
  
  #normalisation by library depth
  names(condition) <- colnames(meG[,-1])
  os <- colSums(mexp[,-1])/1e08
  #os <- os/mean(os)
  os <- total_counts[names(os)]/1e07
  os <- os/mean(os)
  
  etib <- mexp %>% pivot_longer(-pos) %>% dplyr::rename(expr=value) %>% 
    left_join(meG %>% pivot_longer(-pos), by=c("pos","name")) %>% 
    dplyr::rename(G=value) %>% left_join(as_tibble(adar,rownames="name"),by="name") %>% dplyr::rename(adar=value) %>%
    left_join(as_tibble(condition,rownames="name"),by="name") %>% dplyr::rename(condition=value) %>%
    left_join(as_tibble(os,rownames="name"),by="name") %>% dplyr::rename(os=value)
  
  etib <- etib %>% mutate(actual = G/expr)
  etib <- etib %>% group_by(pos,condition) %>% mutate(mean_prop = mean(actual)) %>% ungroup()
  
  #filtering
  prop_supp_pos <- etib %>% group_by(pos,condition) %>% summarise(mp = mean(mean_prop)) %>% 
    pivot_wider(id_cols=pos,names_from=condition,values_from=mp) %>% ungroup() %>%
    dplyr::filter( (contrl>0.005 | treat>0.005) & (contrl<0.99 & treat<0.99)) %>% pull(pos)
  
  rep_supp_pos <- etib %>% group_by(pos,condition) %>% 
    summarise(rep_supp = length(G[G>0])) %>% ungroup() %>% 
    pivot_wider(id_cols=pos,names_from=condition,values_from=rep_supp) %>% dplyr::filter(contrl>2 | treat>2) %>% pull(pos)
  
  expr_supp_pos <- etib %>% group_by(pos,condition) %>%
    summarise(expr = sum(expr)) %>% ungroup() %>%
    pivot_wider(id_cols=pos,names_from=condition,values_from=expr) %>% dplyr::filter(contrl>50 | treat>50) %>% pull(pos)
  
  pos_use <- Reduce(intersect,list(expr_supp_pos,rep_supp_pos,prop_supp_pos))
  
  etib <- etib %>% filter(pos %in% pos_use)
  
  COUNT_DATA <- etib
  return(COUNT_DATA)
}


makeModelData <- function(etibS)
{
  #prepare the data for modelling
  etibS$logADAR <- log(etibS$adar)
  etibS$logEXPR <- log(etibS$expr+1)
  etibS$POS <- as.numeric(factor(etibS$pos,levels=unique(etibS$pos)))
  etibS$COND <- as.numeric(factor(etibS$condition,levels=c("contrl","treat")))
  etibS$SAMP <- as.numeric(factor(etibS$name,levels=unique(etibS$name)))
  etibS$EXPR_BIN <- as.numeric(cut(etibS$logEXPR,5))
  etibS$ADAR_BIN <- as.numeric(cut(etibS$logADAR,3))
  
  ftable(etibS$COND,etibS$ADAR_BIN,etibS$EXPR_BIN)
  table(etibS$EXPR_BIN)
  data_fit <- etibS %>% ungroup() %>% dplyr::select(pos,expr,G,condition,logADAR,ADAR_BIN,EXPR_BIN,COND,POS)
  dim(data_fit)

  toPRED <- etibS %>% group_by(pos) %>% 
    mutate(meanADAR = 2) %>% group_by(pos,condition) %>% 
    mutate(expr = round(mean(expr)), G = round(mean(G)), ADAR_BIN=mean(meanADAR), EXPR_BIN = round(mean(EXPR_BIN)))
  
  toPRED <- toPRED %>% 
    ungroup() %>% 
    dplyr::select(expr,G,condition,ADAR_BIN,EXPR_BIN,COND,POS) %>% distinct()
  toPRED$G <- NA
  
  toPRED2 <- etibS %>% group_by(pos,condition) %>% 
    mutate(meanADAR = 2, 
           expr = round(mean(expr)), 
           G = round(mean(G)), 
           ADAR_BIN=mean(meanADAR), 
           EXPR_BIN = round(mean(EXPR_BIN)))
  
  toPRED2 <- toPRED2 %>% ungroup() %>% 
    dplyr::select(expr,G,condition,ADAR_BIN,EXPR_BIN,COND,POS) %>% 
    distinct()
  toPRED2$G <- NA
  
  data_fit <- data_fit %>% full_join(toPRED) %>% full_join(toPRED2)
  #dim(data_fit)
  return(data_fit)
}

#########################################################################################################
#########################################################################################################
library(INLA)
#make linear combinations in order to test for the significance of trt vs ctrl for each position
getLCs <- function(etibS)
{
  ncombs <- length(unique(etibS$pos))
  ID <- rep(1:ncombs,times=2)
  n.cond <- 2
  B <- matrix(NA,ncol=ncombs*n.cond,nrow=ncombs)
  
  for(i in 1:nrow(B))
  {
    B[i,which(ID==i)] <- c(-1,1)
  }
  
  lc <- INLA::inla.make.lincombs(POS=B)
  return(lc)
}

makeModel <- function(data_fit,lc,mypath)
{
  formula = G ~ 
    f(ADAR_BIN, model = "iid", group=EXPR_BIN) + 
    
  
  formula = G ~ ADAR_BIN*EXPR_BIN + f(POS,model="iid",group = COND)
  
  result = INLA::inla(formula, data = data_fit,
                      family = "binomial",
                      verbose = T,
                      Ntrials = expr,
                      lincomb = lc,
                      num.threads = 20,
                      control.predictor = list(compute = T),
                      control.compute = list(dic = T,config=TRUE,return.marginals.predictor=TRUE),
                      working.directory = paste("./INLA/update/"))
  
  save(result,file=mypath)
  return(result)
}

##########################################################
################RUN FUNCTIONS#############################
##########################################################

data_int <- make_data_int(1)
COUNT_DATA <- makeCountData(data_int)
data_fit <- makeModelData(COUNT_DATA)
lc <- getLCs(COUNT_DATA)
mypath <- "./INLA/update/result_robust_exper_1.Rdat"
result <- makeModel(data_fit,lc,mypath)
#save(result,file=mypath)

data_int <- make_data_int(2)
COUNT_DATA <- makeCountData(data_int)
data_fit <- makeModelData(COUNT_DATA)
lc <- getLCs(COUNT_DATA)
mypath <- "./INLA/update/result_robust_exper_2.Rdat"
result <- makeModel(data_fit,lc,mypath)

#########################################################################################################
#########################################################################################################

collectResults <- function(INT)
{
  mypath <- paste0("./INLA/update/result_robust_exper_",INT,".Rdat")
  
  load(file=mypath)
  
  data_int <- make_data_int(INT)
  COUNT_DATA <- makeCountData(data_int)
  data_fit <- makeModelData(COUNT_DATA)
  dim(data_fit)
  head(data_fit)
  fitted <- result$summary.fitted.values$mean
  length(fitted)
  data_fit$fitted <- exp(fitted)/(1+exp(fitted))
  
  data_fit_subs <- data_fit  %>% filter(is.na(G) & ADAR_BIN==2)
  
  corrected_props <- data_fit_subs %>% pivot_wider(values_from=fitted,id_cols=POS,names_from=condition)
  pos_map <- tapply(data_fit$pos,data_fit$POS,function(x) x[1])
  corrected_props$pos <- pos_map[corrected_props$POS]

  posGR_meta <- data_fit %>% 
    dplyr::select(pos,condition,expr,G) %>% 
    group_by(pos,condition) %>% 
    dplyr::summarise("mean_prop"=mean(G/expr),"mean_expr"=mean(expr),"mean_G"=mean(G)) %>%
    ungroup() %>%
    pivot_wider(names_from=condition,values_from=c(mean_prop,mean_expr,mean_G),id_cols=c(pos))
  
  names_use <- unique(data_fit$pos)
  names_use <- names_use[!is.na(names_use)]
  
  posGR_meta <- as.data.frame(posGR_meta %>% filter(pos %in% names_use))
  row.names(posGR_meta) <- posGR_meta$pos
  
  posGR <- locsGR[names_use]
  
  z <- result$summary.lincomb.derived$mean/result$summary.lincomb.derived$sd
  pval <- 2*pnorm(1-abs(z))
  names(pval) <- pos_map
  plot(z,-log10(pval))
  
  posGR$padj <- NA
  posGR[names(pval)]$padj <- pval
  posGR$padj <- p.adjust(posGR$padj,method="BH")
  
  posGR$prop_trt <- posGR_meta[names(posGR),]$mean_prop_treat
  posGR$prop_ctrl <- posGR_meta[names(posGR),]$mean_prop_contrl
  
  corrected_trt <- corrected_props$treat
  names(corrected_trt) <- corrected_props$pos
  posGR$prop_trt_correc <- corrected_trt[names(posGR)]
  
  corrected_cntrl <- corrected_props$contrl
  names(corrected_cntrl) <- corrected_props$pos
  posGR$prop_ctrl_correc <- corrected_cntrl[names(posGR)]
  
  posGR$expr_trt <-  posGR_meta[names(posGR),]$mean_expr_treat
  posGR$expr_cntrl <-  posGR_meta[names(posGR),]$mean_expr_contrl
  
  strand(posGR[posGR$ref=="A"]) <- "+"
  strand(posGR[posGR$ref=="T"]) <- "-"
  
  posGR$log2FC <- log2(posGR$prop_trt_correc/posGR$prop_ctrl_correc)
  boxplot(posGR$log2FC)
  abline(h=0,lty=2)
  
  plot(posGR$prop_trt_correc~posGR$prop_ctrl_correc,cex=0.5)
  abline(0,1)
  points(posGR[posGR$padj<0.1]$prop_trt_correc~posGR[posGR$padj<0.1]$prop_ctrl_correc,col="red",cex=0.5)

  
  return(posGR)
}

posGR_E1 <- collectResults(1)
posGR_E2 <- collectResults(2)

library(hyperTRIBER)
library(foreach)
library(doParallel)
library(GenomicRanges)
library(doParallel)
library(foreach)

quantvec <- rowMeans(tpm.mat)

posGR_E1 <- addGenes(gtfGR,posGR_E1,20,quant = quantvec,assignStrand = T,geneids = ids)
posGR_E2 <- addGenes(gtfGR,posGR_E2,20,quant = quantvec,assignStrand = T,geneids = ids)

posGR_EXPER <- list(posGR_E1,posGR_E2)

save(posGR_EXPER,file="posGR_EXPER_annoted_update.Rdat")