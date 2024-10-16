
#######################################################################
########################FIGS8 A########################################
#######################################################################

makeSupportPlot <- function(ALBA4_clip_count,expr)
{
  breaks <- c(0,1,10,50,100,1000)
  HT_support <- tapply(names(ALBA4_clip_count),cut(ALBA4_clip_count,breaks=breaks),function(x) mean(x %in% ALBA4_HT))

  bg_sets <- list()
  for(i in 1:1000)
  {
    rset <- getSetMatchedExpression(ALBA4_HT,log(expr_TRIBE+1),n_breaks = 10)
    HT_support_grey <- tapply(rset[sample(length(ALBA4_clip_count),length(ALBA4_clip_count))],cut(ALBA4_clip_count,breaks=breaks),function(x) mean(x %in% ALBA4_HT))
    bg_sets[[i]] <- HT_support_grey
  }
  
  tib <- as_tibble(t(apply(do.call(cbind,bg_sets),1,function(x) c("y"=mean(x),quantile(x,c(0.025,0.975))))),rownames="interval")
  colnames(tib) <- c("interval","y","low","up")
  tib <- tib %>% mutate(true=HT_support) %>% mutate(interval=factor(interval,levels=unique(interval)))
  tib %>% mutate(x=c(1:nrow(tib))) %>% ggplot(aes(interval, y)) + 
    geom_errorbar(aes(ymin = low, ymax = up)) + geom_point(aes(y=true)) +
    ylab("Proportion supported by ALBA4 TRIBE") + 
    xlab("Number of iCLIP sites in gene") + ylim(0,1) +
    theme_bw()
}

ALBA4_clip_count <- tapply(ALBAGR$gene,ALBAGR$gene,length)
ALBA4_clip_count <- ALBA4_clip_count[-1]
sp1 <- makeSupportPlot(ALBA4_clip_count,expr=expr_TRIBE)

ALBA4_clip_count <- tapply(ALBAGR_R$gene,ALBAGR_R$gene,length)
ALBA4_clip_count <- ALBA4_clip_count[-1]
sp2 <- makeSupportPlot(ALBA4_clip_count,expr=expr_TRIBE)

ts <- grid.arrange(sp1,sp2)
ggsave(plot=ts,filename="ALBA4_iCLIP_support_of_genes_TRIBE_CIs_filtered.pdf",width=5,height=6)

#######################################################################
########################FIGS8 B########################################
#######################################################################

ECT2 <- unique(iCLIP[[3]]$gene)
is_ECT2 <- names(ALBA4_clip_count) %in% ECT2

supp_bar_3 <- as_tibble(do.call(rbind,tapply(ALBA4_clip_count,is_ECT2,function(x) c(sum(x<10),sum(x>=10)))),rownames="ECT2_status") %>%
  dplyr::rename("<10 sites"=V1,">=10 sites"=V2) %>% pivot_longer(-ECT2_status) %>%
  ggplot(aes(x=ECT2_status,y=value,fill=name)) + geom_bar(stat="identity",position="dodge") +
  theme_bw() + ylab("Number of genes") + xlab("Also ECT2 gene (ECT2 iCLIP)") + ggtitle("ALBA4 iCLIP gene counts \n by ECT2 status")

ggsave(plot=supp_bar_3,filename="ALBA4_iCLIP_number_of_sites_support_panel.pdf",width=8,height=4)

#######################################################################
########################FIGS8 C########################################
#######################################################################

ALBA24_tab <- tibble("ALBA2"=log(ALBA2_expr+1),"ALBA4"=log(ALBA4_expr+1)) %>% 
  mutate(gene=names(ALBA2_expr)) %>% 
  mutate(ALBA2_target = gene %in% ALBA2_HT) %>%
  mutate(ALBA4_target = gene %in% ALBA4_HT) %>%
  mutate(ALBA2_specific = ALBA2_target & !ALBA4_target) %>%
  mutate(ALBA4_specific = ALBA4_target & !ALBA2_target)

plot_ts_expr <- ALBA24_tab %>% dplyr::filter(ALBA2_specific | ALBA4_specific) %>%
  ggplot(aes(x=ALBA2,y=ALBA4,colour=ALBA2_specific)) + 
  geom_bin2d(bins=200,aes(fill=ALBA2_specific)) + 
  geom_abline(slope=1,intercept = 0,lty=2) + 
  theme_bw() + ylab("log [ALBA4 TPM + 1]") + 
  xlab("log [ALBA2 TPM + 1]") +
  ggtitle("Expression of specific targets")

ggsave(plot=plot_ts_expr,filename="ALBA2_HT_vs_ALBA4_TRIBE_specific_target_expression.pdf",width=5,height=4)

#######################################################################
########################FIGS8 D########################################
#######################################################################

ALBA4_RANDOM <- getSetMatchedExpression(ALBA4,log(expr+1),n_breaks = 10)
ECT2_RANDOM <- getSetMatchedExpression(ECT2,log(expr+1),n_breaks = 10)

mv_real <- makeVenn(list("ALBA_iCLIP"=ALBA4,"ECT2_iCLIP"=ECT2))
mv_alba_iclip <- makeVenn(list("ALBA_iCLIP"=ALBA4,"ECT2_RANDOM"=ECT2_RANDOM))
mv_ect2_iclip <- makeVenn(list("ALBA_iCIP_RANDOM"=ALBA4_RANDOM,"ECT2_iCLIP"=ECT2))
venns_tp <- grid.arrange(grobs=list(mv_real,mv_alba_iclip,mv_ect2_iclip))
ggsave(plot=venns_tp,filename="ALBA4_iCLIP_ECT2_iCLIP_overlaps_with_random_sets.pdf",width=5,height=12)

mv_alba4 <- makeVenn(list("ALBA4_strict"=ALBA4_strict,"ALBA4_permissive"=ALBA4_permissive))
mv_ect2 <- makeVenn(list("ECT23_strict"=ECT2_strict,"ECT23_permissive"=ECT2_permissive))
venns_tp <- grid.arrange(grobs=list(mv_alba4,mv_ect2))
ggsave(plot=venns_tp,filename="ALBA4_ECT2_strict_and_permissive.pdf",width=3,height=6)

#######################################################################
########################FIGS8 E########################################
#######################################################################

ALBA4_expr <- read.table("ALBA4_TRIBE_expression_TPM.csv")
ALBA2_expr <- read.table("ALBA2_HyperTRIBE_expression_TPM.csv")
expr <- rowSums(cbind(as.data.frame(ALBA4_expr[,c(1:5)]),as.data.frame(ALBA2_expr[,c(1:5)])))
expr <- tapply(expr,gsub("\\.[0-9]+","",names(expr)),mean)


ALBA4_permissive <- as.data.frame(genes_ALBA)[,2][genes_ALBA$ALBA4_permissive]
ECT2_permissive <- as.data.frame(genes_ALBA)[,2][genes_ALBA$ECT23]

ALBA4_permissive_RANDOM <- getSetMatchedExpression(ALBA4_permissive,log(expr+1),n_breaks = 10)
ECT2_permissive_RANDOM <- getSetMatchedExpression(ECT2_permissive,log(expr+1),n_breaks = 10)

mv_real <- makeVenn(list("ALBA4_permissive"=ALBA4_permissive,"ECT2_permissive"=ECT2_permissive))
mv_alba_permissive <- makeVenn(list("ALBA_permissive"=ALBA4_permissive,"ECT2_permissive_RANDOM"=ECT2_permissive_RANDOM))
mv_ect2_permissive <- makeVenn(list("ALBA_permissive_RANDOM"=ALBA4_permissive_RANDOM,"ECT2_permissive"=ECT2_permissive))
venns_tp <- grid.arrange(grobs=list(mv_real,mv_alba_permissive,mv_ect2_permissive))

ggsave(plot=venns_tp,filename="ALBA4_permissive_ECT2_permissive_overlaps_with_random_sets.pdf",width=5,height=12)

#######################################################################
########################FIGS8 F########################################
#######################################################################

ALBAGR_R$is_3UTR <- unlist(lapply(strsplit(ALBAGR_R$transcript.types,","),function(x) "3UTR" %in% x))
ALBA4_3UTR <- (unique(ALBAGR_R[ALBAGR_R$is_3UTR]$gene))
ALBA4_not3UTR <- (unique(ALBAGR_R[!ALBAGR_R$is_3UTR]$gene))

mv_real <- makeVenn(list("ALBA4_3UTR"=ALBA4_3UTR,"ECT2_strict"=ECT2_strict))
mv_alba2_random <- makeVenn(list("ALBA4_3UTR"=ALBA4_3UTR,"ECT23_RANDOM"=ECT2_strict_RANDOM))
mv_ect23_random <- makeVenn(list("ALBA4_3UTR_RANDOM"=ALBA4_3UTR_RANDOM,"ECT2_strict"=ECT2_strict))
venns_tp <- grid.arrange(grobs=list(mv_real,mv_alba2_random,mv_ect23_random))
ggsave(plot=venns_tp,filename="ALBA4_3UTR_ECT2_strict_overlaps_with_random_sets.pdf",width=5,height=12)

mv_realA <- makeVenn(list("ALBA4_3UTR"=ALBA4_3UTR,"ECT2_strict"=ECT2_strict))
mv_realB <- makeVenn(list("ALBA4_not3UTR"=ALBA4_not3UTR,"ECT2_strict"=ECT2_strict))
venns_tp <- grid.arrange(grobs=list(mv_realA,mv_realB))
ggsave(plot=venns_tp,filename="ALBA4_3UTR_ECT2_strict_ALBA4_not_3UTR_ECT2_strict.pdf",width=5,height=12)