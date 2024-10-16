load(file="posGR_EXPER_annoted_update.Rdat")
load(file="ALBA_project_SET_LIST.Rdat")

names(posGR_EXPER) <- c("ALBA2 ect234","ECT2 alba1245")

makePLOT <- function(K)
{
  res_tib <- posGR_EXPER[[K]]
  res_tib <- data.frame("pos"=res_tib$names,"trt"=res_tib$prop_trt_correc,"ctrl"=res_tib$prop_ctrl_correc,"gene"=res_tib$gene,"padj"=res_tib$padj,"fc"=res_tib$log2FC,"expr_trt"=res_tib$expr_trt,"expr_cntrl"=res_tib$expr_cntrl)
  row.names(res_tib) <- res_tib$pos
  
  res_tib <- res_tib[!is.na(res_tib$gene),]
  res_tib$sig <- res_tib$padj<0.1 & abs(res_tib$fc)>.5 & (res_tib$trt>0.01 | res_tib$ctrl>0.01)
  
  res_tib$ECT2_targ <- res_tib$gene %in% set_list$ECT2_permissive
  res_tib$ALBA4_targ <- res_tib$gene %in% set_list$ALBA4_3UTR
  res_tib$ECT2_targ_strict <- res_tib$gene %in% set_list$ECT2_strict
  res_tib$ALBA2_targ <- res_tib$gene %in% set_list$ALBA2_permissive
  res_tib$dual_targ <- res_tib$gene %in% set_list$ECT2_ALBA4_dual
  
  p1 <- res_tib %>% 
    mutate(colour = paste(sig,dual_targ,sep="_")) %>%
    mutate(colour=factor(colour,levels=c("FALSE_FALSE","FALSE_TRUE","TRUE_FALSE","TRUE_TRUE"))) %>%
    arrange(colour) %>% 
    ggplot(aes(x=ctrl,y=trt,colour=colour,group=colour)) + 
    geom_point(size=0.33) +
    scale_color_manual(values=c("light grey","dark grey","blue","orange")) +
    geom_abline(slope = 1,intercept = 0,lty=2) +
    scale_y_sqrt() + scale_x_sqrt() +
    ggtitle(names(posGR_EXPER)[K]) + 
    xlab("Editing proportions A/(A+G) in CTRL (corrected)") +
    ylab("Editing proportions A/(A+G) in TRT (corrected)") +
    theme_bw()
  
  #p1
  prop_dual_bound_sig <- tapply(res_tib[res_tib$sig,]$dual_targ,res_tib[res_tib$sig,]$fc>0,mean)
  prop_ECT2_sig <- tapply(res_tib[res_tib$sig,]$ECT2_targ,res_tib[res_tib$sig,]$fc>0,mean)
  prop_ALBA4_sig <- tapply(res_tib[res_tib$sig,]$ALBA4_targ,res_tib[res_tib$sig,]$fc>0,mean)
  prop_ALBA2_sig <- tapply(res_tib[res_tib$sig,]$ALBA2_targ,res_tib[res_tib$sig,]$fc>0,mean)
  
  vlines <- tibble(name=c("dual_TARG","ECT2_TARG","ALBA4_TARG","ALBA2_TARG"),
                   mean=c(mean(res_tib$dual_targ),mean(res_tib$ECT2_targ),mean(res_tib$ALBA4_targ),mean(res_tib$ALBA2_targ)))
  
  p3 <- tibble("type"=c("down","up"),
               "dual_TARG"=prop_dual_bound_sig,
               "ECT2_TARG"=prop_ECT2_sig,
               "ALBA4_TARG"=prop_ALBA4_sig,
               "ALBA2_TARG"=prop_ALBA2_sig) %>% 
    pivot_longer(-type) %>%
    ggplot(aes(x=type,group=name,y=value,fill=name,colour=name)) + 
    geom_bar(stat="identity",show.legend=F) + facet_grid(~name) + 
    ylab("Proportion which are targets") +
    xlab("Change type from HyperTRIBE") + 
    geom_hline(data = vlines,aes(yintercept=mean),colour="black",lty=2) + 
    theme_bw()
  
  degenes <- unique(res_tib$gene)
  FC <-  design_DE$DE[[K]]$log2FoldChange
  names(FC) <- row.names(design_DE$DE[[K]])
  FC_mean <- tapply(FC, gsub("\\.[0-9]","",as.vector(names(FC))),mean)
  res_tib$DE <- FC_mean[res_tib$gene]
  p4 <- res_tib  %>% 
    mutate(FC = fc>0) %>% mutate(status=ifelse(sig,"signif","not_signif")) %>%
    mutate(status = ifelse(sig & (fc>0),"up",status)) %>% 
    mutate(status = ifelse(sig & (fc<0),"down",status)) %>%
    group_by(status,gene) %>% summarise(DE=mean(DE)) %>% drop_na() %>%
    ggplot(aes(y=DE,x=status))  +
    geom_boxplot(outlier.shape = NA,notch=F) + 
    ylim(-1.25,1.25) + ylab("DEA log2 fold-change") +
    theme_bw() + geom_hline(yintercept = 0,lty=2)
  #[res_tib$sig,]
  p5 <- res_tib %>% ggplot(aes(x=as.factor(dual_targ),y=fc,fill=as.factor(dual_targ))) + 
    geom_boxplot() + 
    theme_bw() + geom_hline(yintercept=0,lty=2) + 
    ylab("log2 fold change in editing prop.") 
  
  print(t.test(res_tib[res_tib$sig,]$fc ~ res_tib[res_tib$sig,]$dual_targ))
  
  lay <- rbind(c(1,1,4,4),
               c(1,1,4,4),
               c(1,1,4,4),
               c(2,2,5,5),
               c(2,2,5,5))
  TS <- grid.arrange(grobs = list(p1,p3,p4,p5), layout_matrix = lay)
  
  return(TS)
}


TS <- makePLOT(1)
ggsave(TS,file="INTERACTION_plots_robust_exper1_updated.pdf",width=15,height=10)

TS <- makePLOT(2)
ggsave(TS,file="INTERACTION_plots_robust_exper2_updated.pdf",width=15,height=10)