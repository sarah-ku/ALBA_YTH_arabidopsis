expr <- rowMeans(txi$abundance[,c(1:5,11:15)])

########################################################################
#################BINDING BASES#########################################
########################################################################
seqnames(BSgenome.Athaliana.TAIR.TAIR9) <- c(rep(1:5),"Mt","Pt")
ALBA_seq <- getSeq(BSgenome.Athaliana.TAIR.TAIR9,ALBAGR[which(strand(ALBAGR) %in% c("+","-"))])

seq_tab <- tibble("seq" = c("A","C","G","T"),
                  "all" = table(as.vector(ALBA_seq)),
                  "strong" = table(as.vector(ALBA_seq[names(ALBAGR_R)])))

base_XL <- seq_tab %>% tidyr::pivot_longer(-seq) %>% group_by(name) %>% mutate(value = value/sum(value)) %>%
  ggplot(aes(x=seq,y=value,fill=name)) + 
  geom_bar(stat="identity",position="dodge",colour="black") + 
  theme_bw() + ggtitle("Which base is XL-d, ALBA1 iCLIP")

ggsave(plot=base_XL,filename = "./ALBA_binding_base_all_strong.pdf",width=5,height=6)

########################################################################
#################METAGENE PLOTS#########################################
########################################################################

ALBA_transcripts <- unique(unlist(lapply(base::strsplit(ALBAGR$transcripts,","),function(x) x[1])))
gtfGR_bg <- gtfGR[gtfGR$transcript_id %in% ALBA_transcripts]
gtfGR_bg <- gtfGR_bg[gtfGR_bg$type %in% c("3UTR","CDS","5UTR")]

dens_alba <- getFgBg(feGR = ALBAGR_R,
                     bgGR = sample(gtfGR_bg,50000),
                     gtfGR = gtfGR,
                     expr=expr,
                     dp=2)



dens_alba_ect2_strict <- getFgBg(feGR = ALBAGR[ALBAGR$gene %in% set_list$ECT2_strict],
                     bgGR = sample(gtfGR_bg,5000),
                     gtfGR = gtfGR,
                     expr=expr,
                     dp=2)

dens_alba_ect2_permissive <- getFgBg(feGR = ALBAGR[ALBAGR$gene %in% set_list$ECT2_permiss],
                                 bgGR = sample(gtfGR_bg,5000),
                                 gtfGR = gtfGR,
                                 expr=expr,
                                 dp=2)

dens_alba_not_ect2 <- getFgBg(feGR = ALBAGR[!(ALBAGR$gene %in% set_list$ECT2_permiss)],
                                     bgGR = sample(gtfGR_bg,5000),
                                     gtfGR = gtfGR,
                                     expr=expr,
                                     dp=2)

dens_alba_all <- getFgBg(feGR = ALBAGR,
                     bgGR = sample(gtfGR_bg,50000),
                     gtfGR = gtfGR,
                     expr=expr,
                     dp=2)

dens_list <- list("ALBA_strong" = dens_alba,"ALBA_all"=dens_alba_all)
dlp <- plotDensList2(dens_list,smooth = F)
dlp
ggsave(plot = dlp , filename  ="./METAGENE_ALBA_SCALED.pdf",width=6,height=4)


dens_alba_XL <- getFgBg(feGR = ALBAGR_R_STR,
                         bgGR = sample(gtfGR_bg,5000),
                         gtfGR = gtfGR,
                         expr=expr,
                         dp=2)

ECT2_transcripts <- unique(unlist(lapply(base::strsplit(iCLIP[[3]]$transcripts,","),function(x) x[1])))
dens_ect2 <- getFgBg(feGR = iCLIP[[3]],
                     bgGR = sample(gtfGR[gtfGR$transcript_id %in% ECT2_transcripts],42429),
                     gtfGR = gtfGR,
                     expr = expr,
                     scale = T,
                     dp=2)

nano_transcripts <- unique(unlist(lapply(base::strsplit(nanoGR$transcripts,","),function(x) x[1])))
dens_nano <- getFgBg(feGR = nanoGR,
                     bgGR = sample(gtfGR[gtfGR$transcript_id %in% nano_transcripts],50000),
                     gtfGR = gtfGR,
                     expr = expr,
                     dp=2)

dens_list <- list("ALBA_all" = dens_alba_all, "ALBA_strong" = dens_alba)
dlp <- plotDensList2(dens_list,smooth = FALSE)

ggsave(plot = dlp , filename ="./METAGENE_strength_score_ALBA_SCALED.pdf",width=6,height=4)

dens_list <- list("ALBA" = dens_alba , "ECT2" = dens_ect2,"m6A"=dens_nano)
dlp <- plotDensList2(dens_list,smooth = FALSE)
dlp
ggsave(plot = dlp , filename  ="./METAGENE_m6A_ECT2_ALBA_SCALED.pdf",width=6,height=4)

dens_list <- list("ALBA_not_ect2" = dens_alba_not_ect2, "ALBA_ECT2_strict" = dens_alba_ect2_strict, "ALBA_ECT2_perm" = dens_alba_ect2_permissive)
dlp <- plotDensList2(dens_list,smooth = TRUE)
dlp
ggsave(plot = dlp , filename  ="./METAGENE_ALBA_ECT2_status_SCALED.pdf",width=6,height=4)


########################################################################
#################ALBA NUMBERS#############################################
########################################################################

pgenes <- unique(nanoGR$gene)
nanoGR_BG <- returnDens(nanoGR,gtfGR[gtfGR$gene_id %in% pgenes],getBACKGROUND=T)$bg
mp_alba_m6a <- plot_ol(nanoGR,ALBAGR_R,norm=T,myn = 250)
mp_alba_m6a_MATCH <- plot_ol(nanoGR_BG,ALBAGR_R,norm=T,myn=250)

mp_ect_m6a <- plot_ol(nanoGR,iCLIP[[3]],norm=T,myn = 250)
mp_ect_m6a_MATCH <- plot_ol(nanoGR_BG,iCLIP[[3]],norm=T,myn=250)

restib <- tibble("coord"=-250:250,"m6a"=mp_alba_m6a,"m6a_match"=mp_alba_m6a_MATCH)
dlpA <- restib %>% tidyr::pivot_longer(-coord) %>%
  ggplot(aes(x=coord,y=value,group=name,colour=name)) + 
  geom_line() + theme_minimal() + geom_vline(xintercept = 0,lty=2) + 
  xlab("m6A site") + ylab("Number of ALBA sites per 1000 m6A sites") 
dlpA


restib <- tibble("coord"=-250:250,"m6a"=mp_ect_m6a,"m6a_match"=mp_ect_m6a_MATCH)
dlpB <- restib %>% tidyr::pivot_longer(-coord) %>%
  ggplot(aes(x=coord,y=value,group=name,colour=name)) + 
  geom_line() + theme_minimal() + geom_vline(xintercept = 0,lty=2) + 
  xlab("m6A site") + ylab("Number of ECT2 sites per 1000 m6A sites") 
dlpB

TP <- grid.arrange(dlpA,dlpB,ncol=2)

ggsave(TP,file="/projects/renlab/data/projects/projects_with_PB/ALBA/plots/ALBA_ECT2_around_m6A_with_BG.pdf",width=12,height=6)


ECT_BG_ST <- nanoGR[(nanoGR$gene %in% set_list$ECT2_strict)]
ECT_BG_PER <- nanoGR[(nanoGR$gene %in% set_list$ECT2_permiss)]
ECT_BG_NO <- nanoGR[!(nanoGR$gene %in% set_list$ECT2_permiss)]
ECT_BG_BG <- returnDens(nanoGR,gtfGR[(gtfGR$gene_id %in% set_list$ECT2_permiss)],getBACKGROUND=T)$bg

ECT_BG_NO <- ECT_BG_NO[sample(1:length(ECT_BG_NO),length(ECT_BG_PER),replace=T)]
ECT_BG_ST <- ECT_BG_ST[sample(1:length(ECT_BG_ST),length(ECT_BG_PER),replace=T)]
ECT_BG_BG <- ECT_BG_BG[sample(1:length(ECT_BG_BG),length(ECT_BG_PER),replace=T)]

mp_alba_m6a <- plot_ol(nanoGR,ALBAGR_R,norm=T,myn = 250)
mp_alba_m6a_ST <- plot_ol(ECT_BG_ST,ALBAGR_R,norm=T,myn=250)
mp_alba_m6a_PER <- plot_ol(ECT_BG_PER,ALBAGR_R,norm=T,myn=250)
mp_alba_m6a_NO <- plot_ol(ECT_BG_NO,ALBAGR_R,norm=T,myn=250)
mp_alba_m6a_BG <- plot_ol(ECT_BG_BG,ALBAGR_R,norm=T,myn=250)

restib <- tibble("coord"=-250:250,"ect_strict"=mp_alba_m6a_ST,"ect_perm"=mp_alba_m6a_PER,"ect_no"=mp_alba_m6a_NO,"m6a_bg"=mp_alba_m6a_BG)
dlp1 <- restib %>% tidyr::pivot_longer(-coord) %>%
  ggplot(aes(x=coord,y=value,group=name,colour=name)) + 
  geom_line() + theme_minimal() + geom_vline(xintercept = 0,lty=2) + 
  xlab("m6A site") + ylab("Number of ALBA sites per 1000 m6A sites") 
dlp1


##########################################################
#############FIG5E########################################
##########################################################

mp_ect_m6a <- plot_ol(nanoGR,iCLIP[[3]],norm=T,myn = 250)
mp_ect_m6a_ST <- plot_ol(ECT_BG_ST,iCLIP[[3]],norm=T,myn=250)
mp_ect_m6a_PER <- plot_ol(ECT_BG_PER,iCLIP[[3]],norm=T,myn=250)
mp_ect_m6a_NO <- plot_ol(ECT_BG_NO,iCLIP[[3]],norm=T,myn=250)
mp_ect_m6a_BG <- plot_ol(ECT_BG_BG,iCLIP[[3]],norm=T,myn=250)

restib <- tibble("coord"=-250:250,"ect_strict"=mp_ect_m6a_ST,"ect_perm"=mp_ect_m6a_PER,"ect_no"=mp_ect_m6a_NO,"m6a_bg"=mp_ect_m6a_BG)
dlp2 <- restib %>% tidyr::pivot_longer(-coord) %>%
  ggplot(aes(x=coord,y=value,group=name,colour=name)) + 
  geom_line() + theme_minimal() + geom_vline(xintercept = 0,lty=2) + 
  xlab("m6A site") + ylab("Number of ECT2 sites per 1000 m6A sites") 
dlp2

TP <- grid.arrange(dlp1,dlp2,ncol=2)
ggsave(TP,file="./ALBA_ECT2_around_m6A_site_groups.pdf",width=12,height=6)


##########################################################
#############FIG5E########################################
##########################################################

ECT_BG <- returnDens(iCLIP[[3]],gtfGR[(gtfGR$gene_id %in% set_list$ECT2_strict)],getBACKGROUND=T)$bg

mp_ect_alba <- plot_ol(iCLIP[[3]],ALBAGR_R,norm=T,myn=250)
mp_ect_alba_bg <- plot_ol(ECT_BG,ALBAGR_R,norm=T,myn=250)

restib <- tibble("coord"=-250:250,"ect_alba"=mp_ect_alba,"ect_bg_alba"=mp_ect_alba_bg)
dlp2 <- restib %>% tidyr::pivot_longer(-coord) %>%
  ggplot(aes(x=coord,y=value,group=name,colour=name)) + 
  geom_line() + theme_minimal() + geom_vline(xintercept = 0,lty=2) + 
  xlab("ECT2 XL site") + ylab("Number of ALBA4 sites per 1000 m6A sites") 
dlp2

ggsave(dlp2,file="./ALBA_around_ECT2_with_bg.pdf",width=6,height=6)
