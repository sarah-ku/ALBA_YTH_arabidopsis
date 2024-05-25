#######################################################################
########################FIGS9########################################
#######################################################################


  
#######################################################################
########################FIGS9 A########################################
#######################################################################

ALBA4_RANDOM <- getSetMatchedExpression(ALBA4,log(expr+1),n_breaks = 10)
ALBA4_HT_RANDOM <- getSetMatchedExpression(ALBA4_HT,log(expr+1),n_breaks = 10)

mv_real <- makeVenn(list("ALBA_iCLIP"=ALBA4,"ALBA_HT"=ALBA4_HT))
mv_iclip <- makeVenn(list("ALBA_iCLIP"=ALBA4,"ALBA_HT_RANDOM"=ALBA4_HT_RANDOM))
mv_ht <- makeVenn(list("ALBA_iCIP_RANDOM"=ALBA4_RANDOM,"ALBA_HT"=ALBA4_HT))
venns_tp <- grid.arrange(grobs=list(mv_real,mv_iclip,mv_ht))

ggsave(plot=venns_tp,filename="ALBA4_iCLIP_ALBA4_TRIBE_overlaps_with_random_sets_filtered.pdf",width=5,height=12)

#######################################################################
########################FIGS9 B########################################
#######################################################################

ALBA4_RANDOM <- getSetMatchedExpression(ALBA4,log(expr+1),n_breaks = 10)
ECT2_RANDOM <- getSetMatchedExpression(ECT2,log(expr+1),n_breaks = 10)

mv_real <- makeVenn(list("ALBA_iCLIP"=ALBA4,"ECT2_iCLIP"=ECT2))
mv_alba_iclip <- makeVenn(list("ALBA_iCLIP"=ALBA4,"ECT2_RANDOM"=ECT2_RANDOM))
mv_ect2_iclip <- makeVenn(list("ALBA_iCIP_RANDOM"=ALBA4_RANDOM,"ECT2_iCLIP"=ECT2))
venns_tp <- grid.arrange(grobs=list(mv_real,mv_alba_iclip,mv_ect2_iclip))

ggsave(plot=venns_tp,filename="ALBA4_iCLIP_ECT2_iCLIP_overlaps_with_random_sets.pdf",width=5,height=12)

#######################################################################
########################FIGS9 C########################################
#######################################################################

ALBA4_strict <- as.data.frame(genes_ALBA)[,2][genes_ALBA$ALBA4_strict]
ECT2_strict <- as.data.frame(genes_ALBA)[,2][genes_ALBA$ECT2_strict]

ALBA4_strict_RANDOM <- getSetMatchedExpression(ALBA4_strict,log(expr+1),n_breaks = 10)
ECT2_strict_RANDOM <- getSetMatchedExpression(ECT2_strict,log(expr+1),n_breaks = 10)

mv_real <- makeVenn(list("ALBA4_strict"=ALBA4_strict,"ECT2_strict"=ECT2_strict))
mv_alba_strict <- makeVenn(list("ALBA_strict"=ALBA4_strict,"ECT2_strict_RANDOM"=ECT2_strict_RANDOM))
mv_ect2_strict <- makeVenn(list("ALBA_strict_RANDOM"=ALBA4_strict_RANDOM,"ECT2_strict"=ECT2_strict))
venns_tp <- grid.arrange(grobs=list(mv_real,mv_alba_strict,mv_ect2_strict))

ggsave(plot=venns_tp,filename="ALBA4_strict_ECT2_strict_overlaps_with_random_sets.pdf",width=5,height=12)

#######################################################################
########################FIGS9 D########################################
#######################################################################

ALBA2_permissive <- as.data.frame(genes_ALBA)[,2][genes_ALBA$ALBA2_permissive]
ECT2_ECT3_shoots <- as.data.frame(genes_ALBA)[,2][genes_ALBA$ECT23_shoots]

ALBA2_RANDOM <- getSetMatchedExpression(ALBA2_permissive,log(expr_ALBA2+1),n_breaks = 10)
ECT23_shoots_RANDOM <- getSetMatchedExpression(ECT2_ECT3_shoots,log(expr_ALBA2+1),n_breaks = 10)

mv_real <- makeVenn(list("ALBA2"=ALBA2_permissive,"ECT23"=ECT2_ECT3_shoots))
mv_alba2_random <- makeVenn(list("ALBA2"=ALBA2_permissive,"ECT23_RANDOM"=ECT23_shoots_RANDOM))
mv_ect23_random <- makeVenn(list("ALBA2_RANDOM"=ALBA2_RANDOM,"ECT23"=ECT2_ECT3_shoots))
venns_tp <- grid.arrange(grobs=list(mv_real,mv_alba2_random,mv_ect23_random))

ggsave(plot=venns_tp,filename="ALBA2_ECT23_shoots_overlaps_with_random_sets.pdf",width=5,height=12)