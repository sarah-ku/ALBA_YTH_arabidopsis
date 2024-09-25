load(file = "design_with_results.Rdat")
ect2_inter <- design_full$posGR_filtered_NO_ADAR[[4]]
filt <- ect2_inter$padj<0.1 & (ect2_inter$prop>0.01 | ect2_inter$prop_ctrl>0.01) & ect2_inter$tags_treat>0 & ect2_inter$tags_control>0
summary(filt)
ect2_inter <- ect2_inter[filt]

plot(ect2_inter$prop,ect2_inter$prop_ctrl)
abline(0,1)

load(file = "design_with_results_TWO_SAMPLE.Rdat")
two_samp_ECT2 <- design$posGR_filtered_NO_ADAR[[4]]
filt <- two_samp_ECT2$padj<0.1 & (two_samp_ECT2$prop>0.01 | two_samp_ECT2$prop_ctrl>0.01) & two_samp_ECT2$tags_treat>0 & two_samp_ECT2$tags_control>0
summary(filt)
two_samp_ECT2 <- two_samp_ECT2[filt]


ect2_int <- read_csv2(file="ECT2_INTERACTION_supplementary_csv.csv")
length(ect2_int$seqnames)
ect2_genes <- (unique(ect2_int$gene))
length(ect2_genes[!(ect2_genes=="")])
inla_pos <- as.data.frame(ect2_int)[,1]
inla_genes <- as.data.frame(ect2_int$name)[,1]

a <- makeVenn(list("all"=names(ect2_inter),"inla"=inla_pos,"two_samp"=names(two_samp_ECT2)))
b <- makeVenn(list("all"=unique(ect2_inter$name),"inla"=unique(inla_genes),"two_samp"=unique(two_samp_ECT2$name)))
TS <- grid.arrange(a,b)
ggsave(filename = "INLA_VS_ALL_VS_TWO.pdf",plot = TS,height=8,width=3)
