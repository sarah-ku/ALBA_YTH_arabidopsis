load(file="./m6A_site_network_set.Rdat")
load(file="./parker19_micGR.Rdat")
load(file="./parker19_nanoGR.Rdat")
load(file="./wang24_sacsGR.Rdat")

ol <- plot_ol(m6AGR[m6AGR$experiment=="micGR"],sacsGR)
ol <- plot_ol(m6AGR[m6AGR$experiment=="nanoGR"],sacsGR)

ol_mic <- plot_ol(micGR,sacsGR,norm = T)
ol_mic_corr <- plot_ol(m6AGR[m6AGR$experiment=="micGR"],sacsGR,norm = T)
ol_nano <- plot_ol(m6AGR[m6AGR$experiment=="nanoGR"],sacsGR,norm = T)

TS <- tibble(ind=-100:100,ol_mic,ol_mic_corr,ol_nano) %>%
  pivot_longer(-ind) %>%
  ggplot(aes(x=ind,y=value)) + 
  geom_line() + 
  facet_grid(~name) + 
  theme_bw()

ggsave(TS,file="./SAC_seq_support_MIC_corrected_m6A_sites.pdf",width=8,height=4)
