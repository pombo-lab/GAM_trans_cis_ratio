#plot of change of trans-cis ratios stratified by ESC nucleolus association and domain melting status 

library(tidyverse)
library(ggpubr)
library(ggnewscale)

#set working directory to the directory code is in, to make sure that files are found
#for Rstudio:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#for command line R:
#setwd(getSrcDirectory()[1])

#load data
melting_scores <- read_tsv('../data/master_table_long_genes_withClustIdent.tsv.gz') %>% 
  dplyr::mutate(group_vta_R1 = case_when(meltingScore_DN_R1 > 5 & nad_assoc50 == TRUE ~ 'melt+nad',
                                         meltingScore_DN_R1 > 5 & nad_assoc50 == FALSE ~ 'melt+NOnad',
                                         meltingScore_DN_R1 <= 5 & nad_assoc50 == TRUE ~ 'NOmelt+nad',
                                         meltingScore_DN_R1 <= 5 & nad_assoc50 == FALSE ~ 'NOmelt+NOnad'),
                group_vta_R2 = case_when(meltingScore_DN_R2 > 5 & nad_assoc50 == TRUE ~ 'melt+nad',
                                         meltingScore_DN_R2 > 5 & nad_assoc50 == FALSE ~ 'melt+NOnad',
                                         meltingScore_DN_R2 <= 5 & nad_assoc50 == TRUE ~ 'NOmelt+nad',
                                         meltingScore_DN_R2 <= 5 & nad_assoc50 == FALSE ~ 'NOmelt+NOnad'),
                group_PGN_R1 = case_when(meltingScore_PGN_R1 > 5 & nad_assoc50 == TRUE ~ 'melt+nad',
                                         meltingScore_PGN_R1 > 5 & nad_assoc50 == FALSE ~ 'melt+NOnad',
                                         meltingScore_PGN_R1 <= 5 & nad_assoc50 == TRUE ~ 'NOmelt+nad',
                                         meltingScore_PGN_R1 <= 5 & nad_assoc50 == FALSE ~ 'NOmelt+NOnad'),
                group_PGN_R2 = case_when(meltingScore_PGN_R2 > 5 & nad_assoc50 == TRUE ~ 'melt+nad',
                                         meltingScore_PGN_R2 > 5 & nad_assoc50 == FALSE ~ 'melt+NOnad',
                                         meltingScore_PGN_R2 <= 5 & nad_assoc50 == TRUE ~ 'NOmelt+nad',
                                         meltingScore_PGN_R2 <= 5 & nad_assoc50 == FALSE ~ 'NOmelt+NOnad'),
                group_olg = case_when(meltingScore_OLGs > 5 & nad_assoc50 == TRUE ~ 'melt+nad',
                                      meltingScore_OLGs > 5 & nad_assoc50 == FALSE ~ 'melt+NOnad',
                                      meltingScore_OLGs <= 5 & nad_assoc50 == TRUE ~ 'NOmelt+nad',
                                      meltingScore_OLGs <= 5 & nad_assoc50 == FALSE ~ 'NOmelt+NOnad'))

#Figure ED 6i: correlation of trans cis ratios in replicates
melting_scores %>% ggplot()+
  geom_point(aes(x=median_of_mean_trans_cis_pgn_R1, y=median_of_mean_trans_cis_pgn_R2))+
  geom_smooth(aes(x=median_of_mean_trans_cis_pgn_R1, y=median_of_mean_trans_cis_pgn_R2),method = 'lm')
cor.test(melting_scores$median_of_mean_trans_cis_pgn_R1, melting_scores$median_of_mean_trans_cis_pgn_R2)

melting_scores %>% ggplot()+
  geom_point(aes(x=median_of_mean_trans_cis_dn_R1, y=median_of_mean_trans_cis_dn_R2))+
  geom_smooth(aes(x=median_of_mean_trans_cis_dn_R1, y=median_of_mean_trans_cis_dn_R2),method = 'lm')
cor.test(melting_scores$median_of_mean_trans_cis_dn_R1, melting_scores$median_of_mean_trans_cis_dn_R2)

#Figure ED6 l+m: violin plots of trans-cis ratio and NAD association
plot_violin <- function(celltype= c('oligo', 'vtaR1', 'vtaR2', 'ca1R1', 'ca1R2')){
  if(celltype == 'oligo'){
    color = '#800080'
    grouping = 'group_olg'
    data_melt = melting_scores %>% 
      dplyr::filter(meltingScore_OLGs > 5) %>% 
      pivot_longer(cols=c(median_of_mean_trans_cis_esc, median_of_mean_trans_cis_olig)) %>% 
      dplyr::mutate(name = factor(name, levels=c('median_of_mean_trans_cis_esc', 'median_of_mean_trans_cis_olig'))) 
    data_nonmelt = melting_scores %>% 
      dplyr::filter(meltingScore_OLGs <= 5) %>% 
      pivot_longer(cols=c(median_of_mean_trans_cis_esc, median_of_mean_trans_cis_olig)) %>% 
      dplyr::mutate(name = factor(name, levels=c('median_of_mean_trans_cis_esc', 'median_of_mean_trans_cis_olig'))) 
  }
  if(celltype == 'vtaR1'){
    color = '#259A37'
    grouping = 'group_vta_R1'
    data_melt = melting_scores %>% 
      dplyr::filter(meltingScore_DN_R1 > 5) %>% 
      pivot_longer(cols=c(median_of_mean_trans_cis_esc, median_of_mean_trans_cis_dn_R1)) %>% 
      dplyr::mutate(name = factor(name, levels=c('median_of_mean_trans_cis_esc', 'median_of_mean_trans_cis_dn_R1'))) 
    data_nonmelt = melting_scores %>% 
      dplyr::filter(meltingScore_DN_R1 <= 5) %>% 
      pivot_longer(cols=c(median_of_mean_trans_cis_esc, median_of_mean_trans_cis_dn_R1)) %>% 
      dplyr::mutate(name = factor(name, levels=c('median_of_mean_trans_cis_esc', 'median_of_mean_trans_cis_dn_R1'))) 
  }
  if(celltype == 'vtaR2'){
    color = '#259A37'
    grouping = 'group_vta_R2'
    data_melt = melting_scores %>% 
      dplyr::filter(meltingScore_DN_R2 > 5) %>% 
      pivot_longer(cols=c(median_of_mean_trans_cis_esc, median_of_mean_trans_cis_dn_R2)) %>% 
      dplyr::mutate(name = factor(name, levels=c('median_of_mean_trans_cis_esc', 'median_of_mean_trans_cis_dn_R2'))) 
    data_nonmelt = melting_scores %>% 
      dplyr::filter(meltingScore_DN_R2 <= 5) %>% 
      pivot_longer(cols=c(median_of_mean_trans_cis_esc, median_of_mean_trans_cis_dn_R2)) %>% 
      dplyr::mutate(name = factor(name, levels=c('median_of_mean_trans_cis_esc', 'median_of_mean_trans_cis_dn_R2'))) 
  }
  if(celltype == 'ca1R1'){
    color = '#6367DC'
    grouping = 'group_PGN_R1'
    data_melt = melting_scores %>% 
      dplyr::filter(meltingScore_PGN_R1 > 5) %>% 
      pivot_longer(cols=c(median_of_mean_trans_cis_esc, median_of_mean_trans_cis_pgn_R1)) %>% 
      dplyr::mutate(name = factor(name, levels=c('median_of_mean_trans_cis_esc', 'median_of_mean_trans_cis_pgn_R1'))) 
    data_nonmelt = melting_scores %>% 
      dplyr::filter(meltingScore_PGN_R1 <= 5) %>% 
      pivot_longer(cols=c(median_of_mean_trans_cis_esc, median_of_mean_trans_cis_pgn_R1)) %>% 
      dplyr::mutate(name = factor(name, levels=c('median_of_mean_trans_cis_esc', 'median_of_mean_trans_cis_pgn_R1'))) 
  }
  if(celltype == 'ca1R2'){
    color = '#6367DC'
    grouping = 'group_PGN_R2'
    data_melt = melting_scores %>% 
      dplyr::filter(meltingScore_PGN_R2 > 5) %>% 
      pivot_longer(cols=c(median_of_mean_trans_cis_esc, median_of_mean_trans_cis_pgn_R2)) %>% 
      dplyr::mutate(name = factor(name, levels=c('median_of_mean_trans_cis_esc', 'median_of_mean_trans_cis_pgn_R2'))) 
    data_nonmelt = melting_scores %>% 
      dplyr::filter(meltingScore_PGN_R2 <= 5) %>% 
      pivot_longer(cols=c(median_of_mean_trans_cis_esc, median_of_mean_trans_cis_pgn_R2)) %>% 
      dplyr::mutate(name = factor(name, levels=c('median_of_mean_trans_cis_esc', 'median_of_mean_trans_cis_pgn_R2'))) 
  }
  #plotting: 
  meltplot <- data_melt %>% ggplot()+
    geom_violin(aes(x=!!sym(grouping), y=value, color= name)) +
    geom_boxplot(aes(x=!!sym(grouping), y=value, color= name), width = 0.4, outlier.shape=NA, position=position_dodge(0.9))+
    scale_color_manual(values=c('grey50', 'grey50'))+
    new_scale_color() +
    geom_jitter(aes(x=!!sym(grouping), y=value, color= name), size = 0.8,  position=position_jitterdodge(dodge.width=0.9))+
    scale_color_manual(values=c('darkorange', color))+
    scale_y_continuous(limits=c(0,100))+
    ylab('trans-cis ratio [chrom percentile]') +
    theme(legend.position = 'none')
  meltplot
  nonmeltplot <- data_nonmelt %>% ggplot()+
    geom_violin(aes(x=!!sym(grouping), y=value, color= name)) +
    geom_boxplot(aes(x=!!sym(grouping), y=value, color= name), width = 0.4, outlier.shape=NA, position=position_dodge(0.9))+
    scale_color_manual(values=c('grey50', 'grey50'))+
    new_scale_color() +
    geom_jitter(aes(x=!!sym(grouping), y=value, color= name), size = 0.8,  position=position_jitterdodge(dodge.width=0.9))+
    scale_color_manual(values=c('darkorange', 'grey70'))+
    scale_y_continuous(limits=c(0,100))+
    theme(legend.position = 'none',axis.title.y = element_blank()) 
  nonmeltplot
  toreturn <- ggarrange(meltplot, nonmeltplot, nrow=1)
  return(toreturn)
}
plot_violin(celltype = 'ca1R2') # can be one of: c('oligo', 'vtaR1', 'vtaR2', 'ca1R1', 'ca1R2')

o1 <- wilcox.test((melting_scores %>% dplyr::filter(group_olg == 'melt+nad'))$median_of_mean_trans_cis_esc,
                  (melting_scores %>% dplyr::filter(group_olg == 'melt+nad'))$median_of_mean_trans_cis_olig)
o2 <- wilcox.test((melting_scores %>% dplyr::filter(group_olg == 'melt+NOnad'))$median_of_mean_trans_cis_esc,
                  (melting_scores %>% dplyr::filter(group_olg == 'melt+NOnad'))$median_of_mean_trans_cis_olig)
o3 <- wilcox.test((melting_scores %>% dplyr::filter(group_olg == 'NOmelt+nad'))$median_of_mean_trans_cis_esc,
                  (melting_scores %>% dplyr::filter(group_olg == 'NOmelt+nad'))$median_of_mean_trans_cis_olig)
o4 <- wilcox.test((melting_scores %>% dplyr::filter(group_olg == 'NOmelt+NOnad'))$median_of_mean_trans_cis_esc,
                  (melting_scores %>% dplyr::filter(group_olg == 'NOmelt+NOnad'))$median_of_mean_trans_cis_olig)

v11 <- wilcox.test((melting_scores %>% dplyr::filter(group_vta_R1 == 'melt+nad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_vta_R1 == 'melt+nad'))$median_of_mean_trans_cis_dn_R1)
v21 <- wilcox.test((melting_scores %>% dplyr::filter(group_vta_R1 == 'melt+NOnad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_vta_R1 == 'melt+NOnad'))$median_of_mean_trans_cis_dn_R1)
v31 <- wilcox.test((melting_scores %>% dplyr::filter(group_vta_R1 == 'NOmelt+nad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_vta_R1 == 'NOmelt+nad'))$median_of_mean_trans_cis_dn_R1)
v41 <- wilcox.test((melting_scores %>% dplyr::filter(group_vta_R1 == 'NOmelt+NOnad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_vta_R1 == 'NOmelt+NOnad'))$median_of_mean_trans_cis_dn_R1)

v12 <- wilcox.test((melting_scores %>% dplyr::filter(group_vta_R2 == 'melt+nad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_vta_R2 == 'melt+nad'))$median_of_mean_trans_cis_dn_R2)
v22 <- wilcox.test((melting_scores %>% dplyr::filter(group_vta_R2 == 'melt+NOnad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_vta_R2 == 'melt+NOnad'))$median_of_mean_trans_cis_dn_R2)
v32 <- wilcox.test((melting_scores %>% dplyr::filter(group_vta_R2 == 'NOmelt+nad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_vta_R2 == 'NOmelt+nad'))$median_of_mean_trans_cis_dn_R2)
v42 <- wilcox.test((melting_scores %>% dplyr::filter(group_vta_R2 == 'NOmelt+NOnad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_vta_R2 == 'NOmelt+NOnad'))$median_of_mean_trans_cis_dn_R2)

p11 <- wilcox.test((melting_scores %>% dplyr::filter(group_PGN_R1 == 'melt+nad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_PGN_R1 == 'melt+nad'))$median_of_mean_trans_cis_pgn_R1)
p21 <- wilcox.test((melting_scores %>% dplyr::filter(group_PGN_R1 == 'melt+NOnad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_PGN_R1 == 'melt+NOnad'))$median_of_mean_trans_cis_pgn_R1)
p31 <- wilcox.test((melting_scores %>% dplyr::filter(group_PGN_R1 == 'NOmelt+nad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_PGN_R1 == 'NOmelt+nad'))$median_of_mean_trans_cis_pgn_R1)
p41 <- wilcox.test((melting_scores %>% dplyr::filter(group_PGN_R1 == 'NOmelt+NOnad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_PGN_R1 == 'NOmelt+NOnad'))$median_of_mean_trans_cis_pgn_R1)

p12 <- wilcox.test((melting_scores %>% dplyr::filter(group_PGN_R2 == 'melt+nad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_PGN_R2 == 'melt+nad'))$median_of_mean_trans_cis_pgn_R2)
p22 <- wilcox.test((melting_scores %>% dplyr::filter(group_PGN_R2 == 'melt+NOnad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_PGN_R2 == 'melt+NOnad'))$median_of_mean_trans_cis_pgn_R2)
p32 <- wilcox.test((melting_scores %>% dplyr::filter(group_PGN_R2 == 'NOmelt+nad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_PGN_R2 == 'NOmelt+nad'))$median_of_mean_trans_cis_pgn_R2)
p42 <- wilcox.test((melting_scores %>% dplyr::filter(group_PGN_R2 == 'NOmelt+NOnad'))$median_of_mean_trans_cis_esc,
                   (melting_scores %>% dplyr::filter(group_PGN_R2 == 'NOmelt+NOnad'))$median_of_mean_trans_cis_pgn_R2)

test_results =tibble(comparison=c('melt+nad', 'melt+Nonad', 'NOmelt+nad', 'NOmelt+NOnad'),
                     olg = c(o1$p.value, o2$p.value, o3$p.value, o4$p.value),
                     dn_R1 = c(v11$p.value, v21$p.value, v31$p.value, v41$p.value),
                     dn_R2 = c(v12$p.value, v22$p.value, v32$p.value, v42$p.value),
                     pgn_R1 = c(p11$p.value, p21$p.value, p31$p.value, p41$p.value),
                     pgn_R2 = c(p12$p.value, p22$p.value, p32$p.value, p42$p.value))
