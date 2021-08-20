## plot the trans cis ratio of long genes in ESC and the 3 brain cell-types with pretty raincloud plots. Stratify by melting and non-melting genes 

library(tidyverse)
library(lemon)
library(ggpubr)
`%notin%` <- Negate(`%in%`)

#set working directory to the directory code is in, to make sure that files are found
#for Rstudio:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#for command line R:
#setwd(getSrcDirectory()[1])

#read in trans-cis ratio per gene
npmi_ratio_longgenes <- read_tsv('../data/trans_cis_ratio_longgenes.tsv.gz')
#read in melting scores
melting_scores <- read_tsv('../data/melting_scores.tsv.gz') %>% dplyr::select(gene_id, starts_with('melting_in'))
npmi_ratio_longgenes_melting_score <- left_join(npmi_ratio_longgenes, melting_scores)

#add jitter to the dot plot manually:
set.seed(77)
pj <- position_jitterdodge(jitter.width=0.2, seed=9,
                           jitter.height = 0,
                           dodge.width = 0.1)

plot_transcis <- function(celltype=c('OLG', 'DN_R1', 'DN_R2', 'PGN_R1', 'PGN_R2'), melting_status=TRUE){
  if(celltype == 'OLG'){
    color_brain <- if_else(isTRUE(melting_status),'#800080', 'grey70')
    melting_genes <- npmi_ratio_longgenes_melting_score %>% dplyr::filter(melting_in_OLG == melting_status) # select all melting genes
    melting_genes_long <- melting_genes %>% pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_olig)) %>% #for density plots
      dplyr::mutate(name= factor(name, levels=c('trans_cis_ratio_esc', 'trans_cis_ratio_olig')))
    melting_genes_tohighlight <- melting_genes %>% # slect genes that shall be labeled and indicated with dark grey lines
      dplyr::filter(gene_id %in% c( 'Magi2','Csmd1', 'Rbfox1', 'Lrfn5')) %>% #top_n(n=10, wt = ratiochange_oligo) #highlighted in Fig2g: Csmd1, Pde4d, Rbfox1, Dab1, Robo2 #bottom to top transitions: 'Dab1', 'Ank3' #top to bottom:  'Fam172a', 'Lrfn5'
      pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_olig)) %>% 
      dplyr::mutate(name = case_when(name == 'trans_cis_ratio_esc' ~ 'ESC', name == 'trans_cis_ratio_olig' ~ 'OLG'),
                    name = factor(name, levels=c('ESC', 'OLG')))
    melting_genes_background <- melting_genes %>% dplyr::filter(gene_id %notin% melting_genes_tohighlight$gene_id) %>%  # select all other genes
      pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_olig)) %>% 
      dplyr::mutate(name = case_when(name == 'trans_cis_ratio_esc' ~ 'ESC', name == 'trans_cis_ratio_olig' ~ 'OLG'),
                    name = factor(name, levels=c('ESC', 'OLG')))
    median_lines <- tibble(celltype=c('ESC', 'OLG'), median =c(median(melting_genes$trans_cis_ratio_esc), median(melting_genes$trans_cis_ratio_olig)),
                           color = c('darkorange', color_brain))
  }  
  if(celltype == 'DN_R1'){
    color_brain <- if_else(isTRUE(melting_status),'#259A37', 'grey70')
    melting_genes <- npmi_ratio_longgenes_melting_score %>% dplyr::filter(melting_in_DN_R1 == melting_status) # select all melting genes
    melting_genes_long <- melting_genes %>% pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_dn_R1)) %>%  #for density plots
      dplyr::mutate(name= factor(name, levels=c('trans_cis_ratio_esc', 'trans_cis_ratio_dn_R1')))
    melting_genes_tohighlight <- melting_genes %>% # slect genes that shall be labeled and indicated with dark grey lines
      dplyr::filter(gene_id %in% c('Dscam', 'Nrxn3', 'Cdk14', 'Rbfox3')) %>% #top_n(n=10, wt = ratiochange_dn_R1) #'Fgf14', 'Auts2', 
      pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_dn_R1)) %>% 
      dplyr::mutate(name = case_when(name == 'trans_cis_ratio_esc' ~ 'ESC', name == 'trans_cis_ratio_dn_R1' ~ 'DN_R1'),
                    name = factor(name, levels=c('ESC', 'DN_R1')))
    melting_genes_background <- melting_genes %>% dplyr::filter(gene_id %notin% melting_genes_tohighlight$gene_id) %>%  # select all other genes
      pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_dn_R1)) %>% 
      dplyr::mutate(name = case_when(name == 'trans_cis_ratio_esc' ~ 'ESC', name == 'trans_cis_ratio_dn_R1' ~ 'DN_R1'),
                    name = factor(name, levels=c('ESC', 'DN_R1')))
    median_lines <- tibble(celltype=c('ESC', 'DN_R1'), median =c(median(melting_genes$trans_cis_ratio_esc), median(melting_genes$trans_cis_ratio_dn_R1)),
                           color = c('darkorange', color_brain))
  }
  if(celltype == 'DN_R2'){
    color_brain <- if_else(isTRUE(melting_status),'#259A37', 'grey70')
    melting_genes <- npmi_ratio_longgenes_melting_score %>% dplyr::filter(melting_in_DN_R2 == melting_status) # select all melting genes
    melting_genes_long <- melting_genes %>% pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_dn_R2)) %>%  #for density plots
      dplyr::mutate(name= factor(name, levels=c('trans_cis_ratio_esc', 'trans_cis_ratio_dn_R2')))
    melting_genes_tohighlight <- melting_genes %>% # slect genes that shall be labeled and indicated with dark grey lines
      dplyr::filter(gene_id %in% c('Dscam', 'Nrxn3', 'Cdk14', 'Rbfox3'))  %>%  #top_n(n=10, wt = ratiochange_dn_R2) %>% 
      pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_dn_R2)) %>% 
      dplyr::mutate(name = case_when(name == 'trans_cis_ratio_esc' ~ 'ESC', name == 'trans_cis_ratio_dn_R2' ~ 'DN_R2'),
                    name = factor(name, levels=c('ESC', 'DN_R2')))
    melting_genes_background <- melting_genes %>% dplyr::filter(gene_id %notin% melting_genes_tohighlight$gene_id) %>%  # select all other genes
      pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_dn_R2)) %>% 
      dplyr::mutate(name = case_when(name == 'trans_cis_ratio_esc' ~ 'ESC', name == 'trans_cis_ratio_dn_R2' ~ 'DN_R2'),
                    name = factor(name, levels=c('ESC', 'DN_R2')))
    median_lines <- tibble(celltype=c('ESC', 'DN_R2'), median =c(median(melting_genes$trans_cis_ratio_esc), median(melting_genes$trans_cis_ratio_dn_R2)),
                           color = c('darkorange', color_brain))
  }
  if(celltype == 'PGN_R1'){
    color_brain <- if_else(isTRUE(melting_status),'#6367DC', 'grey70')
    melting_genes <- npmi_ratio_longgenes_melting_score %>% dplyr::filter(melting_in_PGN_R1 == melting_status) # select all melting genes
    melting_genes_long <- melting_genes %>% pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_pgn_R1)) %>%  #for density plots
      dplyr::mutate(name= factor(name, levels=c('trans_cis_ratio_esc', 'trans_cis_ratio_pgn_R1')))
    melting_genes_tohighlight <- melting_genes %>% # slect genes that shall be labeled and indicated with dark grey lines
      dplyr::filter(gene_id %in% c('Rbfox1', 'Grik2', 'Fggy', 'Lingo2')) %>% #top_n(n=10, wt = ratiochange_pgn_R1) %>% #'Fgf14','Asic2',
      pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_pgn_R1)) %>% 
      dplyr::mutate(name = case_when(name == 'trans_cis_ratio_esc' ~ 'ESC', name == 'trans_cis_ratio_pgn_R1' ~ 'PGN_R1'),
                    name = factor(name, levels=c('ESC', 'PGN_R1')))
    melting_genes_background <- melting_genes %>% dplyr::filter(gene_id %notin% melting_genes_tohighlight$gene_id) %>%  # select all other genes
      pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_pgn_R1)) %>% 
      dplyr::mutate(name = case_when(name == 'trans_cis_ratio_esc' ~ 'ESC', name == 'trans_cis_ratio_pgn_R1' ~ 'PGN_R1'),
                    name = factor(name, levels=c('ESC', 'PGN_R1')))
    median_lines <- tibble(celltype=c('ESC', 'PGN_R1'), median =c(median(melting_genes$trans_cis_ratio_esc), median(melting_genes$trans_cis_ratio_pgn_R1)),
                           color = c('darkorange', color_brain))
  }
  if(celltype == 'PGN_R2'){
    color_brain <- if_else(isTRUE(melting_status),'#6367DC', 'grey70')
    melting_genes <- npmi_ratio_longgenes_melting_score %>% dplyr::filter(melting_in_PGN_R1 == melting_status) # select all melting genes
    melting_genes_long <- melting_genes %>% pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_pgn_R2)) %>%  #for density plots
      dplyr::mutate(name= factor(name, levels=c('trans_cis_ratio_esc', 'trans_cis_ratio_pgn_R2')))
    melting_genes_tohighlight <- melting_genes %>% # slect genes that shall be labeled and indicated with dark grey lines
      dplyr::filter(gene_id %in% c('Rbfox1', 'Grik2', 'Fggy', 'Lingo2')) %>% #top_n(n=10, wt = ratiochange_pgn_R2) %>% 
      pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_pgn_R2)) %>% 
      dplyr::mutate(name = case_when(name == 'trans_cis_ratio_esc' ~ 'ESC', name == 'trans_cis_ratio_pgn_R2' ~ 'PGN_R2'),
                    name = factor(name, levels=c('ESC', 'PGN_R2')))
    melting_genes_background <- melting_genes %>% dplyr::filter(gene_id %notin% melting_genes_tohighlight$gene_id) %>%  # select all other genes
      pivot_longer(cols=c(trans_cis_ratio_esc, trans_cis_ratio_pgn_R2)) %>% 
      dplyr::mutate(name = case_when(name == 'trans_cis_ratio_esc' ~ 'ESC', name == 'trans_cis_ratio_pgn_R2' ~ 'PGN_R2'),
                    name = factor(name, levels=c('ESC', 'PGN_R2')))
    median_lines <- tibble(celltype=c('ESC', 'PGN_R2'), median =c(median(melting_genes$trans_cis_ratio_esc), median(melting_genes$trans_cis_ratio_pgn_R2)),
                           color = c('darkorange', color_brain))
  }
  ###PLOTTING:
  a<- ggplot() +
    lemon::geom_pointpath(data= melting_genes_background, aes(x=name, y=value, group=interaction(gene_id), colour=name), linecolor = 'grey90', distance =10, alpha=0.7, size=3, linesize=0.4, 
                          position = pj)+
    lemon::geom_pointpath(data= melting_genes_tohighlight, aes(x=name, y=value, group=interaction(gene_id), colour=name), linecolor = 'grey75', size=3 , linesize=0.4, #distance = 10, 
                          position = pj)+
    scale_color_manual(values=median_lines$color)+
    ylim(c(0,100))+
    ylab('trans-cis ratio [chromosome percentile]')+
    scale_x_discrete(labels=c(median_lines$celltype), expand = c(0.1, 0.1))
#uncomment for labeling of individual genes. But then the function only works for melting genes  
  #    geom_label_repel(data= melting_genes_tohighlight %>% dplyr::filter(name=='ESC'),aes(x=name, y=value, label=gene_id),  nudge_x = -0.2, segment.size = 0) + 
  #    geom_label_repel(data= melting_genes_tohighlight %>% dplyr::filter(name!='ESC'),aes(x=name, y=value, label=gene_id),  nudge_x = 0.2, segment.size = 0) 
  #a
  b <- ggplot() +
    geom_density(data= melting_genes_long , aes(x=value, fill=name), alpha=0.7)+
    xlim(c(0,100))+
    scale_y_continuous(breaks =c(0,0.01))+
    scale_fill_manual(values=median_lines$color)+
    geom_vline(xintercept = median_lines$median, linetype='twodash', size=2, color=median_lines$color)+
    coord_flip()
  #  b
  ab_comb <- ggarrange(a+theme(legend.position = 'none', axis.title=element_blank(), axis.text = element_blank(), plot.margin = unit(c(0,4,0,0), 'lines')),
                       b+theme(axis.text.y = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), legend.position = 'none'), # plot.margin = unit(c(0,8,0,0), 'lines')
                       ncol=2, widths = c(1,0.2), align='h')
  # ab_comb
  return(ab_comb)
}  
plot_transcis(celltype = 'OLG', melting_status = TRUE)
plot_transcis(celltype = 'OLG', melting_status = FALSE)
plot_transcis('DN_R1', melting_status = TRUE)
plot_transcis('DN_R1', melting_status = FALSE)
plot_transcis('DN_R2', melting_status = TRUE)
plot_transcis('DN_R2', melting_status = FALSE)
plot_transcis('PGN_R1', melting_status = TRUE)
plot_transcis('PGN_R1', melting_status = FALSE)
plot_transcis('PGN_R2', melting_status = TRUE)
plot_transcis('PGN_R2', melting_status = FALSE)

