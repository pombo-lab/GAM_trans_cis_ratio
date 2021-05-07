## calculation of GAM NPMI trans-cis ratio per 250 kb bin

library(tidyverse)


#read in list of long genes (GCF_000001635.26_GRCm38.p6_genomic.gtf)
gene_list <- read_tsv('../data/long_genes.tsv.gz')

#need 13 GB of memory to read in the file. Currently not possible to store data on github

#npmi_long_fourcelltypes <- read_tsv('../output/npmi_esc_dn_pgn_olig_withRepl_wdfcleaned_long.tsv') %>% 
#  dplyr::mutate(A_chrom = factor(A_chrom, levels=str_sort(unique(A_chrom), numeric = TRUE)))

#per bin: summarise all NPMI values in cis and trans with mean(). Tried median() as well and results are very similar
#npmi_trans_cis_perbin <- npmi_long_fourcelltypes %>% 
#  group_by(A_chrom, A_start, contact_type) %>% 
#  summarise(mean_npmi_esc=mean(npmi_esc, na.rm=T),
#            mean_npmi_pgn_R1=mean(npmi_pgn_R1, na.rm=T),
#            mean_npmi_pgn_R2=mean(npmi_pgn_R2, na.rm=T),
#            mean_npmi_dn_R1=mean(npmi_dn_R1, na.rm=T),
#            mean_npmi_dn_R2=mean(npmi_dn_R2, na.rm=T),
#            mean_npmi_olig=mean(npmi_olig, na.rm=T))

#read in mean trans and cis npmi values per bin
npmi_trans_cis_perbin <- read_tsv('../data/npmi_trans_cis_per_bin.tsv.gz') %>% 
  dplyr::mutate(A_chrom = factor(A_chrom, levels=str_sort(unique(A_chrom), numeric = TRUE)))

#plot trans and cis contacts per chromosome in the different cell types
npmi_trans_cis_perbin %>% 
  pivot_longer(cols = c(mean_npmi_esc, mean_npmi_pgn_R1, mean_npmi_pgn_R2, mean_npmi_dn_R1, mean_npmi_dn_R2, mean_npmi_olig), names_prefix='mean_npmi_', names_to='cell_type', values_to='npmi') %>% 
  dplyr::mutate(cell_type=factor(cell_type, levels=c('esc', 'dn_R1', 'dn_R2', 'olig','pgn_R1', 'pgn_R2'))) %>% 
  ggplot() + 
  geom_boxplot(aes(x= contact_type, y=npmi, color=A_chrom), position='dodge') + 
  facet_wrap(~cell_type)+
  scale_color_viridis_d() 

##calculate the trans-cis ratio 
npmi_trans_cis_perbin_ratio <- npmi_trans_cis_perbin %>% 
  pivot_wider(names_from=contact_type, values_from=c(mean_npmi_esc, mean_npmi_pgn_R1, mean_npmi_pgn_R2, mean_npmi_dn_R1, mean_npmi_dn_R2, mean_npmi_olig)) %>% 
  dplyr::mutate(trans_cis_ratio_mean_esc = mean_npmi_esc_trans/mean_npmi_esc_cis,
                trans_cis_ratio_mean_pgn_R1 = mean_npmi_pgn_R1_trans/mean_npmi_pgn_R1_cis,
                trans_cis_ratio_mean_pgn_R2 = mean_npmi_pgn_R1_trans/mean_npmi_pgn_R2_cis,
                trans_cis_ratio_mean_dn_R1 = mean_npmi_dn_R1_trans/mean_npmi_dn_R1_cis,
                trans_cis_ratio_mean_dn_R2 = mean_npmi_dn_R2_trans/mean_npmi_dn_R2_cis,
                trans_cis_ratio_mean_olig = mean_npmi_olig_trans/mean_npmi_olig_cis) %>% 
  group_by(A_chrom) %>% 
  dplyr::mutate(trans_cis_ratio_mean_esc_chrompercentile = ntile(trans_cis_ratio_mean_esc, 100),
                trans_cis_ratio_mean_pgn_R1_chrompercentile = ntile(trans_cis_ratio_mean_pgn_R1, 100),
                trans_cis_ratio_mean_pgn_R2_chrompercentile = ntile(trans_cis_ratio_mean_pgn_R2, 100),
                trans_cis_ratio_mean_dn_R1_chrompercentile = ntile(trans_cis_ratio_mean_dn_R1, 100),
                trans_cis_ratio_mean_dn_R2_chrompercentile = ntile(trans_cis_ratio_mean_dn_R2, 100),
                trans_cis_ratio_mean_olig_chrompercentile = ntile(trans_cis_ratio_mean_olig, 100))

#prepare tibble that will be used in the for loop
npmi_ratio_longgenes <- tibble(gene_id=NA,
                               trans_cis_ratio_esc=NA,
                               trans_cis_ratio_pgn_R1=NA, trans_cis_ratio_pgn_R2=NA,
                               trans_cis_ratio_dn_R1=NA, trans_cis_ratio_dn_R2=NA,
                               trans_cis_ratio_olig=NA)

#calculate the median() trans-cis ratio of long genes. Tried mean() as well and results are very similar. 
for(i in gene_list$gene_id){
  #where exactly is the gene of interest:
  query_entry <- gene_list %>% dplyr::filter(gene_id ==i)
  query_chrom <- as.character(query_entry$chrom)
  query_startbin <- query_entry$start_bin
  query_endbin <- query_entry$end_bin
  
  #from the trans-cis table, select the rows that match the coordinates of the gene
  trans_cis_ratios_entry <- npmi_trans_cis_perbin_ratio %>% dplyr::filter(A_chrom == query_chrom,
                                                                          A_start >= query_startbin,
                                                                          A_start < query_endbin)
  
  #calculate the median() trans-cis ratio of all the gene bins in each of the 6 samples
  medians_of_mean <- matrixStats::colMedians(as.matrix(trans_cis_ratios_entry %>% ungroup() %>% dplyr::select(starts_with('trans_cis_ratio_mean') & ends_with('chrompercentile'))))
  
  #if statement makes sure that the order of the samples is correct before writing results to output tibble
  if(paste(colnames(trans_cis_ratios_entry %>% ungroup() %>% dplyr::select(starts_with('trans_cis_ratio_mean') & ends_with('chrompercentile'))), collapse = " ") == 
     'trans_cis_ratio_mean_esc_chrompercentile trans_cis_ratio_mean_pgn_R1_chrompercentile trans_cis_ratio_mean_pgn_R2_chrompercentile trans_cis_ratio_mean_dn_R1_chrompercentile trans_cis_ratio_mean_dn_R2_chrompercentile trans_cis_ratio_mean_olig_chrompercentile'){
    
    #write median trans-cis ratio per gene to tibble
    npmi_ratio_longgenes <- add_row(npmi_ratio_longgenes, 
                                    gene_id=i,
                                    trans_cis_ratio_esc=medians_of_mean[1],
                                    trans_cis_ratio_pgn_R1=medians_of_mean[2], trans_cis_ratio_pgn_R2=medians_of_mean[3], 
                                    trans_cis_ratio_dn_R1=medians_of_mean[4], trans_cis_ratio_dn_R2=medians_of_mean[5], 
                                    trans_cis_ratio_olig=medians_of_mean[6])
  }
}

#write_tsv(npmi_ratio_longgenes %>% drop_na(), file = '../data/cis_trans_ratio_longgenes.tsv.gz')

#plot correlation of replicates:
npmi_ratio_longgenes %>% 
  ggplot()+
  geom_point(aes(x=trans_cis_ratio_dn_R1, y=trans_cis_ratio_dn_R2))

npmi_ratio_longgenes %>% 
  ggplot()+
  geom_point(aes(x=trans_cis_ratio_pgn_R1, y=trans_cis_ratio_pgn_R2))

