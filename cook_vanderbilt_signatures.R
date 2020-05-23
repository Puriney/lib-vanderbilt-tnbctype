suppressPackageStartupMessages(library(magrittr))
"
#--------------------------
# Goal:Create the Vanderbilt signatures in GMT format.
#
# What is GMT format? Visit https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
#--------------------------
# Author: Yun Yan (yun.yan@uth.tmc.edu)
#--------------------------
" -> doc_help
timestamp()
suppressPackageStartupMessages({
  library(tidyverse); library(readr); library(fs); library(limma)
})


setwd('.')
saveto_gmt_plus <- 'vanderbilt_tnbc_subtypes_gene_signatures_plus.gmt'
saveto_gmt_minus <- 'vanderbilt_tnbc_subtypes_gene_signatures_minus.gmt'

files_deg_direction <- list(
  'BL1' = 'BL1.deg.txt',
  'BL2' = 'BL2.deg.txt',
  'IM'  = 'IM.deg.txt',
  'LAR' = 'LAR.deg.txt',
  'M'   = 'M.deg.txt',
  'MSL' = 'MSL.deg.txt'
)
files_deg_qval <- list(
  'BL1' = 'BL1.qval.txt',
  'BL2' = 'BL2.qval.txt',
  'IM'  = 'IM.qval.txt',
  'LAR' = 'LAR.qval.txt',
  'M'   = 'M.qval.txt',
  'MSL' = 'MSL.qval.txt'
)

signature_master_names <- names(files_deg_direction)

if (file.exists(saveto_gmt_plus)){file_delete(saveto_gmt_plus)}
if (file.exists(saveto_gmt_minus)){file_delete(saveto_gmt_minus)}

all_genes_vanderbilt <- c()
for (i in seq_along(signature_master_names)){
  cat('## Processing', signature_master_names[i], '------ \n')
  deg_direction <- suppressMessages(read_tsv(files_deg_direction[[i]], col_names = F))
  deg_qval <- suppressMessages(read_tsv(files_deg_qval[[i]], col_names = F))
  deg_sig <- deg_qval %>% dplyr::filter(X3 < 0.05) %>% pull(X1)

  deg_up <- deg_direction %>% dplyr::filter(X2 == 'UP') %>% pull(X1)
  deg_dn <- deg_direction %>% dplyr::filter(X2 == 'DOWN') %>% pull(X1)
  stopifnot(length(deg_up) + length(deg_dn) == nrow(deg_direction))

  deg_up_sig <- intersect(deg_sig, deg_up)
  deg_dn_sig <- intersect(deg_sig, deg_dn)

  cat(nrow(deg_direction) - (length(deg_up_sig) + length(deg_dn_sig)),
      'genes did not pass the qval selection so excluded.\n')

  cat('UP:', length(deg_up_sig), 'genes.\n')
  cat('DN:', length(deg_dn_sig), 'genes.\n')
  if (!is_empty(intersect(deg_up_sig, deg_dn_sig))){
    cat(c('[warn] The original author must make mistakes',
              'because ', length(intersect(deg_up_sig, deg_dn_sig)),
              'genes are reported up-and down-regulated in the same time:'))
    cat(intersect(deg_up_sig, deg_dn_sig), '\n')
  }

  # C11orf2 has two individual gene names IFT46, VPS51.
  deg_up_sig <- unlist(lapply(deg_up_sig, function(g){
    G <- limma::alias2Symbol(g, species = "Hs")
    if (is_empty(G)){return(g)} ## character(0)
    return(G)
  }))
  deg_dn_sig <- unlist(lapply(deg_dn_sig, function(g){
    G <- limma::alias2Symbol(g, species = "Hs")
    if (is_empty(G)){return(g)}
    return(G)
  }))
  deg_up_sig <- unique(as.character(deg_up_sig))
  deg_dn_sig <- unique(as.character(deg_dn_sig))
  cat('UP:', length(deg_up_sig), 'genes with the official gene symbols.\n')
  cat('DN:', length(deg_dn_sig), 'genes with the official gene symbols.\n')

  text_up <- paste(sprintf('%s_plus', signature_master_names[i]),
                   sprintf('%s significantly up regulated', signature_master_names[i]),
                   paste(deg_up_sig, collapse = '\t'), sep='\t')
  text_dn <- paste(sprintf('%s_minus', signature_master_names[i]),
                   sprintf('%s significantly down regulated', signature_master_names[i]),
                   paste(deg_dn_sig, collapse = '\t'), sep='\t')
  write_lines(x = text_up, path = saveto_gmt_plus, append = T)
  write_lines(x = text_dn, path = saveto_gmt_minus, append = T)
  cat('\n')
  all_genes_vanderbilt <- c(all_genes_vanderbilt, deg_up_sig, deg_dn_sig)
}

all_genes_vanderbilt <- sort(unique(all_genes_vanderbilt))
cat(length(all_genes_vanderbilt), 'genes are involved in the vanderbilt signatures.\n')
## 2188; this is equivalent to the gene names listed in the supplement excel file
## highly indicating I am accurately repeating the author's operation.
write_lines(all_genes_vanderbilt, 'vanderbilt_genes_2188.txt')

# genes_2188 <- read_lines('genes_2188.txt')
# str(intersect(all_genes_vanderbilt, genes_2188))
# str(setdiff(all_genes_vanderbilt, genes_2188))
# str(setdiff(genes_2188, all_genes_vanderbilt))
