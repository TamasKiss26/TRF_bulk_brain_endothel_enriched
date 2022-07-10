library(tximport)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(plotly)
library(DESeq2)
library(org.Mm.eg.db)
library(gprofiler2)
library(ggpubr)
library(fgsea)
library(msigdbr)
library(GSVA)
library(mouse430a2.db)


# https://dnatech.genomecenter.ucdavis.edu/faqs/when-should-i-trim-my-illumina-reads-and-how-should-i-do-it/

## transcript of gene database
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
k <- AnnotationDbi::keys(x = txdb, keytype = 'TXNAME')
tx2gene <- AnnotationDbi::select(x = txdb, keys = k, columns = 'GENEID', keytype = 'TXNAME')
tx2gene <- tx2gene[!base::is.na(tx2gene$GENEID),]


###---
sample_selected <- c(1, 3:5, 7:13, 15:19)
my_threshold <- 1.5
##----


## read files
abundance_files <- base::paste0('./data/quant/S', sample_selected, '/abundance.tsv')
base::names(abundance_files) <- base::paste0('sample', sample_selected)
txi.kallisto.tsv <- tximport::tximport(
  abundance_files, 
  type = 'kallisto', 
  tx2gene = tx2gene, 
  ignoreAfterBar = TRUE,
  countsFromAbundance = 'lengthScaledTPM'
  )

sample_info <- base::data.frame(
  sample = base::paste0('sample', sample_selected),
  group = c(
    base::rep('young', 7),
    base::rep('aged', 6),
    base::rep('agedTRF', 6)
  )[sample_selected]
)
base::rownames(sample_info) <- sample_info$sample


## PCA
pca.res <- stats::prcomp(x = txi.kallisto.tsv$counts)

pca_plot <- ggpubr::ggarrange(
  plotlist = base::list(
    plot_1_2 = pca.res$rotation %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column('sample') %>%
      dplyr::left_join(., sample_info, by = 'sample') %>%
      ggplot2::ggplot(data =., mapping = ggplot2::aes(x = PC1, y = PC2, label = sample, color = group)) +
      geom_point() +
      ggplot2::geom_text(nudge_y = -0.02) +
      ggplot2::theme(legend.position = 'none'),
    plot_1_3 = pca.res$rotation %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column('sample') %>%
      dplyr::left_join(., sample_info, by = 'sample') %>%
      ggplot2::ggplot(data =., mapping = ggplot2::aes(x = PC1, y = PC3, label = sample, color = group)) +
      geom_point() +
      ggplot2::geom_text(nudge_y = -0.02) +
      ggplot2::theme(legend.position = 'none'),
    plot_2_3 = pca.res$rotation %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column('sample') %>%
      dplyr::left_join(., sample_info, by = 'sample') %>%
      ggplot2::ggplot(data =., mapping = ggplot2::aes(x = PC2, y = PC3, label = sample, color = group)) +
      geom_point() +
      ggplot2::geom_text(nudge_y = -0.02) +
      ggplot2::theme(legend.position = c(.8, .2))
  ),
  nrow = 1
)

pca_plot
  
ggplot2::ggsave(
  filename = 'pca_plot.tiff',
  plot = pca_plot,
  device = 'tiff',
  width = 15,
  height = 5,
  dpi ='retina',
  bg = 'white'
)


## heat map
detected_genes <- txi.kallisto.tsv$counts %>% 
  base::as.data.frame() %>%
  dplyr::mutate(sample_all = base::rowSums(.)) %>%
  dplyr::filter(sample_all != 0) %>%
  base::rownames()

hm_all_genes <- txi.kallisto.tsv$counts[detected_genes,] %>% 
  pheatmap::pheatmap(
    show_rownames = F,
    show_colnames = F,
    cluster_cols = T,
    scale = 'row',
    annotation_col = dplyr::select(sample_info, 2)
  )

ggplot2::ggsave(
  filename = 'hm_all_genes.tiff',
  plot = hm_all_genes,
  device = 'tiff',
  width = 5,
  height = 5,
  dpi ='retina',
  bg = 'white'
)


## DESeq2
dds <- DESeq2::DESeqDataSetFromTximport(
  txi = txi.kallisto.tsv,
  colData = sample_info,
  design= ~ group
)
dds <- DESeq2::DESeq(dds)

de_aging <- DESeq2::results(object = dds, contrast = c('group', 'aged', 'young'))
de_trf <- DESeq2::results(object = dds, contrast = c('group', 'agedTRF', 'aged'))


## combine DESeq2 tables
de_aging <- de_aging %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column('EntrezID') %>%
  dplyr::filter(baseMean != 0)

de_trf <- de_trf  %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column('EntrezID') %>%
  dplyr::filter(baseMean != 0)

de_comb <- dplyr::inner_join(x = de_aging, y = de_trf, by = 'EntrezID', suffix = c('_aging', '_trf'))
de_comb <- dplyr::inner_join(x = tibble::rownames_to_column(base::as.data.frame(txi.kallisto.tsv$counts), var = 'EntrezID'), y = de_comb, by = 'EntrezID')


## definition of de genes
de_comb <- de_comb %>% dplyr::mutate(
  'DE_aging' = base::ifelse(
    test = (pvalue_aging < .05 & (log2FoldChange_aging > base::log(my_threshold,2) | log2FoldChange_aging < -base::log(my_threshold,2))),
    yes = 'DE',
    no = 'notDE'
  ),
  'DE_trf' = base::ifelse(
    test = (pvalue_trf < .05 & (log2FoldChange_trf > base::log(my_threshold,2) | log2FoldChange_trf < -base::log(my_threshold,2))),
    yes = 'DE',
    no = 'notDE'
  ),
  'DE_discordant' = base::ifelse(
    test = (
      (DE_aging == 'DE' | DE_trf == 'DE') &
        (
          (log2FoldChange_aging > base::log(my_threshold,2) & log2FoldChange_trf < -base::log(my_threshold,2)) |
            (log2FoldChange_aging < -base::log(my_threshold,2) & log2FoldChange_trf > base::log(my_threshold,2))
        )
      ),
    yes = 'Discordant',
    no = 'notDiscordant'
  )
) 


## save data
readr::write_csv(x = de_comb, file = './data/DE_results.csv')


## heat map de genes
vector_de_genes <- de_comb %>% dplyr::filter((DE_trf == 'DE' | DE_aging == 'DE')) %>% dplyr::pull(EntrezID)
hm_de_genes <- txi.kallisto.tsv$counts[vector_de_genes,] %>% 
  pheatmap::pheatmap(
    show_rownames = F,
    show_colnames = F,
    cluster_cols = T,
    scale = 'row',
    annotation_col = dplyr::select(sample_info, 2)
  )

ggplot2::ggsave(
  filename = 'hm_de_genes.tiff',
  plot = hm_de_genes,
  device = 'tiff',
  width = 5,
  height = 5,
  dpi ='retina',
  bg = 'white'
)


## venn diagram 
VennDiagram::venn.diagram(
  x = base::list(
    'aging effect' = de_comb %>% dplyr::filter(DE_aging == 'DE') %>% dplyr::pull(EntrezID),
    'TRF effect' = de_comb %>% dplyr::filter(DE_trf == 'DE') %>% dplyr::pull(EntrezID),
    'discordant genes' = de_comb %>% dplyr::filter(DE_discordant == 'Discordant') %>% dplyr::pull(EntrezID)
  ), 
  filename = 'de_venn.tiff'
)


## save de gene info
de_aging_genes <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = de_comb %>% dplyr::filter(DE_aging == 'DE') %>% dplyr::pull(EntrezID),
  columns = c('ENTREZID', 'SYMBOL', 'GENENAME'),
  keytype = 'ENTREZID'
)
readr::write_csv(x = de_aging_genes, file = './data/de_genes/de_aging_genes.csv')

de_trf_genes <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = de_comb %>% dplyr::filter(DE_trf == 'DE') %>% dplyr::pull(EntrezID),
  columns = c('ENTREZID', 'SYMBOL', 'GENENAME'),
  keytype = 'ENTREZID'
)
readr::write_csv(x = de_trf_genes, file = './data/de_genes/de_trf_genes.csv')

de_discordant_genes <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = de_comb  %>% dplyr::filter(DE_discordant == 'Discordant') %>% dplyr::pull(EntrezID),
  columns = c('ENTREZID', 'SYMBOL', 'GENENAME'),
  keytype = 'ENTREZID'
)
readr::write_csv(x = de_discordant_genes, file = './data/de_genes/de_discordant_genes.csv')


## save de gene annotation
de_aging_annot <- gprofiler2::gost(
  query =  de_comb %>% dplyr::filter(DE_aging == 'DE') %>% dplyr::pull(EntrezID),
  organism = 'mmusculus',
  ordered_query = T,
  evcodes = T,
  sources = 'GO'
)
de_aging_annot <- de_aging_annot$result
readr::write_csv(x = de_aging_annot, file = './data/de_genes/de_aging_annot.csv')

de_trf_annot <- gprofiler2::gost(
  query =  de_comb %>% dplyr::filter(DE_trf == 'DE') %>% dplyr::pull(EntrezID),
  organism = 'mmusculus',
  ordered_query = T,
  evcodes = T,
  sources = 'GO'
)
de_trf_annot <- de_trf_annot$result
readr::write_csv(x = de_trf_annot, file = './data/de_genes/de_trf_annot.csv')

de_discordant_annot <- gprofiler2::gost(
  query =  de_comb %>% dplyr::filter(DE_discordant == 'Discordant') %>% dplyr::pull(EntrezID),
  organism = 'mmusculus',
  ordered_query = T,
  evcodes = T,
  sources = 'GO'
)
de_discordant_annot <- de_discordant_annot$result
readr::write_csv(x = de_discordant_annot, file = './data/de_genes/de_discordant_annot.csv')


## discordant genes
sp_discordant <- de_comb %>%
  dplyr::filter(DE_aging == 'DE' | DE_trf == 'DE') %>%
  dplyr::mutate(
    DE_both = dplyr::case_when(
      DE_aging == 'DE'    & DE_trf == 'notDE' ~ 'aging DE',
      DE_aging == 'notDE' & DE_trf == 'DE'    ~ 'trf DE',
      DE_aging == 'DE'    & DE_trf == 'DE'    ~ 'both DE'
    ),
    DE_discordant = forcats::fct_rev(base::as.factor(DE_discordant))
  ) %>%
  ggplot2::ggplot(data = ., mapping = ggplot2::aes(x = log2FoldChange_aging, y = log2FoldChange_trf, alpha = DE_discordant, color = DE_both)) +
  ggplot2::geom_point() +
  ggplot2::geom_hline(yintercept = 0) +
  ggplot2::geom_vline(xintercept = 0) +
  ggplot2::geom_hline(yintercept = c(base::log(my_threshold,2), -base::log(my_threshold,2)), color = 'grey') +
  ggplot2::geom_vline(xintercept = c(base::log(my_threshold,2), -base::log(my_threshold,2)), color = 'grey') +
  ggplot2::xlim(-6, 6) +
  ggplot2::ylim(-6, 6) +
  ggplot2::theme(
    legend.title = ggplot2::element_blank()
  ) +
  ggplot2::coord_fixed(ratio = 1)

sp_discordant

ggplot2::ggsave(
  filename = 'sp_discordant.tiff',
  plot = sp_discordant,
  device = 'tiff',
  width = 5,
  height = 5,
  dpi ='retina',
  bg = 'white'
)


## GSEA
gsdb <- msigdbr::msigdbr(
  species = 'Mus musculus'
)

# autophagy
mysets_names_autophagy <- c(
  'GOBP_AUTOPHAGY_OF_MITOCHONDRION', 
  'GOBP_AUTOPHAGY_OF_NUCLEUS',
  'GOBP_CHAPERONE_MEDIATED_AUTOPHAGY',
  'GOBP_LYSOSOMAL_MICROAUTOPHAGY',
  'GOBP_MACROAUTOPHAGY',
  'GOBP_NEGATIVE_REGULATION_OF_AUTOPHAGY',
  'GOBP_NEGATIVE_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION',
  'GOBP_NEGATIVE_REGULATION_OF_MACROAUTOPHAGY',
  'GOBP_PIECEMEAL_MICROAUTOPHAGY_OF_THE_NUCLEUS',
  'GOBP_POSITIVE_REGULATION_OF_AUTOPHAGY',
  'GOBP_POSITIVE_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION',
  'GOBP_POSITIVE_REGULATION_OF_MACROAUTOPHAGY',
  'GOBP_REGULATION_OF_AUTOPHAGY',
  'GOBP_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION',
  'GOBP_REGULATION_OF_CHAPERONE_MEDIATED_AUTOPHAGY',
  'GOBP_REGULATION_OF_MACROAUTOPHAGY',
  'GOBP_SELECTIVE_AUTOPHAGY',
  'KEGG_REGULATION_OF_AUTOPHAGY',
  'REACTOME_AUTOPHAGY',
  'REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY',
  'REACTOME_LATE_ENDOSOMAL_MICROAUTOPHAGY',
  'REACTOME_SELECTIVE_AUTOPHAGY',
  'GOBP_MITOPHAGY',
  'GOBP_NEGATIVE_REGULATION_OF_MITOPHAGY',
  'GOBP_POSITIVE_REGULATION_OF_MITOPHAGY',
  'GOBP_REGULATION_OF_MITOPHAGY',
  'REACTOME_MITOPHAGY',
  'REACTOME_PINK1_PRKN_MEDIATED_MITOPHAGY',
  'REACTOME_RECEPTOR_MEDIATED_MITOPHAGY'
)
mysets_autophagy <- purrr::map(
  .x = mysets_names_autophagy,
  .f = function(x = mysets_names_autophagy){gsdb %>% dplyr::filter(gs_name == x) %>% dplyr::pull(entrez_gene)}
)
base::names(mysets_autophagy) <- mysets_names_autophagy

gsea_res_autophagy_aging <- fgsea::fgsea(
  pathways = mysets_autophagy, 
  stats    = de_comb %>% dplyr::select(EntrezID, stat_aging) %>% dplyr::pull(name = EntrezID),
  minSize  = 15,
  maxSize  = 500
)

gsea_res_autophagy_trf <- fgsea::fgsea(
  pathways = mysets_autophagy, 
  stats    = de_comb %>% dplyr::select(EntrezID, stat_trf) %>% dplyr::pull(name = EntrezID),
  minSize  = 15,
  maxSize  = 500
)

readr::write_csv(x = gsea_res_autophagy_aging, file = './data/de_genes_GSEA/gsea_res_autophagy_aging.csv')
readr::write_csv(x = gsea_res_autophagy_trf, file = './data/de_genes_GSEA/gsea_res_autophagy_trf.csv')


# mitochondrion
mysets_names_mitochondrion <- c(
  'ELECTRON_TRANSPORT_GO_0006118',
  'GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT',
  'GOBP_ELECTRON_TRANSPORT_CHAIN',
  'GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN',
  'GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE',
  'GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_UBIQUINOL_TO_CYTOCHROME_C',
  'GOBP_REGULATION_OF_ELECTRON_TRANSFER_ACTIVITY',
  'GOBP_REGULATION_OF_MITOCHONDRIAL_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT',
  'GOBP_REGULATION_OF_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE',
  'GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN',
  'GOMF_ELECTRON_TRANSFER_ACTIVITY',
  'REACTOME_RESPIRATORY_ELECTRON_TRANSPORT',
  'REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS',
  'REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT',
  'WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA',
  'GOBP_ANTEROGRADE_AXONAL_TRANSPORT_OF_MITOCHONDRION',
  'GOBP_AUTOPHAGY_OF_MITOCHONDRION',
  'GOBP_AXONAL_TRANSPORT_OF_MITOCHONDRION',
  'GOBP_CALCIUM_IMPORT_INTO_THE_MITOCHONDRION',                                               
  'GOBP_ESTABLISHMENT_OF_MITOCHONDRION_LOCALIZATION',
  'GOBP_MAINTENANCE_OF_PROTEIN_LOCATION_IN_MITOCHONDRION',
  'GOBP_MITOCHONDRION_DISTRIBUTION',
  'GOBP_MITOCHONDRION_ENDOPLASMIC_RETICULUM_MEMBRANE_TETHERING',
  'GOBP_MITOCHONDRION_LOCALIZATION',
  'GOBP_MITOCHONDRION_MORPHOGENESIS',
  'GOBP_MITOCHONDRION_ORGANIZATION',
  'GOBP_NEGATIVE_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION',
  'GOBP_NEGATIVE_REGULATION_OF_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_MITOCHONDRION',
  'GOBP_NEGATIVE_REGULATION_OF_MITOCHONDRION_ORGANIZATION',
  'GOBP_NEGATIVE_REGULATION_OF_PROTEIN_TARGETING_TO_MITOCHONDRION',
  'GOBP_POSITIVE_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION',
  'GOBP_POSITIVE_REGULATION_OF_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_MITOCHONDRION',
  'GOBP_POSITIVE_REGULATION_OF_MITOCHONDRION_ORGANIZATION',
  'GOBP_PROTEIN_LOCALIZATION_TO_MITOCHONDRION',
  'GOBP_PROTEIN_PROCESSING_INVOLVED_IN_PROTEIN_TARGETING_TO_MITOCHONDRION',
  'GOBP_PROTEIN_TARGETING_TO_MITOCHONDRION',
  'GOBP_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION',
  'GOBP_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION_IN_RESPONSE_TO_MITOCHONDRIAL_DEPOLARIZATION',
  'GOBP_REGULATION_OF_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_MITOCHONDRION',
  'GOBP_REGULATION_OF_MITOCHONDRION_ORGANIZATION',
  'GOBP_RNA_IMPORT_INTO_MITOCHONDRION',
  'GOCC_MITOCHONDRION',
  'GOMF_MITOCHONDRION_TARGETING_SEQUENCE_BINDING',
  'HP_ABNORMALITY_OF_THE_MITOCHONDRION',
  'REACTOME_RRNA_MODIFICATION_IN_THE_MITOCHONDRION',
  'REACTOME_RRNA_PROCESSING_IN_THE_MITOCHONDRION',
  'REACTOME_TRNA_MODIFICATION_IN_THE_MITOCHONDRION',
  'REACTOME_TRNA_PROCESSING_IN_THE_MITOCHONDRION'
)
mysets_mitochondrion <- purrr::map(
  .x = mysets_names_mitochondrion,
  .f = function(x = mysets_names_mitochondrion){gsdb %>% dplyr::filter(gs_name == x) %>% dplyr::pull(entrez_gene)}
)
base::names(mysets_mitochondrion) <- mysets_names_mitochondrion

gsea_res_mitochondrion_aging <- fgsea::fgsea(
  pathways = mysets_mitochondrion, 
  stats    = de_comb %>% dplyr::select(EntrezID, stat_aging) %>% dplyr::pull(name = EntrezID),
  minSize  = 15,
  maxSize  = 500
)

gsea_res_mitochondrion_trf <- fgsea::fgsea(
  pathways = mysets_mitochondrion, 
  stats    = de_comb %>% dplyr::select(EntrezID, stat_trf) %>% dplyr::pull(name = EntrezID),
  minSize  = 15,
  maxSize  = 500
)

readr::write_csv(x = gsea_res_mitochondrion_aging, file = './data/de_genes_GSEA/gsea_res_mitochondrion_aging.csv')
readr::write_csv(x = gsea_res_mitochondrion_trf, file = './data/de_genes_GSEA/gsea_res_mitochondrion_trf.csv')


# to all gene sets
mysets_names_all <- gsdb %>% dplyr::filter(gs_cat == 'C5') %>% dplyr::pull(gs_name) %>% base::unique()
mysets_all <- purrr::map(
  .x = mysets_names_all,
  .f = function(x = mysets_names_all){gsdb %>% dplyr::filter(gs_name == x) %>% dplyr::pull(entrez_gene)}
)
base::names(mysets_all) <- mysets_names_all

gsea_res_all_aging <- fgsea::fgsea(
  pathways = mysets_all, 
  stats    = de_comb %>% dplyr::select(EntrezID, stat_aging) %>% dplyr::pull(name = EntrezID),
  minSize  = 15,
  maxSize  = 500
)

gsea_res_all_trf <- fgsea::fgsea(
  pathways = mysets_all, 
  stats    = de_comb %>% dplyr::select(EntrezID, stat_trf) %>% dplyr::pull(name = EntrezID),
  minSize  = 15,
  maxSize  = 500
)

readr::write_csv(x = gsea_res_all_aging, file = './data/de_genes_GSEA/gsea_res_all_aging.csv')
readr::write_csv(x = gsea_res_all_trf, file = './data/de_genes_GSEA/gsea_res_all_trf.csv')


## GSVA
# autophagy
gsva_res_autophagy <- GSVA::gsva(
  expr = txi.kallisto.tsv$counts,
  gset.idx.list = mysets_autophagy
)

hm_gsva_res_autophagy <- gsva_res_autophagy %>% 
  pheatmap::pheatmap(
    show_rownames = T,
    show_colnames = F,
    cluster_cols = F,
    scale = 'row',
    annotation_col = dplyr::select(sample_info, 2),
    gaps_col = c(5,11),
    cellwidth = 10
  )

ggplot2::ggsave(
  filename = 'hm_gsva_res_autophagy.tiff',
  plot = hm_gsva_res_autophagy,
  device = 'tiff',
  width = 13,
  height = 5,
  dpi ='retina',
  bg = 'white'
)


# mitochondrion
gsva_res_mitochondrion <- GSVA::gsva(
  expr = txi.kallisto.tsv$counts,
  gset.idx.list = mysets_mitochondrion
)

hm_gsva_res_mitochondrion <- gsva_res_mitochondrion %>% 
  pheatmap::pheatmap(
    show_rownames = T,
    show_colnames = F,
    cluster_cols = F,
    scale = 'row',
    annotation_col = dplyr::select(sample_info, 2),
    gaps_col = c(5,11),
    cellwidth = 10
  )

ggplot2::ggsave(
  filename = 'hm_gsva_res_mitochondrion.tiff',
  plot = hm_gsva_res_mitochondrion,
  device = 'tiff',
  width = 18,
  height = 7,
  dpi ='retina',
  bg = 'white'
)


## SIRT1 regulated genes
sirt1_brain_genes <- readr::read_delim(file = './data/external/GSE28790.top.table.tsv', delim = '\t')

sirt1_brain_genes_id2ensembl <- AnnotationDbi::select(
  x = mouse430a2.db,
  keys = sirt1_brain_genes$ID,
  columns = c('PROBEID', 'ENTREZID'),
  keytype = 'PROBEID'
)
sirt1_brain_genes_id2ensembl <- sirt1_brain_genes_id2ensembl %>% dplyr::filter(!is.na(ENTREZID))

sirt1_brain_genes <- base::list(
  sirt1_brain_genes_UP = sirt1_brain_genes %>% dplyr::filter(adj.P.Val < .05 & logFC > base::log(1.5,2)) %>% dplyr::pull(ID),
  sirt1_brain_genes_DOWN = sirt1_brain_genes %>% dplyr::filter(adj.P.Val < .05 & logFC < -base::log(1.5,2)) %>% dplyr::pull(ID)
)

sirt1_brain_genes <- purrr::map(
  .x = sirt1_brain_genes,
  .f = function(x){sirt1_brain_genes_id2ensembl %>% dplyr::filter(PROBEID %in% x) %>% dplyr::pull(ENTREZID)}
)

gsea_res_sirt1_brain_aging <- fgsea::fgsea(
  pathways = sirt1_brain_genes, 
  stats    = de_comb %>% dplyr::select(EntrezID, stat_aging) %>% dplyr::pull(name = EntrezID),
  minSize  = 15,
  maxSize  = 500
)

gsea_res_sirt1_brain_trf <- fgsea::fgsea(
  pathways = sirt1_brain_genes, 
  stats    = de_comb %>% dplyr::select(EntrezID, stat_trf) %>% dplyr::pull(name = EntrezID),
  minSize  = 15,
  maxSize  = 500
)

readr::write_csv(x = gsea_res_sirt1_brain_aging, file = './data/de_genes_GSEA/gsea_res_sirt1_brain_aging.csv')
readr::write_csv(x = gsea_res_sirt1_brain_trf, file = './data/de_genes_GSEA/gsea_res_sirt1_brain_trf.csv')

gsva_res_sirt1_brain <- GSVA::gsva(
  expr = txi.kallisto.tsv$counts,
  gset.idx.list = sirt1_brain_genes
)

hm_gsva_res_sirt1_brain <- gsva_res_sirt1_brain %>% 
  pheatmap::pheatmap(
    show_rownames = T,
    show_colnames = F,
    cluster_cols = F,
    scale = 'row',
    annotation_col = dplyr::select(sample_info, 2),
    gaps_col = c(5,11)
  )

ggplot2::ggsave(
  filename = 'hm_gsva_res_sirt1_brain.tiff',
  plot = hm_gsva_res_sirt1_brain,
  device = 'tiff',
  width = 10,
  height = 5,
  dpi ='retina',
  bg = 'white'
)


## simple heat maps of gene sets
# mitochondrion
vector_cc_mitochondrion <- gsdb %>% dplyr::filter(gs_name == 'GOCC_MITOCHONDRION') %>% dplyr::pull(entrez_gene) %>% as.character()
hm_cc_mitochondrion_genes <- de_comb %>%
  dplyr::filter(EntrezID %in% vector_cc_mitochondrion) %>%
  dplyr::select(2:17) %>%
  base::as.matrix() %>%
  pheatmap::pheatmap(
    show_rownames = F,
    show_colnames = F,
    cluster_cols = F,
    cluster_rows = T,
    scale = 'row',
    annotation_col = dplyr::select(sample_info, 2),
    gaps_col = c(5,11)
  )

hm_cc_mitochondrion_genes

ggplot2::ggsave(
  filename = 'hm_cc_mitochondrion_genes.tiff',
  plot = hm_cc_mitochondrion_genes,
  device = 'tiff',
  width = 5,
  height = 5,
  dpi ='retina',
  bg = 'white'
)

# autophagy - mytophagy
vector_wp_autophagy <- gsdb %>% dplyr::filter(gs_name == 'WP_AUTOPHAGY') %>% dplyr::pull(entrez_gene) %>% as.character()
vector_bp_autophagy <- gsdb %>% dplyr::filter(gs_name == 'GOBP_REGULATION_OF_AUTOPHAGY') %>% dplyr::pull(entrez_gene) %>% as.character()
vector_kegg_autophagy <- gsdb %>% dplyr::filter(gs_name == 'KEGG_REGULATION_OF_AUTOPHAGY') %>% dplyr::pull(entrez_gene) %>% as.character()
vector_bp_mitophagy <- gsdb %>% dplyr::filter(gs_name == 'GOBP_MITOPHAGY') %>% dplyr::pull(entrez_gene) %>% as.character()
vector_bp_selective_autophagy <- gsdb %>% dplyr::filter(gs_name == 'GOBP_SELECTIVE_AUTOPHAGY') %>% dplyr::pull(entrez_gene) %>% as.character()

vector_autophagy_mitophagy <- base::unique(c(vector_wp_autophagy, vector_bp_autophagy, vector_bp_mitophagy, vector_kegg_autophagy, vector_bp_selective_autophagy))

hm_autophagy_mitophagy_genes <- de_comb %>%
  dplyr::filter(EntrezID %in% vector_autophagy_mitophagy) %>%
  dplyr::select(2:17) %>%
  base::as.matrix() %>%
  pheatmap::pheatmap(
    show_rownames = F,
    show_colnames = F,
    cluster_cols = F,
    cluster_rows = T,
    scale = 'row',
    annotation_col = dplyr::select(sample_info, 2),
    gaps_col = c(5,11)
  )

hm_autophagy_mitophagy_genes

ggplot2::ggsave(
  filename = 'hm_autophagy_mitophagy_genes.tiff',
  plot = hm_autophagy_mitophagy_genes,
  device = 'tiff',
  width = 5,
  height = 5,
  dpi ='retina',
  bg = 'white'
)

