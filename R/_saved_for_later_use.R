pca_df <- pca.res$rotation %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column('sample') %>%
  dplyr::left_join(., sample_info, by = 'sample')

plotly::plot_ly(
  pca_df,
  type = 'scatter3d',
  showlegend = F,
  mode = "markers+text",
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~group,
  symbol = ~group,
  colors = c('#e60000', '#009900','#0C4B8E'),
  symbols = c('diamond', 'star-triangle-up', 'square'),
  text = ~row.names(pca_df)
)