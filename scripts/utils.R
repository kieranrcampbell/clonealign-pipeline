library(ggplot2)
library(ggrepel)
library(glue)

make_volcano <- function(coef, fit, sce_de, ca) {
  all_clones <- sort(unique(ca$clone))
  modelled_clones <- colnames(fit$coefficients)[-1]
  base_clone <- setdiff(all_clones, modelled_clones)
  
  comparison <- glue("Clone {modelled_clones[coef-1]} vs clone {base_clone}")
  
  qvals <- p.adjust(fit$p.value[,coef], method = 'BH')
  
  df_limma <- data_frame(comparison = comparison,
                         log2foldchange = fit$coefficients[,coef], 
                         pval = fit$p.value[,coef],
                         qval = qvals,
                         ensembl_id = rownames(sce_de),
                         gene_symbol = rowData(sce_de)$Symbol,
                         is_significant = qval < 0.05) 
  
  sig_cols <- c("FALSE" = "grey80", "TRUE" = "darkred")
  theme_set(cowplot::theme_cowplot(font_size = 11))
  
  df_text <- top_n(df_limma, 40, abs(log2foldchange))
  
  plot <- ggplot(df_limma, aes(x = log2foldchange, y = -log10(qval), color = is_significant)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = sig_cols, name = "Significantly differentially expressed") +
    geom_text_repel(data = df_text, aes(label = gene_symbol), color = 'black', size = 3) +
    labs(x = expression(log[2]~"(fold change)"),
         y = expression(-log[10]~"(q-value)"),
         title = "Differential expression on non-CNV genes",
         subtitle = comparison) +
    theme(legend.position = "bottom",
          legend.box.background = element_rect(linetype = 1, size = .1),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8))
  
  list(
    plot = plot,
    df = df_limma
  )
}


