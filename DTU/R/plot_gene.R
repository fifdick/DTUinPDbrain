#' Build ggplot obj and dataframe 
#' 
#' Build the ggplot object, to plot data extracted with DTU::plot_genes_cohort()i
#' see plot_genes_cohort to see which inputs the function gets 
#' @param all_vals df created in DTU::plot_genes_cohort with fitted and observed fractions and gene level counts 
#' @param p_gene_id char string identifying gene which should be plotted  
#' @param genes_toplot gene annotation from which to extract only the transcript annotation specific to the gene in p_gene_id
#' @param out
#' @param plot
#' @param tool 
#' @param only_nom_sig 
#' @param selected_samples
#' @param want_jitter
#' @param conditions
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @import ggplot2 
plot_gene <- function(all_vals, p_gene_id, genes_toplot, out = "./", plot = T, tool, only_nom_sig = TRUE,
		      selected_samples = NULL, want_jitter = TRUE, conditions=c("Control", "Case")) {	
 # retrieve gene name from gene info dataframe
 gene_name <- unique(subset(genes_toplot, gene_id == p_gene_id)$gene_name)
 # factors (cant remember why I did that)
 all_vals$gene_id <- as.factor(all_vals$gene_id)
 all_vals$sample_id <- as.factor(all_vals$sample_id)
 #split into lists of dfs, one df for each gene
 all_vals <- split(all_vals, all_vals$gene_id)
 # select all transcripts of that gene
 df <- all_vals[[p_gene_id]] %>% dplyr::group_by(.data$sample_id)
 if (is.null(df)) {
    print("Gene:")
    print(p_gene_id)
    print("was not found in one of the count lists")
    print("using only the counts that have been found")
 } 
 debug_d <- df
 # get tx biotype info for strip labels
 gene_info <- get_gene_info(unique(df$transcript_id), tx = TRUE)
 df  <- dplyr::left_join(df, gene_info[, c("gene_id", "transcript_id", "tx_biotype")], by = c("gene_id", "tx_id"))
 df %<>% dplyr::mutate(tx_id_label = dplyr::case_when(
  grepl("ENSG", .data$transcript_id) ~ paste(.data$transcript_id, "Gene_expression", "(count scale)", sep = "\n"),
  .data$tx_biotype == "nonsense_mediated_decay" ~ paste(.data$transcript_id, "nonsense", "mediated_decay", sep = "\n"),
  TRUE ~ paste(.data$transcript_id, .data$tx_biotype, sep = "\n"))) %>% 
  dplyr::group_by(.data$transcript_id, .data$condition) 
 
 # function to create custom jitter such that observed and fitted values are visually separated in the plot
 jitter_ud <- function(x, amount = 0.15) {
  jitter <- stats::runif(length(x), -amount, amount)
  x <- sapply(seq(1, length(x)), function(i) {
   if (x[i] == 1 && jitter[i] > 0) {
    if (x[i] + jitter[i] > 1.2) {
     return(x[i] + jitter[i])
    } else { 
     return(x[i] + jitter[i] + 0.2)
    }
   }
   if (x[i] == 1 && jitter[i] < 0) {
    if (x[i] + jitter[i] < 0.8)
     return(x[i] + jitter[i])
    else
     return(x[i] + jitter[i] - 0.2)
   }
   if (x[i] == 0 && jitter[i] > 0) {
    if (x[i] + jitter[i] > 0.2)
     return(x[i] + jitter[i])
    else
     return(x[i] + jitter[i] + 0.2)
   }
   if (x[i] == 0 && jitter[i] < 0) {
    if (x[i] + jitter[i] < -0.2)
     return(x[i] + jitter[i])
    else
     return(x[i] + jitter[i] - 0.2)
   }
  })
  return(x)
 }
 
 selective_jitter <- function(x, g) {	
  x <-as.numeric(x) - 1
  x[which(g == "1")] <- jitter_ud(x[which(g == "1")])
   return(x)	
 }

 # create helper col to identify which are fitted which are observed values
 # create column x defined by jitter function
 # create helper function to summerize significance
 df %<>% dplyr::mutate(fitted = ifelse(grepl("observed", .data$countType), "0", "1")) %>% 
  dplyr::mutate(x = selective_jitter(.data$condition, .data$fitted)) %>% 
  dplyr::mutate(boxCol = dplyr::case_when(
   .data$nomSig == 1 && .data$sig == 0 ~ 1,
   .data$nomSig == 0 && .data$sig == 0 ~ 0,
   .data$nomSig == 1 && .data$sig == 1 ~ 2,
   TRUE ~ NA )) #a transcript is not nominally significant but its adj. pvalue is?!
    #(After correction significance is read off the main_df and the tx_pvalueStageR column)

 the_plot <- function(vals, want_jitter = TRUE, tool, conditions) {
  if (want_jitter == TRUE) {	
   p <- ggplot(subset(vals, countType == fitted_tx), aes(x = as.numeric(condition) - 1, y = frac, color = as.factor(boxCol))) +
    geom_boxplot(aes(group = condition), lwd = 0.2, outlier.size = 0, coef = 0) + 
    geom_point(data = subset(vals, countType == "fitted_tx" | countType == "observed_tx"),
     aes(x = x, y = frac, shape = fitted), size = .9, alpha = 0.7) +
    scale_x_continuous(breaks = seq(0, 1), labels = c("0" = conditions[1], "1" = conditions[2])) +
    scale_color_manual(name = "DTU", label = c("1" = "p-value < 0.05", "2" = "adj. p-value < 0.05",
     "0" = "not significant"), values = c("1" = "deepskyblue3", "2" = "darkblue", "0"= "black")) + 
    scale_shape_manual(name = "Data point:", values = c(1, 4))  
  } else {
  # add additional factor to the significance tag column. Used for colouring with ggplot: observed values will be lightgrey and werent assessed for significance (obviously)
  vals %<>% dplyr::mutate(boxCol = ifelse(.data$fitted == 1, boxCol, 3))
  p <- ggplot(vals, aes(x = condition, y = frac, fill = fitted)) +
   geom_boxplot(aes(x = condition, col = as.factor(boxCol)), lwd = 0.4, outlier.shape = NA) +
   geom_point(position = position_jitterdodge(jitter.width = 0.3), aes(group = countType, col = as.factor(boxCol)), size = 0.7) +
   scale_fill_manual(name = "Data point", labels = c("0" = "observed", "1" = "fitted"), values = c("0" = "grey89", "1"="white")) +
   scale_color_manual(name = "DTU", label = c("1" = "p-value < 0.05", "2" = "adj. p-value < 0.05", "0" = "not significant", "3" = "not tested"),
		      values = c("1" = "orange", "2" = "red", "0" = "black", "3" = "grey89")) +
   scale_x_discrete(labels = c("0" = conditions[1], "1" = conditions[2])) +
   guides(shape = F)
  }
  p <- p + 
   facet_wrap(. ~ tx_id_label, scales = "free_y") + 
   ggtitle(paste0(gene_name, " (", tool, ")")) + 
   theme_bw() + 
   theme(plot.title = element_text(lineheight = 0.8, face = "bold"), axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.2)) + 
   theme(strip.background = element_rect(fill = "black")) + 
   theme(strip.text = element_text(color = "white", face = "bold")) + 
   labs(x = "Ensembl transcript (v75)", y = "Transcript abundance", color = "nominal significance") +
   theme(axis.title.x = element_blank()) +
   theme(axis.text.x = element_text(angle = 0)) +
   geom_text(aes(label = name), hjust = 0, vjust = 0, col = "red")
  return(p)
 }
 make_plot <- function(vals, tool, want_jitter = TRUE) {
  vals %<>% dplyr::mutate(name = ifelse(.data$sample_id %in% selected_samples, as.character(.data$sample_id), ""))
  if (only_nom_sig == TRUE) {
   vals <- subset(vals, nomSig == 1)
   p <- the_plot(vals, want_jitter = want_jitter, tool = tool, conditions = conditions)
  } else {
  p <- the_plot(vals, want_jitter = want_jitter, tool = tool, conditions = conditions)		
  }
  return(p)
 }
 my.seed <- sample(1:10000, 1)
 set.seed(my.seed)
 tool <- ifelse(tool == "dex", "DEXSeq", "DRIMSeq")
 count_plot <- make_plot(df, tool = tool, want_jitter = want_jitter)
 grDevices::pdf(paste0(out, gene_name, "_", my.seed, "_counts_", tool, ".pdf"))
 grDevices::dev.off()
 if (plot == T) {
  print(count_plot)
 }       
 result <- list(data = df, plot = count_plot, debug = debug_d)
 return(result)   
}









