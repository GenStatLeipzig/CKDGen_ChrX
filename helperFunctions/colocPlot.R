colocPlot <- function(x,
                      title = "") {
  
  # Plot Co-Localisation Analysis Results
  # 
  # Args:
  #   x: data.frame to plot;
  #      first column with gene names
  #      one gene per row, one column per tissue;
  #      each cell consists of max(PP(H3), PP(H4)), in case of H3 with "-" sign
  #   title: Title for plot
  # 
  # Returns:
  #   ggplot2;
  #   tissues on x axis;
  #   genes on y axis;
  #   legend at right hand side, red for H3, blue for H4, white neither nor, grey for 'NA'
  
  gene.names             <- as.character(x[, 1])
  tissue.names           <- names(x)[-1]
  coloc.matrix           <- as.matrix(x[, -1])
  rownames(coloc.matrix) <- gene.names
  colnames(coloc.matrix) <- tissue.names
  coloc.matrix.t         <- t(coloc.matrix)
  
  ggcp <- ggcorrplot(coloc.matrix.t,
                     title = title) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = "grey",
                         limits = c(-1, 1),
                         name = "PP",
                         breaks = seq(-1, 1, by = 0.25)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  print(ggcp)
}

# Example:
  # .libPaths("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/amanMRO/")
  # library("ggcorrplot")
  # gene.names   <- toupper(letters)[c(1 : 10)]
  # tissue.names <- paste0("T", c(1 : 9))
  # test         <- data.frame(Trait = gene.names,
  #                            T1 = runif(5, -1, 1),
  #                            T2 = runif(5, -0.75, 0.75),
  #                            T3 = runif(5, -0.5, 0.5),
  #                            T4 = runif(5, -1, 0.75),
  #                            T5 = runif(5, -1, 0.5),
  #                            T6 = runif(5, -0.75, 0.5),
  #                            T7 = runif(5, -0.75, 1),
  #                            T8 = runif(5, -0.5, 1),
  #                            T9 = runif(5, -0.5, 0.75))
  # colocPlot(test, title = "Test")









