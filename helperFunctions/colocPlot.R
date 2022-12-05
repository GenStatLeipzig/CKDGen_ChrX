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
  #   legend at right hand side, red for H3, blue for H4, white neither nor, grey for 'NA' values
  
  # Create matrix with posterior probabilities only (including transposition)
  gene.names             <- as.character(x[, 1])
  tissue.names           <- names(x)[-1]
  coloc.matrix           <- as.matrix(x[, -1])
  rownames(coloc.matrix) <- gene.names
  colnames(coloc.matrix) <- tissue.names
  coloc.matrix.t         <- t(coloc.matrix)
  
  # Create matrix indicating which cells consists of 'NA' values
  coloc.matrix.t.na                   <- is.na(coloc.matrix.t)
  coloc.matrix.t.na.numeric           <- apply(coloc.matrix.t.na, 2, as.numeric)
  rownames(coloc.matrix.t.na.numeric) <- tissue.names
  colnames(coloc.matrix.t.na.numeric) <- gene.names
  
  # Substitute 'NA' values in original data by '0'
  coloc.matrix.t.nato0                              <- coloc.matrix.t
  coloc.matrix.t.nato0[is.na(coloc.matrix.t.nato0)] <- 0
  
  # Plotting
  ggcp <- ggcorrplot(coloc.matrix.t.nato0,
                     title = title,
                     p.mat = coloc.matrix.t.na.numeric,
                     sig.level = 0.5,
                     insig = "pch") +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = "grey",
                         limits = c(-1, 1),
                         name = "PP",
                         breaks = seq(-1, 1, by = 0.25)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  print(ggcp)
}

# Example:
  # # R package 'ggcorrplot' on aman
  # .libPaths("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/amanMRO/")
  # library("ggcorrplot")
  # 
  # # Data without 'NA' values
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
  # 
  # # Data with 'NA' values
  # test.na        <- test
  # test.na[8, 3]  <- NA
  # test.na[10, 5] <- NA
  # colocPlot(test.na, title = "Test with 'NA' Values")









