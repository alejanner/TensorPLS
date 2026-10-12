#' This function creates Venn diagrams to visualize overlaps between
#' different sets of features (genes, metabolites, or other biological entities).
#' Apply this function to see for example commons features between different VIPs objects can be useful 
#' @param feature_lists A list of character vectors representing different feature sets.
#'   Each vector should contain unique identifiers (e.g., gene names, metabolite IDs).
#'   Minimum 2 sets, maximum 5 sets for optimal visualization.
#' @param colors A character vector of colors for the Venn diagram circles.
#'   Default is c("orange", "blue", "red", "green", "purple").
#'
#' @return Displays the Venn diagram and invisibly returns the plot grobs.
#'
#' @examples
#' # Example with gene sets
#' genes_set1 <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS")
#' genes_set2 <- c("TP53", "EGFR", "PIK3CA", "AKT1", "PTEN")
#' genes_set3 <- c("BRCA1", "TP53", "CDKN2A", "RB1", "APC")
#' 
#' gene_lists <- list("Oncogenes" = genes_set1, 
#'                    "Tumor_Suppressors" = genes_set2,
#'                    "DNA_Repair" = genes_set3)
#' 
#' # Basic usage
#' draw_venn_diagram(gene_lists)
#' 
#' # With custom colors
#' draw_venn_diagram(gene_lists, colors = c("#FF6B6B", "#4ECDC4", "#45B7D1"))
#'
#' @import VennDiagram
#' @import grid
#' @export
draw_venn_diagram <- function(feature_lists, colors = c("orange", "blue", "red", "green",
  "purple")) {
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    stop("VennDiagram package is required. Please install it with: install.packages('VennDiagram')")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("grid package is required. Please install it with: install.packages('grid')")
  }
  if (!is.list(feature_lists)) {
    stop("feature_lists must be a list of character vectors")
  }
  if (length(feature_lists) < 2) {
    stop("At least 2 feature sets are required")
  }
  if (length(feature_lists) > 5) {
    warning("More than 5 sets may result in unclear visualization. Consider reducing the number of sets.")
  }
  if (!all(sapply(feature_lists, is.character))) {
    stop("All elements in feature_lists must be character vectors")
  }
  feature_lists <- feature_lists[sapply(feature_lists, length) > 0]
  if (length(feature_lists) < 2) {
    stop("After removing empty sets, at least 2 non-empty sets are required")
  }
  if (is.null(names(feature_lists))) {
    names(feature_lists) <- paste("Set", seq_along(feature_lists))
  }
  if (length(colors) < length(feature_lists)) {
    colors <- rep(colors, length.out = length(feature_lists))
    warning("Not enough colors provided. Colors have been recycled.")
  }
  selected_colors <- colors[1:length(feature_lists)]
  n_sets <- length(feature_lists)
  if (n_sets == 2) {
    venn_grobs <- VennDiagram::venn.diagram(x = feature_lists,
      category.names = names(feature_lists),
      filename = NULL, lwd = 2, lty = 'blank', fill = selected_colors[1:2],
      cex = 1.5, fontface = "bold", fontfamily = "sans", cat.cex = 1.2,
      cat.fontface = "bold", cat.default.pos = "outer", cat.pos = c(-27, 27),
      cat.dist = c(0.055, 0.055), cat.fontfamily = "sans", alpha = 0.5,
      col = "white")
  } else
    if (n_sets == 3) {
    venn_grobs <- VennDiagram::venn.diagram(x = feature_lists,
      category.names = names(feature_lists),
      filename = NULL, lwd = 2, lty = 'blank', fill = selected_colors[1:3],
      cex = 1.5, fontface = "bold", fontfamily = "sans", cat.cex = 1.2,
      cat.fontface = "bold", cat.default.pos = "outer", cat.pos = c(-40, 40,
        180), cat.dist = c(0.055, 0.055, 0.055), cat.fontfamily = "sans",
      alpha = 0.5, col = "white")
  } else
    if (n_sets == 4) {
    venn_grobs <- VennDiagram::venn.diagram(x = feature_lists,
      category.names = names(feature_lists),
      filename = NULL, lwd = 2, lty = 'blank', fill = selected_colors[1:4],
      cex = 1.5, fontface = "bold", fontfamily = "sans", cat.cex = 1.2,
      cat.fontface = "bold", cat.pos = c(0, 0, 0, 0), cat.dist = c(0.22,
        0.22, 0.11, 0.11), cat.fontfamily = "sans", alpha = 0.5,
      col = "white")
  } else
    if (n_sets == 5) {
    venn_grobs <- VennDiagram::venn.diagram(x = feature_lists,
      category.names = names(feature_lists),
      filename = NULL, lwd = 2, lty = 'blank', fill = selected_colors[1:5],
      cex = 1.5, fontface = "bold", fontfamily = "sans", cat.cex = 1.2,
      cat.fontface = "bold", cat.fontfamily = "sans", alpha = 0.5, col = "white")
  }
  grid::grid.newpage()
  grid::grid.draw(venn_grobs)
  return(invisible(venn_grobs))
}
