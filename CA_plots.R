##
## Plotting functions for the scRNA-seq cell annotation pipeline (3a-3f).
## plot_dimplot()/plot_featureplot()/plot_resolution_grid() live in CL_plots.R;
## source that file alongside this one.
##
library(ggplot2)
library(Seurat)
library(patchwork)

# ---------------------------------------------------------------
# Marker gene dot plot per cluster/group - size
# ---------------------------------------------------------------
plot_marker_dotplot <- function(object, marker_groups, results_path,
                                 filename = "marker_dotplot.png",
                                 group_by = "seurat_clusters",
                                 title = "Marker expression by group") {
  markers <- unique(unlist(marker_groups))
  markers <- markers[markers %in% rownames(object)]  # skip genes not found in the object

  # Map each gene -> "(lineage)_(gene)" label
  gene_label <- sapply(markers, function(g) {
    grp <- names(marker_groups)[sapply(marker_groups, function(x) g %in% x)][1]
    paste0(grp, "_", g)
  })
  gene_label <- gene_label[markers]

  p <- DotPlot(object, features = markers, group.by = group_by) +
    coord_flip() +
    scale_color_gradient(low = "blue", high = "red") +
    scale_x_discrete(labels = gene_label) +
    labs(title = title) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) # rotate labels and align to the top

  ggsave(paste0(results_path, filename), p,
         width = 10, height = 10, dpi = 300, limitsize = FALSE,
         bg = "white")
}