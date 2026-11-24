##
## Colour palettes. Fixed colours for the labels used so far, so
## figures of different runs are comparable; celltype_palette() adds stable colours for any new label.
##

pam50_subtype_colors <- c("Lum" = "#00BA38", "TNBC" = "#619CFF", "Her2+" = "#F8766D", "Unknown" = "grey70")
scsubtype_colors     <- c("Basal_SC" = "#619CFF", "Her2E_SC" = "#F8766D", "LumA_SC" = "lightgreen", "LumB_SC" = "darkgreen")
igg_level_colors     <- c("High" = "#A11212", "Medium" = "#F5A623", "Low" = "#F2E205")
pcr_colors           <- c("No" = "white", "Yes" = "black")
tls_colors           <- c("0" = "#D9D9D9", "1a" = "#9013FE", "1b" = "#C77DD6", "2a" = "#F9B8D0", "2b" = "#FF00FF")

lineage_colors <- c("Tumor" = "#875cb8", "Leukocytes" = "#d6790e", "Stromal" = "#0a6807",
                    "Normal_Epithelium" = "#E7298A", "Noise" = "grey70")

# Neighbouring populations on the UMAP get contrasting colours
celltype_colors <- c(
  # Leukocytes
  "TCell_cyto" = "#1F78B4", "TCell_naive" = "#A6CEE3", "Treg" = "#6A3D9A", "NK" = "#41B6C4",
  "BCell" = "#E31A1C", "PlasmaBlast" = "#FB9A99", "TAM" = "#FF7F00",
  "cDC" = "#FDBF6F", "actDC" = "#B15928", "pDC" = "#CAB2D6", "Mast" = "#E7298A", "Proliferative" = "#252525",
  # Stromal
  "mCAF" = "#33A02C", "ecm_myCAF" = "#006D2C", "iCAF" = "#B2DF8A", "Pericyte" = "#7F7F00",
  "Adipocyte" = "#FFD92F", "Endothelial" = "#9E0142", "Lymphatic" = "#F46D43",
  "Myoepithelial" = "#BF812D", "CAF_Mixed" = "#C8D9B4",
  # Tumour
  "Lum_Her2" = "#C51B7D", "Basal" = "#253494", "Cycling" = "#41B6C4", "Interferon" = "#FEC44F",
  "Tumor_Stressed" = "#D6604D", "LowRNA" = "#BDBDBD", "Tumor_Ciliated" = "#80CDC1",
  # Normal epithelium
  "Normal_LuminalProgenitor" = "#F781BF",
  # Lineage-level labels (before sub-annotation) and noise/contamination populations
  "Tumor" = "#875cb8", "Leukocytes" = "#d6790e", "Stromal" = "#0a6807", "Normal_Epithelium" = "#E7298A",
  "Noise" = "grey70", "TAM_Noise" = "#D9C8B4", "CAF_Noise" = "#C8D9B4", "Tumor_Noise" = "#C8B4D9"
)

# ---------------------------------------------------------------
# Named colour vector for a set of labels: fixed colours where known,
# otherwise a hue derived from the label text (same label -> same colour in every figure).
# ---------------------------------------------------------------
celltype_palette <- function(labels) {
  labels <- sort(unique(as.character(labels)))
  out <- setNames(unname(celltype_colors[labels]), labels)
  new <- labels[is.na(out)]
  if (length(new)) {
    hue <- vapply(new, function(l) (sum(utf8ToInt(l) * seq_along(utf8ToInt(l))) %% 360), 0)
    out[new] <- grDevices::hcl(h = hue, c = 65, l = 55)
  }
  setNames(unname(out), labels)
}
