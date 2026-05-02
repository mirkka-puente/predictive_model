# Instalar paquetes (solo la primera vez)
#install.packages("readxl")
#install.packages("tidyverse")
#install.packages("limma")
#install.packages("pheatmap")
#install.packages("ggplot2")

# Cargar librerías
library(readxl)
library(tidyverse)
library(dplyr)
library(limma)
library(pheatmap)
library(ggplot2)

#--------------------------------------------------
# 1. Leer archivo 
#--------------------------------------------------
archivo <- read.csv("S3_FPKM.csv", row.names = 1, check.names = FALSE)
head(archivo)
str(archivo)

# Cambio de columnas y filas
archivo <- archivo[,-c(19,20)]
dt <- as.data.frame(t(archivo))
str(dt)

#--------------------------------------------------
# 2. Asignar tratamiento y hpi
#--------------------------------------------------
treatment <- c(
  "mock", "avr", "vir",   # Sample_1,  2,  3  — 1hpi
  "mock", "avr", "vir",   # Sample_7,  8,  9  — 6hpi
  "mock", "avr", "vir",   # Sample_10, 11, 12 — 12hpi
  "mock", "avr", "vir",   # Sample_13, 14, 15 — 1hpi
  "mock", "avr", "vir",   # Sample_19, 20, 21 — 6hpi
  "mock", "avr", "vir"    # Sample_22, 23, 24 — 12hpi
)
dt$tratamiento <- treatment

hpi <- c(
  1,  1,  1,    # Sample_1,  2,  3
  6,  6,  6,    # Sample_7,  8,  9
  12, 12, 12,   # Sample_10, 11, 12
  1,  1,  1,    # Sample_13, 14, 15
  6,  6,  6,    # Sample_19, 20, 21
  12, 12, 12    # Sample_22, 23, 24
)
dt$hour <- hpi

# Eliminar mock
gene_set <- dt %>% filter(tratamiento != "mock")

# Separar metadata y expresión
treatment_hours <- gene_set[, c("tratamiento", "hour")]
gene_set_num    <- gene_set[, !(colnames(gene_set) %in% c("tratamiento", "hour"))]

#--------------------------------------------------
# 3. Curación de datos
#--------------------------------------------------
any(is.na(gene_set_num))

# Filtro 1: eliminar genes con cero expresión
keep <- colSums(gene_set_num > 0) > 0
gene_set_filtered <- gene_set_num[, keep]

# Filtro 2: gen debe tener FPKM > 1 en al menos
# 3 muestras de avr O 3 muestras de vir
avr_idx <- which(treatment_hours$tratamiento == "avr")
vir_idx <- which(treatment_hours$tratamiento == "vir")
expressed_avr <- colSums(gene_set_filtered[avr_idx, ] > 1) >= 3
expressed_vir <- colSums(gene_set_filtered[vir_idx, ] > 1) >= 3
keep2 <- expressed_avr | expressed_vir

cat("Genes antes del filtro: ", ncol(gene_set_filtered), "\n")
cat("Genes después del filtro:", sum(keep2), "\n")

gene_set_filtered <- gene_set_filtered[, keep2]

# Transformación logarítmica
gene_set_log <- log2(gene_set_filtered + 1)

#--------------------------------------------------
# 4. Chequeo
#--------------------------------------------------
cat("Muestras en expresión: ", nrow(gene_set_log), "\n")
cat("Muestras en metadata:  ", nrow(treatment_hours), "\n")
cat("Labels coinciden:      ", all(rownames(gene_set_log) == rownames(treatment_hours)), "\n")

#==================================================
# PASO 1 — PCA
#==================================================
pca    <- prcomp(gene_set_log, scale. = TRUE)
var_exp <- summary(pca)$importance[2, ] * 100

pca_df <- data.frame(
  PC1         = pca$x[, 1],
  PC2         = pca$x[, 2],
  tratamiento = treatment_hours$tratamiento,
  hpi         = factor(treatment_hours$hour)
)

ggplot(pca_df, aes(x = PC1, y = PC2, 
                   color = tratamiento, 
                   shape = hpi,
                   label = rownames(pca_df))) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text(vjust = -0.8, size = 3) +          # etiquetas de muestra
  scale_color_manual(values = c(avr = "#5DCAA5", vir = "#D85A30")) +
  labs(
    title = "PCA — avr vs vir",
    x     = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y     = paste0("PC2 (", round(var_exp[2], 1), "%)"),
    color = "Tratamiento",
    shape = "hpi"
  ) +
  theme_classic(base_size = 13)

#==================================================
# PASO 2 — HEATMAP
#==================================================
# Top 500 genes más variables
gene_var <- apply(gene_set_log, 2, var)
top500   <- names(sort(gene_var, decreasing = TRUE)[1:500])
mat      <- t(gene_set_log[, top500])   # pheatmap: genes como filas

# Anotación de columnas
annotation_col <- data.frame(
  Treatment = treatment_hours$tratamiento,
  hpi       = factor(treatment_hours$hour)
)
rownames(annotation_col) <- rownames(gene_set_log)

ann_colors <- list(
  Treatment = c(avr = "#5DCAA5", vir = "#D85A30"),
  hpi       = c("1" = "#EEEDFE", "6" = "#7F77DD", "12" = "#3C3489")
)

pheatmap(
  mat,
  annotation_col           = annotation_col,
  annotation_colors        = ann_colors,
  show_rownames            = FALSE,
  show_colnames            = TRUE,
  scale                    = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  color                    = colorRampPalette(c("#185FA5", "white", "#D85A30"))(100),
  legend_breaks            = c(-2, -1, 0, 1, 2),
  legend_labels            = c("-2", "-1", "0", "1", "2\nRow Z-score"),
  fontsize_col             = 9,
  border_color             = NA,
  main                     = "Top 500 genes variables — avr vs vir"
)

#==================================================
# PASO 3 — LIMMA
#==================================================
treatment_factor <- factor(treatment_hours$tratamiento, levels = c("avr", "vir"))
design <- model.matrix(~ treatment_factor)
colnames(design) <- c("Intercept", "vir_vs_avr")

# Fit con trend=TRUE (recomendado para FPKM log-transformado)
fit <- lmFit(t(gene_set_log), design)
fit <- eBayes(fit, trend = TRUE)

# Resultados
res <- topTable(fit, coef = "vir_vs_avr",
                number = Inf, adjust.method = "BH")

# Diagnóstico
hist(res$P.Value, breaks = 50, col = "#5DCAA5",
     main = "Distribución de p-valores (avr vs vir)",
     xlab = "P-value")

# Top 10 genes
cat("\n=== TOP 10 GENES ===\n")
print(head(res[order(res$P.Value),
               c("logFC", "AveExpr", "P.Value", "adj.P.Val")], 10))

# Filtrar genes DE
sig <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ]

cat("\n=== RESULTADOS FINALES ===\n")
cat("Genes DE totales: ", nrow(sig), "\n")
cat("Más altos en vir: ", nrow(sig[sig$logFC > 1,  ]), "\n")
cat("Más altos en avr: ", nrow(sig[sig$logFC < -1, ]), "\n")


