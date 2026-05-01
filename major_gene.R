# Instalar paquetes
#install.packages("readxl")
#install.packages("tidyverse")

# Cargar librerías
library(readxl)
library(tidyverse)
library(dplyr)

#--------------------------------------------------
# 1. Leer archivo 
#--------------------------------------------------
archivo <- read.csv("S3_FPKM.csv", row.names = 1, check.names = FALSE)
head(archivo)
str(archivo)

# Cambio de columnas y filas
archivo <- archivo[,-c(19,20)]
dt <- as.data.frame(t(archivo))

# Revisar estructura
str(dt)

#--------------------------------------------------
# 2. Asignar hpi y tratamiento
#--------------------------------------------------

# Manualmente las columnas de tratamiento y hpi
#tratamiento
treatment <- c(
  "mock",    # para Sample_1
  "avr",    # para Sample_2
  "vir",    # para Sample_3
  "mock",    # para Sample_7
  "avr",    # para Sample_8
  "vir",    # para Sample_9
  "mock",    # para Sample_10
  "avr",    # para Sample_11
  "vir",    # para Sample_12
  "mock",    # para Sample_13
  "avr",    # para Sample_14
  "vir",    # para Sample_15
  "mock",    # para Sample_19
  "avr",    # para Sample_20
  "vir",    # para Sample_21
  "mock",    # para Sample_22
  "avr",    # para Sample_23
  "vir"    # para Sample_24
)

dt$tratamiento <- treatment

#Hours post inoculation
hpi <- c(
  1,    # Sample_1
  1,    # Sample_2
  1,    # Sample_3
  6,    # Sample_7
  6,   # Sample_8
  6,   # Sample_9
  12,    # Sample_10
  12,    # Sample_11
  12,   # Sample_12
  1,    # Sample_13
  1,    # Sample_14
  1,   # Sample_15
  6,   # Sample_19
  6,    # Sample_20
  6,    # Sample_21
  12,   # Sample_22
  12,   # Sample_23
  12   # Sample_24
)
dt$hour <- hpi

# Eliminar los mock
gene_set <- dt %>% filter(tratamiento != "mock")

# Dividir mi data en numerica y no numerica
treatment_hours <- gene_set[, c("tratamiento", "hour")]
gene_set_num <- gene_set[, !(colnames(gene_set) %in% c("tratamiento", "hour"))]


#--------------------------------------------------
# 2. Curacion de Data 
#--------------------------------------------------
# Hay algun NA en el dataset?
any(is.na(gene_set_num))

# Eliminar las columnas con cero expresion genica
keep <- colSums(gene_set_num > 0) > 0
gene_set_filtered <- gene_set_num[, keep]


# Transformacion logaritmica
gene_set_log <- log2(gene_set_filtered + 1)

#--------------------------------------------------
# 3. PCA
#--------------------------------------------------

library(ggplot2)

# PCA
pca <- prcomp(gene_set_log, scale. = TRUE)

# Armar dataframe con resultados + metadatos
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  tratamiento = treatment_hours$tratamiento,
  hour        = factor(treatment_hours$hour)
)

# Varianza explicada
var_exp <- summary(pca)$importance[2, ] * 100

# Plot components
ggplot(pca_df, aes(x = PC1, y = PC2, color = tratamiento, shape = hour)) +
  geom_point(size = 4, alpha = 0.85) +
  labs(
    x     = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y     = paste0("PC2 (", round(var_exp[2], 1), "%)"),
    color = "Treatment",
    shape = "hpi"
  ) +
  theme_classic(base_size = 10)

#--------------------------------------------------
# 4. Exploracion de datos
#--------------------------------------------------
#install.packages("pheatmap")
library(pheatmap)

# Top 500 most variable genes
gene_var <- apply(gene_set_log, 2, var)
top500 <- names(sort(gene_var, decreasing = TRUE)[1:500])
mat <- t(gene_set_log[, top500])  # genes deben ser filas

# Anotacion para el top del heatmap
annotation_col <- data.frame(
  Treatment = treatment_hours$tratamiento,
  hpi       = factor(treatment_hours$hour)
)
rownames(annotation_col) <- rownames(gene_set_log)

# Colors for annotation
ann_colors <- list(
  Treatment = c(avr = "#5DCAA5", vir = "#D85A30"),
  hpi       = c("1" = "#EEEDFE", "6" = "#7F77DD", "12" = "#3C3489")
)

# Plot
pheatmap(
  mat,
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  scale             = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  color             = colorRampPalette(c("#185FA5", "white", "#D85A30"))(100),
  fontsize_col      = 9,
  border_color      = NA,
  main              = "Top 500 variable genes — avr vs vir",
  legend_breaks     = c(-2, -1, 0, 1, 2),
  legend_labels     = c("Row Z-score \n -2", "-1", "0", "1", "2")
)

#--------------------------------------------------
# 5. Limma model 
#--------------------------------------------------
#install.packages("BiocManager")
BiocManager::install("limma")
library(limma)

# Cambiar las hhoras a factores
treatment_hours$hour <- factor(treatment_hours$hour, levels = c("1", "6", "12"))

# Matrix
design <- model.matrix(~ tratamiento * hour, data = treatment_hours)
colnames(design)

# Modelo Limma 
fit <- limma::lmFit(t(gene_set_log), design)
fit <- limma::eBayes(fit)

# --- Comparacion entre tratamientos ---

# 1. avr vs vir (sin tiempo)
res_overall <- topTable(fit, coef = "tratamientovir", 
                        number = Inf, adjust.method = "BH")

# 2. Interaccion: does avr/vir difference change at 6 hpi vs 1 hpi?
res_int6  <- topTable(fit, coef = "tratamientovir:hour6",  
                      number = Inf, adjust.method = "BH")

# 3. Interaction: does avr/vir difference change at 12 hpi vs 1 hpi?
res_int12 <- topTable(fit, coef = "tratamientovir:hour12", 
                      number = Inf, adjust.method = "BH")

# --- Filter significant genes (adj.P < 0.05, |logFC| > 1) ---
sig_overall <- res_overall[res_overall$adj.P.Val < 0.05 & abs(res_overall$logFC) > 1, ]
sig_int6    <- res_int6[res_int6$adj.P.Val < 0.05 & abs(res_int6$logFC) > 1, ]
sig_int12   <- res_int12[res_int12$adj.P.Val < 0.05 & abs(res_int12$logFC) > 1, ]

# Quick summary
cat("Overall DE genes:", nrow(sig_overall), "\n")
cat("Interaction at 6 hpi:", nrow(sig_int6), "\n")
cat("Interaction at 12 hpi:", nrow(sig_int12), "\n")

# --- DIAGNOSIS ---

# 1. How many genes pass BEFORE the logFC filter?
sig_ponly <- res_overall[res_overall$adj.P.Val < 0.05, ]
cat("Genes with adj.P < 0.05 (no logFC filter):", nrow(sig_ponly), "\n")

# 2. What if you relax to raw p-value?
sig_raw <- res_overall[res_overall$P.Value < 0.05, ]
cat("Genes with raw P < 0.05:", nrow(sig_raw), "\n")

# 3. What does the p-value distribution look like?
hist(res_overall$P.Value, breaks = 50, 
     main = "P-value distribution", xlab = "P-value",
     col = "#7F77DD")
# A good experiment shows enrichment near 0 (left spike)
# A flat/uniform histogram = no signal or model problem

# 4. Check your design matrix — does it look right?
print(design)
# You should see 18 rows, columns for tratamientovir, hour6, hour12, interactions

# 5. Check sample sizes per group
table(treatment_hours$tratamiento, treatment_hours$hour)

# How many degrees of freedom does your model have?
fit$df.residual
# If you see values of 6 or less — that's your problem

