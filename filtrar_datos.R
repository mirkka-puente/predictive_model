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
#--------------------------------------------------
# FIX: modelo que controla por tiempo
#--------------------------------------------------

# hour es factor
treatment_hours$hour <- factor(treatment_hours$hour, levels = c("1", "6", "12"))
treatment_factor     <- factor(treatment_hours$tratamiento, levels = c("avr", "vir"))

# Diseño aditivo: controla por hora, compara avr vs vir dentro de cada hora
design <- model.matrix(~ hour + treatment_factor, data = treatment_hours)
colnames(design)

print(design)  # debe tener 12 filas, 4 columnas

# Fit
fit <- lmFit(t(gene_set_log), design)
fit <- eBayes(fit, trend = TRUE)

# Extraer SOLO el efecto de tratamiento (controlando por hora)
res <- topTable(fit, coef = "treatment_factorvir",
                number = Inf, adjust.method = "BH")

# Diagnóstico
par(mfrow = c(1, 2))
hist(res$P.Value, breaks = 50, col = "#5DCAA5",
     main = "P-values (modelo con hora)",
     xlab = "P-value")
hist(res$logFC, breaks = 50, col = "#7F77DD",
     main = "logFC distribution",
     xlab = "logFC")
par(mfrow = c(1, 1))

#--------------------------------------------------
# RESULTADOS FINALES — umbral ajustado para n pequeño
#--------------------------------------------------

# P.Value < 0.05 sin corrección (exploratorio)
sig_raw <- res[res$P.Value < 0.05 & abs(res$logFC) > 1, ]
cat("Genes con raw P < 0.05 y |logFC| > 1: ", nrow(sig_raw), "\n")

cat("\n--- Con raw P < 0.05 y |logFC| > 1 ---\n")
cat("Más altos en vir: ", nrow(sig_raw[sig_raw$logFC > 1,  ]), "\n")
cat("Más altos en avr: ", nrow(sig_raw[sig_raw$logFC < -1, ]), "\n")

#--------------------------------------------------
# paso 4. VOLCANO PLOT 
#--------------------------------------------------

# Instalar
#install.packages("ggrepel")
#install.packages("caret")
#install.packages("randomForest")

library(ggrepel)
library(caret)
library(randomForest)

res$significancia <- "No significativo"
res$significancia[res$P.Value < 0.05 & res$logFC >  1] <- "Alto en vir"
res$significancia[res$P.Value < 0.05 & res$logFC < -1] <- "Alto en avr"

# Etiquetas para top 15 genes
res$gen     <- rownames(res)
top_labels  <- res[order(res$P.Value), ][1:15, ]

ggplot(res, aes(x = logFC, y = -log10(P.Value), color = significancia)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "No significativo" = "grey70",
    "Alto en vir"      = "#D85A30",
    "Alto en avr"      = "#5DCAA5"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data         = top_labels,
    aes(label    = gen),
    size         = 3,
    color        = "black",
    max.overlaps = 15
  ) +
  labs(
    title = "Volcano plot — avr vs vir (controlado por hpi)",
    x     = "log2 Fold Change (vir vs avr)",
    y     = "-log10(P-value)",
    color = NULL
  ) +
  theme_classic(base_size = 12)

#--------------------------------------------------
# paso 5. MODELO PREDICTIVO
#--------------------------------------------------

# Feature matrix: 185 genes DE (differencially expressed) como predictores
de_genes <- rownames(sig_raw)
X <- gene_set_log[, de_genes]
y <- factor(treatment_hours$tratamiento, levels = c("avr", "vir"))

cat("Dimensiones X:", nrow(X), "muestras x", ncol(X), "genes\n")
cat("Labels y:", as.character(y), "\n")

# LOOCV — única opción válida con n=12
ctrl <- trainControl(
  method          = "LOOCV",
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# --- Random Forest ---
set.seed(42)
rf_model <- train(
  x         = X,
  y         = y,
  method    = "rf",
  metric    = "ROC",
  trControl = ctrl
)

# --- SVM ---
set.seed(42)
svm_model <- train(
  x         = X,
  y         = y,
  method    = "svmRadial",
  metric    = "ROC",
  trControl = ctrl
)

# --- Lasso ---
set.seed(42)
lasso_model <- train(
  x         = X,
  y         = y,
  method    = "glmnet",
  metric    = "ROC",
  trControl = ctrl,
  tuneGrid  = expand.grid(
    alpha  = 1,
    lambda = seq(0.001, 0.1, length = 20)
  )
)

#--------------------------------------------------
# COMPARAR MODELOS
#--------------------------------------------------
resultados <- resamples(list(
  RandomForest = rf_model,
  SVM          = svm_model,
  Lasso        = lasso_model
))

cat("\n=== COMPARACIÓN DE MODELOS ===\n")
summary(resultados)

# Visualizar comparación
bwplot(resultados, metric = "ROC",
       main = "Comparación de modelos — AUC-ROC (LOOCV)")

#--------------------------------------------------
# FEATURE IMPORTANCE — top genes predictivos
#--------------------------------------------------
importancia <- varImp(rf_model)$importance
importancia$gen <- rownames(importancia)
importancia <- importancia[order(importancia$Overall, decreasing = TRUE), ]

cat("\n=== TOP 20 GENES MÁS PREDICTIVOS ===\n")
print(head(importancia, 20))

# Graficar top 20
ggplot(head(importancia, 20),
       aes(x = reorder(gen, Overall), y = Overall)) +
  geom_col(fill = "#5DCAA5", alpha = 0.85) +
  coord_flip() +
  labs(
    title = "Top 20 genes predictivos — Random Forest",
    x     = NULL,
    y     = "Importancia"
  ) +
  theme_classic(base_size = 12)

#--------------------------------------------------
# MÉTRICAS FINALES DEL MEJOR MODELO
#--------------------------------------------------
cat("\n=== MÉTRICAS RANDOM FOREST ===\n")
cat("AUC-ROC:     ", max(rf_model$results$ROC), "\n")
cat("Sensitivity: ", max(rf_model$results$Sens), "\n")
cat("Specificity: ", max(rf_model$results$Spec), "\n")

# Guardar resultados
write.csv(importancia,       "feature_importance_RF.csv")
write.csv(sig_raw,           "limma_DE_185genes.csv")
write.csv(res,               "limma_todos_los_genes.csv")

cat("\n¡Pipeline completo! Archivos guardados.\n")




