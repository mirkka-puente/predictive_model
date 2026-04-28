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
  theme_classic(base_size = 13)

