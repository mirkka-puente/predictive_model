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
# 3. Chequeo — numero de columnas
#--------------------------------------------------
cat("Numero de muestras:", nrow(gene_set_log), "\n")
cat("Columnas en tratamientos y horas: ", length(treatment_hours), "\n")

#--------------------------------------------------
# 4. limma — avr vs vir
#--------------------------------------------------
#install.packages("BiocManager")
#BiocManager::install("limma")
library(limma)


# El "bloque" es el experimento biológico repetido en el tiempo
# Tienes 2 réplicas biológicas × 3 tiempos × 2 tratamientos = 12 muestras
# Las réplicas biológicas son el bloque

bloque <- factor(c(
  1, 1,   # avr/vir rep1 - 1hpi
  1, 1,   # avr/vir rep1 - 6hpi
  1, 1,   # avr/vir rep1 - 12hpi
  2, 2,   # avr/vir rep2 - 1hpi
  2, 2,   # avr/vir rep2 - 6hpi
  2, 2    # avr/vir rep2 - 12hpi
))

# Diseño simple avr vs vir
treatment_factor <- factor(treatment_hours$tratamiento, levels = c("avr", "vir"))
design <- model.matrix(~ treatment_factor)
colnames(design) <- c("Intercept", "vir_vs_avr")

print(design)

# Paso 1: estimar correlación entre réplicas del mismo bloque
corfit <- duplicateCorrelation(t(gene_set_log), design, block = bloque)
cat("Correlación entre bloques:", corfit$consensus.correlation, "\n")
# Valor entre 0.1–0.9 = bien. Negativo o >0.99 = revisar bloques

# Paso 2: fit con correlación
fit <- lmFit(t(gene_set_log), design,
             block = bloque,
             correlation = corfit$consensus.correlation)
fit <- eBayes(fit)

# Paso 3: resultados
res <- topTable(fit, coef = "vir_vs_avr",
                number = Inf, adjust.method = "BH")

# Verificar que adj.P.Val ya no está clavado
print(head(res[order(res$P.Value), 
               c("logFC", "AveExpr", "P.Value", "adj.P.Val")], 10))

# Histograma — ahora debe tener spike cerca de 0
hist(res$P.Value, breaks = 50, col = "#5DCAA5",
     main = "P-values con duplicateCorrelation",
     xlab = "P-value")

# Filtrar genes DE
sig <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ]
cat("Genes DE totales:  ", nrow(sig), "\n")
cat("Más altos en vir:  ", nrow(sig[sig$logFC > 1,  ]), "\n")
cat("Más altos en avr:  ", nrow(sig[sig$logFC < -1, ]), "\n")


