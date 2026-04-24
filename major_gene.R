# Instalar paquetes
install.packages("readxl")
install.packages("tidyverse")

# Cargar librerías
library(readxl)
library(tidyverse)

#--------------------------------------------------
# 1. Leer archivo Excel
#--------------------------------------------------
file_path <- "S3_FPKM.xlsx"

# Ver hojas disponibles
excel_sheets(file_path)

# Leer hoja específica (ajusta el nombre si es diferente)
genes_raw <- read_excel(file_path, sheet = "IQ.OWLS Genes")

# Ver primeras filas
head(genes_raw)

# Revisar estructura
str(genes_raw)

#--------------------------------------------------
# 2. Separar IDs de genes
#--------------------------------------------------

# Asumiendo que la primera columna contiene IDs tipo AT1G01010
gene_ids <- genes_raw[[1]]

# Remover primera columna para quedarte solo con expresión
expr_matrix <- genes_raw[, -1]

# Convertir a matriz numérica
expr_matrix <- as.data.frame(expr_matrix)

expr_matrix[] <- lapply(expr_matrix, as.numeric)

# Convertir a matrix
expr_matrix <- as.matrix(expr_matrix)

# Asignar nombres de filas
rownames(expr_matrix) <- gene_ids

# Revisar dimensiones
dim(expr_matrix)