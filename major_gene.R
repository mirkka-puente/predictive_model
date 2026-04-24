# Instalar paquetes
#install.packages("readxl")
#install.packages("tidyverse")

# Cargar librerías
library(readxl)
library(tidyverse)
library(dplyr)

#--------------------------------------------------
# 1. Leer archivo Excel
#--------------------------------------------------
archivo <- "S3_FPKM.xlsx"

# Ver hojas disponibles
excel_sheets(archivo)

# Leer hoja específica
dataset <- read_excel(archivo, sheet = "IQ.OWLS Genes")
dt <- as.data.frame(t(dataset))
colnames(dt) <- as.character(dt[1, ])
dt <- dt[-c(1,20,21),]


# Ver primeras filas
head(dt)

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




