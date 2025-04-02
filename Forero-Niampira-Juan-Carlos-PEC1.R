# creacion directorios
dir.create("datos")
dir.create("resultados")


#cargar datos

metaDatos <- read.csv("datos/DataInfo_S013.csv", row.names = 1)
datos <- read.csv("datos/DataValues_S013.csv", row.names = 1)
class(datos)

infoSamples <- datos[, 1:5]
if (ncol(datos) == 695) datos <- datos[, -c(1:5)]
if (nrow(metaDatos) == 695) metaDatos <- metaDatos[-c(1:5), ]

#creación contenedor

if (!require(BiocManager)) install.packages("BiocManager")
if (!require(SummarizedExperiment)) BiocManager.install("SummarizedExperiment")

library(dplyr)
# Convertir las columnas 2 (SURGERY), 4 (GENDER) y 5 (Group) a factores
infoSamples <- infoSamples %>%
  mutate(SURGERY = as.factor(SURGERY), GENDER = as.factor(GENDER), Group = as.factor(Group))
# Crear el nuevo nombre para las filas
infoSamples <- infoSamples %>%
  mutate(RowName = paste0(toupper(substr(SURGERY, 1, 1)), "_", Group, "_", SUBJECTS))
rownames(infoSamples) <- infoSamples$RowName
infoSamples <- infoSamples %>%
  select(-RowName)
head(infoSamples)

#conversión matriz

matDatos <- as.matrix(datos)
if (sum(rownames(matDatos) != infoSamples$SUBJECTS) == 0) rownames(matDatos) <- rownames(infoSamples)

mySE <- SummarizedExperiment(assays = list(rawValues = matDatos), rowData = infoSamples,
                             colData = metaDatos)
show(mySE)

#análisis exploratorio

data <- rowData(mySE)[, c("SURGERY", "AGE", "GENDER", "Group")]
# Convertir a data.frame para facilidad de manipulación
data <- as.data.frame(data)
# Calcular tablas de frecuencias para las columnas 1 (SURGERY), 3 (GENDER) y 4
# (Group)
(freq_surgery <- table(data$SURGERY))

(freq_gender <- table(data$GENDER))
(freq_group <- table(data$Group))
(freq_surgery_Group <- table(data$SURGERY, data$Group))
(freq_surgery_Group_Gender <- table(data$SURGERY, data$Group, data$GENDER))
(summary_age <- summary(data$AGE))


# identificación y clasificación de los valores faltantes
if (!require(naniar)) install.packages("naniar")
library(naniar)
# Visualización de valores faltantes
dataf <- as.data.frame(t(assay(mySE, "rawValues")))
# rownames(dataf) <- rownames(rawVals) colnames(dataf) <- colnames(rawVals)
library(ggplot2)
# Visualización etiquetas
vis_miss(dataf) +
  theme_dark() +
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, color = "black"))

n_miss(dataf)
n_complete(dataf)
prop_miss(dataf)
prop_complete(dataf)
pct_miss(dataf)
pct_complete(dataf)

rawValues <- assay(mySE, "rawValues")
rawNoMissings <- rawValues
rawNoMissings[is.na(rawNoMissings)] <- 1
sum(is.na(rawNoMissings)
    
assays(mySE)$rawNoMissings <- rawNoMissings
show(mySE)
    
# Distribución muestras

boxplot(assay(mySE, "rawNoMissings"), col = rainbow(ncol(mySE), alpha = 0.6), ylab = "Concentraciones", cex.lab = 0.7, horizontal = FALSE, cex.axis = 0.3, las = 2, main = "Distribucion de valores por metabolito", cex.main = 0.7, border = "blue")
boxplot(t(assay(mySE, "rawNoMissings")), col=heat.colors(ncol(assay(mySE))), xlab = "Muestra", ylab = "Concentraciones", cex.lab = 0.8, horizontal = FALSE, cex.axis = 0.6, las = 2, main = "Distribucion de los valores", cex.main = 0.8, border = "blue")

#trasformación dtos

logNoMissings <- log(assay(mySE, "rawNoMissings"))
assays(mySE)$logNoMissings <- logNoMissings
show(mySE)

boxplot(assay(mySE, "logNoMissings"), ylab = "logConcentraciones", cex.lab = 0.7,
        horizontal = FALSE, cex.axis = 0.3, las = 1, main = "Distribucion de los logvalores por metabolito",
        cex.main = 0.8, border = "blue")

boxplot(t(assay(mySE, "logNoMissings")), col=heat.colors(ncol(assay(mySE))), xlab = "Muestra", ylab = "log-Concentraciones",
        cex.lab = 0.7, horizontal = FALSE, cex.axis = 0.5, las = 1, main = "Distribucion de los log", cex.main = 0.7)

#heatmap

if (!require(pheatmap)) install.packages("pheatmap")
library(pheatmap)
# Heatmap sin clustering
pheatmap(t(assay(mySE, "rawNoMissings")), cluster_rows = FALSE, cluster_cols = FALSE, main = "Heatmap sin Clustering", fontsize_row = 4, fontsize_col = 7, color = colorRampPalette(c("white", "red"))(200))

# Heatmap sin clustering
pheatmap(t(assay(mySE, "rawValues")), cluster_rows = FALSE, cluster_cols = TRUE, main = "Heatmap sin Clustering", fontsize_row = 4, fontsize_col = 7, color = colorRampPalette(c("white", "red"))(200))

#PCA

if (!require(PCAtools)) BiocManager::install("PCAtools", dep = TRUE)

library(PCAtools)
matriz <- t(assay(mySE, "rawNoMissings"))
pclas <- pca(matriz, metadata = rowData(mySE), removeVar = 0.1)
class(pclas)

screeplot(pclas, axisLabSize = 10, titleLabSize = 1)
biplot(pclas, showLoadings = FALSE, labSize = 2, pointSize = 2, axisLabSize = 6, title = "PCA datos" , sizeLoadingsNames = 2, colby = "SURGERY", hline = 0, vline = 0, legendPosition = "right")
biplot(pclas, showLoadings = FALSE, labSize = 2, pointSize = 2, axisLabSize = 6, ellipse = TRUE,
       title = "PCA de los datos", sizeLoadingsNames = 3, colby = "SURGERY", hline = 0,
       vline = 0, legendPosition = "right")

biplot(pclas, showLoadings = TRUE, labSize = 2, pointSize = 1, axisLabSize = 6, title = "PCA datos", sizeLoadingsNames = 2, colby = "SURGERY", lab = NULL, hline = 0, vline = 0, legendPosition = "right")


    