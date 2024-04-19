library(NOISeq)
library(GO.db)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(gplots)
library(plotly)
library(topGO)
library(Rgraphviz)
library(limma)
library(gProfileR)
library(biomaRt)
library(EnhancedVolcano)
library(pheatmap)

#FILTRAR DESeq2 PARA LOG2(FC)

filter <- read.delim("deseq2_final.tabular", header=FALSE)
names(filter) <- c("Gene_ID","base mean","log2(FC)","StdErr","Wald-Stats","p-value","P-adj")

datos_ordenados <-  filter[order(-filter$`log2(FC)`),]
datos_ordenados <- datos_ordenados %>%
  filter(`base mean` > 200)

sobreexpresadas_control <- head(datos_ordenados,20)
sobreexpresadas_control

datos_filtrados_sobre <- subset(datos_ordenados,`log2(FC)` > 0.5)
datos_filtrados_sobre
datos_filtrados_infra <- datos_ordenados %>%
  filter(`log2(FC)` < -0.5) 
datos_filtrados_infra



#PASAR GENEID A SYMBOL

gene_ids <- as.character(datos_filtrados_sobre$Gene_ID) 
symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID")

datos_filtrados_sobre$symbol <- symbols
datos_filtrados_sobre$Geneid = NULL

symbol_index <- which(names(datos_filtrados_sobre) == "symbol")

datos_filtrados_sobre <- datos_filtrados_sobre[, c(symbol_index, setdiff(1:ncol(datos_filtrados_sobre), symbol_index))]
datos_filtrados_sobre



#GENERAR FICHERO NOISeq

file_names <- c("count_13.tabular", "count_16.tabular", "count_15.tabular", "counts_12.tabular", "counts_09.tabular", "count_11.tabular")

for (file in file_names) {
  table_data <- read.table(file, header = TRUE, row.names = NULL)
  assign(sub(".tabular", "", file), table_data)
}

unified_table_control <- merge(merge(count_11, counts_12, by = "Geneid"), counts_09, by = "Geneid")
unified_table_muta <- merge(merge(count_15, count_16, by = "Geneid"), count_13, by = "Geneid")

unified_table_control$Media_Expresion_control <- rowSums(unified_table_control[, -1], na.rm = TRUE) / 3
unified_table_muta$Media_Expresion_muta <- rowSums(unified_table_muta[, -1], na.rm = TRUE) / 3

medias_expresion <- merge(unified_table_control, unified_table_muta, by = "Geneid")
medias_expresion$X11.bam= NULL
medias_expresion$X12.bam= NULL
medias_expresion$X13.bam= NULL
medias_expresion$X15.bam= NULL
medias_expresion$X16.bam= NULL
medias_expresion$X09.bam= NULL
medias_expresion
medias_expr_ordenado <-  medias_expresion[order(-medias_expresion$Media_Expresion_control),]



#ANÁLISIS NOISeq

datos<- read.csv('medias_expr_ord.csv', header = TRUE, sep=",",stringsAsFactors = FALSE)

colnames(datos)<-(c("control","mutante","symbol"))
rownames(datos) <- datos$`symbol`
datos_reducidos <- as.data.frame(datos[1:2])


factores = data.frame(Muestras= c("control","mutante"))
misdatos <- NOISeq::readData(data = datos_reducidos,  factors = factores)
head(assayData(misdatos)$exprs)

datos.tidy <- datos_reducidos %>%  rownames_to_column("gene") 

genes.tidy <- pivot_longer(datos.tidy,cols = - gene, names_to = "sample", values_to = "counts")  %>% 
  group_by(gene)  %>% 
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() 

X11()
p_mean_scaled <- ggplot(genes.tidy, aes(x = log10(gene_average), y = log10(gene_stdev))) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  labs(x = "Conteo medio por gen (log10 scale)",
       y = "Desviación típica de los conteos por gen (log10 scale)") +
  ggtitle("Relación entre la media y la desviación estándar del conteo por gen\n(sin estabilización de varianza)")
p_mean_scaled


comparacion=c("control","mutante")
mynoiseq1 = noiseq(misdatos, conditions=comparacion,k = 0.5, norm="n", factor = "Muestras",nss = 5, lc = 1,pnr = 0.2, v = 0.02, replicates = "technical")
mynoiseq1.deg = degenes(mynoiseq1, q = 0.9, M = NULL)
mynoiseq1.deg1 = degenes(mynoiseq1, q = 0.9, M = "up")
mynoiseq1.deg2 = degenes(mynoiseq1, q = 0.9, M = "down")

X11()
DE.plot <- DE.plot(mynoiseq1, q = 0.9, graphic = "expr", log.scale = TRUE)
DE.plot

X11()
MD.plot <- DE.plot(mynoiseq1, q = 0.9, graphic = "MD")
MD.plot

prob <- 1-mynoiseq1@results[[1]][,5]
datosVolcano <-cbind(as.data.frame(mynoiseq1@results[[1]][,3]), prob)
row.names(datosVolcano)<- row.names(mynoiseq1@results[[1]])
colnames(datosVolcano)<-c("LogFC","prob")

EnhancedVolcano(datosVolcano,
                lab = rownames(datosVolcano),
                x = "LogFC",
                y = 'prob',
                xlim = c(-6, 6),
                ylim = c(0, 5.5), 
                pCutoff = 0.2,
                pointSize = 1.0,
                FCcutoff = 1,
                cutoffLineCol = 'black',
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8
)



# RECUPERAR BIOTIPO  DE BIOMART
genes_diff <- row.names(mynoiseq1.deg)

gene_diff_convert <- bitr(genes_diff, fromType = "SYMBOL", 
                          toType = c("ENSEMBL"),
                          OrgDb = org.Hs.eg.db)
Homomart = useMart(biomart = "ENSEMBL_MART_ENSEMBL")
homoData <- useDataset(mart = Homomart, dataset = "hsapiens_gene_ensembl")



# RECUPERAR LAS ANOTACIONES GO
homo_data <-getBM(attributes = c("ensembl_gene_id","gene_biotype"), 
                  filters="ensembl_gene_id", 
                  values=gene_diff_convert$ENSEMBL, 
                  mart=homoData)
View(homo_data)
