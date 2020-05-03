# CNV Plotting Functions
library(shiny)
library(data.table)
library(DT)
library(ComplexHeatmap)
library(circlize)
library(plyr)
options(shiny.maxRequestSize=1000*1024^2)


# Read peaks
readSeg <- function(file.path, chromCoords){
  seg <- data.frame(fread(file.path, stringsAsFactors = FALSE, sep="\t"))
  seg <- data.frame(Sample=seg[,1], chromosome=gsub("chr", "", seg[,2]), start=as.numeric(seg[,3]), end=as.numeric(seg[,4]), n_probes=seg[,5], segment_mean=seg[,6])
  seg$start <- as.numeric(lapply(1:nrow(seg), function(i) subset(chromCoords, chr==seg$chr[i])$start + seg$start[i]))
  seg$end <- as.numeric(lapply(1:nrow(seg), function(i) subset(chromCoords, chr==seg$chr[i])$start + seg$end[i]))
  seg[complete.cases(seg),]
}


# Transfer the values from data frame into the matrix
populateCNVMatrix <- function(chromCoords, coordMat, segDf){
  coordMat <- subset(coordMat, chr %in% segDf$chromosome)
  chromCoords <- subset(chromCoords, chr %in% segDf$chromosome)
  samples <- unique(segDf$Sample)
  cnvMat <- matrix(rep(rep(0,nrow(coordMat)),length(samples)),ncol=length(samples))
  colnames(cnvMat) <- samples
  populate.binned.data <- function(start, end, column, value){ cnvMat[(start/100000):(end/100000), column] <<- value}
  for(i in samples) { apply(subset(segDf, Sample == i)[,c(3, 4, 6)], 1 , function(j) populate.binned.data(j[1],j[2],i,j[3]))}
  row.labels <- unlist(lapply(chromCoords$chr, function(i) rep(i, nrow(subset(coordMat, chr==i)))))
  row.names(cnvMat) <- row.labels
  cnvMat
}

makeHeatmapAnnotations <- function(cnvMat, chromCoords, coordMat, metadata, column_anno){
  coordMat <- subset(coordMat, chr %in% row.names(cnvMat))
  chromCoords <- subset(chromCoords, chr %in% row.names(cnvMat))
  set.seed(1)
  if(!is.null(metadata)){
    metadata <- subset(metadata, Sample %in% colnames(cnvMat))
    metadata <- metadata[order(metadata[,column_anno]),]
    }
  row.labels <- row.names(cnvMat)
  rowAnnoDf <- data.frame(chr=row.labels)
  chromColorList <- list(chr=chromCoords$color)
  names(chromColorList[[1]]) <- chromCoords$chr
  chromLabels <- unlist(lapply(chromCoords$chr, function(i) c(rep("", floor(nrow(subset(coordMat, chr==i))/2)), i, rep("", ceiling(nrow(subset(coordMat, chr==i))/2)-1))))
  hmtx <- HeatmapAnnotation(text = row_anno_text(chromLabels), which = "row", annotation_width = max_text_width(chromLabels), annotation_height = max_text_height(chromLabels))
  hma <- HeatmapAnnotation(df = rowAnnoDf, col = chromColorList, which = "row", simple_anno_size= unit(2, "mm"), show_annotation_name = FALSE, show_legend = FALSE)
  annos <- list()
  annos$hmtx <- hmtx
  annos$hma <- hma
  if(column_anno == "None" || is.null(column_anno)){
    return(annos)
  } else {
    conditionDf <- data.frame(metadata[,column_anno], row.names=metadata$Sample)
    names(conditionDf) <- column_anno
    columnAnno <- HeatmapAnnotation(df = conditionDf, which = "column", show_annotation_name = FALSE, name = column_anno)
    annos$columnAnno <- columnAnno
    annos$colOrder <- metadata$Sample
    return(annos)
  }
}
# ComplexHeatmap
plotCNVHeatmap <- function(cnvMat, annos){
  if(is.null(annos$colOrder)){
    cnvMat <- cnvMat
  } else {
    cnvMat <- cnvMat[,as.character(annos$colOrder)]
  }
  hm <- Heatmap(cnvMat, bottom_annotation = annos$columnAnno,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                show_row_names = FALSE,
                col = colorRamp2(c(-1.5, 0, 1.5), c("blue3", "white", "red3")),
                name = "Copy Number Change")
  ht_list <- annos$hmtx + annos$hma + hm
  draw(ht_list)  
}

readMeta <- function(metadataFilePath, sep=sep){
  meta <- read.csv(metadataFilePath, sep=sep)
  names(meta)[1] <- "Sample"
  meta
}
