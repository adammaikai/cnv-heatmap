hg19.genome.sizes <- "hg19.chrom.sizes.txt"
hg38.genome.sizes <- "hg38.chrom.sizes.txt"
mm10.genome.sizes <- "mm10.chrom.sizes.txt"

# Chromosome sizes, for binning and plotting coordinates
binGenome <- function(chromSizesPath){
  chromCoords <- read.csv(chromSizesPath, sep="\t", stringsAsFactors = FALSE, header = FALSE, col.names = c("chr", "size"))
  chromCoords$chr <- gsub("chr", "", chromCoords$chr)
  chromCoords$start <- cumsum(as.numeric(c(1,chromCoords$size[1:nrow(chromCoords)-1])))
  chromCoords$cumsum <- cumsum(as.numeric(chromCoords$size))
  chromCoords$color <- ifelse(c(1:length(chromCoords$chr))%%2 ==0, yes="black", no="gray")
  chromCoords
}

chromCoordsHg19 <- binGenome(hg19.genome.sizes)
chromCoordsHg38 <- binGenome(hg38.genome.sizes)
chromCoordsMm10 <- binGenome(mm10.genome.sizes)
save(chromCoordsHg19, chromCoordsHg38, chromCoordsMm10, file="CNVPlottingPrep.RData")
