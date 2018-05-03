rm(list=ls())
hg19.info <- 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes'
hg19 <- readLines(hg19.info)
dir.create('./data/hg19/',showWarnings=FALSE)
writeLines(
	hg19[1:24],
	'./data/hg19/hg19-chrom-sizes.txt')
