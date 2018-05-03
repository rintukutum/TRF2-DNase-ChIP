rm(list=ls())
metafile.encode <- 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/files.txt'
gz.files <- list.files(
	path = './data/DNase/',
	pattern = '.gz',
	full.names = TRUE
	)
meta.info <- readLines(metafile.encode)
get.fileInfo <- function(x){
	x.parts <- strsplit(x, split='\t')[[1]]
	x.filename <- x.parts[1]
	x.otherinfo <- gsub(' ', '', strsplit(x.parts[2], split=';')[[1]])
	
	out <- sapply(x.otherinfo,function(x){strsplit(x,split='\\=')[[1]][2]})
	names(out) <- sapply(x.otherinfo,function(x){strsplit(x,split='\\=')[[1]][1]})
	out <- c(out,x.filename)
	names(out)[length(out)] <- 'file.name'
	return(out)
}
out.path = './data/recode-DNase/'
dir.create(out.path,showWarnings=FALSE)
library('progress')
pb <- progress_bar$new(
      format = paste0("Processing "," [:bar] :percent in :elapsed"),
      total = length(meta.info),
      clear = FALSE, width= 60)
for(i in 1:length(meta.info)){
	file.info <- get.fileInfo(meta.info[i])
	idx.file <- grep(file.info['file.name'],gz.files)
	dat <- read.table(gz.files[idx.file],
		header=FALSE,
		stringsAsFactors=FALSE)
	out.filename <- paste(
		out.path,
		file.info['cell'],
		'-',
		file.info['treatment'],
		'.bed',
		sep='')
	pb$tick()
	write.table(
		dat,
		out.filename,
		quote = FALSE,
		row.names = FALSE,
		col.names = FALSE
	)
}
