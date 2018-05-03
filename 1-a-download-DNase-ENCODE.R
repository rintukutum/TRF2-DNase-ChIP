rm(list=ls())
metafile.encode <- 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/files.txt'
download.file(url=metafile.encode,dest='./data/metafile-encode.txt',quiet=TRUE)
encode.url <- 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/'
####
#meta.data <- read.table('./data/metafile-encode.txt',header=FALSE,sep='\t',stringsAsFactors=FALSE)
meta.data <- read.table(metafile.encode,header=FALSE,sep='\t',stringsAsFactors=FALSE)
dir.create('./data/DNase/',showWarnings=FALSE)
for(i in 1:nrow(meta.data)){
	url2download <- paste(encode.url,meta.data[i,1],sep='')
	outfile <- paste('./data/DNase/',meta.data[i,1],sep='')
	download.file(
		url = url2download,
		dest = outfile,
		quiet = TRUE
	)
	message(paste(i,': ', outfile,' done!',sep=''))
}
