##################
### Author: Rintu Kutum
##################
############### format dist_json
rm(list=ls())
output.json <- list.files(path='./data/Histones/dist_jsons/',pattern='.json', full.names=TRUE)

cutStr <- function(x, mark){
	file.name <- gsub('StdAln.json','',gsub('wgEncode','',strsplit(x,split='\\/')[[1]][6]))
	file.name <- gsub('UcdAln.json','',file.name)
	file.name <- gsub(mark,'',file.name)
	file.name <- paste(mark,'_',paste(rev(strsplit(file.name,split='Histone')[[1]]),collapse = '_'),'.json',sep='')
	file.name <- gsub('b_','_',file.name)
	return(file.name)
}

marks <- c('H3k4me1', 'H3k4me2', 'H3k4me3', 'H3k9ac', 'H3k27ac',
	'H4k20me1','H3k9me1', 'H3k9me3', 'H3k27me3','H2az',
	'H3k36me3', 'H3k79me2')

dir.create('./data/Histones/dist_json_formatted',showWarnings=FALSE)

for(i in 1:length(marks)){
		mark.files <- output.json[grep(marks[i],output.json)]
		mark.format <- list()
		for(j in 1:length(mark.files)){
			out.filename <- as.character(cutStr(mark.files[j],mark = marks[i]))
			file.copy(from = mark.files[j], to = paste('./data/Histones/dist_json_formatted/',out.filename,sep=''))
		}		
}

