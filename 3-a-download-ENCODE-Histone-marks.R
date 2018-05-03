##################
### Author: Rintu Kutum
##################
rm(list=ls())
dir.create('./data/Histones/',showWarnings=FALSE)
cmd <- "curl -l ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/peaks/jan2011/histone_macs/optimal/hub/ > ./data/Histones/temp.txt"
system(cmd)
system("grep '.bb' ./data/Histones/temp.txt > ./data/Histones/encode-big-bed-files.txt")
###
histone.files <- readLines('./data/Histones/encode-big-bed-files.txt')
url <- 'http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/peaks/jan2011/histone_macs/optimal/hub/'
dir.create('./data/Histones/bb/',showWarnings=FALSE)
for(i in 1:length(histone.files)){
	url.file <- paste(url,histone.files[i],sep='')
	output.file <- paste('./data/Histones/bb/',histone.files[i],sep='')
	download.file(
		url = url.file,
		dest = output.file,
		quiet = TRUE
	)
	message(paste(i,'. ',output.file,sep=''))
}
dir.create('./soft',showWarnings=FALSE)
download.file(
	url = 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed',
	dest = './soft/bigBedToBed')
system('chmod 755 ./soft/bigBedToBed')
rm(list=ls())
histone.files <- list.files(
	path='./data/Histones/bb/',
	pattern='.bb',
	full.names=TRUE
)
bedFileName <- function(x,path='./data/Histones/bed/'){
	paste(path,'/',gsub('.bb','',strsplit(x,split='\\/')[[1]][6]),'.bed',sep='')
}
dir.create('./data/Histones/bed/',showWarnings=FALSE)
for(i in 1:length(histone.files)){
	output <- bedFileName(histone.files[i])
	cmd <- paste("./soft/./bigBedToBed ./",histone.files[i], ' ',output,sep='')
	system(cmd)
}

############### format bed files
rm(list=ls())
bed.files <- list.files(path='./data/Histones/bed/',pattern='.bed', full.names=TRUE)
cutStr <- function(x, mark){
	file.name <- gsub('StdAln.bed','',gsub('wgEncode','',strsplit(x,split='\\/')[[1]][6]))
	file.name <- gsub('UcdAln.bed','',file.name)
	file.name <- gsub(mark,'',file.name)
	file.name <- paste(mark,'_',paste(rev(strsplit(file.name,split='Histone')[[1]]),collapse = '_'),'.bed',sep='')
	file.name <- gsub('b_','_',file.name)
	return(file.name)
}

marks <- c('H3k4me1', 'H3k4me2', 'H3k4me3', 'H3k9ac', 'H3k27ac',
	'H4k20me1','H3k9me1', 'H3k9me3', 'H3k27me3','H2az',
	'H3k36me3', 'H3k79me2')

dir.create('./data/Histones/bed_formatted/',showWarnings=FALSE)
for(i in 1:length(marks)){
		mark.files <- bed.files[grep(marks[i],bed.files)]
		for(j in 1:length(mark.files)){
			out.filename <- as.character(cutStr(mark.files[j],mark = marks[i]))
			file.copy(from = mark.files[j], to = paste('./data/Histones/bed_formatted/',out.filename,sep=''))
		}		
}
