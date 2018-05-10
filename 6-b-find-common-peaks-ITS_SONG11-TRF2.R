#######################
# author : Rintu Kutum
#######################
rm(list=ls())
source('./util-funcs.R')
its_bed <- read.table(
	'./data/TRF2/ITS.bed',
	header=FALSE,
	sep='\t',
	stringsAsFactors=FALSE)
song11_bed <- read.table(
	'./data/TRF2/SONG11.bed',
	header=FALSE,
	sep='\t',
	stringsAsFactors=FALSE)
trf2_bed <- read.table(
	'./data/TRF2/all_common_peaks_pankaj_sorted.bed',
	header=FALSE,
	sep='\t',
	stringsAsFactors=FALSE
)
dir.create('./data/intersect-ITS-SONG11',showWarnings=FALSE)
its_common <- findCommonPeaks.TRF2(its_bed,trf2_bed)
song11_common <- findCommonPeaks.TRF2(song11_bed,trf2_bed)

write.csv(
	its_common,
	'./data/intersect-ITS-SONG11/ITS-TRF2-intersect.csv',
	row.names=FALSE
	)

write.csv(
	song11_common,
	'./data/intersect-ITS-SONG11/SONG11-TRF2-intersect.csv',
	row.names=FALSE
	)
#-----------

