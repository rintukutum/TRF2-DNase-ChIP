#############
# author : Rintu Kutum
#############
rm(list=ls())
source('./util-funcs.R')
trf2_bed <- read.table(
	'./data/TRF2/all_common_peaks_pankaj_sorted.bed',
	header=FALSE,
	sep='\t',
	stringsAsFactors=FALSE
)
######
flank_bed <- getFlankingRegions(
	bed = trf2_bed,
	bp_span = 25
)
write.table(
	flank_bed,
	'./data/TRF2/flank_50bp_TRF2_summit.bed',
	col.names=FALSE,
	sep='\t',
	row.names=FALSE,
	quote=FALSE
)
rm(list=ls())
source('./util-funcs.R')
trf2_flank <- read_bed('./data/TRF2/flank_50bp_TRF2_summit.bed')
	########## 1 mismatch
output.df <- getMotifMaps(
	motif = 'TTAGGG',
	bed = trf2_flank,
	seq = trf2_flank$V4,
	max.mismatch = 1
)
dir.create('./data/TRF2/motif-TTAGGG/',showWarnings=FALSE)
write.csv(
output.df,
'./data/TRF2/motif-TTAGGG/motif-TTAGGG-1-mismatch-TRF2-summit-peak.csv',
row.names=FALSE
)
########## perfect match
output.df <- getMotifMaps(
	motif = 'TTAGGG',
	bed = trf2_flank,
	seq = trf2_flank$seq,
	max.mismatch = 0
)
write.csv(
output.df,
'./data/TRF2/motif-TTAGGG/motif-TTAGGG-NO-mismatch-TRF2-summit-peak.csv',
row.names=FALSE
)
######### get the tags 

################# 
rm(list=ls())
motif_mapped <- read.csv(
'./data/TRF2/motif-TTAGGG/motif-TTAGGG-1-mismatch-TRF2-summit-peak.csv',
stringsAsFactors=FALSE
)
coord_hit <- motif_mapped
idx <- coord_hit$hit >= 4
coord_hit_idx <- plyr::dlply(coord_hit[idx,],'idx')
source('./util-funcs.R')
hits.status <- lapply(coord_hit_idx,getConsHits)
hits.status.df <- na.omit(plyr::ldply(hits.status))[,-1]


