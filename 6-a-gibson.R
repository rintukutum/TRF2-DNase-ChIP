#######################
# author : Rintu Kutum
#######################
rm(list=ls())
source('./util-funcs.R')
mono_bed <- read_bed(
	'./data/Gibson-TRF2/gibson-TRF2-monoclonal-peak.bed'
)
poly_bed <- read_bed(
	'./data/Gibson-TRF2/gibson-TRF2-ployclonal-peak.bed'
)
trf2_bed <- read_bed(
	'./data/TRF2/all_common_peaks_pankaj_sorted.bed'
)
dir.create('./data/intersect-gibson',showWarnings=FALSE)
mono_common <- findCommonPeaks.TRF2(
	query = mono_bed,
	search = trf2_bed
)
poly_common <- findCommonPeaks.TRF2(
	query  = poly_bed,
	search = trf2_bed
)

write.csv(
	mono_common,
	'./data/intersect-gibson/gibson_mono_intersect.csv',
	row.names=FALSE
	)

write.csv(
	poly_common,
	'./data/intersect-gibson/gibson_poly_intersect.csv',
	row.names=FALSE
	)
#-----------

