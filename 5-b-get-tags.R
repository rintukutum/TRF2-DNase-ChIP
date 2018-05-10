rm(list=ls())
source('./util-funcs.R')
bed.file <- './data/TRF2/flank_50bp_TRF2_summit.bed'
bed.dat <- read_bed(bed.file)
write.table(
	bed.dat[,1:3],
	'./data/TRF2/flank_50bp_TRF2_coord.bed',
	sep = '\t',
	col.names = FALSE,
	row.names = FALSE,
	quote = FALSE
)
#---------#
rm(list=ls())
source('./util-funcs.R')
bed.file <- './data/TRF2/flank_50bp_TRF2_coord.bed'
rep1.bam <- '../DNase-ENCODE-common-peaks/sorted_rep1_dedup.bam'
rep2.bam <- '../DNase-ENCODE-common-peaks/sorted_rep2_dedup.bam'
####
# rep1 coverage
getCoverageSamtools(
	bed.file = bed.file,
	bam.file = rep1.bam,
	output.name = './data/TRF2/flank_rep1.bed'
)
####
# rep2 coverage
getCoverageSamtools(
	bed.file = bed.file,
	bam.file = rep2.bam,
	output.name = './data/TRF2/flank_rep2.bed'
)
#############
rm(list=ls())
source('./util-funcs.R')
rep1_reads <- read_bed('./data/TRF2/flank_rep1.bed',)
rep2_reads <- read_bed('./data/TRF2/flank_rep2.bed')
reads_mean <- ((rep1_reads[,'V4'])/23 + (rep2_reads[,'V4'])/24)/2
#############
# flank 50bp reads
flank_50bp_reads <- data.frame(rep1_reads[,1:3],mean=reads_mean)
write.table(
	x = flank_50bp_reads,
	file = './data/TRF2/flank_50bp_TRF2_reads.bed',
	col.names = FALSE,
	row.names = FALSE,
	sep = '\t',
	quote = FALSE
)

