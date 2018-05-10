#############
# author : Rintu Kutum
#############
# Libraries
library('Biostrings')
library('GenomicRanges')
### progress bars
getPB <-function(msg='', n){
	pb <- progress::progress_bar$new(
		format = paste("  ",msg," [:bar] :percent eta: :eta", sep =''),
		clear = FALSE,
		width = 60,
		total = n
	)
	return(pb)
}
### create the flanking regions from the mid peaks
getFlankingRegions <- function(
	bed,
	bp_span
){
	flank_bed <- bed
	p.bar <- getPB(
		msg = 'Processing flanking regions',
		n  = nrow(bed)
	)
	seq <- c()
	for(i in 1:nrow(bed)){
		x <- bed[i,]
		p.bar$tick()
		chr <- as.character(x[1])
		start <- as.numeric(x[2])
		end <- as.numeric(x[3])
		mp <- round((end + start)/2)
		s <- mp - bp_span
		flank_bed[i,2] <- s
		e <- mp + bp_span
		flank_bed[i,3] <- e
		seq[i] <- as.character(
			getSeq(
				BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
				names=chr,
				start=s,
				end=e
			)
		)
	}
	flank_bed$seq <- seq
	return(flank_bed)
}
#------------
# get the count of motifs mapped to a given sequence
getCountQuery2Flank <- function(
	motif,
	seq,
	max.mismatch
){
	mapped <- c()
	p.bar <- getPB(
		msg = 'Counting motifs',
		n  = length(seq)
	)
	for(i in 1:length(seq)){
		p.bar$tick()
		mapped[i] <- countPattern(
			motif,
			seq[i],
			max.mismatch = max.mismatch
		)
	}
	return(mapped)
}
#-----------
# get the coordinates where motif is present in the sequence
# 
getMotifMaps <- function(
	motif,
	bed,
	seq,
	max.mismatch = 1
){
	map.count <- getCountQuery2Flank(
		motif = motif,
		seq = seq,
		max.mismatch = max.mismatch
	)	
	idx <- which(map.count != 0)
	output <- list()
	p.bar <- getPB(
		msg = 'Mapping motifs',
		n  = length(idx)
	)
	for(i in 1:length(idx)){
		p.bar$tick()
		hit <- matchPattern(
			motif,
			seq[idx[i]],
			max.mismatch=1
		)
		output[[i]] <- cbind(
			data.frame(hit@ranges),
			query = motif,
			flank.seq = as.character(hit@subject),
			chr = bed[idx[i],1],
			flank.start = bed[idx[i],2],
			flank.end = bed[idx[i],3],
			hit = length(hit@ranges),
			idx = idx[i]
		)
	}
	output.df <- plyr::ldply(output)
	return(output.df)
}
#-----------
# get consensus 4 motifs
getConsHits <- function(
	coord,
	hits = 4
){
	cons <- list()
	z <- 1
	cons.status <-c()
	for(i in 1:nrow(coord)){
		if(i == 1){
			cons.status[i] <- FALSE
			e = coord[i,'end']+1
			score <- 0
			cons.idx <- c()
		}else{
			s = coord[i,'start']
			cons.status[i] <- FALSE
			if(e == s){
				# score
				score <- score + 1
				cons.idx[score] <- i
				if(score == hits){
					cons.status[i] <- TRUE
					cons[[z]] <- cons.idx
					z <- z + 1
				}
			}else{
				cons.status[i] <- FALSE
				score <- 0
				cons.idx <- c()
				cons.idx[score] <- i
			}
			e = coord[i,'end']+1
		}
	}
	if(any(cons.status)){
		cons.coord <- coord[unlist(cons),]
	}else{
		cons.coord <- coord[1,]
		cons.coord[1,] <- NA
	}
	return(cons.coord)
}
#---------------
# if start is greater than end, then swap
checkSEswap <- function(s,e){
	if(s > e){
		sd <- e
		ed <- s
		# swap
		s <- sd
		e <- ed
	}
	return(c(s=s,e=e))
}
#--------------
# find common peaks against our TRF2 peaks
findCommonPeaks.TRF2 <- function(
	query,
	search
){
	# format
	colnames(query) <- colnames(search) <- c('chr','start','end')
	# sorting
	for(i in 1:nrow(query)){
		s <- query[i,2]
		e <- query[i,3]
		x <- checkSEswap(s,e)
		query[i,2:3] <- x
	}
	query.s <- plyr::ddply(query, 'chr', function(x){x[order(x$start),]})
	search.s <- plyr::ddply(search, 'chr', function(x){x[order(x$start),]})
	q.gr <- makeGRangesFromDataFrame(query.s)
	s.gr <- makeGRangesFromDataFrame(search.s)
	sq_overlap <- findOverlaps(s.gr,q.gr)
	sq_count <- countOverlaps(s.gr,q.gr)
	#
	q.coord <- q.gr[sq_overlap@to,]
	q.ov <- data.frame(q.coord)[,1:4]
	colnames(q.ov)[1] <- 'chr'
	colnames(q.ov)[-1] <- paste(colnames(q.ov)[-1],'-q',sep='')
	s.coord <- s.gr[sq_overlap@from,]
	s.ov <- data.frame(s.coord)[,2:4]
	colnames(s.ov) <- paste(colnames(s.ov),'-s',sep='')
	if(nrow(s.ov) == 0){
		s.ov[1,] <- NA
		q.ov[1,] <- NA
	}
	output <- data.frame(
		q.ov,
		s.ov
	)
	return(output)
}
#########
# include only knowm 
read_bed <- function(
	bed.file,
	delim = '\t',
	header = FALSE
){
	bed <- read.table(
		file = bed.file,
		header = header,
		sep = delim,
		stringsAsFactors = FALSE
		)
	idx.chr1_24 <- grep('chr[1-9|X|Y]',bed[,1])
	return(bed[idx.chr1_24,])
}
