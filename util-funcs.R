#############
# Libraries
library('Biostrings')
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
