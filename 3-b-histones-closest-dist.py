import pybedtools as pbt
import sys
import numpy as np
import json
import os

def make_simple_bed(pbt_file):
    bed = []
    pbt_open = open(pbt_file, 'r');
    for i in pbt_open:
        j = str(i).split('\t')[:3]
        bed.append("\t".join([j[0], j[1], j[2]]))
        
    bed = pbt.BedTool("\n".join(bed), from_string=True)
    return bed

def get_closest_dist_array(bed1, bed2, window_size, bin_size):
    dist = []
    closest = bed1.closest(bed2)
    for i in closest:
        j = str(i).split('\t')
        mid_mark = int(j[4]) + (int(j[5]) - int(j[4])) / 2
        d = (int(j[1]) - mid_mark)
        if abs(d) < window_size:
            d = d + window_size
            dist.append(d / bin_size)
	y = np.zeros(2 * window_size / bin_size)
    for i in dist:
        y[i] += 1
    return list(y)			

def get_peak_mids(bed_file):
    mid_peaks = []
    h = open(bed_file, 'r')
    for line in h:
		j = line.rstrip('\n').split('\t')
		mid = int(j[1]) + ((int(j[2]) - int(j[1])) / 2)
		mid_peaks.append("\t".join([str(j[0]), str(mid), str(mid + 1)]))
    return mid_peaks

bed_file = './data/TRF2/all_common_peaks_pankaj_sorted.bed'
peaks = get_peak_mids(bed_file)
peaks_bed = pbt.BedTool("\n".join(peaks), from_string=True).sort()
#-------------
# all histone files
files = os.listdir('./data/Histones/bed/')
histone_files = []
for f in files:
	histone_files.append(''.join(['./data/Histones/bed/',f]))


# create the directory
if not os.path.exists('./data/Histones/dist_jsons'):
	os.makedirs('./data/Histones/dist_jsons')
seed = 1
for h in histone_files:
	histones = make_simple_bed(h)
	window = 10000
	binsize = 200
	closest_dist = get_closest_dist_array(peaks_bed, histones, window, binsize)
	shuffle_dist = []
	print(''.join(['Histone: ',str(seed)]))
	for i in range(100):
		print "\r%d" % int(i+1),
		shuffled_histones = histones.shuffle(genome = 'hg19',chrom = True,seed = seed + i)
		shs = shuffled_histones.sort()
		shuffle_dist.append(get_closest_dist_array(peaks_bed, shs, window, binsize))
	seed = seed + 1
	dist_out = h.replace('.bed','.json').replace('bed','dist_jsons')
	with open(dist_out, 'w') as OUT:
		json.dump({'closest_dist': closest_dist, 'shuffle_dist': shuffle_dist}, OUT, indent=2)
