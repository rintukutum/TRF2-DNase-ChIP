##################
### original author: Parashar Dhapola
### modified by Rintu Kutum
##################
import numpy as np
import os
import glob
import sys
## all dnase bed files
dnase_bed = sorted(glob.glob('./data/recode-DNase/*.bed'))

chrom_indices = {}
chrom_lens = []
bin_size = 1000
with open('./data/hg19/hg19-chrom-sizes.txt') as h:
    for n, l in enumerate(h):
        c = l.rstrip('\n').split('\t')
        chrom_indices[c[0]] = n
        chrom_lens.append(int(c[1]))
chrom_array = [np.zeros(int(x / bin_size) + 1) for x in chrom_lens]

############## enrichment
for bed in dnase_bed:
	with open(bed) as h:
		for n, l in enumerate(h):
			print "\r%d\t%s" % (n,bed),
			sys.stdout.flush()
			c = l.rstrip('\n').split(' ')
			s = int(c[1]) / bin_size
			e = int(c[2]) / bin_size
			if s == e:
				e = e + 1
			chrom_array[chrom_indices[c[0]]][s:e] += int(c[6])  # signal value column
		print "\n"

# create the directory
if not os.path.exists('./data/DES/'):
	os.makedirs('./data/DES/')

for n, i in enumerate(chrom_array):
    for k, v in chrom_indices.items():
        if v == n:
            fn = './data/DES/%s.npy' % k
            print "Saving to file %s" % fn
            np.save(fn, i)
###############            
