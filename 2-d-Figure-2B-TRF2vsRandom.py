##################
### original author: Parashar Dhapola
### modified by Rintu Kutum
##################
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import ndimage
import math
#########################
chrom_indices = {}
chrom_lens = []
bin_size = 1000
with open('./data/hg19/hg19-chrom-sizes.txt') as h:
    for n, l in enumerate(h):
        c = l.rstrip('\n').split('\t')
        chrom_indices[c[0]] = n
        chrom_lens.append(int(c[1]))

##########
## create chrom for TRF2
trf2_chrom_array_c = [np.zeros(int(x / bin_size) + 1) for x in chrom_lens]
h = open('./data/TRF2/all_common_peaks_pankaj_sorted.bed','r')
for l in h:    
	c = l.rstrip('\n').split('\t')
	s = int(c[1]) // 1000
	e = int(c[2]) // 1000
	if s == e:
		e+=1
	trf2_chrom_array_c[chrom_indices[c[0]]][s:e] += 1

if not os.path.exists('./data/TRF2/trf2_binned_array_count/'):
	os.makedirs('./data/TRF2/trf2_binned_array_count/')
for n, i in enumerate(trf2_chrom_array_c):
    for k, v in chrom_indices.items():
        if v == n:
            fn1 = './data/TRF2/trf2_binned_array_count/%s.npy' % k
            print "Saving to file %s" % fn1
            sys.stdout.flush()
            np.save(fn1, i)

def get_binned_array(a, b):
    remainder = len(a) % b
    if remainder > 0:
        binned = a[:-remainder]
    else:
        binned = a.copy()
    binned = binned.reshape(len(a) // b, b).mean(axis=1)
    if remainder > 0:
        return np.append(binned, a[-remainder:].mean())
    return binned   

all_chrom_subtracted = np.array([])
points = []
bin_size = 400
np.random.seed(1261)
for chrom in [x for x in range(1,23)]+['X', 'Y']:
    #print '\r%s' % chrom, end='',
    print '\r%s' % chrom
    trf2_peaks = np.load('./data/TRF2/trf2_binned_array_count/chr%s.npy' % chrom)
    trf2_peaks = get_binned_array(trf2_peaks, bin_size)
    dnase = np.load('./data/DES/chr%s.npy' % chrom)
    dnase = get_binned_array(dnase, bin_size)
    dnase = np.log2(dnase+1)
    ####################
    specific = dnase.copy()
    specific[trf2_peaks == 0] = 0
    ####################
    chrom_len = len(trf2_peaks)
    random_peaks = []
    for i in range(1):
        got_nums = []
        rand_p= np.zeros(chrom_len)
        for j in range(len(np.nonzero(trf2_peaks)[0])):
            while True:
                rand_num = np.random.randint(chrom_len)
                if rand_num not in got_nums:
                    got_nums.append(rand_num)
                    break
            rand_p[rand_num] = 1
        random_peaks.append(rand_p)
    background = dnase.copy()
    idx = random_peaks[0] == 0
    background[idx] = 0
    
    if len(all_chrom_subtracted) == 0:
        all_chrom_subtracted = specific.copy()
        all_chrom_background = background.copy()
    else:
        all_chrom_subtracted = np.hstack((all_chrom_subtracted, specific))
        all_chrom_background = np.hstack((all_chrom_background, background))
    points.append(len(dnase))

ref_subtract = all_chrom_subtracted - all_chrom_background

#############################
#############################
#########################
breaks_pos = np.cumsum([0.0] + points)
breaks_pos = np.radians((breaks_pos/breaks_pos[-1])*360)
xticks = []
for i in range(len(breaks_pos)-1):
    xticks.append(((breaks_pos[i+1]-breaks_pos[i])/2)+breaks_pos[i])
polar_pos = np.radians(np.linspace(0,360,len(ref_subtract)))
##################
pos = ref_subtract.copy()
pos[pos < 0] = 0
pos[pos > 10] = 10
neg = ref_subtract.copy()
neg[neg > 0] = 0
neg[neg < -10] = -10

#########################
fig = plt.figure(figsize=(7, 7))
gs = mpl.gridspec.GridSpec(6, 6, hspace=0.5)
ax = fig.add_subplot(gs[:, :], projection='polar')
#######
smoothen_factor = 5
a = ndimage.gaussian_filter1d(pos, smoothen_factor)
b = ndimage.gaussian_filter1d(neg, smoothen_factor)
#######
track_pos = 7
ax.fill_between(polar_pos, track_pos, a + track_pos,
                color='crimson', alpha=0.7, edgecolor='none', label='DES in TRF2\npeaks only')
ax.fill_between(polar_pos, track_pos, b + track_pos,
                color='grey', alpha=0.7, edgecolor='none', label='DES in random\npeaks only')
for i in breaks_pos:
    ax.axvline(i, ymin=0.4, ymax=0.9, ls='-.', c='k', lw=1)
for i,j in zip([x for x in range(1,23)]+['X', 'Y'], xticks):
    rotation = math.degrees(j)
    if 270 > rotation > 90:
        rotation-=180
    ax.text(j, track_pos + 3, 'chr'+str(i), fontsize=10, color='grey', alpha=0.7,
            horizontalalignment='center', verticalalignment='center',  rotation=rotation)
ax.set_xticks(xticks)
ax.xaxis.grid(False)
ax.yaxis.grid(False)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_ylim((0,11))
ax.spines['polar'].set_visible(False)
ax.legend(loc='center', frameon=False, fontsize=11)
plt.savefig('./figures/Figure-2B-dnase_TRF2_vs_Random.png', dpi=300)

