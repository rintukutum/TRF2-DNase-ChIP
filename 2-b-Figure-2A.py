##################
### original author: Parashar Dhapola
### modified by Rintu Kutum
##################
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.stats import sem
import scipy
import os

#######################################
binned_trf2_peak_spans = {}
with open('./data/TRF2/all_common_peaks_pankaj_sorted.bed') as h:
    for n, l in enumerate(h):
        c = l.rstrip('\n').split('\t')
        s = int(c[1]) // 1000
        e = int(c[2]) // 1000
        if s == e:
            e+=1
        if c[0] not in binned_trf2_peak_spans:
            binned_trf2_peak_spans[c[0]] = []
        binned_trf2_peak_spans[c[0]].append((s,e))

binned_random_spans = {}
np.random.seed(1261)
for chrom in binned_trf2_peak_spans:
    min_val = min(binned_trf2_peak_spans[chrom][0])
    max_val = max(binned_trf2_peak_spans[chrom])[1]
    binned_random_spans[chrom] = []
    for i in binned_trf2_peak_spans[chrom]:
        rand_val = np.random.randint(min_val, max_val)
        binned_random_spans[chrom].append((rand_val, rand_val+i[1]-i[0])) # keeping span sizes equivalent

dnase_enrichement_in_trf2_peaks = []
dnase_enrichement_in_random_peaks = []
dnase_trf2_peaks_chr = {}
dnase_random_peaks_chr = {}
chroms = ['chr%d' % x for x in range(1,23)] + ['chrX', 'chrY']
for chrom in chroms:
    dnase_signal_array = np.load('./data/DES/%s.npy' % chrom)
    dnase_trf2_peaks_chr[chrom] = []
    dnase_random_peaks_chr[chrom] = []
    for i in range(len(binned_trf2_peak_spans[chrom])):
        trf2_span = binned_trf2_peak_spans[chrom][i]
        random_span = binned_random_spans[chrom][i]
        dnase_enrichement_in_trf2_peaks.append(sum(dnase_signal_array[trf2_span[0]:trf2_span[1]]))
        dnase_enrichement_in_random_peaks.append(sum(dnase_signal_array[random_span[0]:random_span[1]]))
        ##
        dnase_trf2_peaks_chr[chrom].append(sum(dnase_signal_array[trf2_span[0]:trf2_span[1]]))
        dnase_random_peaks_chr[chrom].append(sum(dnase_signal_array[random_span[0]:random_span[1]]))
dnase_enrichement_in_trf2_peaks = np.array(dnase_enrichement_in_trf2_peaks)
dnase_enrichement_in_random_peaks = np.array(dnase_enrichement_in_random_peaks)


##################
# create the directory
if not os.path.exists('./figures'):
	os.makedirs('./figures')
#############
fig = plt.figure(figsize=(5,3))
gs = mpl.gridspec.GridSpec(100,100)
fontsize=10

ax = fig.add_subplot(gs[:88, :])
a = np.log2(dnase_enrichement_in_trf2_peaks+1)
b = np.log2(dnase_enrichement_in_random_peaks+1)
ax.hist(a, bins=100, color='crimson', alpha=0.6,
        edgecolor='none', label='TRF2 peaks')
ax.hist(b, bins=100, color='darkgrey', alpha=0.9,
        edgecolor='none', label='Random peaks')
ax.set_ylabel('Frequency of peaks', fontsize=fontsize)
ax.set_yticklabels([0, 100, 200, 300, 400, 500, 600], fontsize=fontsize)
ax.set_xticklabels([x for x in range(0,18,2)], fontsize=fontsize)
ax.legend(frameon=False, loc=(0.23, 0.75), fontsize=fontsize)
ax.set_xlim((0,17))
ax.set_ylim((0,600))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax = fig.add_subplot(gs[97:, :])
sns.heatmap([[x for x in range(0,100)], [x for x in range(0,100)]], ax=ax, cbar=False, cmap='Greys')
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('Log2 (DNase enrichment signal (DES))', fontsize=fontsize)

ax = fig.add_subplot(gs[5:55, 75:98])
a = np.log2(dnase_enrichement_in_trf2_peaks+1)
b = np.log2(dnase_enrichement_in_random_peaks+1)
pval = 'p = %.2E' % (scipy.stats.mannwhitneyu(a,b))[1]
boxprops = dict(linestyle='-', linewidth=2, color='black')
medianprops = dict(linestyle='-', linewidth=2, color='crimson')
whiskerprops = dict(linestyle='--', linewidth=1.5, color='black')
ax.boxplot([[], [], a, [], [], b, [], []], sym='+', widths=2,
           boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops)
ax.set_xticks([3,6])
ax.set_xticklabels(['TRF2', 'Random'], fontsize=fontsize,
                   rotation=40, horizontalalignment='center')
ax.set_yticks([x for x in range(0,25,5)])
ax.set_yticklabels([x for x in range(0,25,5)], fontsize=fontsize)
ax.set_ylabel('Log2 (DES)', fontsize=fontsize)
ax.text(0.05, 1, pval.replace('E', 'e'), transform=ax.transAxes, fontsize=fontsize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.savefig('./figures/Figure-2A-dnase_histogram-final.png', dpi=300)
####
out = []
header = ['chr','trf2_DES','random_DES']
out.append('\t'.join(header))
pvals = []
p_header = ['chr','p-value']
pvals.append('\t'.join(p_header))
for chrom in dnase_trf2_peaks_chr:	
	a = np.array(dnase_trf2_peaks_chr[chrom])
	b = np.array(dnase_random_peaks_chr[chrom])
	pval = scipy.stats.mannwhitneyu(a,b)[1]
	pvals.append('\t'.join([chrom,str(pval)]))
	status = pval < 0.05
	print '\r%s\t%.2E\t%s' % (chrom,pval,str(status))
	sys.stdout.flush()
	chr_out = []
	for i in range(len(a)):
		val = '\t'.join([chrom,str(a[i]),str(b[i])])
		chr_out.append(val)
	chr_val = '\n'.join(chr_out)
	out.append(chr_val)
output = '\n'.join(out)
outfile = open('./data/TRF2/Table-2A-figure.txt','w')
outfile.write(output)
outfile.close()

out_pval = '\n'.join(pvals)
out_pval_file = open('./data/TRF2/Table-2A-p-values.txt','w')
out_pval_file.write(out_pval)
out_pval_file.close()




