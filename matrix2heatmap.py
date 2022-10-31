#!/usr/bin/env python3
import time
print("start: ", time.asctime(time.localtime(time.time())))
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
print("load package: ", time.asctime(time.localtime(time.time())))

# parse input options
parser = argparse.ArgumentParser(prog="matrix2heatmap.py",
	description="read the HiC link matrix of HiC-Pro output and draw contact heatmap among genomic bins in PDF format",
	epilog="Author: Sen Wang, wangsen1993@163.com, 2021/9/2")
parser.add_argument("bed", type=str, help="positions of genomic bins in BED format")
parser.add_argument("matrix", type=str, help="HiC link matrix among bins")
parser.add_argument("compress", type=int, default=1, help="number of 100-Kb bins to plot as a whole, used for ploting heatmap of large genomes (default 1)")
args = parser.parse_args(sys.argv[1:])
binsize = 100000 * args.compress

# read BED file and link bin id and contig id
bins = []
bin = []
ctgs = []
ctg = ""
with open(args.bed, "r") as f:
	for line in f:
		t = line.strip().split("\t")
		ctg = t[0] if ctg == "" else ctg
		bin.append(t[3])
		if t[0] != ctg:
			bin.pop()
			bins.append(bin)
			ctgs.append(ctg)
			ctg = t[0]
			bin = [t[3]]

bins.append(bin)
ctgs.append(ctg)
bin = []
bin2ctg = {}
step = 0
bins2 = {}
for i in range(0, len(ctgs)):
	bin = bins[i]
	ctg = ctgs[i]
	if len(bin) <= args.compress:
		bin2ctg[step] = ctg
		for j in bin:
			bins2[j] = step
		step += 1
	else:
		num = 0
		for j in bin:
			num += 1
			bins2[j] = step
			if num % args.compress == 0:
				step += 1
		if len(bin) % args.compress == 0:
			bin2ctg[step - 1] = ctg
		else:
			bin2ctg[step] = ctg
			step += 1

print("read BED: ", time.asctime(time.localtime(time.time())))
# read HiC contact matrix
bin_list = list(set(sorted(bins2.values())))
binmx = [[0 for i in bin_list] for j in bin_list]
with open(args.matrix, "r") as f:
	for line in f:
		t = line.strip().split("\t")
		if t[0] in bins2 and t[1] in bins2:
			binmx[bins2[t[0]]][bins2[t[1]]] += float(t[2])
			binmx[bins2[t[1]]][bins2[t[0]]] += float(t[2])

print("read matrix: ", time.asctime(time.localtime(time.time())))
# draw heatmap
binmx = np.asarray(binmx, dtype=float)
size = len(bin_list) / 72;
if size < 5:
	size = 5

plt.figure(1, figsize=(size * 1.365, size * 1.1))
plt.imshow(np.log2(binmx+1), cmap="YlOrRd", aspect="auto", interpolation="none", origin="lower")
plt.colorbar(orientation="vertical", shrink=0.5)
plt.title(str(binsize) + " per bin")
pos = 0
for i in sorted(bin2ctg.keys()):
	plt.axvline(x=int(i+1) - 0.5, ls="-", lw=0.1, color="black")
	plt.axhline(y=int(i+1) - 0.5, ls="-", lw=0.1, color="black")
	plt.text((pos + int(i+1)) / 2, 0, bin2ctg[i] + " ", fontsize=8, horizontalalignment="center", verticalalignment="top", rotation="vertical")
	plt.text(0, (pos + int(i+1)) / 2, bin2ctg[i] + " ", fontsize=8, horizontalalignment="right", verticalalignment="center")
	pos = i

plt.xticks([])
plt.yticks([])
ax = plt.gca()
ax.spines['right'].set_linewidth(0.1)
ax.spines['top'].set_linewidth(0.1)
ax.spines['left'].set_linewidth(0.1)
ax.spines['bottom'].set_linewidth(0.1)
print("draw heatmap: ", time.asctime(time.localtime(time.time())))
plt.savefig(args.bed + ".pdf", bbox_inches="tight")
print("save figure: ", time.asctime(time.localtime(time.time())))
