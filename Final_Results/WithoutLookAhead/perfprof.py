#!/usr/bin/env python2

from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import sys

from optparse import OptionParser

# parameters
defLW = 1  # default line width
defMS = 2  # default marker size
marker = 'o'
colors = ['r', 'b', 'y', 'g', 'm', 'c', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']


class CmdLineParser(object):
	def __init__(self):
		self.parser = OptionParser(usage='usage: python2 perfprof.py [options] cvsfile.csv outputfile.pdf')
		# default options
		self.parser.add_option("-D", "--delimiter", dest="delimiter", default=None, help="delimiter for input files")
		self.parser.add_option("-M", "--maxratio", dest="maxratio", default=4, type=float, help="maxratio for perf. profile")
		self.parser.add_option("-S", "--shift", dest="shift", default=0, type=float, help="shift for data")
		self.parser.add_option("-L", "--logplot", dest="logplot", action="store_true", default=False, help="log scale for x")
		self.parser.add_option("-T", "--timelimit", dest="timelimit", default=1e99, type=float, help="time limit for runs")
		self.parser.add_option("-P", "--plot-title", dest="plottitle", default=None, help="plot title")
		self.parser.add_option("-X", "--x-label", dest="xlabel", default='Cost Ratio', help="x axis label")
		self.parser.add_option("-B", "--bw", dest="bw", action="store_true", default=False, help="plot B/W")

	def addOption(self, *args, **kwargs):
		self.parser.add_option(*args, **kwargs)

	def parseArgs(self):
		(options, args) = self.parser.parse_args()
		options.input = args[0]
		options.output = args[1]
		return options


def readTable(fp, delimiter):
	"""
	read a CSV file with performance profile specification
	the format is as follows:
	ncols algo1 algo2 ...
	nome_istanza tempo(algo1) tempo(algo2) ...
	...
	"""
	firstline = fp.readline().strip().split(delimiter)
	ncols = int(firstline[0])
	assert(ncols <= len(colors))
	cnames = firstline[1:]
	rnames = []
	rows = []
	for row in fp:
		row = row.strip().split(delimiter)
		rnames.append(row[0])
		rdata = np.empty(ncols)
		for j in range(ncols):
			rdata[j] = float(row[j + 1])
		rows.append(rdata)
	data = np.array(rows)
	return (rnames, cnames, data)


def main():
	parser = CmdLineParser()
	opt = parser.parseArgs()
	print(opt)
	# read data
	rnames, cnames, data = readTable(open(opt.input, 'r'), opt.delimiter)
	nrows, ncols = data.shape
	# add shift
	#lo shift funziona in modo da aumentare tutti i tempi di una quantita' costante per tutti i sample in modo da non far pesare troppo computazioni che sono brevi in ogni caso.
	data = data + opt.shift
	
	# compute ratios
	# trova il minimo di tutte le righe
	minima = data.min(axis=1)
	ratio = np.zeros((nrows,ncols))

	# calcola il ratio in ogni riga 
	for j in range(ncols):
		ratio[:, j] = data[:, j] / minima
	
	# compute maxratio
	
	if opt.maxratio == -1:
		opt.maxratio = ratio.max()
	
	# any time >= timelimit will count as maxratio + bigM (so that it does not show up in plots)
	# penalizza quando timelimit scade in modo da non entrare mai nel rapporto considerato
	for i in range(nrows):
		for j in range(ncols):
			if data[i,j] >= opt.timelimit:
				ratio[i,j] = opt.maxratio + 1e6
	# sort ratios
	ratio.sort(axis=0)
	# plot first
	# evenly spaced values, uno per riga
	# corrisponde effettivamente all'asse y va da 0 a 1 e rappresenta la percentuale di istanze che sono in un ratio in base al valore delle ascisse 
	y = np.arange(nrows, dtype=np.float64) / nrows
	for j in range(ncols):
		options = dict(label=cnames[j],
				linewidth=defLW, linestyle='steps-post',
				marker=marker, markeredgewidth=defLW, markersize=defMS)
		# plt.step(ratio[:,j], y, label=cnames[j], linewidth=defLW, marker=markers[j], markersize=defMS)
		if opt.bw:
			options['markerfacecolor'] = 'w'
			options['markeredgecolor'] = 'k'
			options['color'] = 'k'
		else:
			options['color'] = colors[j]
		if opt.logplot:
			plt.semilogx(ratio[:, j], y, **options)
		else:
			plt.plot(ratio[:, j], y, **options)
	plt.axis([1, opt.maxratio, 0, 1])
	plt.legend(loc='lower right')
	if opt.plottitle is not None:
		plt.title(opt.plottitle)
	plt.xlabel(opt.xlabel)
	plt.savefig(opt.output)

if __name__ == '__main__':
	main()
