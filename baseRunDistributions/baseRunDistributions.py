#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: baseRunDistributions.py      -i FASTA [-c INT] [--include_n] [-h|--help]

    Options:
        -h --help                           show this
        -i, --infile <FASTA>                FASTA file.
        -c, --run_min <INT>                 Minimum length of nucleotide run to count [default: 10]
        --include_n                         Include N's in count
"""

from __future__ import division
from docopt import docopt
import re

import sys
from itertools import groupby
from collections import Counter
from os.path import basename, isfile, abspath, join

from numpy import linspace, std, median
import matplotlib.colors as colors
import matplotlib as mat
mat.use('agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
plt.style.use('ggplot')

def progress(iteration, steps, max_value):
    if int(iteration == max_value):
        sys.stdout.write('\r')
        print "[PROGRESS]\t: %d%%" % (100)
    elif int(iteration) % int(steps) == 0:
        sys.stdout.write('\r')
        print "[PROGRESS]\t: %d%%" % (float(int(iteration)/int(max_value))*100),
        sys.stdout.flush()
    else:
        pass

def read_fasta(infile):
    if not isfile(infile):
         sys.exit("[ERROR]\t- %s is not a file")
    with open(infile) as fh:
        header, seqs = '', []
        for l in fh:
            if l[0] == '>':
                if (header):
                    yield header, ''.join(seqs)
                header, seqs = l[1:-1].split()[0], [] # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        yield header, ''.join(seqs)

def count_base_runs(infile, run_min, include_n_flag, bases):
    base_run_count = {base : [] for base in bases}
    sequences = {}
    print "[STATUS]\t- Parsing %s." % (infile)
    for header, seq in read_fasta(infile):
        sequences[header] = seq
    print "[STATUS]\t- Getting counts from sequences."
    iteration = 0
    sequences_count = len(sequences)
    steps = int(sequences_count/100)
    for header, seq in sequences.items():
        iteration += 1
        for base, group in groupby(seq):
            if base in base_run_count:
                length_group_list = len(list(group))
                if length_group_list >= run_min:
                    base_run_count[base].append(length_group_list)
        progress(iteration, 100, sequences_count)
    base_run_counter = {base : Counter(base_run_count[base]) for base in base_run_count.keys()}
    if not (base_run_count):
        sys.exit("[STATUS]\t- No runs found.")
    return base_run_counter

def plot_base_run_count(base_run_counter, out_plot_f):
    f, ax = plt.subplots(figsize=(24, 12))
    idx = 1
    colours = {'G' : u'#608000', 'C' : u'#ff9900', 'A' : u'#e60000', 'T': u'#9999ff', 'N': u'#9ca4a5'}
    #colours = plt.cm.Set1(linspace(0, 1, len(base_run_count)))

    for base in base_run_counter:
        counter = base_run_counter[base]
        idx += 1
        x_values = []
        y_values = []
        for run_length, count in sorted(counter.items()):
            x_values.append(run_length)
            y_values.append(count)
        ax.plot(x_values, y_values, label=base, marker='o', linestyle = '--', linewidth = (5 - (idx*0.5)), markersize=(15 - idx), alpha=1.0, color =colours[base],  markerfacecolor=colours[base])

    plt.xlabel("Length of run of bases")
    plt.ylabel("Occurence")
    plt.yscale('log')
    plt.xlim(run_min-1, )
    plt.ylim(0.9, )
    plt.legend(loc='upper right', numpoints=1, fontsize = 20 , ncol=1, frameon=True)
    ax.grid(True, which='both')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.f'))
    plt.tight_layout()
    plt.savefig(out_plot_f, format='png')
    print "plotting %s" % out_plot_f
    plt.close()

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    print args
    infile = args['--infile']
    run_min = int(args['--run_min'])
    include_n_flag = args['--include_n']
    bases = ['A', 'G', 'C', 'T']
    if include_n_flag:
        bases.append('N')
    out_plot_f = "%s.base_run_distribution.png" % (basename(infile))
    out_results_f = "%s.base_run_results.txt" % (basename(infile))
    base_run_count = count_base_runs(infile, run_min, include_n_flag, bases)
    plot_base_run_count(base_run_count, out_plot_f)
