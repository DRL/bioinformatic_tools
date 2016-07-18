#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: plot_cov.py      -i cov... [-h|--help]

    Options:
        -h --help                           show this
        -i, --infile <COV...>               COV file.
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

def plot(infiles):
    f, ax = plt.subplots(figsize=(24, 12))
    idx = 1
    out_plot_f = "plot.png"
    colours = {'G' : u'#608000', 'C' : u'#ff9900', 'A' : u'#e60000', 'T': u'#9999ff', 'N': u'#9ca4a5'}
    #colours = plt.cm.Set1(linspace(0, 1, len(base_run_count)))
    for infile in infiles:
        print infile
        idx = 0
        x_values = []
        cov_values = []
        std_values = []
        with open(infile) as fh:
            for line in fh:
                temp = line.rstrip("\n").split()
                cov = float(temp[2])
                std = float(temp[3])
                x_values.append(idx)
                cov_values.append(cov)
                std_values.append(std)
                idx += 1
        ax.plot(x_values, cov_values, label = "cov." + infile)
        ax.plot(x_values, std_values, label = "std." + infile)

    plt.xlabel("")
    plt.ylabel("")
    #plt.xlim(run_min-1, )
    #plt.ylim(0.9, )
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
    infiles = args['--infile']
    plot(infiles)
