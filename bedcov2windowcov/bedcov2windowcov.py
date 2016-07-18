#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: bedcov2windowcov.py      -g <BED> -c <GENOMECOV> -b <BLAST>... -v <VCF>
                                [--window_cov <INT>] [--window_var <INT>] [-h|--help]

    Options:
        -h --help                           show this
        -g, --bedfile <BED>                 Bedfile of assembly
        -c, --genomecov <GENOMECOV>         Genomecov file
        -b, --blast <BLAST>...              Blast file(s)
        -v, --vcf <VCF>                     VCF file
        --window_cov <INT>                  Window size for coverage calculation [default: 5000]
        --window_var <INT>                  Window size for variant calculation [default: 5000]
"""

from __future__ import division
import sys
import os
from docopt import docopt
from os.path import basename, isfile, abspath, join
#from scipy import signal

from numpy import linspace, std, mean
import matplotlib.colors as colors
import matplotlib as mat
mat.use('agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.patches as patches
plt.style.use('ggplot')

def progress(iteration, steps, max_value):
    if int(steps) >= 1:
        if int(iteration == max_value):
            sys.stdout.write('\r')
            print "[PROGRESS]\t: %d%%" % (100)
        elif int(iteration) % int(steps) == 0:
            sys.stdout.write('\r')
            print "[PROGRESS]\t: %d%%" % (float(int(iteration)/int(max_value))*100),
            sys.stdout.flush()
        else:
            pass

def read_file(infile):
    with open(infile) as fh:
        for line in fh:
            if not line.startswith("#"):
                yield line.rstrip("\n").split("\t")

class DataObj():
    def __init__(self):
        self.contigs_by_contigID = {}
        self.contig_order = []

    def parse_bed(self, bedfile):
        print "[STATUS] - Parsing %s" % bedfile
        for field in read_file(bedfile):
            contigID = field[0]
            contig_start = int(field[1])
            contig_length = int(field[2])
            self.contigs_by_contigID[contigID] = ContigObj(contigID, contig_length)
            print contigID, self.contigs_by_contigID[contigID].length
            self.contig_order.append(contigID)

    def parse_cov(self, covfile):
        print "[STATUS] - Parsing %s" % covfile
        for field in read_file(covfile):
            contigID = field[0]
            cov = int(field[2])
            #if not contigID in self.contigs_by_contigID:
            #    self.contigs_by_contigID[contigID] = ContigObj(contigID)
            #    self.contig_order.append(contigID)
            self.contigs_by_contigID[contigID].add_cov(cov)
        #self.contigs_by_contigID[contigID].length = len(self.contigs_by_contigID[contigID].covs)

    def parse_blast(self, blastfile):
        print "[STATUS] - Parsing %s" % blastfile
        for field in read_file(blastfile):
            q_cov = float(field[14])
            evalue = float(field[10])
            queryID = field[0]
            subjectID = field[1]
            q_start = int(field[6]) -1
            q_end = int(field[7]) - 1
            s_start = int(field[8]) - 1
            s_end = int(field[9]) - 1
            if q_cov >= 80.0 and evalue <= 0.0:
                blastHitObj = BlastHitObj(blastfile, queryID, subjectID, q_start, q_end, s_start, s_end)
                self.contigs_by_contigID[subjectID].add_blasthit(blastHitObj)

    def parse_vcf(self, vcffile):
        for field in read_file(vcffile):
            contigID = field[0]
            contig_pos = int(field[1])
            self.contigs_by_contigID[contigID].add_variant(contig_pos)

    def yield_contigs(self):
        for contigID in self.contig_order:
            yield self.contigs_by_contigID[contigID]

class BlastHitObj():
    def __init__(self, blastfile, queryID, subjectID, q_start, q_end, s_start, s_end):
        self.blastfile = basename(blastfile)
        self.queryID = queryID
        self.subjectID = subjectID
        self.q_start = q_start
        self.q_end = q_end
        self.s_start = s_start
        self.s_end  = s_end

class ContigObj():
    def __init__(self, contigID, contig_length):
        self.contigID = contigID
        self.length = contig_length

        self.covs = []
        self.covs_mean = []

        self.variants = self.length * [0]
        self.variants_count = []

        self.blast_hits = {}

    def add_blasthit(self, blastHitObj):
        blast_file = blastHitObj.blastfile
        subjectID = blastHitObj.subjectID
        queryID = blastHitObj.queryID
        if subjectID == self.contigID:
            if not blastfile in self.blast_hits:
                self.blast_hits[blastfile] = {}
            if not queryID in self.blast_hits[blastfile]:
                self.blast_hits[blastfile][queryID] = []
            self.blast_hits[blastfile][queryID].append(blastHitObj)

    def add_variant(self, pos):
        self.variants[pos] = 1

    def add_cov(self, cov):
        self.covs.append(int(cov))

    def calculate_covs_mean(self, window):
        print "[STATUS] - Calculating mean coverage in window"
        w = int(window/2)
        start = 0
        stop = w
        window_values = []
        steps = self.length/1000
        positions = self.length
        for idx, cov in enumerate(self.covs):
            if idx < w:
                window_values = self.covs[start:stop]
                stop += 1
            elif idx >= w and idx < (self.length - w):
                start = idx - w
                stop = idx + w + 1
                window_values = self.covs[start:stop]
            elif idx >= (self.length - w):
                start = idx - w + 1
                stop = self.length + 1
                window_values = self.covs[start:stop]
            else:
                sys.exit("[ERROR1] - WTF?")
            self.covs_mean.append(sum(window_values)/len(window_values))
            progress(idx+1, steps, positions)
        if not len(self.covs_mean) == self.length:
            sys.exit("[ERROR2] - WTF?")

    def calculate_variant_count(self, window):
        print "[STATUS] - Calculating variant count"
        w = int(window/2)
        start = 0
        stop = w
        window_values = []
        values = self.variants
        steps = self.length/1000
        positions = self.length
        for idx, cov in enumerate(self.covs):
            if idx < w:
                window_values = values[start:stop]
                stop += 1
            elif idx >= w and idx < (self.length - w):
                start = idx - w
                stop = idx + w + 1
                window_values = values[start:stop]
            elif idx >= (self.length - w):
                start = idx - w + 1
                stop = self.length + 1
                window_values = values[start:stop]
            else:
                sys.exit("[ERROR1] - WTF?")
            self.variants_count.append(sum(window_values))
            progress(idx+1, steps, positions)
        if not len(self.variants_count) == self.length:
            sys.exit("[ERROR2] - WTF?")

    def plot(self):
        #f, ax = plt.subplots(figsize=(24, 12))
        out_plot_f = "%s.png" % self.contigID
        x_values = [x for x in range(self.length)]
        #covs_mean_values = [x+10 if x < 200.0 else 200.0 for x in self.covs_mean]
        blast_hits = len(self.blast_hits)
        fig, axarray = plt.subplots(blast_hits + 1, sharex=True, figsize=(40,40))
        xmin = 0 - 20
        xmax = self.length + 20
        covs_mean_values = [y if y <= MAX_COV else MAX_COV for y in self.covs_mean]
        variant_count_values = [y for y in self.variants_count]
        if (blast_hits):
            for ax in axarray:
                ax.set_xlim(xmin, xmax)
            axarray[-1].set_ylim(-10, 210)
            axarray[-1].set_title('Coverage and Variants', fontsize= 24)
            axarray[-1].plot(x_values, covs_mean_values, color = "royalblue", linewidth=4.0, label = "Mean coverage (w=%s)" % (window_cov))
            axarray[-1].plot(x_values, variant_count_values, color = "tomato", linewidth=4.0, label = "Variant count (w=%s)" % (window_var))
            axarray[-1].legend(loc='upper right', numpoints=1, fontsize = 24 , ncol=1, frameon=True)
        else:
            axarray.set_ylim(-10, 210)
            axarray.set_title('Coverage and Variants', fontsize= 24)
            axarray.plot(x_values, covs_mean_values, color = "royalblue", linewidth=4.0, label = "Mean coverage (w=%s)" % (window_cov))
            axarray.plot(x_values, variant_count_values, color = "tomato", linewidth=4.0, label = "Variant count (w=%s)" % (window_var))
            axarray.legend(loc='upper right', numpoints=1, fontsize = 24 , ncol=1, frameon=True)

        if (blast_hits):
            for blast_file_idx, blastfile in enumerate(self.blast_hits):
                num_queries = len(self.blast_hits[blastfile])
                y_positions = linspace(0.1, 0.9, num_queries)
                if num_queries == 1:
                    y_positions = [0.5]
                elif num_queries == 2:
                    y_positions = [0.33, 0.66]
                else:
                    pass
                y_positions = linspace(0.1, 0.9, num_queries)
                y_height = 1/(num_queries*2)
                #print y_height
                y_mins = [y-y_height/2 for y in y_positions]
                #print y_mins
                y_labels = []
                axarray[blast_file_idx].yaxis.grid(False, which='both')
                axarray[blast_file_idx].set_yticks(y_positions)
                axarray[blast_file_idx].set_title("HSPs from %s" % blastfile, fontsize = 24)
                for hit_idx, queryID in enumerate(self.blast_hits[blastfile]):
                    y_pos = y_positions[hit_idx]
                    y_min = y_mins[hit_idx]
                    axarray[blast_file_idx].axhline(y=y_pos, zorder=1, color='white', label=queryID)
                    y_labels.append(queryID)
                    for query_idx, blastHitObj in enumerate(self.blast_hits[blastfile][queryID]):
                        xmin = blastHitObj.s_start
                        xmax = blastHitObj.s_end
                        x_width = xmax - xmin
                        #print xmin, y_min, x_width, y_height
                        patch = patches.Rectangle((xmin, y_min), x_width, y_height, facecolor="royalblue", zorder=10)
                        axarray[blast_file_idx].add_patch(patch)
                    labels = [item.get_text() for item in axarray[blast_file_idx].get_yticklabels()]
                axarray[blast_file_idx].set_yticklabels(y_labels)
        #axarray[2].grid(True, which='both')
        plt.xlabel("Position")

        #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.f'))
        #plt.tight_layout()
        fig.savefig(out_plot_f, format='png')
        print "plotting %s" % out_plot_f
        plt.close()

    def __str__(self):
        string = ''
        for idx in range(self.length):
            string += "%s\t%s\t%s\t%s\n" % (self.contigID, idx, self.covs_mean[idx], self.covs_std[idx])
        return string


if __name__ == "__main__":

    __version__ = 0.1
    args = docopt(__doc__)
    #print args
    bedfile = args['--bedfile']
    covfile = args['--genomecov']
    blastfiles = args['--blast']
    vcffile = args['--vcf']
    window_cov = int(args['--window_cov'])
    window_var = int(args['--window_var'])

    MAX_COV = 200.0

    dataObj = DataObj()
    dataObj.parse_bed(bedfile)
    dataObj.parse_cov(covfile)
    dataObj.parse_vcf(vcffile)
    for blastfile in blastfiles:
        dataObj.parse_blast(blastfile)

    idx = 1
    contig_count = len(dataObj.contig_order)
    steps = contig_count/1000

    for contigObj in dataObj.yield_contigs():
        print "[STATUS] - contig %s (%s out of %s contigs)" % (contigObj.contigID, idx, contig_count)
        if (contigObj.length):
            contigObj.calculate_covs_mean(window_cov)
            contigObj.calculate_variant_count(window_var)
            contigObj.plot()
        #progress(idx, steps, contig_count)
        idx += 1
