#!/usr/bin/env python
# encoding: utf-8
"""
bed_from_genbank.py

grab the gene records from a genbank file (edit for other record types).

- requires:  biopython

"""

from Bio import SeqIO
import sys

def gb2bed(infile):
    out_fh = open("%s.bed" % (infile), 'w')
    for record in SeqIO.parse(open(infile, "rU"), "genbank") :
        print record.__dict__
        chrom = record.name
        for feature in record.features:
            print feature.__dict__
            if feature.type == 'gene':
                start = feature.location.start.position
                stop = feature.location.end.position
                name = feature.qualifiers['label'][0]
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = "%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start, stop, name, 0, strand)
                out_fh.write(bed_line)
    out_fh.close()


if __name__ == '__main__':

    infile = ''
    try:
        infile = sys.argv[1]
    except:
        sys.exit("gb2bed.py INFILE")
    gb2bed(infile)
