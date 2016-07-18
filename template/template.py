#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from os.path import basename, isfile, abspath, join
import sys

def test():
    pass

if __name__ == "__main__":

    SCRIPT = basename(__file__)

    infile, value, outfile = '', 0,''
    try:
        infile = sys.argv[1]
        value = sys.argv[2]
        outfile = sys.argv[3]
    except:
        sys.exit("USAGE: ./%s %s %s %s" % (SCRIPT, "IN", "INTEGER", "OUT"))

    test()
