#!/usr/bin/env python

import pandas as pd
import sys
import math

from optparse import OptionParser

def jc(b,do_correction,a):
    deps = 0.0001
    if abs(b) < deps:
        return 0
    loc = 1 - (4.0 * b / 3)

    if (loc <= 0):
        return 5.0 # max dist
    if do_correction:
        if a < deps:
            a = deps
        ret = a * (0.75 * (loc**(-1.0/a) - 1))
    else:
        ret = -0.75 * (math.log(loc))
    return ret

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-d", "--distances", dest="dist_fp",
                      help="path to the table of observed distances", metavar="FILE")
    parser.add_option("-a", "--alpha", dest="alpha", default=1.0,
                      help="gamma value", metavar="FLOAT")
    parser.add_option("-g", "--gamma", dest="do_correction", action='store_true',
                      help="do the correction for the rate heterogeneity, modeled with gamma distribution")

    (options, args) = parser.parse_args()
    dist_fp = options.dist_fp
    alpha = float(options.alpha)
    do_correction = options.do_correction

    df = pd.read_csv(dist_fp, sep="\s+", header = 0, index_col = 0)
    ndf = df.applymap(lambda x: jc(x, do_correction, alpha))
    ndf.to_csv(sys.stdout, sep='\t')

