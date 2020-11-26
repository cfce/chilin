#!/usr/bin/env python

"""
The Stringent version, use to collect high quality proper paired reads, should be even number.
based on pair_reads_statistics_recollect_stringent.py to filter SAM files
add standards: exclude 97 and 145 those are the same coordinates,
               add insert size

"""

import sys
import linecache

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(1)



def flag2table(sam):
    """
    convert SAM to FLAG => Freq
    """

    map_table = {}
    map_table["flag"] = []
    map_table["mapq"] = []
## open SAM and generate flag => frequency, mapq => freq
    n = 0
    ##pairs_table = [] ## only scan file once
    pairs_table_index = []
    pair1 = []
    pair2 = []
    list_index = 0
    fout = open(sam+"_PE.sam", "w")
    proper = 0
## only scan files once
    with open(sam, "r") as inf:
        for line in inf:
            if line.startswith("@"):
                n = 0
                fout.write(line)
                continue
            n += 1
            lines = line.strip().split()
            ## [(pair1_flag, pair1_quality, pair2..)]
            if n%2 == 1:
                # line[1]: flag, line[4]: mapq, line[8]: insert size
                pair1 = (lines[1], lines[4], lines[8], line)
            elif n%2 == 0:
                pair2 = (lines[1], lines[4], lines[8], line)
                if pass_criteria([pair1, pair2], 1000):
                    proper += 1
                    fout.write(pair1[3])
                    fout.write(pair2[3])
    fout.close()

    print "\t".join(["high quality proper pairs"])
    print "\t".join([str(proper), str(n)])

def pass_criteria(l, insert):
    """
    high quality paired mapped reads should be
    add insert size standards
    flag: 2, to kick out weird (no 4, no 8)
    both have mapping quality >= 1
    input: list, [[(pair1_flag, pair1_quality), (pair2)],  [another...]]
    """
    ## l[0] pair1_flag, pair2_quality
    ## l[1] pair1_flag, pair2_quality
    ## stringent pass standard
    return all([int(l[0][1]) >= 1, int(l[1][1]) >= 1, 0 < abs(int(l[0][2])) < insert, 0 < abs(int(l[1][2])) < insert,bit_flag(int(l[0][0])), bit_flag(int(l[1][0]))])

def bit_flag(x):
    proper_pair = x&2
    query_unmapped = x&4
    mate_unmapped = x&8
    return proper_pair and (not query_unmapped) and (not mate_unmapped)

def str_flag(f):
    """
    convert bit to string
    """
    flag = \
{1:     "paired read",
2:      "proper pair",
4:        "query unmapped",
8:        "mate unmapped",
16:       "strand of the query (1 -> reverse)",
32:       "strand of the mate",
64:        "first read in pair",
128:        "second read in pair",
256:        "alignment is not primary",
512:        "does not pass quality check",
1024:       " PCR or optical duplicate",
2048:       " Supplementary Alignment or chimeric alignment", ## new bwa mem flags
}
    return flag[f]


def main():
    sam = sys.argv[1]
    flag2table(sam)

if __name__ == "__main__":
    main()
