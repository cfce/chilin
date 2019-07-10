#!/usr/bin/env python

"""
sampling sam pair end to a specific count
from 1M, 5M, to 15, 20, 25, 30M
usage: sampling_pe_sam.py sam
use samtools view -h to generate sam with header
use shell to get interface
"""

import sys
from random import sample
from random import seed
import os

seed(999)

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(1)

class Counter(object):
    def __init__(self, number, wanted, header, pair):
        self.pair = True if pair == 'True' else False
        self.header = header
        if self.pair:
            self.number = number/2
            self.wanted = wanted/2
        else:
            self.number = number
            self.wanted = wanted

    def __call__(self):
        r = range(self.number)
        if len(r) < self.wanted:
            raise Exception
        self.rand_num = sample(r, self.wanted)
        if self.pair:
            rand_num = range(self.header) + [ 2*n+self.header for n in self.rand_num ] + [ 2*n+1+self.header for n in self.rand_num ]
            for z in sorted(rand_num):
                yield z
        else:
            rand_num = range(self.header) + [ i + self.header for i in self.rand_num ]
            for i in sorted(rand_num):
                yield i

def sampling(sam, number,format, pair=True, name_order=True):
    import os
    from random import sample
    n = os.popen("samtools view -Sc %s " % sam)
    h = os.popen("samtools view -SH %s | wc -l |cut -d' ' -f1 " % sam)
    count = int(n.read().strip())
    header = int(h.read().strip())
    print(count,header)
    n.close()
    h.close()

    c = Counter(count, int(number), header,pair)
    k = 0
    index = c()
    ix = index.next()
    samf = open(sam)
    out = open(sam + ".%s"%number, "w")

    while 1:
        reads = has_next(samf)
        if not reads:
            break
        if ( k == ix ):
            out.write(reads)
            ix = has_next(index)
            if not ix:
                break
        k += 1

def has_next(iterable):
    """ pick element from generator
    """
    try:
        return next(iterable)
    except StopIteration:
        return None
    finally:
        del iterable


def main():
    sampling(sys.argv[1], sys.argv[2], "sam", sys.argv[3])

if __name__ == "__main__":
    sampling(sys.argv[1], sys.argv[2], "sam", sys.argv[3])
