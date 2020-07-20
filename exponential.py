#
# Count number of linear pieces for an exponential example
#

import sys
import time

from fractions import Fraction

import explicit

assert sys.hexversion >= 0x3000000

# --------------------------------------------------------------------

def count_pieces(d):
  t1 = time.perf_counter()
  R = explicit.compute_R(d)
  t2 = time.perf_counter()
  sys.stderr.write("Running time for %d points: %g secs\n" % (len(d), t2-t1))
  for r in R[2:]:
    print(len(r))
  
# --------------------------------------------------------------------

if __name__ == "__main__":
  #a = 2 * 65536 + 16536
  #b = 65536 + 32768 + 8192 + 2048 + 512
  #d = [ None, 11569, 49184, 65536, 98304, 109056, 145408, 146432, 147456, 32768 ]
  #d = [ None, 11569, 32768 + 16384 + 32, 65536, 65536 + 32768, b, a - 2048, a - 1024, a, 32768 ]
  d = [ None, 0x2d31, 0xc020, 0x10000, 0x18000, 0x1aa00, 0x24000 - 0x800, 0x24000 - 0x400, 0x24000, 0x8000 ]
  d = list(map(lambda x: Fraction(x) if x is not None else None, d))
  count_pieces(d)

# --------------------------------------------------------------------
