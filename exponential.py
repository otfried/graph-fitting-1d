#
# Count number of linear pieces for different distributions
#

import sys
import time

import explicit

from fractions import Fraction

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
  d = [ None, Fraction(11569, 65536), Fraction(1537, 2048), Fraction(1), Fraction(3, 2),
        Fraction(213, 128), Fraction(71, 32), Fraction(143, 64), Fraction(9, 4), Fraction(1, 2) ]
  p = count_pieces(d)

# --------------------------------------------------------------------
