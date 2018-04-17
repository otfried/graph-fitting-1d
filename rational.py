#
# Experiments with rational arithmetic
#

import sys
import time

import qp
import explicit
from fitutils import *

from fractions import Fraction

assert sys.hexversion >= 0x3000000

# --------------------------------------------------------------------

def fixup_weights(w):
  """Make sure that weights actually satisfy the constraints."""
  n = len(w)
  if w[1] < 1:
    w[1] = 1
  if w[n-1] < 1:
    w[n-1] = 1
  for i in range(3, n-1):
    if w[i] + w[i-1] < 1:
      w[i] = 1 - w[i-1]

def rationalize_vector(d):
  for i in range(1, len(d)):
    d[i] = Fraction(d[i])

def qp_and_rational(n):
  d = make_distances(n, 50)
  print(d)
  wq, t = qp.min_weights(d)
  fixup_weights(wq)
  rationalize_vector(d)
  we = explicit.min_weights(d)
  #print(wq)
  #print(we)
  print(compare_weights(wq, we))
  qq = quality(d, wq)
  qe = quality(d, we)
  print(qq)
  print(qe)
  print(float(qe))

def rational(n):
  d = make_distances(n, 50)
  rationalize_vector(d)
  t1 = time.perf_counter()
  we = explicit.min_weights(d)
  t2 = time.perf_counter()
  #print(we)
  qe = quality(d, we)
  print(qe)
  print(float(qe))
  sys.stderr.write("Rational arithmetic on %d points: %g secs\n" %
                   (n, t2 - t1))

# --------------------------------------------------------------------

if __name__ == "__main__":
  #qp_and_rational(100)
  rational(3200)

# --------------------------------------------------------------------
