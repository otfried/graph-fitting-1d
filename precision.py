#
# Experiments with high-precision
#

import sys
import time

import qp
import explicit
from fitutils import *

from mpmath import mpf, mp

assert sys.hexversion >= 0x3000000

# --------------------------------------------------------------------

def fixup_weights(w):
  """Make sure that weights actually satisfy the constraints.
The QP solver uses a tolerance of 1e-6."""
  n = len(w)
  if w[1] < 1:
    w[1] = 1
  if w[n-1] < 1:
    w[n-1] = 1
  for i in range(3, n-1):
    if w[i] + w[i-1] < 1:
      w[i] = 1 - w[i-1]

def mpf_vector(d):
  for i in range(1, len(d)):
    d[i] = mpf(d[i])

# --------------------------------------------------------------------

if __name__ == "__main__":
  mp.dps = 50
  d = make_distances(30, 50)
  print(d)
  wq, t = qp.min_weights(d)
  fixup_weights(wq)
  w1 = explicit.min_weights(d)
  mpf_vector(d)
  we = explicit.min_weights(d)
  #print(wq)
  #print(w1)
  #print(we)
  print(compare_weights(wq, w1))
  print(compare_weights(wq, we))
  print(compare_weights(w1, we))
  qq = quality(d, wq)
  q1 = quality(d, w1)
  qe = quality(d, we)
  print(qq)
  print(q1)
  print(qe)

# --------------------------------------------------------------------
