#
# Utilities for 1D graph fitting
#

import sys
import random, time

assert sys.hexversion >= 0x3000000

# --------------------------------------------------------------------

def compare_weights(w1, w2):
  """Return L_infty norm of difference between two weight vectors."""
  maxDiff = 0
  for i in range(1, len(w1)):
    maxDiff = max(maxDiff, abs(w1[i] - w2[i]))
  return maxDiff

# --------------------------------------------------------------------

def make_distances(n, max_weight = 10000):
  """Make a distance vector for n points by generating integer
interpoint distances uniformly at random from {1, ... , max_weight}."""
  distances = [1] * n
  distances[0] = None
  for i in range(1, n):
    distances[i] = random.randrange(1, max_weight + 1)
  return distances

# --------------------------------------------------------------------

def make_gauss(n):
  """Make a distance vector for n points by generating interpoint
distances using Gaussian distribution."""
  distances = [1] * n
  distances[0] = None
  for i in range(1, n):
    d = random.gauss(100.0, 30.0)
    if d < 1.0:
      d = 1.0
    distances[i] = d
  return distances

# --------------------------------------------------------------------

def quality(w, d):
  """Compute quality Q for given weights and distance vector."""
  n1 = len(d) - 1
  q = (d[1] * w[1])**2 + (d[n1] * w[n1])**2
  for j in range(2, n1+1):
    q += (d[j] * w[j] - d[j-1] * w[j-1])**2
  return q

# --------------------------------------------------------------------
