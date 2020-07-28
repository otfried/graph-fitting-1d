#
# Solve 1D graph fitting problem using explicit representation
# of piecewise linear functions
#

from pwlf import Piece, PieceWiseLinear

import sys
assert sys.hexversion >= 0x3000000

# --------------------------------------------------------------------

# compute R[i+1] from R[i]
def next_R(Ri, Rn, d, i):
  xi = 2 * d[i] * d[i+1]
  wio = Ri.finv(0)
  if wio >= 1:
    # case 1: wio > 1
    #print("i = %d, wio = %g" % (i, wio))
    j = Ri.piece_inv_at(0)
    for p in Ri.pieces[j:]:
      nbeta = max(0, p.fbeta() / xi)
      Rn.add(Piece(4*d[i+1]**2 - xi**2/p.a, xi * p.b / p.a, nbeta))
  else:
    # case 2: wio < 1
    wstar, jstar = Ri.find_diamond(xi)
    j1 = Ri.piece_at(1)
    #print("i = %d, w* = %g, j* = %d, j1 = %d" % (i, wstar, jstar, j1))
    # Need pieces for when w_i = 1 - w_i+1
    nbeta = 0
    j = j1
    while j >= jstar:
      p = Ri.pieces[j]
      Rn.add(Piece(p.a + 4 * (d[i] + d[i+1]) * d[i+1], -p.a - p.b - xi, nbeta))
      nbeta = 1 - p.beta
      j -= 1
    for p in Ri.pieces[jstar:]:
      nbeta = max(p.f(wstar), p.fbeta()) / xi
      Rn.add(Piece(4*d[i+1]**2 - xi**2/p.a, xi * p.b / p.a, nbeta))
  
# compute the array of R's
def compute_R(d):
  R = [ None, None ]
  n1 = len(d) - 1
  R.append(PieceWiseLinear())
  R[2].add(Piece(4 * d[2]**2, -2*d[1]*d[2]))
  R[2].add(Piece(3 * d[2]**2, 0, 2*d[1]/d[2]))
  #print("R[2] = %s" % repr(R[2]))
  for i in range(2, n1):
    R.append(PieceWiseLinear())
    Ri = R[i]
    Rn = R[i+1]
    next_R(Ri, Rn, d, i)
    #print("R[%d] = %s" % (i+1, repr(Rn)))
  return R

def compute_weights(R, d):
  n1 = len(d) - 1
  w = [None] + ([0] * n1)
  w[n1] = max(R[n1].finv(0), 1)
  j = n1 - 1
  while j > 1:
    ww = R[j].finv(2 * d[j] * d[j+1] * w[j+1])
    if ww + w[j+1] < 1:
      ww = 1 - w[j+1]
    w[j] = ww
    j -= 1
  w[1] = max(w[2] * d[2] / d[1] / 2, 1)
  return w

# --------------------------------------------------------------------

def min_weights(d):
  """Compute optimal weights for the given interpoint distances,
using explicit representation of PWLF.

Weights are returned as a list starting with None, so w1 is at index 1."""
  R = compute_R(d)
  w = compute_weights(R, d)
  return w
    
# --------------------------------------------------------------------

if __name__ == "__main__":
  if len(sys.argv) >= 2:
    dfile = open(sys.argv[1], "r")
    distances = [ None ]
    for line in dfile.readlines():
      distances.append(float(line))
  else:
    distances = [ None, 7, 3, 29, 30, 12, 13, 17, 5, 8, 12, 29, 3, 5]
  w = min_weights(distances)
  for weight in w[1:]:
    sys.stdout.write("%g\n" % weight)

# --------------------------------------------------------------------
