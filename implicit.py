#
# Solve 1D graph fitting problem using implicit representation
# of functions R_i
#

from pwlf import Piece

import sys
assert sys.hexversion >= 0x3000000

# --------------------------------------------------------------------

class Matrix():
  def __init__(self, xx=1, xy=0, yx=0, yy=1, x0=0, y0=0):
    self.xx = xx
    self.xy = xy
    self.yx = yx
    self.yy = yy
    self.x0 = x0
    self.y0 = y0

  def __mul__(self, rhs):
    "Matrix multiplication."
    return Matrix(self.xx * rhs.xx + self.xy * rhs.yx,
                  self.xx * rhs.xy + self.xy * rhs.yy,
                  self.yx * rhs.xx + self.yy * rhs.yx,
                  self.yx * rhs.xy + self.yy * rhs.yy,
                  self.xx * rhs.x0 + self.xy * rhs.y0 + self.x0,
                  self.yx * rhs.x0 + self.yy * rhs.y0 + self.y0)
    
  def __call__(self, v):
    "Matrix * vector multiplication."
    x, y = v
    x1 = self.xx * x + self.xy * y + self.x0
    y1 = self.yx * x + self.yy * y + self.y0
    return (x1, y1)

  def __str__(self):
    return "[%g %g %g %g %g %g]" % (self.xx, self.xy, self.yx, self.yy,
                                    self.x0, self.y0)

# --------------------------------------------------------------------  

class Node():
  """A Node is the representation of one R_i."""
  def __init__(self, type):
    assert type in ["leaf", "pass", "split"]
    self.type = type

  def __str__(self):
    if self.type == "leaf":
      return "Leaf(%g, %g): %s / %s" % (self.x, self.y, self.left, self.right)
    if self.type == "pass":
      return "Pass"
    if self.type == "split":
      return "Split(%g, %g) %s %s" % (self.x, self.y, self.Mleft, self.Mright)

# --------------------------------------------------------------------  

class Derivatives():
  """Store the sequence of R_i implicitly."""
  def __init__(self, d):
    self.d = d
    L = Node("leaf")
    L.left = Piece(4 * d[2]**2, -2*d[1]*d[2])
    L.right = Piece(3 * d[2]**2, 0)
    L.x = 2*d[1]/d[2]
    L.y = L.right.f(L.x)
    self.R = [ None, None, L ]

  def __getitem__(self, i):
    return self.R[i]

  def side(self, flip, x, x0):
    if flip:
      return x > x0
    else:
      return x < x0

  def f(self, i, x):
    M = Matrix()
    flip = False
    while True:
      n = self[i]
      if n.type == "pass":
        M = M * n.matrix
        i -= 1
      elif n.type == "split":
        x0, y0 = M((n.x, n.y))
        if self.side(flip, x, x0):
          M = M * n.Mleft
          flip = not flip
        else:
          M = M * n.Mright
        i -= 1
      elif n.type == "leaf":
        x0, y0 = M((n.x, n.y))
        if self.side(flip, x, x0):
          lf = n.left
          flip = not flip
        else:
          lf = n.right
        a1 = M.xx + M.xy * lf.a
        b1 = M.xy * lf.b + M.x0
        a2 = M.yx + M.yy * lf.a
        b2 = M.yy * lf.b + M.y0
        p = Piece(a2/a1, b2 - a2*b1/a1)
        return p.f(x)

  def finv(self, i, y):
    M = Matrix()
    flip = False
    while True:
      n = self[i]
      if n.type == "pass":
        M = M * n.matrix
        i -= 1
      elif n.type == "split":
        x0, y0 = M((n.x, n.y))
        if self.side(flip, y, y0):
          M = M * n.Mleft
          flip = not flip
        else:
          M = M * n.Mright
        i -= 1
      elif n.type == "leaf":
        x0, y0 = M((n.x, n.y))
        if self.side(flip, y, y0):
          lf = n.left
          flip = not flip
        else:
          lf = n.right
        a1 = M.xx + M.xy * lf.a
        b1 = M.xy * lf.b + M.x0
        a2 = M.yx + M.yy * lf.a
        b2 = M.yy * lf.b + M.y0
        p = Piece(a2/a1, b2 - a2*b1/a1)
        return p.finv(y)

  # finds x s.t. x + Ri(x)/xi = 1
  def find_diamond(self, i, xi):
    M = Matrix()
    flip = False
    while True:
      n = self[i]
      if n.type == "pass":
        M = M * n.matrix
        i -= 1
      elif n.type == "split":
        x0, y0 = M((n.x, n.y))
        if self.side(flip, 1.0, x0 + y0/xi):
          M = M * n.Mleft
          flip = not flip
        else:
          M = M * n.Mright
        i -= 1
      elif n.type == "leaf":
        x0, y0 = M((n.x, n.y))
        if self.side(flip, 1.0, x0 + y0/xi):
          lf = n.left
          flip = not flip
        else:
          lf = n.right
        a1 = M.xx + M.xy * lf.a
        b1 = M.xy * lf.b + M.x0
        a2 = M.yx + M.yy * lf.a
        b2 = M.yy * lf.b + M.y0
        p = Piece(a2/a1, b2 - a2*b1/a1)
        return (xi - p.b) / (xi + p.a)

  def next_R(self, i):
    assert i+1 == len(self.R)
    d = self.d
    xi = 2 * d[i] * d[i+1]
    wio = self.finv(i, 0.0)
    if wio >= 1.0:
      # case 1: wio > 1
      n = Node("pass")
      n.matrix = Matrix(0, 1/xi, -xi, 4 * d[i+1]**2 / xi)
    else:
      # case 2: wio < 1
      wdia = self.find_diamond(i, xi)
      n = Node("split")
      n.x = 1 - wdia
      n.y = 4 * d[i+1]**2 * n.x - xi * wdia
      n.Mleft = Matrix(-1, 0, -(2*xi + 4 * d[i+1]**2), -1,
                       1, xi + 4 * d[i+1]**2)
      n.Mright = Matrix(0, 1/xi, -xi, 4 * d[i+1]**2 / xi)
    self.R.append(n)

# --------------------------------------------------------------------    

def compute_weights(R, d):
  n1 = len(d) - 1
  w = [None] + ([0.0] * n1)
  w[n1] = max(R.finv(n1, 0.0), 1.0)
  j = n1 - 1
  while j > 1:
    ww = R.finv(j, 2 * d[j] * d[j+1] * w[j+1])
    if ww + w[j+1] < 1.0:
      ww = 1.0 - w[j+1]
    w[j] = ww
    j -= 1
  w[1] = max(w[2] * d[2] / d[1] / 2.0, 1.0)
  return w

# --------------------------------------------------------------------    

def min_weights(d):
  """Compute optimal weights for the given interpoint distances,
using implicit representation of the R_i.

Weights are returned as a list starting with None, so w1 is at index 1."""
  R = Derivatives(d)
  for i in range(3, len(d)):
    R.next_R(i-1)
  w = compute_weights(R, d)
  return w
    
# --------------------------------------------------------------------

if __name__ == "__main__":
  distances = [ None, 7, 3, 29, 30, 12, 13, 17, 5, 8, 12, 29, 3, 5]
  w = min_weights(distances)
  print(w)

# --------------------------------------------------------------------
