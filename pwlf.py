#
# A class representing a piecewise linear function
#

import sys
assert sys.hexversion >= 0x3000000

# --------------------------------------------------------------------

class Piece(object):
  """A linear function y = a * x + b.
  beta is the breakpoint where this linear piece starts."""
  def __init__(self, a, b, beta = 0):
    self.a = a
    self.b = b
    self.beta = beta
  def __repr__(self):
    return "Piece(%s, %s, %s)" % (self.a, self.b, self.beta)
  def fbeta(self):
    """Evaluate function at beta."""
    return self.a * self.beta + self.b
  def f(self, x):
    """Evaluate function at x."""
    return self.a * x + self.b
  def finv(self, y):
    """Evaluate inverse function at y."""
    return (y - self.b) / self.a

class PieceWiseLinear(object):
  """A piecewise linear function (PWLF), where each piece is a Piece
  object."""
  def __init__(self):
    self.pieces = []

  def __len__(self):
    return len(self.pieces)

  def add(self, piece):
    """Add a linear piece to the PWLF."""
    self.pieces.append(piece)

  def __repr__(self):
    return "[" + ", ".join(map(lambda x: repr(x), self.pieces)) + "]"

  def piece_at(self, x):
    """Return the index of the piece spanning x-coordinate x."""
    i = len(self.pieces) - 1
    while i > 0 and self.pieces[i].beta > x:
      i -= 1
    return i
    
  def piece_inv_at(self, y):
    """Return the index of the piece spanning y-coordinate y."""
    i = len(self.pieces) - 1
    while i > 0 and self.pieces[i].fbeta() > y:
      i -= 1
    return i
    
  def f(self, x):
    """Evaluate the PWLF at x-coordinate x."""
    return self.pieces[self.piece_at(x)].f(x)
  
  def finv(self, y):
    """Evaluate the inverse of the PWLF at y-coordinate y."""
    i = self.piece_inv_at(y)
    return self.pieces[i].finv(y)

  def find_diamond(self, xi):
    """Computes x such that x + f(x)/xi = 1.
    The method returns a pair (x, i), where i is the index of the piece 
    containing x."""
    i = len(self.pieces) - 1
    while True:
      p = self.pieces[i]
      if p.fbeta()/xi + p.beta <= 1:
        # found the piece!
        return (xi - p.b) / (xi + p.a), i
      i -= 1

# --------------------------------------------------------------------
