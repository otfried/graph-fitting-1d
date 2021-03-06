#
# Solve 1D graph fitting problem using quadratic programming
#

import picos as pic

import time
import sys
assert sys.hexversion >= 0x3000000

# --------------------------------------------------------------------

def make_problem(dist):
  n = len(dist) # number of points
  prob = pic.Problem()
  # create variables w (for weight)
  w = [ None ] 
  for i in range(1, n):
    w.append(prob.add_variable('w[%d]' % i, 1))
  # create variables v
  v = [ None ]
  for i in range(1, n+1):
    v.append(prob.add_variable('v[%d]' % i, 1))
  # add constraints on weight and computation of v[i]
  prob.add_constraint(w[1] >= 1)
  prob.add_constraint(v[1] == w[1] * dist[1])
  prob.add_constraint(w[n-1] >= 1)
  prob.add_constraint(v[n] == -w[n-1] * dist[n-1])
  for i in range(2, n):
    prob.add_constraint(v[i] == w[i] * dist[i] - w[i-1] * dist[i-1])
    prob.add_constraint(w[i-1] + w[i] >= 1)
  squares = [v[i]**2 for i in range(1, n+1)]
  prob.set_objective('min', sum(squares))
  return prob, w, v

# --------------------------------------------------------------------

def min_weights(distances, solver='cplex'):
  """Compute optimal weights for the given interpoint distances
using a quadratic solver through PICOS.

Weights are returned as a list starting with None, so w1 is at index 1."""
  t0 = time.perf_counter()
  prob, ws, vs = make_problem(distances)
  t1 = time.perf_counter()
  prob.solve(solver=solver, verbose=True, tol=1e-9)
  t2 = time.perf_counter()  
  assert prob.status == 'optimal'
  weights = [ None ]
  for w in ws[1:]:
    weights.append(w.value)
  return weights, t2 - t1, t1 - t0
    
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
