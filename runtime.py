#
# Measure running time of different algorithms
#

import sys
import time

import qp
import explicit
import implicit

from fitutils import *

assert sys.hexversion >= 0x3000000

# --------------------------------------------------------------------

def call_one(d, f, descr):
  t1 = time.perf_counter()
  w = f(d)
  t2 = time.perf_counter()  
  sys.stderr.write("Calling %s on %d points: %g secs\n" %
                   (descr, len(d), t2-t1))
  return w, t2-t1
    
def call_three(n):
  sys.stderr.write("%d points:\n" % n)
  d = make_distances(n, 10000)
  w1, t1, t0 = qp.min_weights(d)
  sys.stderr.write("Calling QP on %d points: %g secs (%g secs for modelling)\n" % (n, t1, t0))
  w2, t2 = call_one(d, explicit.min_weights, "Explicit PWLF")  
  w3, t3 = call_one(d, implicit.min_weights, "Implicit Ri")
  dw1 = compare_weights(w1, w2)
  dw2 = compare_weights(w2, w3)
  sys.stderr.write("Delta W = %g, %g\n" % (dw1, dw2))
  return t1, t2, t3
  
# --------------------------------------------------------------------

if __name__ == "__main__":
  tim = [ 100, 200, 400, 800, 1600, 3200]
  timings = {}
  for n in tim:
    timings[n] = call_three(n)
  out = open("timings.tex", "w")
  out.write("""\\begin{tabular}{|c|ccc|}
\\hline
$n$ & QP & Explicit & Implicit \\\\
\\hline
""")
  for n in tim:
    t1, t2, t3 = timings[n]
    out.write("%d & %.3g & %.3g & %.3g \\\\\n" % (n, t1, t2, t3))
  out.write("""\\hline
\\end{tabular}
""")

# --------------------------------------------------------------------
