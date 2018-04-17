#
# Count number of linear pieces for different distributions
#

import sys
import time

import explicit
from fitutils import *

assert sys.hexversion >= 0x3000000

# --------------------------------------------------------------------

def count_pieces(d):
  t1 = time.perf_counter()
  R = explicit.compute_R(d)
  t2 = time.perf_counter()
  sys.stderr.write("Running time for %d points: %g secs\n" % (len(d), t2-t1))
  return max(map(lambda r : len(r), R[2:]))
  
# --------------------------------------------------------------------

if __name__ == "__main__":
  tpieces = {}
  mpieces = {}
  tim = [ 100, 1000, 10000, 100000 ]
  distr = [ 'u50', 'u10000', 'gauss' ]
  repeat = 1000
  for n in tim:
    for dis in distr:
      tpieces[n,dis] = 0
      mpieces[n,dis] = 0
  for i in range(repeat):
    for n in tim:
      for dis in distr:
        if dis == 'u50':
          d = make_distances(n, 50)
        elif dis == 'u10000':
          d = make_distances(n, 10000)
        else:
          d = make_gauss(n)
        p = count_pieces(d)
        tpieces[n,dis] += p
        if p > mpieces[n,dis]:
          mpieces[n,dis] = p
  out = open("pieces.tex", "w")
  out.write("""\\begin{tabular}{|c||cc|cc|cc|}
\\hline
  & \\multicolumn{2}{c|}{small uniform} &
  \\multicolumn{2}{c|}{large uniform} &
  \\multicolumn{2}{c|}{Gaussian} \\\\
$n$ & avg & max & avg & max & avg & max \\\\
\\hline
""")
  for n in tim:
    out.write("%d & %g & %d & %g & %d & %g & %d \\\\\n" %
              (n, tpieces[n,'u50']/repeat, mpieces[n,'u50'],
               tpieces[n,'u10000']/repeat, mpieces[n,'u10000'],
               tpieces[n,'gauss']/repeat, mpieces[n,'gauss']))
  out.write("""\\hline
\\end{tabular}
""")

# --------------------------------------------------------------------
