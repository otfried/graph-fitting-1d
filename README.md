# Fitting a graph to one-dimensional data

This code accompanies the paper *Fitting a graph to one-dimensional
data* by Siu-Wing Cheng, Otfried Cheong, and Taegyoung Lee.

It is written in Python 3.

- For reference, `qp.py` solves the problem using the quadratic solver
  in CVXOPT.  You will need the `cvxopt` and `picos` packages.

- `explicit.py` solves the problem by storing intermediate solutions
   explicitely as piecewise linear functions (implemented in
   `pwlf.py`).

- `implicit.py` solves the problem using and implicit representation
   of the intermediate solutions.  This guarantees quadratic running
   time, but is in practice slower than the explicit representation
   (at least for random inputs).

- `runtime.py` compares the running time of the three methods.

- `pieces.py` counts the number of pieces in the intermediate solutions.

- `rational.py` solves the problem exactly using rational arithmetic.

- `precision.py` solves the problem using arbitrary-precision
   arithmetic.  You will need to install the `mpmath` package.