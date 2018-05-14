# F2PY_NS_2D
A python code for 2D Incompressible flow interfaced with Fortran through f2py

Poisson solver - Gauss/Siedel or Multigrid V Cycle

"sh CleanandRun.sh" in terminal for visualizing solution and time taken

May need modification of CleanandRun.sh according to Python version (this example uses 3.6 for 64 bit arch),
compiler for Fortran is gfortran

Change flag in init_domain() function to switch between Python and Fortran execution of time integration/Gauss-Siedel
