rm Fortran_functions.so
rm "Fortran_functions.cpython*"
rm "Multigrid_solver.cpython*"
f2py -c --fcompiler=gfortran -m Fortran_functions Fortran_functions.f95
f2py -c --fcompiler=gfortran -m Multigrid_solver Multigrid_solver.f95
cp -R Fortran_functions.cpython-36m-x86_64-linux-gnu.so Fortran_functions.so
cp -R Multigrid_solver.cpython-36m-x86_64-linux-gnu.so Multigrid_solver.so
python 2D_NS_f2py.py
