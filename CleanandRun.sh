rm Fortran_functions.so
rm Fortran_functions.cpython-36m-x86_64-linux-gnu
f2py -c --fcompiler=gfortran -m Fortran_functions Fortran_functions.f95
cp -R Fortran_functions.cpython-36m-x86_64-linux-gnu.so Fortran_functions.so
python 2D_NS_f2py.py
