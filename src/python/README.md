## Building kfold python extension library

1. `f2py run.f90 -m kfold -h kfold.pyf --overwrite-signature`
2. `make -f Makefile`
3. `f2py --quiet -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -lgfortran -c kfold.pyf *.o`

### Notes

* First f2py invocation generates the .pyf file.
* Only need to run for the subroutines that you want to access from python.
* Compile all fortran source code (generate object files).
* Make sure to use the -fPIC flag in the compiler options to generate the correct code for f2py.
* Second f2py invocation links files together and generates the python wrappers.
* Make sure to specify the gfortran library.
* `--quiet` is a nice option.