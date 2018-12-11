#!/bin/bash
f2py run.f90 -m kfold -h kfold.pyf --overwrite-signature
make -f Makefile
f2py --quiet -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -lgfortran -c kfold.pyf *.o
rm -f *.o *.mod