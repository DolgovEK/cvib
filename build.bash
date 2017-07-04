#!/bin/bash

VERNO=03
MARCH=native

# now included in system startup script
source /opt/intel/mkl/bin/mklvars.sh intel64

cd ${VERNO}
gcc -Wno-format -DMKL_ILP64 -std=c99 -O2 -g -march=${MARCH} -I${MKLROOT}/include -c *.c

cd vector
#gcc -std=c99 -O3 -march=sandybridge -mavx -m64 -mlong-double-128 -ffast-math -ftree-vectorize -fopt-info-vec -I${MKLROOT}/include -c *.c
gcc -Wno-format -DMKL_ILP64 -std=c99 -O3 -g -march=${MARCH} -ffast-math -fopenmp -fopt-info-vec=vec.rep -mavx -I${MKLROOT}/include -c *.c

cd ../fortran
gfortran -O2 -g -march=sandybridge -c *.f90

cd ../../
gcc -Wno-format -DMKL_ILP64 -std=c99 -O2 -g -march=${MARCH} -fopenmp -I${MKLROOT}/include -o cvib.x ${VERNO}/*.o ${VERNO}/vector/*.o ${VERNO}/fortran/*.o -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -lgfortran

