#!/bin/bash

VERNO=02

# now included in system startup script
#source /opt/intel/mkl/bin/mklvars.sh intel64

cd ${VERNO}
gcc -std=c11 -O2 -m64 -mlong-double-128 -I${MKLROOT}/include -c *.c
cd ../

gcc -std=c11 -O2 -m64 -mlong-double-128 -I${MKLROOT}/include -o cvib.x ${VERNO}/*.o  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lquadmath -lpthread -lm -ldl

