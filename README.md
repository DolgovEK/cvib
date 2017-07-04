# cvib
This is README for CVIB version 0.3, July 2017

CVIB is the program for molecular vibration calculations, written in C99.

Contents
1. Introduction
2. Program capabilities
3. Installation
3.1. Prerequisites
3.2. Compilation
4. Input file description
5. Useful hints

1.Introduction

CVIB is created in intention to extend capabilities of the existing VIB SCF / MP2 / CI programs incorporated into the quantum chemistry (QC) packages. CVIB is a separate program, (in perspective) capable to import energy surface data from different QC programs. It's written in pure C in the attempt for the most effective utilisation of the modern computer hardware. CVIB is written for the high-end workstations, not clusters. The critical code is parallelized and vectorized, it wants lots of local RAM (for CAS-CI jobs); but no MPI , no distributed RAM, etc, is neither used nor planned. 

2.Program capabilities

Current version of CVIB is writtwn for simple kinetic energy operator -1/2 d2/dx^2, and evenly distributed gaussians basis set. This implies the use of natural vibrational coordinates - very much like the most of existing VIB SCF / CI programs up to date. 

Version history of program capabilities (inclusive):

Version 0.3 adds vibrational CAS-CI and vibrational SCF followed by MP2, all based on QFF potential (both imported from GAMESS(US) and entered manually). CAS-VCI being the most tough job is vectorized (SIMD) and parallelized (OpenMP), the rest is using parallelism from Intel MKL only.

Version 0.2 consisits of simple model 1D tasks - harmonic oscillator, morse oscillator, as well as general 1D anharmonic solver, which can use polynom and pointwise 1D potentials as input. 1D pointwise potential can be interpolated both locally and globally. The program has GAMESS(US) output parser via friend=gamess option.

3.Installation

Prerequisites

Program intensively uses many functions of Intel MKL set of libraries. MKL is absolutely essential, and should be installed separately. It is available at no charge from www.intel.com (non-commercial license). 

Any OpenMP C99 capable compiler with automatic vectorisation should be able to do the main job. Program is tested mainly with GCC 5.x (Ubuntu). For timing statisctics we use gfortran (because, for example, wall clock timer from standard ISO C library has precision of 1 second only), and this code most likely will not work with other Fortran compilers. From the other side, gfrotran is readily available for wide rande of platforms (for GCC and gfortran look at gcc.gnu.org). Finally, for those ones who still prefer pure C, there is an option to manually comment calls to timestamp_() funnction (or even rewrite it in C). It will not influence the computation results, only the timing statisctics.

Compilation

For compilation simply run the script build.bash. It is written for gcc and gfortran compilers, with -march=native option, so the compiler will take the best vector instructions available; but it is assumed that AVX-capable processor is used. We recommend not to change the optimisation options other than -march and -mavx. The build.bash script searches for MKL shell variables script mklvars.sh in /opt/intel/mkl/bin directory which is the default location for modern MKL. If you have MKL installed in the non-default directory, please change this path to represent the actual location of mklvars.sh, or make sure mklvars.sh is already source'd into your shell prior to the compilation. 

4.Input file description 5.Useful hints

For the time being, please refer to the example input files, they are self-explanatory.
