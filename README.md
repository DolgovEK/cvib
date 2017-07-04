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

1. Introduction

CVIB is created in intention to extend capabilities of the existing VIB SCF / MP2 / CI programs incorporated in the quantum chemistry (QC) packages. CVIB is a separate program, (in perspective) capable to import energy surface data form different QC programs.  It's written on pure C in the attempt of the most effective utilisation of modern computer hardware. CVIB is written for high-end workstations, not clusters. The critical code is parallelized and vectorized, wants lots of local RAM (for CAS-CI jobs), but no MPI , no distributed RAM, etc, is neither used nor planned. 

2. Program capabilities

Current version of CVIB use simple kinetic energy operator -1/2 d2/dx^2, and distributed gaussians basis set, distributed evenly. it implies the use of natural vibrational coordinates - very much like the most of existing VIBSCF / CI programs up to date. 

Version history of program capabilities (inclusive):

Version 0.3 adds vibrational CAS-CI and vibrational SCF followed by MP2, all based on QFF potential (both imported from GAMESS(US) and entered manually). CAS-VCI being the most tough job is vectorized (SIMD) and parallelized (OpenMP), the rest is using parallelism from Intel MKL only.

Version 0.2 consisits of simple model 1D tasks - harmonic oscillator, morse oscillator, as well as general 1D anharmonic solver, which can use polynom and pointwise 1D potentials as input. 1D pointwise potential can be interpolated both locally and globally. The program has GAMESS(US) output parser via friend=gamess option.
