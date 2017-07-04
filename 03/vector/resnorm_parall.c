/*
 * resnorm_parall.c - parallel and vectorized implementation of (G*hwfn, G*hwfn)
 * marix-vector product is formed on the fly 
 *
 * Version 0.3, July 2017
 * 
 * Copyright 2017 Evgeny Dolgov <dolgov.eugeny@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */


#include <stdio.h>
#include <mkl.h>
#include <math.h>

static inline void vec_setzero(const int blocksize, double * restrict vec)
{

vec      = __builtin_assume_aligned(     vec,64);

for (int i=0;i<blocksize;i++)
	vec[i] = 0.0;
	
return;
}

static inline void vec_muladd(const int blocksize, double * restrict vec, const double * restrict mat, const double dtemp)
{

vec      = __builtin_assume_aligned(     vec,64);
mat      = __builtin_assume_aligned(     mat,32); // mat lead dim is min 4 blocks of blocksize!

		for(int vecno=0;vecno<blocksize;vecno++)
			vec[vecno] += mat[vecno]*dtemp;

return;
}

static inline void vec_squareadd(const int blocksize, double * restrict vec, const double * restrict vtemp)
{

vec      = __builtin_assume_aligned(     vec,64);
vtemp    = __builtin_assume_aligned(   vtemp,64); // mat lead dim is min 4 blocks of blocksize!

	for(int vecno=0;vecno<blocksize;vecno++)
		vec[vecno] += vtemp[vecno]*vtemp[vecno]; 

return;
}

static inline void vec_sqrt(const int blocksize, double * restrict vec)
{

vec      = __builtin_assume_aligned(     vec,64);

	for(int vecno=0;vecno<blocksize;vecno++)
		vec[vecno] = sqrt(vec[vecno]); 

return;
}

/******************************************************************************/

int resnorm_parall(const int blocksize,const int iend,const lapack_int horder, const double * restrict hwfn_ci, const double * restrict mat, double * restrict wfn_norm)

/******************************************************************************/
{

//double vtemp[blocksize];

double * vtemp = mkl_malloc( sizeof(double[blocksize]),64);
if (vtemp==NULL){
	printf("resnorm_parall: Error while allocating memory for vtemp \n");
	exit (1);
}

vtemp    = __builtin_assume_aligned(   vtemp,64);
mat      = __builtin_assume_aligned(     mat,64);
wfn_norm = __builtin_assume_aligned(wfn_norm,64);
hwfn_ci  = __builtin_assume_aligned( hwfn_ci,64);

vec_setzero(blocksize,wfn_norm);
//static void vec_setzero(const int blocksize, double * restrict vec)
//for (int i=0;i<blocksize;i++)
//	wfn_norm[i] = 0.0;


#pragma omp parallel for
for (lapack_int ii=0;ii<horder;ii++)
{

	vec_setzero(blocksize,vtemp);
//	for (int i=0;i<blocksize;i++)
//		vtemp[i] = 0.0;

	for (int i=0;i<blocksize;i++)
		vec_muladd(blocksize,vtemp,&mat[i*iend],hwfn_ci[i*horder + ii]);
//static inline void vec_muladd(const int blocksize, double * restrict vec, double * restrict mat, const double dtemp)
//		for(int vecno=0;vecno<blocksize;vecno++)
//			vtemp[vecno] += mat[i*iend + vecno]*hwfn_ci[i*horder + ii];

	vec_squareadd(blocksize,wfn_norm,vtemp);
//static inline void vec_squareadd(const int blocksize, double * restrict vec, const double * restrict vtemp, const double dtemp)
//	for(int vecno=0;vecno<blocksize;vecno++)
//		wfn_norm[vecno] += vtemp[vecno]*vtemp[vecno]; 

}

mkl_free(vtemp);

vec_sqrt(blocksize, wfn_norm);
//static inline void vec_sqrt(const int blocksize, double * restrict vec)
//for(int vecno=0;vecno<blocksize;vecno++)
//	wfn_norm[vecno] = sqrt(wfn_norm[vecno]); 

	return 0;
}

