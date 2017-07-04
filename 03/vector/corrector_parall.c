/*
 * corrector_parall.c - parallel vectorized implementation of corrector in Davidson CI
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
 * 
 */


#include <stdio.h>
#include <mkl.h>

/******************************************************************************/

 static inline void vec_corr(const int i, const lapack_int horder, double * restrict wfn_ci, const double * restrict hdiag, const double dtemp)

/******************************************************************************/
{
	lapack_int ii;

	wfn_ci = __builtin_assume_aligned( wfn_ci,64);  
	hdiag  = __builtin_assume_aligned(  hdiag,64);  
	
	for (ii = 0; ii< horder; ii++)
		wfn_ci[ii] /= (hdiag[ii] - dtemp);

	return;
}

/******************************************************************************/

int corrector_parall(const int blsz, const lapack_int horder, double * restrict wfn_ci, const double * restrict hdiag, const double * restrict levels_ci)

/******************************************************************************/
{
// wfn = res / (H(i,i)- Ei)
	wfn_ci     = __builtin_assume_aligned(   wfn_ci,64);  
	hdiag      = __builtin_assume_aligned(    hdiag,64);  
	levels_ci  = __builtin_assume_aligned(levels_ci,64);  

#pragma omp parallel
{
#pragma omp for schedule(static,1)
	for (int i = 0; i<blsz; i++)
	for (lapack_int ii = 0; ii< horder; ii++)
		wfn_ci[(blsz+i)*horder + ii] /= (hdiag[ii] - levels_ci[i]);
//		vec_corr(i,horder,&wfn_ci[(blsz+i)*horder],hdiag,levels_ci[i]);
//static void vec_corr(const int i, const lapack_int horder, double * restrict wfn_ci, const double * restrict hdiag, const double dtemp)
}	
	return 0;
}

