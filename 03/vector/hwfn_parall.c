/*
 * hwfn_parall.c - parallel and vectorized implementation of H*wfn (H formed on the fly)
 * 
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


#include <mkl.h>
#include <string.h>
#include <stdio.h>

static int         icheck    (const int nmodes, unsigned char csf[][nmodes],const lapack_int ii,const lapack_int jj,int *diff,int *ndiff);
static void inline vec_muladd(const int blocksize, double * restrict wfn_ci, double * restrict hwfn_ci,const double dtemp);
static void inline mat_muladd(const int blocksize, double * restrict wfn_ci, double * restrict hwfn_ci,double * restrict hblock);
static void inline mat_add   (const int blocksize, const int max_excit, const double * restrict int_ijlm, double * restrict hblock);
static void inline mat_copy  (const int blocksize, const int max_excit, const double * restrict int_ijlm, double * restrict hblock);
static void inline diag_add  (const int blocksize, const double dtemp, double * restrict hblock);


/******************************************************************************/

void hwfn_parall(const int blsz, const int iend, const lapack_int horder, const int max_excit, 
				const int nmodes, const int npairs, double wfn_ci[][horder], double hwfn_ci[][horder],
				double int_ijlm[][max_excit][max_excit][max_excit][max_excit],
				const int blocksize, unsigned char csf[][nmodes],const double * restrict hdiag)

/******************************************************************************/
/*
 * Purpose
 * To calculate hwfn=H*wfn, forming H blockwise on the fly from int_ijlm and csf
 * this is Open MP parallelized version of hwfn
 * also using fast memcpy();
 * ATTENTION: the most critical code for computational performance sits here
 */
{
lapack_int ii,jj;
int i,j,k,l,m;//,mode;
//int blocksize = excit[nmodes-1];
const int nblocks = horder/blocksize;

hdiag = __builtin_assume_aligned(hdiag,64);
int_ijlm = __builtin_assume_aligned(int_ijlm,64);
wfn_ci = __builtin_assume_aligned(wfn_ci,64);
hwfn_ci = __builtin_assume_aligned(hwfn_ci,64);

// printf("\n in hwfn horder = %d\n",horder);
	
// zero out hwfn
/* for (i=0;i<blsz;i++)
	for(ii=0;ii<horder;ii++)
		hwfn_ci[i][ii] = 0.0; // */

memset(hwfn_ci,0,sizeof(double[blsz][horder]));


/////////////////////////////////////////////////////////////////////////
// fill Hamoiltonial
// using OpenMP parallel

#pragma omp parallel private(i,j,k,l,m,ii,jj)
{

// Allocate space for block of h 
double (*hblock)[blocksize] = mkl_malloc( sizeof(double[blocksize][blocksize]),64);
if (hblock==NULL){
	printf(" Error while allocating memory for hblock \n");
	exit (1);
}
hblock = __builtin_assume_aligned(hblock,64);


#pragma omp for schedule(static,1)
for (ii=0;ii<nblocks;ii++)
for (jj=0;jj<nblocks;jj++)
{

	int diff[2];
	int ndiff;
	int ib = ii*blocksize;
	int jb = jj*blocksize;
//	double dtemp = 0.0;

	if (! icheck(nmodes,csf,ib,jb,diff,&ndiff)) // Check that we have max 2 modes in difference
		{
			int pair=0;

			switch (ndiff)
			{
////////////////////////////////////////////////////////////////////////
			case 0: // diagonal

				memset(hblock,0,sizeof(double[blocksize][blocksize]));

// off-diagonal part - 1 mode diff on last mode, others are equal

				for (k=0;k<(nmodes-1);k++)
				{
					pair  = k*nmodes - k * (k + 3) / 2 + (nmodes-1) - 1; 
					l = (int) csf[ib][k];

					mat_add(blocksize, max_excit, &int_ijlm[pair][l][l][0][0], &hblock[0][0]);
//static void inline mat_add(const int blocksize, const int max_excit, const double * restrict int_ijlm, double * restrict hblock)
				}	

// diagonal part of the block
				for (k=0;k<blocksize;k++)
					hblock[k][k] = hdiag[ib+k];


// multiply block
				for (k=0;k<blsz     ;k++)
					mat_muladd(blocksize,&wfn_ci[k][jb],&hwfn_ci[k][ib],&hblock[0][0]);

				break;

////////////////////////////////////////////////////////////////////////
			case 1: // 1 mode diff => diff[1] = last mode
				pair  = diff[0]*nmodes - diff[0] * (diff[0] + 3) / 2 + (nmodes - 1) - 1; 
				l = (int) csf[ib][diff[0]];
				m = (int) csf[jb][diff[0]];

// just copy block of ints - can not switch to memcpy because aray strides can be different

					mat_copy(blocksize, max_excit, &int_ijlm[pair][l][m][0][0], &hblock[0][0]);
//static void inline mat_copy(const int blocksize, const int max_excit; const double * restrict int_ijlm, double * restrict hblock)

// on diagonal, add other integrals				
					for (k=0;k<diff[0];k++)
					{
						pair  = k*nmodes - k * (k + 3) / 2 + diff[0] - 1; 
						j = (int) csf[ib][k];

						diag_add(blocksize, int_ijlm[pair][j][j][l][m], &hblock[0][0]);
//static void inline diag_add(const int blocksize, const double dtemp, double * restrict hblock)
					}	
			
					for (k=diff[0]+1;k<(nmodes-1);k++)
					{
						pair  = diff[0]*nmodes - diff[0] * (diff[0] + 3) / 2 + k - 1; 
						j = (int) csf[ib][k];

						diag_add(blocksize, int_ijlm[pair][l][m][j][j], &hblock[0][0]);
					}

// multiply block
				for (k=0;k<blsz     ;k++)
					mat_muladd(blocksize,&wfn_ci[k][jb],&hwfn_ci[k][ib],&hblock[0][0]);

				break;

////////////////////////////////////////////////////////////////////////
			case 2: // 2 mode diff - simplest code
// find the integral needed
			pair  = diff[0]*nmodes - diff[0] * (diff[0] + 3) / 2 + diff[1] - 1; 
					i = (int) csf[ib][diff[0]];
					j = (int) csf[jb][diff[0]];
					l = (int) csf[ib][diff[1]];
					m = (int) csf[jb][diff[1]];
// only 1 number on diagonal
				double dtemp = int_ijlm[pair][i][j][l][m];
			
				for (k=0;k<blsz;k++)
					vec_muladd(blocksize,&wfn_ci[k][jb],&hwfn_ci[k][ib],dtemp);
//void inline vec_muladd(const int blocksize, double * restrict wfn_ci, double * restrict hwfn_ci,const double dtemp)

				break;
			}
		}
//	else // do nothing
}

mkl_free(hblock);
} // end of parallel

return;	
}

/******************************************************************************/

static int icheck(const int nmodes, unsigned char csf[][nmodes],const lapack_int ii,const lapack_int jj,int *diff,int *ndiff)

/******************************************************************************/{
int i,j;	

j=0;
*ndiff=nmodes;

for	(i=0;i<nmodes;i++)
{
	if (csf[ii][i]!=csf[jj][i])
		{
			if (j>1) return 1; // for speed, break if more than 2 mode diff - can be easily changed to any mode diff (take care of diff[] array size!!)
			diff[j]=i;
			j++;
		}
}

*ndiff=j;

return 0;
	
}

/******************************************************************************/

static void inline vec_muladd(const int blocksize, double * restrict wfn_ci, double * restrict hwfn_ci,const double dtemp)

/******************************************************************************/
{
	int i;
	
	wfn_ci  = __builtin_assume_aligned( wfn_ci,64);
	hwfn_ci = __builtin_assume_aligned(hwfn_ci,64);

	for (i=0;i<blocksize;i++)
					hwfn_ci[i] += wfn_ci[i]*dtemp;

return;
} // */

/******************************************************************************/

static void inline mat_muladd(const int blocksize, double * restrict wfn_ci, double * restrict hwfn_ci,double * restrict hblock)

/******************************************************************************/
{
	int i,j;
	
	wfn_ci  = __builtin_assume_aligned( wfn_ci,64);
	hwfn_ci = __builtin_assume_aligned(hwfn_ci,64);
	hblock  = __builtin_assume_aligned( hblock,64);

	for (i=0;i<blocksize;i++)
	for (j=0;j<blocksize;j++)
		hwfn_ci[i] += wfn_ci[j]*hblock[i*blocksize+j];

return;
} // */

/******************************************************************************/

static void inline mat_add(const int blocksize, const int max_excit, const double * restrict int_ijlm, double * restrict hblock)

/******************************************************************************/
{
	int i,j;
	
	int_ijlm = __builtin_assume_aligned(int_ijlm,64); 
	hblock   = __builtin_assume_aligned(  hblock,64);

					for (i=0;i<(blocksize);i++)
					for (j=0;j<(blocksize);j++)
						hblock[i*blocksize+j] += int_ijlm[i*max_excit+j];
//						hblock[i] += int_ijlm[i];


return;
}

/******************************************************************************/

static void inline diag_add(const int blocksize, const double dtemp, double * restrict hblock)

/******************************************************************************/
{
	int i,j;
	
	hblock   = __builtin_assume_aligned(  hblock,64);

				for (i=0;i<blocksize;i++)
						hblock[i*(blocksize+1)] += dtemp;

return;
}


/******************************************************************************/

static void inline mat_copy(const int blocksize, const int max_excit, const double * restrict int_ijlm, double * restrict hblock)

/******************************************************************************/
{
	int i,j;
	
	int_ijlm = __builtin_assume_aligned(int_ijlm,64);  
	hblock   = __builtin_assume_aligned(  hblock,64);  


				for (i=0;i<blocksize;i++)
					for	(j=0;j<blocksize;j++)
						hblock[i*blocksize + j] = int_ijlm[i*max_excit + j];

return;
}
