/* ci.c - subroutine to do VCI/QFF calculation
 * using data from GAMESS(US) output files
 * Vectorized and parallelized version
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mkl.h>


#define AMU 1822.9					// mass conversion factor from a.u.(me==1) to a.m.u (1/12 mc == 1)
#define C 0.7						// basis broadness - hardcoded, as in case of GAMESS(US)
#define PI 3.14159265358979323846 	// just pi
#define CM1 219474.63				// hartree to cm-1 conversion factor
#define CI_MAXIT 30					// default max no. of iterations for Davidson CI 
#define CI_TOL 1.0E-6				// default residual convergence tolerance for CI
#define CI_MAX_EXC  8				// default maximal excitation in CAS-CI
#define XRANGE_DFLT 6.0				// default value for xrange

struct ints4 
{
	double x;
	double x2;
	double x3;
	double x4;		
};	

// Forward function declaration
int ci_qff(const char *filename, const int nmodes,const int npairs,const int nbas);
static double pot_1d(const int nmodes,const int mode,const double x,const double *fdiag);
static int icheck(const int nmodes, unsigned char csf[][nmodes],const lapack_int ii,const lapack_int jj,int *diff,int *ndiff);
static void printmat4(const int dim1, const int dim2, const double *mat);
static double res_norm(const int vecno, const lapack_int horder,const int blocksize,double *wfn,
					   const lapack_int ld_wfn,double *mat,const int ld_mat);
static void hwfn(const int blsz, const int iend, const lapack_int horder, const int max_excit, 
				const int nmodes, const int npairs, double wfn_ci[][horder], double hwfn_ci[][horder],
				double int_ijlm[][max_excit][max_excit][max_excit][max_excit],
				const int blocksize, unsigned char csf[][nmodes],const double * hdiag);

//* Function declaration, hwfn_parall.c
void hwfn_parall(const int blsz, const int iend, const lapack_int horder, const int max_excit, 
				const int nmodes, const int npairs, double wfn_ci[][horder], double hwfn_ci[][horder],
				double int_ijlm[][max_excit][max_excit][max_excit][max_excit],
				const int blocksize, unsigned char csf[][nmodes],const double hdiag[]);

//* Function declaration, corrector_parall.c
int corrector_parall(const int blsz, const lapack_int horder, double * restrict wfn_ci,
					 const double * restrict hdiag, const double * restrict levels_ci);

//* Function declaration, resnorm_parall.c
int resnorm_parall(const int blocksize,const int iend,const lapack_int horder,
				   const double * restrict hwfn_ci, const double * restrict mat, double * restrict wfn_norm);

//* Function declaration, print_and_punch.c
void printmat1(const int dim,const double mat[][dim]);
void printmat3(const int dim,const double *mat); 

//* Function declaration, parser_ci.c
int parse_excit(const char *filename, const int nmodes, int excit[]); 
int parse_ci(const char *filename, int *blocksize, int *maxit, double *convtol);
int parse_basis(const char * filename,const int nmodes, double * xrange, double * deltax);
int parse_qff(const char * filename,const int nmodes,double fdiag[][4],double fpair[][6]);

//* Function declaration, gamess.c
int parse_gamess_out(const int nmodes, const int npairs, double harmonic[nmodes], double xrange[nmodes],
                     double deltax[nmodes],double fdiag[nmodes][4], double foffdiag[npairs][6]);

//* Function declaration, timestamp.f90
void timestamp_();

/******************************************************************************/

int ci_qff(const char *filename, const int nmodes,const int npairs,const int nbas) 

/******************************************************************************/
{
double xrange[nmodes],deltax[nmodes];	//basis parameters
double harmonic[nmodes];				// harmonic frequencies
int mode;								// mode counter
int it,i,j,k,l,m;						// loop counters
lapack_int horder,ii,jj,kk,ll;			// large counters and vars
int excit[nmodes];						// array of maximal excitations
int max_excit; 							//maxval(excit)

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
// Allocations

double (*fdiag)[4] = mkl_malloc(sizeof(double[nmodes][4]),64);
if (fdiag==NULL){
	printf(" Error while allocating memory for fdiag \n");
	exit (1);
}
// fdiag 
// 0 => Hii
// 1 => Tiii
// 2 => Uiiii
//

double (*fpair)[6] = mkl_malloc(sizeof(double[npairs][6]),64);					//parameters of potential
if (fpair==NULL){
	printf(" Error while allocating memory for fpair \n");
	exit (1);
}
// fpair 
// 0 => Tiij
// 1 => Tjji
// 2 => Uiiij
// 3 => Ujjji
// 4 => Uiijj

////////////////////////////////////////////////////////////////////
// Fill in default values

for(i=0;i<nmodes;i++)
	excit[i]=CI_MAX_EXC;
	
for(i=0;i<nmodes;i++)
	xrange[i]=XRANGE_DFLT;

for(i=0;i<nmodes;i++)
	deltax[i]=2.0*xrange[i]/(double) (nbas-1);

for(i=0;i<nmodes;i++)
{
	fdiag[i][0]=0.5; //=AMU; // Hii   = 0.5
	fdiag[i][1]=0.0; // Tiii  = 0.0
	fdiag[i][2]=0.0; // Uiiii = 0.0
}

for (i=0;i<npairs;i++)
{
	fpair[i][0]=0.0;
	fpair[i][1]=0.0;
	fpair[i][2]=0.0;
	fpair[i][3]=0.0;
	fpair[i][4]=0.0;
}

// 2D Henon-Heyles potential

/*fdiag[1][1]=-sqrt(0.0125)/3.0; 	// T111
fpair[0][0]=sqrt(0.0125);	// T001 */

// 3D Henon-Heyles potential - default fill for the forcefield

	fdiag[0][1]=-0.01*AMU*6.0*sqrt(AMU); // T000  = 0.01
	fdiag[1][1]=-0.01*AMU*6.0*sqrt(AMU); // T111  = 0.01 
	fpair[0][1]=-0.1*AMU*6.0*sqrt(AMU);  // T011= 0.1
  	fpair[2][1]=-0.1*AMU*6.0*sqrt(AMU);  // T122= 0.1... */


// print banner
printf ("\n");
printf (" !=====================================================================================\n");
printf (" !  nD VCI subprogram for polynomial potential                   \n");
printf (" !  The potential is analytic polynom in QFF format   \n");
printf (" !  Integrals over polynom are analytical               \n");
printf ("\n");

time_t  curr_time  = time(NULL);
char (*str) = ctime (&curr_time);

printf(" Subprogram started at %s",str);

timestamp_(); // call Fortran90 timing program

printf("\n Initial values:\n\n");
printf(" nmodes = %3d\n",nmodes);
printf(" npairs = %3d\n",npairs);
printf(" nbas   = %3d\n",nbas);

printf("\n Reading GAMESS output file...");
parse_gamess_out(nmodes,npairs,harmonic,xrange,deltax,fdiag,fpair);
printf(" done \n");

printf("\n Checking input file for data override...");
parse_excit(filename,nmodes,excit);
parse_basis(filename,nmodes,xrange,deltax);
parse_qff  (filename,nmodes,fdiag,fpair);
printf(" done \n");

//printf("\n Basis centers:\n");
//for (i=0;i<nmodes;i++)
//	printf(" Mode = %2d, xrange = %12.6f, deltax = %12.6f\n",i,xrange[i],deltax[i]);

printf("\n Formatted data blocks\n");

printf ("\n {sizing nmodes = %3d, nbas = %3d}\n",nmodes,nbas);

printf("\n {basis  mode        xrange       deltax\n"); // print {basis} group
for (i=0;i<nmodes;i++)
	printf(" %12d %+12.9f %+12.9f\n",i,xrange[i],deltax[i]);
printf(" }\n");

printf("\n {excit mode  max_excit\n");
for (i=0;i<nmodes;i++)
	printf(" %10d %10d\n",i,excit[i]);
printf(" }\n");

printf("\n {qff \n");
for (i=0;i<nmodes;i++)
{
	printf(" Mode = %2d\n",i);
	printf(" Hii   = %+15.8e\n",fdiag[i][0]);
	printf(" Tiii  = %+15.8e\n",fdiag[i][1]);
	printf(" Uiiii = %+15.8e\n",fdiag[i][2]);
}
k=0;
for (i=0;i<(nmodes-1);i++) {
for (j=(i+1);j<nmodes;j++) 
{
	printf(" Mode pair = %2d , %2d; pair # %3d\n",i,j,k);
	printf(" Tiij  = %+15.8e\n",fpair[k][0]);
	printf(" Tjji  = %+15.8e\n",fpair[k][1]);
	printf(" Uiiij = %+15.8e\n",fpair[k][2]);
	printf(" Ujjji = %+15.8e\n",fpair[k][3]);
	printf(" Uiijj = %+15.8e\n",fpair[k][4]);
	k++;
}}
printf(" }\n");

// convert to polynom coefficients - for GAMESS QFF

for (i=0;i<nmodes;i++)
{
	fdiag[i][0]/=AMU*2.0;
	fdiag[i][1]/=AMU*6.0*sqrt(AMU);
	fdiag[i][2]/=AMU*AMU*24.0;
}

for (i=0;i<npairs;i++)
{
	fpair[i][0]/=AMU*2.0*sqrt(AMU);
	fpair[i][1]/=AMU*2.0*sqrt(AMU);
	fpair[i][2]/=AMU*AMU*6.0;
	fpair[i][3]/=AMU*AMU*6.0;
	fpair[i][4]/=AMU*AMU*4.0;
}

//////////////////////////////////////////////////////////////////////////
// ALLOCATE memory for large matrices, in C99 style VLAs - nbas is known only at runtime 

double (*ss)[nbas] = mkl_malloc( sizeof(double[nbas][nbas]),64);
if (ss==NULL){
	printf(" Error while allocating memory for ss \n");
	exit (1);
}

double (*wfn)[nbas][nbas] = mkl_malloc( sizeof(double[nmodes][nbas][nbas]),64);
if (wfn==NULL){
	printf(" Error while allocating memory for wfn \n");
	exit (1);
} // */

double (*ifi)[nbas][nbas][4] = mkl_malloc( sizeof(double[nmodes][nbas][nbas][4]),64);
if (ifi==NULL){
	printf(" Error while allocating memory for ifi \n");
	exit (1);
}

double (*levels_1d)[nbas] = mkl_malloc( sizeof(double[nmodes][nbas]),64);									//energy levels
if (levels_1d==NULL){
	printf(" Error while allocating memory for levels_1d \n");
	exit (1);
}

//////////////////////////////////////////////////////////////////////////
// do 1D calculations

for(mode=0;mode<nmodes;mode++)
{
	printf("\n=========================================\n");
	printf(" Starting calculation for mode = %2d\n\n",mode);

// fill in A,X
	double x[nbas];
	for (i = 0; i < nbas; ++i)
		x[i] = -xrange[mode] + deltax[mode] * (double) i;

	double a = C * C /(deltax[mode] * deltax[mode]);
	double a2 = a * a; // a[i] * a[j];
	double bi = 4.0 * a;		// 2* exponent coefficient

// print x and potential
	printf("\n Potential at basis points :\n");
	printf("\n  No               x            V(x)\n");
	for (i = 0; i < nbas; ++i)
		printf(" %3d    %12.6f    %12.9f\n",i,x[i],pot_1d(nmodes,mode,x[i],&fdiag[mode][0]));

	printf("\n Basis broadness parameter a = %12.6f\n",a);

/*	printf("\n Transformed polynom coeffs: :\n");
		printf(" F2 = %12.5e\n",fdiag[mode][0]);
		printf(" F3 = %12.5e\n",fdiag[mode][1]);
		printf(" F4 = %12.5e\n",fdiag[mode][2]); */

// fill H,S
	for(i=0;i<nbas;i++){
	for(j=0;j<nbas;j++){

		double dtemp = x[i] - x[j];
		struct ints4 iffi;
		
//    Compact code for overlap matrix S for equidistant basis
		dtemp = exp(- a * dtemp * dtemp * 0.5);
		ss[i][j] = dtemp;

/*    Special code for kinetic energy - useful for polynomial expansions - and symmetric
 *    -(dphi/dx)^2 * 0.5= (2 AI*AJ *x^2 - 2AI*AJ *x *(XI+XJ) + 2(AI*AJ*XI*XJ))* phi^2 */
 
		double xij = 0.5 * (x[i] + x[j]); 		// point of expansion

//    Recurrent formula for polynomial integrals
		iffi.x  = xij;                     			  // integral <phi*x  *phi>
		iffi.x2 = iffi.x  * xij + 1.0    / bi; 		  // integral <phi*x^2*phi>
		iffi.x3 = iffi.x2 * xij + iffi.x  / bi * 2.0; // integral <phi*x^3*phi>
		iffi.x4 = iffi.x3 * xij + iffi.x2 / bi * 3.0; // integral <phi*x^4*phi>

		wfn[mode][i][j] = dtemp * (
		      iffi.x2 * fdiag[mode][0]            	// potential, x^2 part
		+     iffi.x3 * fdiag[mode][1]            	// potential, x^3 part
		+     iffi.x4 * fdiag[mode][2]            	// potential, x^4 part
		+ 2.0*iffi.x2 * a2               	// kinetic energy , x^2 part
		- 2.0*iffi.x  * a2 * (x[i] + x[j])	// kinetic energy , x   part 
		+ 2.0         * a2 *  x[i] * x[j]);	// kinetic energy , const part

//    Save integrals for future use
		ifi[mode][i][j][0] = dtemp * iffi.x;
		ifi[mode][i][j][1] = dtemp * iffi.x2;
		ifi[mode][i][j][2] = dtemp * iffi.x3;
		ifi[mode][i][j][3] = dtemp * iffi.x4;
		
	}
	}

lapack_int info;
lapack_int itype=1;
char jobz='V',uplo='U';

info=LAPACKE_dsygv(LAPACK_ROW_MAJOR,itype,jobz,uplo,nbas,&wfn[mode][0][0],nbas,&ss[0][0],nbas,&levels_1d[mode][0]);
/* lapack_int LAPACKE_dsygv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, double* a, lapack_int lda,
                          double* b, lapack_int ldb, double* w ); */

if (info != 0) printf("\n!! error in dsygv call, info = %d\n",info);

// check the phase of wf, correct if needed

if (wfn[mode][nbas/2][0] < 0.0 )
	for(i=0;i<nbas;i++)
		for(j=0;j<nbas;j++)
			wfn[mode][i][j] *= -1.0;

// print results	

printf("\n Eigenvalues (hartree) and transition energies (cm-1) for mode = %2d\n\n",mode);

for(i=0;i<nbas;++i){
	printf("%4d %20.15f %20.5f\n",i,levels_1d[mode][i],(levels_1d[mode][i]-levels_1d[mode][0])*CM1);
}

printf("\n");

printf("\n Eigenvectors for mode = %2d \n",mode);
printmat3(nbas,&wfn[mode][0][0]);

} // end of for(mode=0;mode<nmodes;mode++)

// print timing

curr_time  = time(NULL);
str = ctime (&curr_time);

printf("\n 1D part finished at %s\n",str);

timestamp_();

///////////////////////////////////////////////////////////////////////
// start CI

printf("\n========================================================\n");
printf(" CSF-based Complete Active Space VCI program (CAS-VCI)\n");

// memory requirements
horder=1;
for(i=0;i<nmodes;i++) 
	horder *= excit[i];

size_t size = horder*horder;

max_excit = (int) excit[0];
	for(i=1;i<nmodes;i++) 
		max_excit = (excit[i]>max_excit) ? (int) excit[i] : max_excit;

printf("\n MAX excitation = %3d\n",max_excit);
printf(" Order of the Hamiltonian matrix = %i\n",horder);
printf(" Number of elements in the Hamiltonian matrix = %lu\n",size);

//printf("\n Memory requirements for CI (excluding wavefunction and temp storage): \n");

size_t totsize=0;

size = sizeof(unsigned char[horder][nmodes]);
size_t Kb = size / 1024;
size_t Mb = Kb / 1024;
size_t Gb = Mb / 1024;
Mb -= Gb*1024;
Kb -= Mb*1024 + Gb *1024 * 1024;
printf("\n Memory reqired for CSF sorage = ");
if (Gb > 0)   	{printf ("%6.2f GB\n",(float) Gb+0.001 * (float) Mb);}
else if (Mb > 0){printf ("%6.2f Mb\n",(float) Mb+0.001 * (float) Kb);}
else 		  	{printf ("%3d KB\n",(int) Kb);}
totsize += size;

size = sizeof(double[nmodes][max_excit][max_excit][4]);
Kb = size / 1024;
Mb = Kb / 1024;
Gb = Mb / 1024;
Mb -= Gb*1024;
Kb -= Mb*1024 + Gb *1024 * 1024;
printf(" Memory reqired for half-transformed integrals = ");
if (Gb > 0)   	{printf ("%6.2f GB\n",(float) Gb+0.001 * (float) Mb);}
else if (Mb > 0) {printf ("%6.2f Mb\n",(float) Mb+0.001 * (float) Kb);}
else 		  	{printf ("%3d KB\n",(int) Kb);}
totsize += size;

size = sizeof(double[npairs][max_excit][max_excit][max_excit][max_excit]);
Kb = size / 1024;
Mb = Kb / 1024;
Gb = Mb / 1024;
Mb -= Gb*1024;
Kb -= Mb*1024 + Gb *1024 * 1024;
printf(" Memory reqired for completely transformed integrals = ");
if (Gb > 0)   	{printf ("%6.2f GB\n",(float) Gb+0.001 * (float) Mb);}
else if (Mb > 0) {printf ("%6.2f Mb\n",(float) Mb+0.001 * (float) Kb);}
else 		  	{printf ("%3d KB\n",(int) Kb);}
totsize += size;

size = sizeof(double[horder]);
Kb = size / 1024;
Mb = Kb / 1024;
Gb = Mb / 1024;
Mb -= Gb*1024;
Kb -= Mb*1024 + Gb *1024 * 1024;
printf(" Memory reqired for diagonal Hamiltonian = ");
if (Gb > 0)   	{printf ("%6.2f GB\n",(float) Gb+0.001 * (float) Mb);}
else if (Mb > 0) {printf ("%6.2f Mb\n",(float) Mb+0.001 * (float) Kb);}
else 		  	{printf ("%3d KB\n",(int) Kb);}
totsize += size;

size = sizeof(double[horder])+sizeof(lapack_int[horder]);
Kb = size / 1024;
Mb = Kb / 1024;
Gb = Mb / 1024;
Mb -= Gb*1024;
Kb -= Mb*1024 + Gb *1024 * 1024;
printf(" Memory required for CI energies and indices = ");
if (Gb > 0)   	{printf ("%6.2f GB\n",(float) Gb+0.001 * (float) Mb);}
else if (Mb > 0) {printf ("%6.2f Mb\n",(float) Mb+0.001 * (float) Kb);}
else 		  	{printf ("%3d KB\n",(int) Kb);}
totsize += size;


Kb = totsize / 1024;
Mb = Kb / 1024;
Gb = Mb / 1024;
Mb -= Gb*1024;
Kb -= Mb*1024 + Gb *1024 * 1024;
printf("\n Total memory required = ");
if (Gb > 0)   	{printf ("%6.2f GB\n",(float) Gb+0.001 * (float) Mb);}
else if (Mb > 0){printf ("%6.2f Mb\n",(float) Mb+0.001 * (float) Kb);}
else 		  	{printf ("%3d KB\n",(int) Kb);}


///////////////////////////////////////////////////////////////////////
// Allocate arrays for starting guess
//
// Allocate space for array of levels
double (*levels_ci) = mkl_malloc( sizeof(double[horder]),64);
if (levels_ci==NULL){
	printf(" Error while allocating memory for levels_ci \n");
	exit (1);
}

// Allocate space for array of indices
lapack_int (*index_ci) = malloc( sizeof(lapack_int[horder]));
if (index_ci==NULL){
	printf(" Error while allocating memory for index_ci \n");
	exit (1);
}

///////////////////////////////////////////////////////////////////////
// Allocate arrays for CI
// for CSF we use smallest int type to save memory
unsigned char (*csf)[nmodes] = malloc(sizeof(unsigned char[horder][nmodes]));
if (csf==NULL){
	printf(" Error while allocating memory for csf \n");
	exit (1);
}

// Allocate space for IFIT - transformed integrals
double (*ifit)[max_excit][max_excit][4] = mkl_malloc( sizeof(double[nmodes][max_excit][max_excit][4]),64);
if (ifit==NULL){
	printf(" Error while allocating memory for ifit \n");
	exit (1);
}

// Allocate space for completely transformed integrals
double (*int_ijlm)[max_excit][max_excit][max_excit][max_excit] = mkl_malloc( sizeof(double[npairs][max_excit][max_excit][max_excit][max_excit]),64);
if (int_ijlm==NULL){
	printf(" Error while allocating memory for int_ijlm \n");
	exit (1);
}

// Allocate space for Hdiag (Diagonal Hamiltonian)
double (*hdiag) = mkl_malloc( sizeof(double[horder]),64);
if (hdiag==NULL){
	printf(" Error while allocating memory for hdiag \n");
	exit (1);
}



///////////////////////////////////////////////////////////////////////
// form CSF

for(i=0;i<nmodes;i++) 
	csf[0][i]= (unsigned char) 0;		//first CSF is always ground state
	
for(ii=1;ii<horder;ii++)
{
	for(i=0;i<nmodes;i++) 				//copy previous vector
		csf[ii][i] = csf[ii-1][i];
	
	for(i=(nmodes-1);i>=0;i--) 
	{
		if(csf[ii][i]<(excit[i]-1))
		{
			csf[ii][i] += (unsigned char) 1;
			break;
		}
		else
		{
			csf[ii][i]  = (unsigned char) 0;
		}
	}
		
}

////////////////////////////////////////////////////////////////////////
// Integral transformation code from gamess.c
// TODO: for QFF field we need to transform x, x2 and x3 only

for(mode=0;mode<nmodes;mode++)
	{
	for(k=0;k<excit[mode];k++)
		for(l=0;l<excit[mode];l++)
		{
		ifit[mode][k][l][0]=0;
		ifit[mode][k][l][1]=0;
		ifit[mode][k][l][2]=0;
		ifit[mode][k][l][3]=0;
		
		for (i = 0; i < nbas; i++) {
			for (j = 0; j < nbas; j++) 
			{
			double c2 = wfn[mode][i][k] * wfn[mode][j][l];
			ifit[mode][k][l][0]+=ifi[mode][i][j][0] * c2;
			ifit[mode][k][l][1]+=ifi[mode][i][j][1] * c2;
			ifit[mode][k][l][2]+=ifi[mode][i][j][2] * c2;
			ifit[mode][k][l][3]+=ifi[mode][i][j][3] * c2;
			}}
//		printf(" mode = %2d, k = %2d, l = %2d, ifit.x = %15.10f, ifit.x2 = %15.10f, ifit.x = %15.10f\n",mode,k,l,ifit[mode][k][l][0],ifit[mode][k][l][1],ifit[mode][k][l][2]);
		}
	}

////////////////////////////////////////////////////////////////////////
// Completete the transformation of integrals
k=0;
for(int mode1=0; mode1 < (nmodes-1); mode1++)
	for(int mode2=mode1+1; mode2 < nmodes; mode2++)
	{
		for(i=0;i<excit[mode1];i++)
		for(j=0;j<excit[mode1];j++)

		for(l=0;l<excit[mode2];l++)
		for(m=0;m<excit[mode2];m++)
		{
			int_ijlm[k][i][j][l][m] =
			ifit[mode1][i][j][1] * ifit[mode2][l][m][0] * fpair[k][0] +  //integral  Tiij 
			ifit[mode1][i][j][0] * ifit[mode2][l][m][1] * fpair[k][1] +  //integral  Tjji 
			ifit[mode1][i][j][2] * ifit[mode2][l][m][0] * fpair[k][2] +  //integral  Uiiij 
			ifit[mode1][i][j][0] * ifit[mode2][l][m][2] * fpair[k][3] +  //integral  Ujjji 
			ifit[mode1][i][j][1] * ifit[mode2][l][m][1] * fpair[k][4] ;  //integral  Uiijj 
//			printf("int_ijlm[%2d][%1d][%1d][%1d][%1d] = %+10.5e\n",k,i,j,l,m,int_ijlm[k][i][j][l][m]);
		}
//		printf(" mode1 %2d, mode2 %2d, k %2d, excit1 %2d, excit2 %2d\n",mode1,mode2,k,excit[mode1],excit[mode2]);
		k++;
	}


curr_time  = time(NULL);
str = ctime (&curr_time);

printf("\n Integrals transformed at %s\n",str);

timestamp_();

/////////////////////////////////////////////////////////////////////////
// Build Hdiag

for (ii=0;ii<horder;ii++)
{
// zero out the element
	hdiag[ii] = 0.0;

// Sum of 1D energies
	for (mode=0;mode<nmodes;mode++)
		hdiag[ii] += levels_1d[mode][csf[ii][mode]];
// sum of pairs
	int pair=0;
	for(i=0;i<(nmodes-1);i++)
		for(j=i+1;j<nmodes;j++)
		{
			l = (int) csf[ii][i];
			m = (int) csf[ii][j];
			hdiag[ii] += int_ijlm[pair][l][l][m][m];
			pair++;
		}
}
	
////////////////////////////////////////////////////////////////////////
// Choose starting vectors
//

///////////////////////////////////////////////////////////////////////
// fill in levels and index and sort and print

for(ii=0;ii<horder;ii++)
{
	levels_ci[ii]=hdiag[ii];
	index_ci[ii]=ii;
}

// sort and find min. no. of vectors 
k=0;
for(ii=1;ii<horder;ii++)
{
	for(jj=horder-1;jj>ii;jj--)
	{
		if(levels_ci[jj]<levels_ci[jj-1])
		{
			double dtemp = levels_ci[jj];
			levels_ci[jj] = levels_ci[jj-1];
			levels_ci[jj-1] = dtemp;
			lapack_int itemp = index_ci[jj];
			index_ci[jj] = index_ci[jj-1];
			index_ci[jj-1] = itemp;
		}
	}
	i = csf[index_ci[ii]][0];
	for (j=1;j<nmodes;j++)
		i += csf[index_ci[ii]][j];
	if (i == 1)
		k++;
	if (k==nmodes)
	{
		printf("\n Highest fundamental frequency detected, level number = %3d\n\n",ii);
		break;
	}
}

printf(" Recommended minimal block size for Davidson algorithm = %2d\n\n",ii+1);

printf(" Chosen starting vectors: \n\n");

for(jj=0;jj<=ii;jj++)
{
	printf(" %2d    CSF #%6d    ",jj,index_ci[jj]);
	for(i=0;i<nmodes;i++)
		printf("%2d ",csf[index_ci[jj]][i]);
	printf("    Diag. E = %10.5e h, or %10.2f cm-1\n",levels_ci[jj],CM1*(levels_ci[jj]-levels_ci[0]));	
}

////////////////////////////////////////////////////////////////////////
// Start Block Davidson
//

// print banner
printf("\n=========================================\n");
printf(" BLAS3-based Block Davidson algorhitm\n");

int blocksize = ii+1;
int nblocks = 4;

int maxit = CI_MAXIT;
double convtol = CI_TOL;

printf("\n Checking input file for data override...");
parse_ci(filename,&blocksize,&maxit,&convtol);
printf(" done \n");

int iend = nblocks*blocksize;

printf("\n {ci blocksize = %3d, maxit = %3d, convtol = %8.2e}\n",blocksize,maxit,convtol);

// calculate memory
size = 2 * sizeof(double[horder][iend])+sizeof(double[iend][iend])+sizeof(double[blocksize]);
Kb = size / 1024;
Mb = Kb / 1024;
Gb = Mb / 1024;
Mb -= Gb*1024;
Kb -= Mb*1024 + Gb *1024 * 1024;
printf("\n Memory required for wavefunction and temp matrices = ");
if (Gb > 0)   	{printf ("%6.2f GB\n",(float) Gb+0.001 * (float) Mb);}
else if (Mb > 0) {printf ("%6.2f Mb\n",(float) Mb+0.001 * (float) Kb);}
else 		  	{printf ("%3d KB\n",(int) Kb);}

///////////////////////////////////////////////////////////////////////
// Allocate arrays for Block Davidson
//
// Allocate space for wfn
double (*wfn_ci)[horder] = mkl_malloc( sizeof(double[iend][horder]),64);
if (wfn_ci==NULL){
	printf(" Error while allocating memory for wfn_ci \n");
	exit (1);
}

// Allocate space for H*wfn
double (*hwfn_ci)[horder] = mkl_malloc( sizeof(double[iend][horder]),64);
if (hwfn_ci==NULL){
	printf(" Error while allocating memory for hwfn_ci \n");
	exit (1);
}

// Allocate space for temp matrix
double (*mat)[iend] = mkl_malloc( sizeof(double[iend][iend]),64);
if (mat==NULL){
	printf(" Error while allocating memory for mat \n");
	exit (1);
}

// Allocate space for temp vector 
double * wfn_norm = mkl_malloc( sizeof(double[blocksize]),64);
if (wfn_norm==NULL){
	printf("Error while allocating memory for wfn_norm \n");
	exit (1);
}

// initialise arrays
memset(wfn_ci,0,sizeof(double[iend][horder]));

for (i=0;i<blocksize;i++)
	wfn_ci[i][index_ci[i]] = 1.0;

// index_ci is not needed anymore
free(index_ci);

///////////////////////////////////////////////////////////////////////
// Do Davidson
//
// We use here 3 microiteration block - Davidson
// Results are printed only at 3rd microiteration before iteration restart
//

for (it=0;it<maxit;it++)
{
// Effective H in mat (mat = (H*wfn)T * wfn)
//hwfn(blocksize,iend,horder,max_excit,nmodes,npairs,wfn_ci,hwfn_ci,int_ijlm,excit[nmodes-1],csf,hdiag);
hwfn_parall(blocksize,iend,horder,max_excit,nmodes,npairs,wfn_ci,hwfn_ci,int_ijlm,excit[nmodes-1],csf,hdiag);
/*static void hwfn(const int blsz, const int iend, const lapack_int horder, const int max_excit, 
				const int nmodes, const int npairs, double wfn_ci[][horder], double hwfn_ci[][horder],
				double int_ijlm[][max_excit][max_excit][max_excit][max_excit],
				const int blocksize, unsigned char csf[][nmodes], double hdiag[]) */

cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                 blocksize, blocksize, horder, 1.0, &wfn_ci[0][0],
                 horder, &hwfn_ci[0][0], horder, 0.0, &mat[0][0], iend);
/*void cblas_dgemm(const  CBLAS_LAYOUT Layout, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                 const MKL_INT K, const double alpha, const double *A,
                 const MKL_INT lda, const double *B, const MKL_INT ldb,
                 const double beta, double *C, const MKL_INT ldc);*/

// Diagonalise mat

lapack_int info;
char jobz='V',uplo='U';

info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,jobz,uplo,blocksize,&mat[0][0],iend,levels_ci);

/*lapack_int LAPACKE_dsyev( int matrix_layout, char jobz, char uplo, lapack_int n,
                          double* a, lapack_int lda, double* w ); */

if (info != 0) printf("\n!! error in dsyev call, info = %d\n",info);

// Calculate residual hwfn = G*(H*wfn-E*wfn)

// hwfn = -e*wfn+hwfn
for (i=0;i<blocksize;i++)
	cblas_daxpy(horder,-levels_ci[i],&wfn_ci[i][0],1,&hwfn_ci[i][0],1);
/*void cblas_daxpy(const MKL_INT N, const double alpha, const double *X,
                 const MKL_INT incX, double *Y, const MKL_INT incY);*/

// wfn [hblock+i] = mat * hwfn []  
cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                 blocksize, horder, blocksize, 1.0, &mat[0][0],
                 iend, &hwfn_ci[0][0], horder, 0.0, &wfn_ci[blocksize][0], horder);

// perform correction and orthonormalisation of residual thus forming new vectors
// wfn = res / (H(i,i)- Ei)
corrector_parall(blocksize,horder,&wfn_ci[0][0],hdiag,levels_ci);
//int corrector_parall(const int blsz, const lapack_int horder, double * restrict wfn_ci, const double * restrict hdiag, const double * restrict levels_ci)



// orthonorm new part of wfn
for(i=0;i<blocksize;i++)
{
	for (j = 0; j < (blocksize+i); j++)
	{
		double dtemp = cblas_ddot (horder, &wfn_ci[j][0], 1, &wfn_ci[blocksize+i][0],1);
		cblas_daxpy(horder,-dtemp,&wfn_ci[j][0],1,&wfn_ci[blocksize+i][0],1);
	}
	double dtemp = cblas_dnrm2(horder,&wfn_ci[blocksize+i][0],1);
	cblas_dscal(horder,1.0/dtemp,&wfn_ci[blocksize+i][0],1);
}

int blsz = 2*blocksize;
// hwnf_ci = h*wfn_ci
//hwfn(blsz,iend,horder,max_excit,nmodes,npairs,wfn_ci,hwfn_ci,int_ijlm,excit[nmodes-1],csf,hdiag);
hwfn_parall(blsz,iend,horder,max_excit,nmodes,npairs,wfn_ci,hwfn_ci,int_ijlm,excit[nmodes-1],csf,hdiag);

// mat = hwfn_ci * (wfn_ci)T = wfn_ci*H*(wfn_ci)T
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                 blsz, blsz, horder, 1.0, &wfn_ci[0][0],
                 horder, &hwfn_ci[0][0], horder, 0.0, &mat[0][0], iend);


// Diagonalise mat: mat -> diag

//jobz='V';uplo='U';

info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,jobz,uplo,blsz,&mat[0][0],iend,levels_ci);
if (info != 0) printf("\n!! error in dsyev call, info = %d\n",info);

// Calculate residual hwfn = G*(H*wfn-E*wfn)

// hwfn = -e*wfn+hwfn
for (i=0;i<blsz;i++)
	cblas_daxpy(horder,-levels_ci[i],&wfn_ci[i][0],1,&hwfn_ci[i][0],1);
/*void cblas_daxpy(const MKL_INT N, const double alpha, const double *X,
                 const MKL_INT incX, double *Y, const MKL_INT incY);*/

// wfn [2*hblock+i] = mat * hwfn []  
cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                 blsz, horder, blsz, 1.0, &mat[0][0],
                 iend, &hwfn_ci[0][0], horder, 0.0, &wfn_ci[blsz][0], horder);

// perform correction and orthonormalisation of residual thus forming new vectors

// wfn = res / (H(i,i)- Ei)
corrector_parall(blsz,horder,&wfn_ci[0][0],hdiag,levels_ci);
//int corrector_parall(const int blsz, const lapack_int horder, double * restrict wfn_ci, const double * restrict hdiag, const double * restrict levels_ci)

// orthonorm new part of wfn
for(i=0;i<blsz;i++)
{
	for (j = 0; j < (blsz+i); j++)
	{
		double dtemp = cblas_ddot (horder, &wfn_ci[j][0], 1, &wfn_ci[blsz+i][0],1);
		cblas_daxpy(horder,-dtemp,&wfn_ci[j][0],1,&wfn_ci[blsz+i][0],1);
	}
	double dtemp = cblas_dnrm2(horder,&wfn_ci[blsz+i][0],1);
	cblas_dscal(horder,1.0/dtemp,&wfn_ci[blsz+i][0],1);
}
/////////////
// microit 3
blsz = 4*blocksize;
// hwnf_ci = h*wfn_ci
//hwfn(blsz,iend,horder,max_excit,nmodes,npairs,wfn_ci,hwfn_ci,int_ijlm,excit[nmodes-1],csf,hdiag);
hwfn_parall(blsz,iend,horder,max_excit,nmodes,npairs,wfn_ci,hwfn_ci,int_ijlm,excit[nmodes-1],csf,hdiag);

// mat = hwfn_ci * (wfn_ci)T = wfn_ci*H*(wfn_ci)T
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                 blsz, blsz, horder, 1.0, &wfn_ci[0][0],
                 horder, &hwfn_ci[0][0], horder, 0.0, &mat[0][0], iend);

// Diagonalise mat: mat -> diag

//jobz='V';uplo='U';

info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,jobz,uplo,blsz,&mat[0][0],iend,levels_ci);
if (info != 0) printf("\n!! error in dsyev call, info = %d\n",info);

// Calculate residual hwfn = G*(H*wfn-E*wfn)

// hwfn = -e*wfn+hwfn - do not need full size here
for (i=0;i<blocksize;i++)
	cblas_daxpy(horder,-levels_ci[i],&wfn_ci[i][0],1,&hwfn_ci[i][0],1);

// Print results	
printf("\n Iteration %3d\n",it);
printf(" Eigenvalues (hartree), transition energies (cm-1), and residual norm\n\n");

int 	iconv[blocksize];

//find norm of mat*hwfn vectors
resnorm_parall(blocksize,iend,horder,&hwfn_ci[0][0],&mat[0][0],wfn_norm);
//int resnorm_parall(const int blocksize,const int iend,const lapack_int horder, const double * restrict hwfn_ci, const double * restrict mat, double * restrict wfn_norm)

for(i=0;i<blocksize;i++)
{
//	wfn_norm[i] = res_norm(i,horder,blsz,&hwfn_ci[0][0],horder,&mat[0][0],iend);
	iconv[i] = (wfn_norm[i] > convtol) ? 0 : 1;
	printf("%3d %20.15f %20.2f %25.10f\n",i,(levels_ci[i]),(levels_ci[i]-levels_ci[0])*CM1,wfn_norm[i]);
}

printf("\n");

// Check convergence and restart
//hwfn = mat * wfn
cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                 blsz, horder, blsz, 1.0, &mat[0][0],
                 iend, &wfn_ci[0][0], horder, 0.0, &hwfn_ci[0][0], horder);
                 
// wfn = hwfn
cblas_dcopy(iend*horder,&hwfn_ci[0][0],1,&wfn_ci[0][0],1);
/*void cblas_dcopy(const MKL_INT N, const double *X, const MKL_INT incX,
                 double *Y, const MKL_INT incY);*/

printf("\n");
j = 0;
for (i=0;i<blocksize;i++)
{
	j += iconv[i];
	if (iconv[i] > 0) printf(" Eigenvalue %2d converged!\n",i);
}	

fflush(stdout);

if (j == blocksize) break;

} // end of main loop for Davidson

/////////////////////////////////////////////////////////////////////////
// Print CSF analysis

for (ii=0;ii<blocksize;ii++)
{
	printf("\n Level #%3d, VCI energy = %12.5f cm-1\n",ii,(levels_ci[ii]-levels_ci[0])*CM1);

	double dtemp = fabs(wfn_ci[ii][0]);
	lapack_int kk=0;
	for (jj=1;jj<horder;jj++)
		if(fabs(wfn_ci[ii][jj])>dtemp)
			{
				kk=jj;
				dtemp=fabs(wfn_ci[ii][jj]);
			}
		
	printf(" Lead CSF #%8d   ",kk);
	for (j=0;j<nmodes;j++)
		printf(" %2d",(int) csf[kk][j]);
	printf("    coeff*2 = %12.8f\n",dtemp*dtemp);

	dtemp = 0.0;
	lapack_int ll=0;
	for (jj=0;jj<horder;jj++)
		if (jj!=kk)
		if(fabs(wfn_ci[ii][jj])>dtemp)
			{
				ll=jj;
				dtemp=fabs(wfn_ci[ii][jj]);
			}
		
	printf("      CSF #%8d   ",ll);
	for (j=0;j<nmodes;j++)
		printf(" %2d",(int) csf[ll][j]);
	printf("    coeff*2 = %12.8f\n",dtemp*dtemp);

	dtemp = 0.0;
	lapack_int mm=0;
	for (jj=0;jj<horder;jj++)
		if (jj!=kk)
		if (jj!=ll)
		if(fabs(wfn_ci[ii][jj])>dtemp)
			{
				mm=jj;
				dtemp=fabs(wfn_ci[ii][jj]);
			}

	printf("      CSF #%8d   ",mm);
	for (j=0;j<nmodes;j++)
		printf(" %2d",(int) csf[mm][j]);
	printf("    coeff*2 = %12.8f\n",dtemp*dtemp);
}


///////////////////////////////////////////////////////////////////////
// End of program - free memory 
mkl_free(wfn_ci); mkl_free(hwfn_ci); mkl_free(mat); mkl_free (wfn_norm);
mkl_free(wfn); mkl_free(ss); free(csf); 
mkl_free(levels_1d); mkl_free(fdiag); mkl_free(fpair);
mkl_free(ifi); mkl_free(ifit); mkl_free(levels_ci);
mkl_free(hdiag); mkl_free(int_ijlm);

curr_time  = time(NULL);
str = ctime (&curr_time);

printf("\n Subprogram finished at %s\n",str);

timestamp_(); // call Fortran90 timing program

return 0;	
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

static double pot_1d(const int nmodes,const int mode,const double x,const double *fdiag)

/******************************************************************************/{
double result;
	
result= x * x * (x * (x * fdiag[2] + fdiag[1]) + fdiag[0]);	
	
return result;	
}

/******************************************************************************/

static void printmat4(const int dim1, const int dim2, const double *mat) 

/******************************************************************************/
/* 
 * Purpose:
 * To print wfn[dim1=iend][dim2=horder] in the transposed form
 * */
{
const int ncols=16; /* how many columns fit in one page */
int i,j,k;
int nblocks,nlast;

nblocks = dim1 / ncols; /* how many blocks of 16 columns to print */
nlast   = dim1 % ncols; /* how many columns remains in last block*/

for(i=0;i<nblocks;++i)
{
	printf("    ");
	for(k=0;k<ncols;++k)
		printf("%13d",k+i*ncols);	
	printf("\n");
	
	for(j=0;j<dim2;++j)
	{
		printf(" %4d",j);
		for(k=0;k<ncols;++k)
			printf(" %+12.5e",mat[j+ dim2*(k + i * ncols)]);	
		printf("\n");
	}	
	printf("\n");
}		

if (nlast>0) 
{
	printf("    ");
	for(k=0;k<nlast;++k)
		printf("%13d",k+nblocks*ncols);	
	
	printf("\n");
	
	for(j=0;j<dim2;++j)
	{
		printf(" %4d",j);
		for(k=0;k<nlast;++k)
			printf(" %+12.5e",mat[j + dim2*(k + nblocks * ncols)]);	
		
		printf("\n");
	} 	
}
return;
}

/******************************************************************************/

static double res_norm(const int vecno, const lapack_int horder,const int blocksize,double *wfn,const lapack_int ld_wfn,double *mat,const int ld_mat)

/******************************************************************************/
/*
 * Purpose
 * Calculate norm of the vector vecno of mat*wfn matrix in a memory-saving fashion
*/
{
double result,dtemp;
int i;
lapack_int ii;

result = 0.0;
for (ii=0;ii<horder;ii++)
{
	dtemp = 0.0;
	for (i=0;i<blocksize;i++)
		dtemp += mat[i*ld_mat + vecno]*wfn[i*ld_wfn + ii];
	result += dtemp * dtemp;
}

return sqrt(result);
}

/******************************************************************************/

static void hwfn(const int blsz, const int iend, const lapack_int horder, const int max_excit, 
				const int nmodes, const int npairs, double wfn_ci[][horder], double hwfn_ci[][horder],
				double int_ijlm[][max_excit][max_excit][max_excit][max_excit],
				const int blocksize, unsigned char csf[][nmodes],const double hdiag[])

/******************************************************************************/
/*
 * Purpose
 * To calculate hwfn=H*wfn, forming H blockwise on the fly from int_ijlm and csf
 * This is the default implementation
 * Parallel and vectorized version is hwfn_parall()
 */
{
lapack_int ii,jj;
int i,j,k,l,m,mode;
//int blocksize = excit[nmodes-1];
const int nblocks = horder/blocksize;

// printf("\n in hwfn horder = %d\n",horder);
	
// Allocate space for block of h 
double (*hblock)[blocksize] = mkl_malloc( sizeof(double[blocksize][blocksize]),64);
if (hblock==NULL){
	printf(" Error while allocating memory for hblock \n");
	exit (1);
}

// zero out hwfn
for (i=0;i<blsz;i++)
	for(ii=0;ii<horder;ii++)
		hwfn_ci[i][ii] = 0.0;


/////////////////////////////////////////////////////////////////////////
// fill Hamoiltonial

for (jj=0;jj<nblocks;jj++)
for (ii=0;ii<nblocks;ii++)
{

	int diff[2];
	int ndiff;
//	double dtemp = 0.0;

	if (! icheck(nmodes,csf,ii*blocksize,jj*blocksize,diff,&ndiff)) // Check that we have max 2 modes in difference
		{
			int pair=0;

			switch (ndiff)
			{
////////////////////////////////////////////////////////////////////////
			case 0: // diagonal
// zero out the block
//				for (i=0;i<blocksize;i++)
//					for	(j=0;j<blocksize;j++)
//						hblock[i][j] = 0.0;

// diagonal part of the block
				for (k=0;k<blocksize;k++)
					hblock[k][k] = hdiag[ii*blocksize+k];
	// Sum of 1D energies
/*				{				
 * for (mode=0;mode<nmodes;mode++)
		hblock[k][k] += levels_1d[mode][csf[ii*blocksize+k][mode]];

	// sum of pairs
					pair=0;
					for(i=0;i<(nmodes-1);i++)
						for(j=i+1;j<nmodes;j++)
						{
							l = (int) csf[ii*blocksize+k][i];
							m = (int) csf[ii*blocksize+k][j];
							hblock[k][k] += int_ijlm[pair][l][l][m][m];
							pair++;
						}
				}
*/
// off-diagonal part - 1 mode diff on last mode, others are equal
				for (i=0;i<(blocksize-1);i++)
				for (j=i+1;j<(blocksize);j++)
				{
					hblock[i][j] = 0.0;

					for (k=0;k<(nmodes-1);k++)
					{
						pair  = k*nmodes - k * (k + 3) / 2 + (nmodes-1) - 1; 
						l = (int) csf[ii*blocksize][k];
						hblock[i][j] += int_ijlm[pair][l][l][i][j];
					}	

					hblock[j][i] = hblock[i][j];
				}
// multiply block
				for (k=0;k<blsz     ;k++)
				for (i=0;i<blocksize;i++)
				for (j=0;j<blocksize;j++)
					hwfn_ci[k][blocksize*ii+i] += wfn_ci[k][blocksize*jj+j]*hblock[i][j];

				break;
////////////////////////////////////////////////////////////////////////
			case 1: // 1 mode diff => diff[1] = last mode
				pair  = diff[0]*nmodes - diff[0] * (diff[0] + 3) / 2 + (nmodes - 1) - 1; 
				l = (int) csf[ii*blocksize][diff[0]];
				m = (int) csf[jj*blocksize][diff[0]];
// just copy block of ints - switch to memcpy?
				for (i=0;i<blocksize;i++)
					for	(j=0;j<blocksize;j++)
						hblock[i][j] = int_ijlm[pair][l][m][i][j]; 
// on diagonal, add other integrals				
				for (i=0;i<blocksize;i++)
				{

					for (k=0;k<diff[0];k++)
					{
						pair  = k*nmodes - k * (k + 3) / 2 + diff[0] - 1; 
						j = (int) csf[ii*blocksize][k];
						hblock[i][i] += int_ijlm[pair][j][j][l][m];
					}	
			
					for (k=diff[0]+1;k<(nmodes-1);k++)
					{
						pair  = diff[0]*nmodes - diff[0] * (diff[0] + 3) / 2 + k - 1; 
						j = (int) csf[ii*blocksize][k];
						hblock[i][i] += int_ijlm[pair][l][m][j][j];
					}
				}	
// multiply block
				for (k=0;k<blsz     ;k++)
				for (i=0;i<blocksize;i++)
				for (j=0;j<blocksize;j++)
					hwfn_ci[k][blocksize*ii+i] += wfn_ci[k][blocksize*jj+j]*hblock[i][j];

				break;
////////////////////////////////////////////////////////////////////////
			case 2: // 2 mode diff - simplest code
// find the integral needed
			pair  = diff[0]*nmodes - diff[0] * (diff[0] + 3) / 2 + diff[1] - 1; 
					i = (int) csf[ii*blocksize][diff[0]];
					j = (int) csf[jj*blocksize][diff[0]];
					l = (int) csf[ii*blocksize][diff[1]];
					m = (int) csf[jj*blocksize][diff[1]];
// only 1 number on diagonal
				double dtemp = int_ijlm[pair][i][j][l][m];
// multiply block
				for (k=0;k<blsz     ;k++)
				for (i=0;i<blocksize;i++)
					hwfn_ci[k][blocksize*ii+i] += wfn_ci[k][blocksize*jj+i]*dtemp;//*hblock[i][j];

				break;
			}
		}
//	else // do nothing
}

mkl_free(hblock);
return;	
}
