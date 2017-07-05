/* vscf_mp2.c - subroutines to do VSCF followed by VMP2 valculation
 * using data read from GAMESS(US) output files
 * prior to VSCF the small survey CAS CI is done to select the states of interest
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


#define AMU 1822.9					// mass from a.u.(me==1) to a.m.u (1/12 mc == 1)
#define C 0.7						// basis broadness - hardcoded in case of GAMSS(US)
#define PI 3.14159265358979323846 	// just pi
#define CM1 219474.63				// hartree to cm-1 conversion factor
#define MAXIT 100					// max no. of iterations for SCF 
#define SCF_CONVTOL 1.0e-10			// SCF convergence criteria in hartree
#define CI_MAX_EXC  4				// default maximal excitation in pre-SCF CI
#define XRANGE_DFLT 6.0				// default value for xrange
#define ENERGY_DFLT 4500.0			// default energy threshold for automatic determination of max excit

struct ints4 
{
	double x;
	double x2;
	double x3;
	double x4;		
};	

// Forward function declaration
static double pot_1d(const int nmodes,const int mode,const double x, double fdiag[nmodes][4]);
static int icheck(const int nmodes,unsigned char *csf,unsigned char *csf_scf,int *diff,int *ndiff);

//* Function declaration, parser_ci.c
int parse_excit(const char *filename, const int nmodes, int excit[]); 
int parse_ci(const char *filename, int *blocksize, int *maxit, double *convtol);
int parse_basis(const char * filename,const int nmodes, double * xrange, double * deltax);
int parse_qff(const char * filename,const int nmodes,double fdiag[][4],double fpair[][6]);

//* Function declaration, gamess.c
int parse_gamess_out(const int nmodes, const int npairs, double harmonic[nmodes], double xrange[nmodes], double deltax[nmodes],double fdiag[nmodes][4], double foffdiag[npairs][6]);

//* Function declaration, print_and_punch.c
void printmat1(const int dim, double mat[][dim]);
void printmat3(const int dim, double *mat); 

//* Function declaration, timestamp.f90
void timestamp_();

/******************************************************************************/

int vscf_qff(const char *filename, const int nmodes,const int npairs,const int nbas) 

/******************************************************************************/
{
double xrange[nmodes],deltax[nmodes];	// basis parameters
double harmonic[nmodes];				// harmonic frequencies
int mode;								// mode counter
int it,i,j,k,l,m;						// loop counters
lapack_int horder,hsize,ii,jj,kk,ll,csf_no;// large counters and vars
int excit[nmodes];						// array of maximal excitations
int max_excit; 							//maxval(excit)

double (*fdiag)[4] = mkl_malloc(sizeof(double[nmodes][4]),16);
if (fdiag==NULL){
	printf(" Error while allocating memory for fdiag \n");
	exit (1);
}
// fdiag 
// 0 => Hii
// 1 => Tiii
// 2 => Uiiii
//

double (*fpair)[6] = mkl_malloc(sizeof(double[npairs][6]),16);	//parameters of potential, aligned to 16 bytes
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


// Fill in default values

for(i=0;i<nmodes;i++)
	excit[i]=CI_MAX_EXC;
	
for(i=0;i<nmodes;i++)
	xrange[i]=XRANGE_DFLT;

for(i=0;i<nmodes;i++)
	deltax[i]=2.0*xrange[i]/(double) (nbas-1);

for(i=0;i<nmodes;i++)
{
	fdiag[i][0]=0.5*AMU*2.0; // Hii   = 0.5
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

// 3D Henon-Heyles potential

	fdiag[0][1]=-0.01*AMU*2.0*sqrt(AMU); // T000  = 0.01
	fdiag[1][1]=-0.01*AMU*2.0*sqrt(AMU); // T111  = 0.01 
	fpair[0][1]=-0.1 *AMU*2.0*sqrt(AMU);// T011= 0.1
  	fpair[2][1]=-0.1 *AMU*2.0*sqrt(AMU);// T122= 0.1... */


// print banner
printf ("\n");
printf (" !=====================================================================================\n");
printf (" !  nD VSCF subprogram for polynomial potential in QFF form                  \n");
printf (" !  The potential is analytic polynom (to be read from GAMESS(US) output file)   \n");
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

// ALLOCATE memory for large matrices, in C99 style VLAs - nbas is known only at runtime 

double (*ss)[nbas] = mkl_malloc( sizeof(double[nbas][nbas]),32);
if (ss==NULL){
	printf(" Error while allocating memory for ss \n");
	exit (1);
}

double (*wfn)[nbas][nbas] = mkl_malloc( sizeof(double[nmodes][nbas][nbas]),32);
if (wfn==NULL){
	printf(" Error while allocating memory for wfn \n");
	exit (1);
} // */

double (*ifi)[nbas][nbas][4] = mkl_malloc( sizeof(double[nmodes][nbas][nbas][4]),32);
if (ifi==NULL){
	printf(" Error while allocating memory for ifi \n");
	exit (1);
}

double (*levels_1d)[nbas] = mkl_malloc( sizeof(double[nmodes][nbas]),32);									//energy levels
if (levels_1d==NULL){
	printf(" Error while allocating memory for levels_1d \n");
	exit (1);
}

// do calculations

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
		printf(" %3d    %12.6f    %12.9f\n",i,x[i],pot_1d(nmodes,mode,x[i],fdiag));

	printf("\n Basis broadness parameter a = %12.6f\n",a);

	printf("\n Transformed polynom coeffs: :\n");
		printf(" F2 = %12.5e\n",fdiag[mode][0]);
		printf(" F3 = %12.5e\n",fdiag[mode][1]);
		printf(" F4 = %12.5e\n",fdiag[mode][2]); 

// fill H,S
	for(i=0;i<nbas;i++){
	for(j=0;j<nbas;j++){

		double dtemp = x[i] - x[j];
		struct ints4 iffi;
		
//    Compact code for overlap matrix S for equidistant basis
		dtemp = exp(- a * dtemp * dtemp * 0.5);
		ss[i][j] = dtemp;

/*  special code for kinetic energy - useful for polynomial expansions - and symmetric
 *    -(dphi/dx)^2 * 0.5= (2 AI*AJ *x^2 - 2AI*AJ *x *(XI+XJ) + 2(AI*AJ*XI*XJ))* phi^2 */
 
		double xij = 0.5 * (x[i] + x[j]); 		// point of expansion

		iffi.x  = xij;                     			// integral <phi*x  *phi>
		iffi.x2 = iffi.x  * xij + 1.0    / bi; 		// integral <phi*x^2*phi>
		iffi.x3 = iffi.x2 * xij + iffi.x  / bi * 2.0; 	// integral <phi*x^3*phi>
		iffi.x4 = iffi.x3 * xij + iffi.x2 / bi * 3.0; 	// integral <phi*x^4*phi>

		wfn[mode][i][j] = dtemp * (
		      iffi.x2 * fdiag[mode][0]            	// potential, x^2 part
		+     iffi.x3 * fdiag[mode][1]            	// potential, x^3 part
		+     iffi.x4 * fdiag[mode][2]            	// potential, x^4 part
		+ 2.0*iffi.x2 * a2               	// kinetic energy , x^2 part
		- 2.0*iffi.x  * a2 * (x[i] + x[j])	// kinetic energy , x   part 
		+ 2.0         * a2 *  x[i] * x[j]);	// kinetic energy , const part

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

// start initial guess
printf("\n=========================================\n");
printf(" Initial guess for MAX EXCIT\n");

printf("\n Automatic determinantion of maximal excitation per mode based on excitation threshold = %8.1f cm-1\n",ENERGY_DFLT);

for(i=0;i<nmodes;i++)
	for(j=1;j<nbas;j++)
		if (fabs(levels_1d[i][j]-levels_1d[i][0])*CM1 < ENERGY_DFLT)
			{excit[i]=j+1;} else break;

printf("\n Recommended excitation levels:\n");
for (i=0;i<nmodes;i++)
	printf("mode = %2d, excit = %2d, 1D energy = %8.1f\n",i,excit[i],fabs(levels_1d[i][excit[i]-1]-levels_1d[i][0])*CM1);

printf("\n Checking input file for excitation override...");
parse_excit(filename,nmodes,excit);
printf(" done \n");

printf("\n Excitation levels after override:\n");

printf("\n {excit mode  max_excit\n");
for (i=0;i<nmodes;i++)
	printf(" %10d %10d\n",i,excit[i]);
printf(" }\n");


// start CI
printf("\n=========================================\n");
printf(" General CSF-based VCI program\n");

// allocate csf descriptor array
horder=1;
for(i=0;i<nmodes;i++) 
	horder *= excit[i];

hsize = horder*horder;

printf("\n Order of the Hamiltonian matrix = %i\n",horder);
printf("\n Number of elementis in the Hamiltonian matrix = %i\n",hsize);

// Allocate array of CSF descriptors
// we use smallest int type to save memory
unsigned char (*csf)[nmodes] = malloc(sizeof(unsigned char[horder][nmodes]));

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


// print CSFs

/* printf("\n Generated CSFs:\n");
for (ii=0;ii<horder;ii++)
	{
		printf("CSF #%3d :   ",ii);
		for (j=0;j<nmodes;j++)
			printf(" %2d",(int) csf[ii][j]);
		double dtemp = 0.0;
		for (mode=0;mode<nmodes;mode++)
			dtemp += levels_1d[mode][csf[ii][mode]];
		printf(" Energy = %20.5f\n",dtemp);
	}

printf("\n");// */

max_excit = (int) excit[0];
	for(i=1;i<nmodes;i++) 
		max_excit = (excit[i]>max_excit) ? (int) excit[i] : max_excit;

printf("\n MAX excitation = %3d\n",max_excit);

// Allocate space for IFIT - transformed integrals
// TODO: memory usage printout
double (*ifit)[max_excit][max_excit][4] = mkl_malloc( sizeof(double[nmodes][max_excit][max_excit][4]),32);
if (ifit==NULL){
	printf(" Error while allocating memory for ifit \n");
	exit (1);
}

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

// Allocate space for H (Hamiltonian)
// TODO: memory usage printout

double (*h)[horder] = mkl_malloc( sizeof(double[horder][horder]),64);
if (h==NULL){
	printf(" Error while allocating memory for h \n");
	exit (1);
}

// Build H matrix

for(ii=0;ii<horder;ii++)
{
	// Diagonal H
	double dtemp = 0.0;
	
	// Sum of 1D energies
	for (mode=0;mode<nmodes;mode++)
		dtemp += levels_1d[mode][csf[ii][mode]];
	
	//Sum of 2D integrals
	k=0;
	for(i=0;i<(nmodes-1);i++)
		for(j=i+1;j<nmodes;j++)
		{
			l = (int) csf[ii][i];
			m = (int) csf[ii][j];
			dtemp +=
			ifit[i][l][l][1] * ifit[j][m][m][0] * fpair[k][0] +  //integral  Tiij 
			ifit[i][l][l][0] * ifit[j][m][m][1] * fpair[k][1] +  //integral  Tjji 
			ifit[i][l][l][2] * ifit[j][m][m][0] * fpair[k][2] +  //integral  Uiiij 
			ifit[i][l][l][0] * ifit[j][m][m][2] * fpair[k][3] +  //integral  Ujjji 
			ifit[i][l][l][1] * ifit[j][m][m][1] * fpair[k][4] ;  //integral  Uiijj */
//			printf("Diagonal H: (i = %2d, l = %2d), (j = %2d, m = %2d), k = %2d, dtemp = %15.10f\n",i,l,j,m,k,dtemp);
			k++;
		}
	
	h[ii][ii] = dtemp;
	
	for(jj=0;jj<ii;jj++)
	{
	// Off-diagonal H
		int diff[2];
		int ndiff;
		dtemp = 0.0;
		if (icheck(nmodes,&csf[ii][0],&csf[jj][0],diff,&ndiff)==0) // Check that we have max 2 modes in difference
		{
			int pair=0;

			i = (int) csf[ii][diff[0]];
			l = (int) csf[jj][diff[0]];

			switch (ndiff)
			{
			case 1: // 1 mode diff - loop over n-1 mode pairs
			
			for (k=0;k<diff[0];k++)
			{
			pair  = k*nmodes - k * (k + 3) / 2 + diff[0] - 1; 
					j = (int) csf[ii][k];
					m = (int) csf[jj][k]; // = j!
			dtemp +=
			ifit[diff[0]][i][l][0] * ifit[k][j][m][1] * fpair[pair][0] +  //integral  Tiij 
			ifit[diff[0]][i][l][1] * ifit[k][j][m][0] * fpair[pair][1] +  //integral  Tjji 
			ifit[diff[0]][i][l][0] * ifit[k][j][m][2] * fpair[pair][2] +  //integral  Uiiij 
			ifit[diff[0]][i][l][2] * ifit[k][j][m][0] * fpair[pair][3] +  //integral  Ujjji 
			ifit[diff[0]][i][l][1] * ifit[k][j][m][1] * fpair[pair][4] ;  //integral  Uiijj */
//			printf ("case 1 : ii = %3d, jj = %3d, ndiff = %1d, diff = %1d, pair = %2d, dtemp = %15.10e\n",ii,jj,ndiff,diff[0],pair,dtemp);
			}	
			
			for (k=diff[0]+1;k<nmodes;k++)
			{
			pair  = diff[0]*nmodes - diff[0] * (diff[0] + 3) / 2 + k - 1; 
					j = (int) csf[ii][k];
					m = (int) csf[jj][k]; // = j!
			dtemp +=
			ifit[diff[0]][i][l][1] * ifit[k][j][m][0] * fpair[pair][0] +  //integral  Tiij 
			ifit[diff[0]][i][l][0] * ifit[k][j][m][1] * fpair[pair][1] +  //integral  Tjji 
			ifit[diff[0]][i][l][2] * ifit[k][j][m][0] * fpair[pair][2] +  //integral  Uiiij 
			ifit[diff[0]][i][l][0] * ifit[k][j][m][2] * fpair[pair][3] +  //integral  Ujjji 
			ifit[diff[0]][i][l][1] * ifit[k][j][m][1] * fpair[pair][4] ;  //integral  Uiijj */
//			printf ("case 1 : ii = %3d, jj = %3d, ndiff = %1d, diff = %1d, pair = %2d, dtemp = %15.10e\n",ii,jj,ndiff,diff[0],pair,dtemp);
			}	
				break;
				
			case 2: // 2 mode diff - only 1 mode pair

			pair  = diff[0]*nmodes - diff[0] * (diff[0] + 3) / 2 + diff[1] - 1; 
					i = (int) csf[ii][diff[0]];
					j = (int) csf[ii][diff[1]];
					l = (int) csf[jj][diff[0]];
					m = (int) csf[jj][diff[1]];
			dtemp +=
			ifit[diff[0]][i][l][1] * ifit[diff[1]][j][m][0] * fpair[pair][0] +  //integral  Tiij 
			ifit[diff[0]][i][l][0] * ifit[diff[1]][j][m][1] * fpair[pair][1] +  //integral  Tjji 
			ifit[diff[0]][i][l][2] * ifit[diff[1]][j][m][0] * fpair[pair][2] +  //integral  Uiiij 
			ifit[diff[0]][i][l][0] * ifit[diff[1]][j][m][2] * fpair[pair][3] +  //integral  Ujjji 
			ifit[diff[0]][i][l][1] * ifit[diff[1]][j][m][1] * fpair[pair][4] ;  //integral  Uiijj */

//			printf ("case 2 : ii = %3d, jj = %3d, ndiff = %1d, diff = %1d  %1d, pair = %2d, dtemp = %15.10e\n",ii,jj,ndiff,diff[0],diff[1],pair, dtemp);
							
				break;
			}	// end of switch 
			h[ii][jj] = dtemp;
			h[jj][ii] = dtemp;
		} else 
		{// end of if */

		h[ii][jj] = 0.0;
		h[jj][ii] = 0.0;
		} 
		
	} // end of H build
}

// print timing

curr_time  = time(NULL);
str = ctime (&curr_time);

printf("\n H matrix built at %s\n",str);

timestamp_();


/*printf("\n Full H matrix \n");
printmat3(horder,&h[0][0]); // */

// Diagonalise H

double (*levels_ci) = mkl_malloc(sizeof(double[horder]),64);

lapack_int info;
char jobz='V',uplo='U';

info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,jobz,uplo,horder,&h[0][0],horder,levels_ci);

/*lapack_int LAPACKE_dsyev( int matrix_layout, char jobz, char uplo, lapack_int n,
                          double* a, lapack_int lda, double* w ); */

if (info != 0) printf("\n!! error in dsyev call, info = %d\n",info);

// Print results	

printf("\n Some eigenvalues (hartree) and transition energies (cm-1) calculated by VCI \n\n");

for(i=0;i<nbas;++i){
	printf("%4d %20.15f %20.5f\n",i,(levels_ci[i]),(levels_ci[i]-levels_ci[0])*CM1);
}

printf("\n");

// print timing

curr_time  = time(NULL);
str = ctime (&curr_time);

printf("\n CI part finished at %s\n",str);

timestamp_();

// Print CSF analysis and choose CSF for SCF

for (ii=1;ii<horder;ii++)
{
	printf("\n Level #%3d, VCI energy = %12.5f cm-1\n",ii,(levels_ci[ii]-levels_ci[0])*CM1);

	double dtemp = fabs(h[0][ii]);
	lapack_int kk=0;
	for (jj=0;jj<horder;jj++)
		if(fabs(h[jj][ii])>dtemp)
			{
				kk=jj;
				dtemp=fabs(h[jj][ii]);
			}
		
	printf(" Lead CSF #%4d   ",kk);
	for (j=0;j<nmodes;j++)
		printf(" %2d",(int) csf[kk][j]);
	printf("    coeff*2 = %12.8f\n",dtemp*dtemp);

	dtemp = fabs(h[0][ii]);
	lapack_int ll=0;
	for (jj=0;jj<horder;jj++)
		if (jj!=kk)
		if(fabs(h[jj][ii])>dtemp)
			{
				ll=jj;
				dtemp=fabs(h[jj][ii]);
			}
		
	printf("      CSF #%4d   ",ll);
	for (j=0;j<nmodes;j++)
		printf(" %2d",(int) csf[ll][j]);
	printf("    coeff*2 = %12.8f\n",dtemp*dtemp);

	dtemp = fabs(h[0][ii]);
	lapack_int mm=0;
	for (jj=0;jj<horder;jj++)
		if (jj!=kk)
		if (jj!=ll)
		if(fabs(h[jj][ii])>dtemp)
			{
				mm=jj;
				dtemp=fabs(h[jj][ii]);
			}

	printf("      CSF #%4d   ",mm);
	for (j=0;j<nmodes;j++)
		printf(" %2d",(int) csf[mm][j]);
	printf("    coeff*2 = %12.8f\n",dtemp*dtemp);

	if (csf[kk][0]==1)
	{
		printf("\n Highest fundamental frequency detected, level number = %3d\n\n",ii);
		break;
	}
}

int max_states =  (int) (ii+1);
unsigned char (*csf_scf)[nmodes] = malloc(sizeof(unsigned char[max_states][nmodes]));

for (j=0;j<nmodes;j++)
	csf_scf[0][j] = 0;

for (i=1;i<max_states;i++)
{
	double dtemp = fabs(h[0][i]);
	lapack_int kk=0;
	for (jj=0;jj<horder;jj++)
		if(fabs(h[jj][i])>dtemp)
			{
				kk=jj;
				dtemp=fabs(h[jj][i]);
			}
		
	printf(" %3d CSF selected for SCF:   ",i);
	for (j=0;j<nmodes;j++)
	{
		csf_scf[i][j]=csf[kk][j];
		printf(" %2d",(int) csf_scf[i][j]);
	}
	printf("\n");
}	

double (*levels_scf) = malloc(sizeof(double[max_states]));
if (levels_scf==NULL){
	printf(" Error while allocating memory for levels_scf \n");
	exit (1);
}

double (*levels_mp2) = malloc(sizeof(double[max_states]));
if (levels_mp2==NULL){
	printf(" Error while allocating memory for levels_mp2 \n");
	exit (1);
}

double (*levels_1d_scf)[nbas] = mkl_malloc( sizeof(double[nmodes][nbas]),32);									//energy levels
if (levels_1d_scf==NULL){
	printf(" Error while allocating memory for levels_1d_scf \n");
	exit (1);
}

double (*wfn_1d)[nbas][nbas] = mkl_malloc( sizeof(double[nmodes][nbas][nbas]),32);
if (wfn_1d==NULL){
	printf(" Error while allocating memory for wfn_1d \n");
	exit (1);
} 

// Backing up WFN from 1D calculation for SCF start
// TODO: swtitch to memcpy
for(i=0;i<nmodes;i++)
	for(j=0;j<nbas;j++)
		for(k=0;k<nbas;k++)
			wfn_1d[i][j][k]=wfn[i][j][k];
			
			
// Do SCF 

for (int state=0;state<max_states;state++)
{


	printf("\n===============================================\n");
	printf(" Starting SCF calculation for state #%3d\n\n",state);

double scf_energy=0;

for(it=0;it<MAXIT;it++)
{

	for(mode=0;mode<nmodes;mode++)
	{
// fill in A,X
	double x[nbas];
	for (i = 0; i < nbas; ++i)
		x[i] = -xrange[mode] + deltax[mode] * (double) i;

	double a = C * C /(deltax[mode] * deltax[mode]);
	double a2 = a * a; // a[i] * a[j];
	double bi = 4.0 * a;		// 2* exponent coefficient

// fill H,S
	for(i=0;i<nbas;i++){
	for(j=0;j<nbas;j++){

		double dtemp = x[i] - x[j];
		struct ints4 iffi;
		
//    Compact code for overlap matrix S for equidistant basis
		dtemp = exp(- a * dtemp * dtemp * 0.5);
		ss[i][j] = dtemp;

/*  special code for kinetic energy - useful for polynomial expansions - and symmetric
 *    -(dphi/dx)^2 * 0.5= (2 AI*AJ *x^2 - 2AI*AJ *x *(XI+XJ) + 2(AI*AJ*XI*XJ))* phi^2 */
 
		double xij = 0.5 * (x[i] + x[j]); 		// point of expansion

		iffi.x  = xij;                     			// integral <phi*x  *phi>
		iffi.x2 = iffi.x  * xij + 1.0    / bi; 		// integral <phi*x^2*phi>
		iffi.x3 = iffi.x2 * xij + iffi.x  / bi * 2.0; 	// integral <phi*x^3*phi>
		iffi.x4 = iffi.x3 * xij + iffi.x2 / bi * 3.0; 	// integral <phi*x^4*phi>

		wfn[mode][i][j] = dtemp * (
		      iffi.x2 * fdiag[mode][0]            	// potential, x^2 part
		+     iffi.x3 * fdiag[mode][1]            	// potential, x^3 part
		+     iffi.x4 * fdiag[mode][2]            	// potential, x^4 part
		+ 2.0*iffi.x2 * a2               	// kinetic energy , x^2 part
		- 2.0*iffi.x  * a2 * (x[i] + x[j])	// kinetic energy , x   part 
		+ 2.0         * a2 *  x[i] * x[j]);	// kinetic energy , const part

// Effective potential

		for (k=0;k<mode;k++)
			{
			int pair  = k*nmodes - k * (k + 3) / 2 + mode - 1; 
//					l = 0;
//					m = 0; 
					l = m = (int) csf_scf[state][k];
			wfn[mode][i][j] += dtemp * (
			iffi.x  * ifit[k][l][m][1] * fpair[pair][0] +  //integral  Tiij 
			iffi.x2 * ifit[k][l][m][0] * fpair[pair][1] +  //integral  Tjji 
			iffi.x  * ifit[k][l][m][2] * fpair[pair][2] +  //integral  Uiiij 
			iffi.x3 * ifit[k][l][m][0] * fpair[pair][3] +  //integral  Ujjji 
			iffi.x2 * ifit[k][l][m][1] * fpair[pair][4]);  //integral  Uiijj */
//			printf ("case 1 : ii = %3d, jj = %3d, ndiff = %1d, diff = %1d, pair = %2d, dtemp = %15.10e\n",ii,jj,ndiff,diff[0],pair,dtemp);
			}	
			
		for (k=mode+1;k<nmodes;k++)
			{
			int pair  = mode*nmodes - mode * (mode + 3) / 2 + k - 1; 
					//l = 0;
					//m = 0;
					l = m = (int) csf_scf[state][k];
			wfn[mode][i][j] += dtemp * (
			iffi.x2 * ifit[k][l][m][0] * fpair[pair][0] +  //integral  Tiij 
			iffi.x  * ifit[k][l][m][1] * fpair[pair][1] +  //integral  Tjji 
			iffi.x3 * ifit[k][l][m][0] * fpair[pair][2] +  //integral  Uiiij 
			iffi.x  * ifit[k][l][m][2] * fpair[pair][3] +  //integral  Ujjji 
			iffi.x2 * ifit[k][l][m][1] * fpair[pair][4]);  //integral  Uiijj */
//			printf ("case 1 : ii = %3d, jj = %3d, ndiff = %1d, diff = %1d, pair = %2d, dtemp = %15.10e\n",ii,jj,ndiff,diff[0],pair,dtemp);
			}	

	}
	}

// diagonalise

	lapack_int info;
	lapack_int itype=1;
	char jobz='V',uplo='U';

	info=LAPACKE_dsygv(LAPACK_ROW_MAJOR,itype,jobz,uplo,nbas,&wfn[mode][0][0],nbas,&ss[0][0],nbas,&levels_1d_scf[mode][0]);
/* lapack_int LAPACKE_dsygv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, double* a, lapack_int lda,
                          double* b, lapack_int ldb, double* w ); */

	if (info != 0) printf("\n!! error in dsygv call, info = %d\n",info);

// check the phase of wf, correct if needed

	if (wfn[mode][nbas/2][0] < 0.0 )
		for(i=0;i<nbas;i++)
		for(j=0;j<nbas;j++)
			wfn[mode][i][j] *= -1.0;

	}// end of for(mode...

	double dtemp1 = 0.0;
	for (mode=0;mode<nmodes;mode++)
		dtemp1 += levels_1d_scf[mode][csf_scf[state][mode]];

	printf ("SCF iteration No %3d, energy = %15.10f, delta-e = %15.10f\n",it,dtemp1,scf_energy-dtemp1);

	if (fabs(scf_energy-dtemp1)<SCF_CONVTOL)
		{
			scf_energy=dtemp1;
			break;
		}

	scf_energy=dtemp1;

// Transform integrals
// TODO - cut transformation for current CSF only 

for(mode=0;mode<nmodes;mode++)
	{
//	for(k=0;k<excit[mode];k++)
//		for(l=0;l<excit[mode];l++)
//		{
		k = l = (int) csf_scf[state][mode];
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
//		}
	}

}// end of for(it...

// Calculate MP1

double mp1_energy = 0.0;
	
	//Sum of 2D integrals
	k=0;
	for(i=0;i<(nmodes-1);i++)
		for(j=i+1;j<nmodes;j++)
		{
//			l = 0;//(int) csf[ii][i];
//			m = 0;//(int) csf[ii][j];
			l = (int) csf_scf[state][i];
			m = (int) csf_scf[state][j];
			mp1_energy +=
			ifit[i][l][l][1] * ifit[j][m][m][0] * fpair[k][0] +  //integral  Tiij 
			ifit[i][l][l][0] * ifit[j][m][m][1] * fpair[k][1] +  //integral  Tjji 
			ifit[i][l][l][2] * ifit[j][m][m][0] * fpair[k][2] +  //integral  Uiiij 
			ifit[i][l][l][0] * ifit[j][m][m][2] * fpair[k][3] +  //integral  Ujjji 
			ifit[i][l][l][1] * ifit[j][m][m][1] * fpair[k][4] ;  //integral  Uiijj */
//			printf("Diagonal H: (i = %2d, l = %2d), (j = %2d, m = %2d), k = %2d, dtemp = %15.10f\n",i,l,j,m,k,dtemp);
			k++;
		}
// dtemp *= (double) (nmodes-1);

printf ("\n Final result : SCF energy = %15.10f, MP1 correction = %15.10f, Total = %15.10f\n",scf_energy,-mp1_energy,scf_energy-mp1_energy);

mp1_energy = scf_energy-mp1_energy;

// Simple CSF-based MP2
/*
double mp2_energy = 0.0;

for(jj=0;jj<horder;jj++)
	{
	// Off-diagonal H
		int diff[2];
		int ndiff;
		if (icheck(nmodes,&csf[jj][0],&csf_scf[state][0],diff,&ndiff)==0) // Check that we have max 2 modes in difference
		{

			if (ndiff==2)
			{
// 2 mode diff - only 1 mode pair

			int pair;
			double dtemp = 0.0;

			pair  = diff[0]*nmodes - diff[0] * (diff[0] + 3) / 2 + diff[1] - 1; 
					i = (int) csf_scf[state][diff[0]];
					j = (int) csf_scf[state][diff[1]];
					l = (int) csf[jj][diff[0]];
					m = (int) csf[jj][diff[1]];
			
			dtemp =
			ifit[diff[0]][i][l][1] * ifit[diff[1]][j][m][0] * fpair[pair][0] +  //integral  Tiij 
			ifit[diff[0]][i][l][0] * ifit[diff[1]][j][m][1] * fpair[pair][1] +  //integral  Tjji 
			ifit[diff[0]][i][l][2] * ifit[diff[1]][j][m][0] * fpair[pair][2] +  //integral  Uiiij 
			ifit[diff[0]][i][l][0] * ifit[diff[1]][j][m][2] * fpair[pair][3] +  //integral  Ujjji 
			ifit[diff[0]][i][l][1] * ifit[diff[1]][j][m][1] * fpair[pair][4] ;  //integral  Uiijj 

//			printf ("case 2 : ii = %3d, jj = %3d, ndiff = %1d, diff = %1d  %1d, pair = %2d, dtemp = %15.10e\n",ii,jj,ndiff,diff[0],diff[1],pair, dtemp);
							
			double dtemp1 = 0.0;
			for (mode=0;mode<nmodes;mode++)
				dtemp1 += levels_1d_scf[mode][csf[jj][mode]];

// TODO: improve this 
// eliminate coinciding CSF and all large contrs in MP2
			//if (fabs(dtemp1 - scf_energy)>1.0e-10) 
			mp2_energy -= dtemp * dtemp / (dtemp1 - scf_energy);  

			}	// end of if (ndiff=
			
		} // end of if 
		
	} // end of for(jj) */

//Alternative fast MP2 code

double alt_mp2=0.0;

k=0;
for(i=0;i<(nmodes-1);i++)
	for(j=i+1;j<nmodes;j++)
	{

		for(l=0;l<nbas;l++) 
		{
		
		if (l==csf_scf[state][i]) continue;
		
		struct ints4 iffi1;

		iffi1.x  = 0.0;
		iffi1.x2 = 0.0;
		iffi1.x3 = 0.0;
		iffi1.x4 = 0.0;

		for (int i1 = 0; i1 < nbas; i1++) 
			for (int j1 = 0; j1 < nbas; j1++) 
			{
			double c2 = wfn[i][i1][l] * wfn[i][j1][csf_scf[state][i]];
			iffi1.x  +=ifi[i][i1][j1][0] * c2;
			iffi1.x2 +=ifi[i][i1][j1][1] * c2;
			iffi1.x3 +=ifi[i][i1][j1][2] * c2;
			iffi1.x4 +=ifi[i][i1][j1][3] * c2;
			}


		for(m=0;m<nbas;m++)
// TODO : choice of l,m for excited states (and replace 0)!!
		{
			
		if (m==csf_scf[state][j]) continue;
			
		struct ints4 iffi2;

		iffi2.x  = 0.0;
		iffi2.x2 = 0.0;
		iffi2.x3 = 0.0;
		iffi2.x4 = 0.0;

		for (int i1 = 0; i1 < nbas; i1++) 
			for (int j1 = 0; j1 < nbas; j1++) 
			{
			double c2 = wfn[j][i1][m] * wfn[j][j1][csf_scf[state][j]];
			iffi2.x  +=ifi[j][i1][j1][0] * c2;
			iffi2.x2 +=ifi[j][i1][j1][1] * c2;
			iffi2.x3 +=ifi[j][i1][j1][2] * c2;
			iffi2.x4 +=ifi[j][i1][j1][3] * c2;
			}

		double dtemp =
		iffi1.x2 * iffi2.x  * fpair[k][0] +  //integral  Tiij 
		iffi1.x  * iffi2.x2 * fpair[k][1] +  //integral  Tjji 
		iffi1.x3 * iffi2.x  * fpair[k][2] +  //integral  Uiiij 
		iffi1.x  * iffi2.x3 * fpair[k][3] +  //integral  Ujjji 
		iffi1.x2 * iffi2.x2 * fpair[k][4] ;  //integral  Uiijj
		
		double denom = levels_1d_scf[i][l]+levels_1d_scf[j][m]-levels_1d_scf[i][csf_scf[state][i]]-levels_1d_scf[j][csf_scf[state][j]];

//		if ((l!=csf_scf[state][i]) && (m!=csf_scf[state][j]))
//		{	
			alt_mp2 -= dtemp * dtemp / denom;
			if (fabs(denom) < 1.0e-5) printf("warning: small denominator in mp2, mode1 = %2d, mode2 = %2d, l =%3d, m=%3d, denom = %10.5e\n",i,j,l,m,denom);
//		}

		}}// end of for l<nbas
	k++;
} // end of for i< nmodes-1

//printf("\n Alternative full basis mp2 = %15.10f, Total energy = %15.10f\n",alt_mp2,mp1_energy+alt_mp2); // */
printf ("\n MP2 correction : SCF+MP1  = %15.10f, MP2 correction = %15.10f, Total = %15.10f\n",mp1_energy,alt_mp2,mp1_energy+alt_mp2);

levels_scf[state]=mp1_energy;
levels_mp2[state]=mp1_energy+alt_mp2;

//Restore WFN

for(i=0;i<nmodes;i++)
	for(j=0;j<nbas;j++)
		for(k=0;k<nbas;k++)
			wfn[i][j][k]=wfn_1d[i][j][k];

// Transform integrals
// TODO - cut transformation for current CSF only 

for(mode=0;mode<nmodes;mode++)
	{
//	for(k=0;k<excit[mode];k++)
//		for(l=0;l<excit[mode];l++)
//		{
		k = l = (int) csf_scf[state][mode];
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
//		}
	}

} // end of for(state

// Print final result

printf("\n Calculated transition energies (relative to ground state, cm-1):\n\n");

for(i=1;i<max_states;i++)
{
	printf(" #%3d, CSF:  ",i);

	double e_harm=0.0;
	double e_diag=0.0;

	for (j=0;j<nmodes;j++)
	{
		e_harm += harmonic[nmodes-j-1] * (double)csf_scf[i][j];
		e_diag += (levels_1d[j][csf_scf[i][j]] - levels_1d[j][0]);
		printf(" %2d",(int) csf_scf[i][j]);
	}

	printf("   harmonic= %10.2f, diagonal = %10.2f, vci = %10.2f, vscf = %10.2f, vmp2 = %10.2f\n",
	e_harm,e_diag*CM1,(levels_ci[i]-levels_ci[0])*CM1,(levels_scf[i]-levels_scf[0])*CM1,(levels_mp2[i]-levels_mp2[0])*CM1);
}

// End of program - free memory

mkl_free(wfn); mkl_free(ss); free(csf); free(csf_scf); 
mkl_free(levels_1d); mkl_free(fdiag); mkl_free(fpair);
mkl_free(ifi); mkl_free(h); mkl_free(ifit); mkl_free(levels_ci);
free(levels_scf); free(levels_mp2);

// print timing

curr_time  = time(NULL);
str = ctime (&curr_time);

printf("\n Subprogram finished at %s\n",str);

timestamp_();

return 0;	
}

static int icheck(const int nmodes,unsigned char *csf,unsigned char *csf_scf,int *diff,int *ndiff)
{
int i,j;	

j=0;
*ndiff=nmodes;

for	(i=0;i<nmodes;i++)
{
	if (csf_scf[i]!=csf[i])
		{
			if (j>1) return 1; // for speed, break if more than 2 mode diff - can be easily changed to any mode diff (take care of diff[] array size!!)
			diff[j]=i;
			j++;
		}
}

*ndiff=j;

return 0;
	
}

static double pot_1d(const int nmodes,const int mode,const double x, double fdiag[nmodes][4])
{
double result;
	
result= x * x * (x * (x * fdiag[mode][2] + fdiag[mode][1]) + fdiag[mode][0]);	
	
return result;	
}

