/*
 * 1d_polynom_points.c - 1D solvers with polynomial and pointwise potential
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
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl.h>
#include <quadmath.h>

#define PI 3.14159265358979323846 	// just pi
#define MAX_POW 32 					// max power of polynomial
static const double C=0.7;     		// basis set broadness

//* Function declaration, print_and_punch.c
void printmat1(const int dim, double mat[][dim]); 

//* Function declaration, parser.c
int parse_points(const char *filename, const int nbas, double pot[]); 
int parse_power(const char *filename, const int npow, double fi[]); 

static inline void dvand (const int n,const double alpha[MAX_POW], double x[MAX_POW] );
static inline double poly_d(const int pow,const double fi[MAX_POW],const double x);

/******************************************************************************/

void one_d(const char *filename, const double xstart, const double deltax, const int nb, const int pow) 

/******************************************************************************/
/* 
 * Purpose
 *   calculate 1D energy levels
 *   kinetic energy is simple diagonal -1/2
 *   potential is given as analytic polynom
 *   basis is DGB centered at equidistant points
 * 
 *   Potential is integrated analytically
 *   No interpolation is done
 *  
 */
{
//double xstart,xend,deltax;    //* basis parameters 
double x[nb],a;             	//* basis parameters - small, keep them in stack 
double vtemp[nb],levels[nb];    //* temporal array and array with energy levels          
double dtemp1,dtemp2,dtemp3,dtemp4;      //* temporal values 
double ifi[MAX_POW];            //* integrals 
// we keep any variable that is connected to array size as lapack_int (mkl_int) - for easy 32/64-bit switch
//int nbas;                     //* amount of basis centers 
//int i,j,k;                    //* loop counters
lapack_int nbas,i,j,k;
//int pow=MAX_POW;  			// pow - actual max power of polynomial to be used
double fi[MAX_POW];				// polynome coefficients
lapack_int info=0;   			// int LAPACK info

// check pow
if (pow>MAX_POW)
{
	printf("\n !!!!Error!!!\n In one_d_polynom: pow is greater than MAX_POW\n");
	printf(" pow=%5d; MAX_POW=%5d \n",pow,MAX_POW);
	exit(1);
}

// tune sizing
nbas = (lapack_int) nb;
if (nbas<pow) 
{
	printf("\n !!!!Error!!!\n In one_d_polynom: pow is greater than nbas (= not enough data for interpolation)\n");
	printf(" pow=%5d; nbas=%5d \n",pow,nbas);
	exit(1);
}


// print banner
printf ("\n");
printf (" !=====================================================================================\n");
printf (" !  1D subprogram for polynomial potential                   \n");
printf (" !  The potential is analytic polynom (read from {polynom} group)   \n");
printf (" !  Integrals over polynom are analytical               \n");
printf ("\n");

// FILL in X and A values for equidistant basis
for (i = 0; i < nbas; ++i)
	x[i] = xstart + deltax * i;

a = C * C /(deltax * deltax);

// FILL default potential
for(i=0;i<pow;i++) fi[i]=0;
fi[2]=0.5;

// READ potential
printf("\n Reading {polynom} group... ");
parse_power(filename,pow,fi); 
printf("done \n\n");


printf("Complete {polynom} group:\n");
printf("{polynom\n");
for(i=0;i<pow;i++)
	printf(" %3i %+20.15e\n",i,fi[i]);
printf("}\n");


// ALLOCATE memory for large S (overlap) and H (hamiltonian) matrices, in C99 style VLAs - nbas is known only at runtime 

double (*s)[nbas] = malloc( sizeof(double[nbas][nbas]));
if (s==NULL){
	printf(" Error while allocating memory for s \n");
	exit (1);
}

double (*h)[nbas] = malloc( sizeof(double[nbas][nbas]));
if (h==NULL){
	printf(" Error while allocating memory for h \n");
	exit (1);
}


//DO CALCULATION

// FILL in S,H 

		dtemp2 = a + a; // a[i] + a[j];
		dtemp3 = a * a; // a[i] * a[j];
	double bi  = 2.0 * dtemp2;		// 2* exponent coefficient

for(i=0;i<nbas;++i){
	for(j=0;j<nbas;++j){
		dtemp1 = x[i] - x[j];
		dtemp4 = a * (x[i] + x[j]);//a[i] * x[i] + a[j] * x[j];
		
/*    Compact code for overlap matrix S
 * S(I,J)=sqrt(sqrt(A(I)*A(J))*TWO/(A(I)+A(J)))*exp(-A(I)*A(J)*(X(I)-X(J))**2/(A(I)+A(J))) */
		s[i][j]=sqrt(2.0 * sqrt(dtemp3) / dtemp2) * 
		exp(- dtemp3 * dtemp1 *dtemp1 / dtemp2);

/*  special code for kinetic energy - useful for polynomial expansions - and symmetric
    -(dphi/dx)^2 * 0.5= (2 AI*AJ *x^2 - 2AI*AJ *x *(XI+XJ) + 2(AI*AJ*XI*XJ))* phi^2 

    original code for power integrals
    AI=A(i)*X(i)+A(j)*X(j)
    BI=A(i)+A(j)
    IFI(1)=AI/BI 
    IFI(2)=ONE/(TWO*BI)+IFI(1)*AI/BI - 
    do k=3,norder
      IFI(k)=IFI(k-2)*dble(k-1)/(TWO*BI)+IFI(k-1)*AI/BI
    end do
*/
		double xij = dtemp4 / dtemp2; 	// point of expansion
        ifi[0] = 1.0;
		ifi[1] = xij;                     // integral <phi*x  *phi>
		ifi[2] = ifi[0] / bi + ifi[1] * xij; // integral <phi*x^2*phi>

		h[i][j] = (    fi[0]
		+     ifi[1] * fi[1]						// potential, const and x part (normally can be ignored)
		+     ifi[2] * fi[2]                     	// potential, x^2 part
		+ 2.0*ifi[2] * dtemp3               	// kinetic energy , x^2 part
		- 2.0*ifi[1] * dtemp3 * (x[i] + x[j])	// kinetic energy , x   part 
		+ 2.0*ifi[0] * dtemp3 *  x[i] * x[j]);	// kinetic energy , const part

		for(k=3;k<pow;++k){
			ifi[k] = ifi[k-1] * xij + ifi[k-2] / bi * (double) (k-1); 	// integral <phi*x^k*phi>
			h[i][j] += fi[k] * ifi[k];								    // potential, x^k part
		} 
		
		h[i][j] *= s[i][j];

	}
}

//printf("\n   S matrix printout \n");
//printmat1(nbas,s);
		
//printf("\n   H matrix printout \n");
//printmat1(nbas,h);
//END OF FILL

// DIAGONALISE and print

//lapack_int info;
lapack_int itype=1;
char jobz='V',uplo='U';

info=LAPACKE_dsygv(LAPACK_ROW_MAJOR,itype,jobz,uplo,nbas,&h[0][0],nbas,&s[0][0],nbas,levels);
/* lapack_int LAPACKE_dsygv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, double* a, lapack_int lda,
                          double* b, lapack_int ldb, double* w ); */

if (info != 0) printf("\n!! error in dsygv call, info = %d\n",info);

printf("\n Eigenvalues found:\n\n");

for(i=0;i<nbas;++i){
	printf("%4d %20.15f\n",i,levels[i]);
}

printf("\n");


printf("\n   Eigenvector printout \n");
printmat1(nbas,h);

// END OF DIAGONALISE

// CREATE the matrix of wavefunction values
// we will save on one matrix (and pay more computation cost), and compute phi(xi) on the fly, in vtemp

for(i=0;i<nbas;++i){
	for(j=0;j<nbas;++j){
		dtemp1=x[i]-x[j];
		vtemp[j]=sqrt(sqrt(2.0 * a / PI)) * exp (- a * dtemp1 *dtemp1);
	}
	for(k=0;k<nbas;++k){
		s[i][k]=0.0;
		for(j=0;j<nbas;++j)
			s[i][k]	+= h[j][k]*vtemp[j];
	}
}		

printf("\n   Wavefunction values at basis centers \n");
printmat1(nbas,s);

// END OF CREATE 
	
// END OF CALCULATE

free(h); free(s);
return;
}



/******************************************************************************/

void one_d_polynom(const char *filename, const double xstart, const double deltax, const int nb, const int pow) 

/******************************************************************************/
/* 
 * Purpose
 *   calculate 1D energy levels
 *   kinetic energy is simple diagonal -1/2
 *   potential is given as equidistant points
 *   basis is DGB centered at potential points
 * 
 *   Potential is interpolated by one polynom centered at x=0, the polynom is integrated analytically
 *   Interpolation is done by Jacobi SVD from Intel MKL
 *  
 */
{
//double xstart,xend,deltax;     //* basis parameters 
double x[nb],a;             	 //* basis parameters - small, keep them in stack 
double pot[nb]; 				 //* potential values and polynom data  
double vtemp[nb],levels[nb];     //* temporal array and array with energy levels          
double dtemp1,dtemp2,dtemp3,dtemp4;      //* temporal values 
double ifi[MAX_POW];             //* integrals 
// we keep any variable that is connected to array size as lapack_int (mkl_int) - for easy 32/64-bit switch
//int nbas;                      //* amount of basis centers 
//int i,j,k;                     //* loop counters
lapack_int nbas,i,j,k;
//int pow=MAX_POW;  					 	// pow - actual max power of polynomial to be used
double fi[MAX_POW];					 	// polynome coefficients
lapack_int info=0,istat[3];   			// int SVD stuff, and pow - actual max power of polynomial to be used
double scale,sva[MAX_POW],stat[7];      // singular values (scale,sva) and other SVD-related stuff
double v[MAX_POW][MAX_POW];             // array for right singular vectors
//__float128 qai,qbi,qfi[MAX_POW],qifi[MAX_POW];  //increased precision variables

// check pow
if (pow>MAX_POW)
{
	printf("\n !!!!Error!!!\n In one_d_polynom: pow is greater than MAX_POW\n");
	printf(" pow=%5d; MAX_POW=%5d \n",pow,MAX_POW);
	exit(1);
}

// tune sizing
nbas = (lapack_int) nb;
if (nbas<pow) 
{
	printf("\n !!!!Error!!!\n In one_d_polynom: pow is greater than nbas (= not enough data for interpolation)\n");
	printf(" pow=%5d; nbas=%5d \n",pow,nbas);
	exit(1);
}


// print banner
printf ("\n");
printf (" !=====================================================================================\n");
printf (" !  1D subprogram for pointwise potential                   \n");
printf (" !  The potential is interpolated by polynomial at x=0 using Jacobi SVD algorhitm    \n");
printf (" !  Integrals over polynomial are analytical               \n");
printf ("\n");

// FILL in X and A values for equidistant basis
for (i = 0; i < nbas; ++i)
	x[i] = xstart + deltax * i;

a = C * C /(deltax * deltax);

// FILL in default (quadratic) potential
for (i = 0; i < nbas; ++i)
	pot[i] = 0.5 * x[i] * x[i];

// Read potential
printf("\n Reading {points} group... ");
parse_points(filename,nb,pot); 
printf("done \n");

// ALLOCATE memory for large S (overlap) and H (hamiltonian) matrices, in C99 style VLAs - nbas is known only at runtime 

double (*s)[nbas] = malloc( sizeof(double[nbas][nbas]));
if (s==NULL){
	printf(" Error while allocating memory for s \n");
	exit (1);
}

double (*h)[nbas] = malloc( sizeof(double[nbas][nbas]));
if (h==NULL){
	printf(" Error while allocating memory for h \n");
	exit (1);
}

// ALLOCATE Vandermonde matrix in double (dvm)

double (*dvm)[MAX_POW] = malloc( sizeof(double[nbas][MAX_POW]));
if (dvm==NULL){
	printf(" Error while allocating memory for h \n");
	exit (1);
}

// ALLOCATE U matrix - left singular vectors

double (*u)[nbas] = malloc( sizeof(double[nbas][nbas]));
if (u==NULL){
	printf(" Error while allocating memory for u \n");
	exit (1);
}

// END OF ALLOCATE 

// Print potential and basis parameters at basis centers

printf("\n   N     X value     A value  Potential(x)\n");
for (i = 0; i < nbas; ++i){
	printf(" %3d  %10.6f  %10.6f  %12.6f\n",i+1,x[i],a,pot[i]);
}
printf("\n Complete {points} group:\n");
printf("{points\n");
for(i=0;i<nbas;i++)
	printf(" %3i %+20.15e\n",i,pot[i]);
printf("}\n");
//END OF PRINT

// BUILD full Vandermonde martix

for(i=0;i<nbas;++i){
	dvm[i][0]=1.0;
	for(j=1;j<MAX_POW;++j)
		dvm[i][j] = dvm[i][j-1] * x[i];
	
}


// printf("\n  Full Vandermonde matrix\n\n");

// printmat2d(nbas,dvm);

// DO SVD
/*lapack_int LAPACKE_dgejsv( int matrix_layout, char joba, char jobu, char jobv,
*                           char jobr, char jobt, char jobp, lapack_int m,
*                           lapack_int n, double* a, lapack_int lda, double* sva,
*                           double* u, lapack_int ldu, double* v, lapack_int ldv,
*                           double* stat, lapack_int* istat ); */

char joba='F';          // set to R to remove small singular values
char jobu='F',jobv='V'; //request full set of vectors
char jobr='N';          // set to R to remove large singular values
char jobt='N',jobp='P';

printf("\n\n CALLING LAPACK \n\n");

info=LAPACKE_dgejsv(LAPACK_ROW_MAJOR,joba,jobu,jobv,jobr,jobt,jobp,
					nbas,pow,&dvm[0][0],MAX_POW,sva,&u[0][0],nbas,&v[0][0],MAX_POW,stat,istat);

if (info != 0) printf("\n!! error in dsygv call, info = %d\n",info);

printf(" LAPACK SVD (dgejsv) completed!\n\n Stat array values = ");
for(i=0;i<7;i++)
	printf("%10.4e",stat[i]);
printf("\n");

printf("Istat array values = ");
for(i=0;i<3;i++)
	printf("%10d",istat[i]);
printf("\n\n");

scale=stat[0]/stat[1];

printf("Singular values found: \n");

for(i=0;i<pow;i++)
	printf("%3d %+15.10e\n",i,sva[i]*scale);
printf("\n");

// END OF SVD

// FIND SOLUTION - transform pot to coefficients
// vtemp=U(T)*pot  
/*for(i=0;i<nbas;i++){
	vtemp[i]=0.0;
	for(j=0;j<nbas;j++)
		vtemp[i]+=u[j][i]*pot[j];
}

// vtemp=Sigma(-1)*vtemp  
for(i=0;i<pow;i++)
    vtemp[i]=vtemp[i]/(sva[i]*scale);

// fi=V*vtemp=V*Sigma(-1)*U(T)*pot  
for(i=0;i<pow;i++){
	fi[i]=0.0;
	for(j=0;j<pow;j++)
		fi[i]+=v[i][j]*vtemp[j];
}
*/

// alternative code - we create single transformation matrix V*Sigma(-1)*U(T) in dvm for future use
// we found it numerically less stable

double (*o)[nbas] = malloc( sizeof(__float128[nbas][nbas]));
if (o==NULL){
	printf(" Error while allocating memory for o \n");
	exit (1);
}

// o(T)= V*Sigma(-1)*U(T)
for(i=0;i<nbas;++i){
	for(j=0;j<pow;++j){
		o[i][j]=(__float128)0.0;
		for(k=0;k<pow;++k){
			o[i][j] += (__float128) v[j][k] * (__float128) u[i][k] / (__float128) (sva[k]*scale);
		}
	}
}

//transform pot to coefficients
//fi=dvm(T)*pot
for(i=0;i<pow;++i){
	fi[i]=0.0;
	for(j=0;j<nbas;++j)
		fi[i]+= (double) (o[j][i]* (__float128) pot[j]);
}
// end of alternative code 

printf("Power coefficients found: \n");

printf("\n  N              Fit   \n");

for(i=0;i<pow;i++) 
	printf("F%2d = %+15.10e\n",i,fi[i]);

printf("\n");


// END OF SOLUTION

// PRINT FIT
printf("Fit accuracy : \n");
printf("\n  N            x    V(given)        V(fit)         Diff\n");

dtemp1=0.0;

for(i=0;i<nbas;i++){
 dtemp2 = poly_d(pow,fi,x[i]);
 dtemp3 = pot[i]-dtemp2;
 dtemp1 += dtemp3*dtemp3;
 printf("%3d  %+10.4e  %+10.4e  %+10.4e  %+10.4e\n",i,x[i],pot[i],dtemp2,dtemp3);
}

printf("Sigma (sqrt(sum((V(given)-V(fit))**2))) = %+10.4e\n\n",dtemp1);

printf("Complete polynom group:\n");
printf("{polynom\n");
for(i=0;i<pow;i++)
	printf(" %3i %+20.15e\n",i,fi[i]);
printf("}\n");

//END OF PRINT FIT

//DO CALCULATION

// FILL in S,H 

		dtemp2 = a + a; // a[i] + a[j];
		dtemp3 = a * a; // a[i] * a[j];
	double bi  = 2.0 * dtemp2;		// 2* exponent coefficient

for(i=0;i<nbas;++i){
	for(j=0;j<nbas;++j){
		dtemp1 = x[i] - x[j];
		dtemp4 = a * (x[i] + x[j]);//a[i] * x[i] + a[j] * x[j];
		
/*    Compact code for overlap matrix S
 * S(I,J)=sqrt(sqrt(A(I)*A(J))*TWO/(A(I)+A(J)))*exp(-A(I)*A(J)*(X(I)-X(J))**2/(A(I)+A(J))) */
		s[i][j]=sqrt(2.0 * sqrt(dtemp3) / dtemp2) * 
		exp(- dtemp3 * dtemp1 *dtemp1 / dtemp2);

/*  special code for kinetic energy - useful for polynomial expansions - and symmetric
    -(dphi/dx)^2 * 0.5= (2 AI*AJ *x^2 - 2AI*AJ *x *(XI+XJ) + 2(AI*AJ*XI*XJ))* phi^2 

    original code for power integrals
    AI=A(i)*X(i)+A(j)*X(j)
    BI=A(i)+A(j)
    IFI(1)=AI/BI 
    IFI(2)=ONE/(TWO*BI)+IFI(1)*AI/BI - 
    do k=3,norder
      IFI(k)=IFI(k-2)*dble(k-1)/(TWO*BI)+IFI(k-1)*AI/BI
    end do
*/
		double xij = dtemp4 / dtemp2; 	// point of expansion
        ifi[0] = 1.0;
		ifi[1] = xij;                     // integral <phi*x  *phi>
		ifi[2] = ifi[0] / bi + ifi[1] * xij; // integral <phi*x^2*phi>

		h[i][j] = (    fi[0]
		+     ifi[1] * fi[1]						// potential, const and x part (normally can be ignored)
		+     ifi[2] * fi[2]                     	// potential, x^2 part
		+ 2.0*ifi[2] * dtemp3               	// kinetic energy , x^2 part
		- 2.0*ifi[1] * dtemp3 * (x[i] + x[j])	// kinetic energy , x   part 
		+ 2.0*ifi[0] * dtemp3 *  x[i] * x[j]);	// kinetic energy , const part

		for(k=3;k<pow;++k){
			ifi[k] = ifi[k-1] * xij + ifi[k-2] / bi * (double) (k-1); 	// integral <phi*x^k*phi>
			h[i][j] += fi[k] * ifi[k];								    // potential, x^k part
		} 
		
		h[i][j] *= s[i][j];

	}
}

//printf("\n   S matrix printout \n");
//printmat1(nbas,s);
		
//printf("\n   H matrix printout \n");
//printmat1(nbas,h);
//END OF FILL


// DIAGONALISE and print

//lapack_int info;
lapack_int itype=1;
char jobz='V',uplo='U';

info=LAPACKE_dsygv(LAPACK_ROW_MAJOR,itype,jobz,uplo,nbas,&h[0][0],nbas,&s[0][0],nbas,levels);
/* lapack_int LAPACKE_dsygv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, double* a, lapack_int lda,
                          double* b, lapack_int ldb, double* w ); */

if (info != 0) printf("\n!! error in dsygv call, info = %d\n",info);

printf("\n Eigenvalues found:\n\n");

for(i=0;i<nbas;++i){
	printf("%4d %20.15f\n",i,levels[i]);
}

printf("\n");


printf("\n   Eigenvector printout \n");
printmat1(nbas,h);

// END OF DIAGONALISE

// CREATE the matrix of wavefunction values
// we will save on one matrix (and pay more computation cost), and compute phi(xi) on the fly, in vtemp

for(i=0;i<nbas;++i){
	for(j=0;j<nbas;++j){
		dtemp1=x[i]-x[j];
		vtemp[j]=sqrt(sqrt(2.0 * a / PI)) * exp (- a * dtemp1 *dtemp1);
	}
	for(k=0;k<nbas;++k){
		s[i][k]=0.0;
		for(j=0;j<nbas;++j)
			s[i][k]	+= h[j][k]*vtemp[j];
	}
}		

printf("\n   Wavefunction values at basis centers \n");
printmat1(nbas,h);

// END OF CREATE 
	
// END OF CALCULATE


free(dvm); free(u);
free(s); free(h);
return;
}

/******************************************************************************/

void one_d_points(const char *filename, const double xstart, const double deltax, const int nb) 

/******************************************************************************/
/* 
 * Purpose
 *   calculate 1D energy levels
 *   kinetic energy is simple diagonal -1/2
 *   potential is given as equidistant points
 *   basis is DGB centered at potential points
 * 
 *   Potential intergals are evaluated by seies of local interpolations using Bjorck-Pereyra algorhitm 
 *   at each basis center, followed by local analytical integration 
 * 
 */
{
//double xstart,xend,deltax;     //* basis parameters 
double x[nb],a;             	 //* basis parameters - small, keep them in stack 
double pot[nb]; 				 //* potential values and polynom data  
double vtemp[nb],levels[nb];     //* temporal array and array with energy levels          
double dtemp1,dtemp2,dtemp3,dtemp4;      //* temporal values 
double ifi;             //* integrals 
// we keep any variable that is connected to array size as lapack_int (mkl_int) - for easy 32/64-bit switch
//int nbas;                      //* amount of basis centers 
//int i,j,k;                     //* loop counters
lapack_int nbas,i,j,k;
int pow=MAX_POW;  					 // pow - actual max power of polynomial to be used
double fi[MAX_POW],alpha[MAX_POW];	 // polynome coefficients
double fii[nb+nb-1];					 // intermediate integrals

// tune sizing
nbas = (lapack_int) nb;
if (nbas<pow) pow=nbas; 

// print banner
printf ("\n");
printf (" !=====================================================================================\n");
printf (" !  1D subprogram for pointwise potential                   \n");
printf (" !  The potential is interpolated by polynomials locally using Bjorck-Pereyra algorhitm    \n");
printf (" !  Integrals on polynomials are analytical               \n");
printf ("\n");

// FILL in X and A values for equidistant basis
for (i = 0; i < nbas; ++i)
	x[i] = xstart + deltax * i;

a = C * C /(deltax * deltax);

// FILL in default (quadratic) potential
for (i = 0; i < nbas; ++i)
	pot[i] = 0.5 * x[i] * x[i];

// Read potential
printf("\n Reading {points} group... ");
parse_points(filename,nb,pot); 
printf("done \n");

// ALLOCATE memory for large S (overlap) and H (hamiltonian) matrices, in C99 style VLAs - nbas is known only at runtime 

double (*s)[nbas] = malloc( sizeof(double[nbas][nbas]));
if (s==NULL){
	printf(" Error while allocating memory for s \n");
	exit (1);
}

double (*h)[nbas] = malloc( sizeof(double[nbas][nbas]));
if (h==NULL){
	printf(" Error while allocating memory for h \n");
	exit (1);
}

// END OF ALLOCATE 

// Print potential and basis parameters at basis centers

printf("\n   N     X value     A value  Potential(x)\n");
for (i = 0; i < nbas; ++i){
	printf(" %3d  %10.6f  %10.6f  %12.6f\n",i+1,x[i],a,pot[i]);
}
printf("\n Complete {points} group:\n");
printf("{points\n");
for(i=0;i<nbas;i++)
	printf(" %3i %+20.15e\n",i,pot[i]);
printf("}\n");
//END OF PRINT

//DO CALCULATION

// FILL in S,H 

printf("\n Calculating S and H matrices\n");

printf("\n For equidistant basis, sigma = %15e, tenth-wide interval = %15e, while deltax = %15e\n",0.5/sqrt(a),4.29193*0.5/sqrt(a),deltax);

dtemp2 = a + a;
dtemp3 = a * a;

// Fill in fij
// Bjorck-Pereyra fit for local derivatives + analytic integration over taylor series

for(i=0;i<8;i++)
{
	double xij = xstart + deltax * 0.5 * (double)i ; // point of expansion
//	printf("\n i = %2i  xij = %+15.10e\n",i,xij);
	
	for(k=0;k<9;k++) 	// fill in data for Bjorck-Pereyra fit
	{
		alpha[k] = x[k] - xij;
		fi[k] = pot[k];
//		printf(" k = %2i  alpha[k] + xij = %15.10e  fi[k] = %15.10e\n",k,alpha[k]+xij,fi[k]);
	}

	dvand(9,alpha,fi); // Do Bjorck-Pereyra fit

	ifi=1.0; 												// first integral x^0 = sij * const
	fii[i] = fi[0];
	
	for(k=2;k<9;k+=2)	// Integrate
	{
		ifi *= 0.5 / dtemp2 * (double)(k-1);
		fii[i] +=fi[k]*ifi;
	} 
}

for(i=8;i<(nbas+nbas-9);i++)
{
	double xij = xstart + deltax * 0.5 * (double)i ; // point of expansion
//	printf("\n i = %2i  xij = %+15.10e\n",i,xij);

	j=i/2-4;

	for(k=0;k<9;k++) 	// fill in data for Bjorck-Pereyra fit
	{
		alpha[k] = x[k+j] - xij;
		fi[k] = pot[k+j];
//		printf(" k = %2i  alpha[k] = %15.10e  fi[k] = %15.10e\n",k,alpha[k]+xij,fi[k]);
	}

	dvand(9,alpha,fi); // Do Bjorck-Pereyra fit

	ifi=1.0; 												// first integral x^0 = sij * const
	fii[i] = fi[0];
	
	for(k=2;k<9;k+=2)	// Integrate
	{
		ifi *= 0.5 / dtemp2 * (double)(k-1);
		fii[i] +=fi[k]*ifi;
	} 
}

j=(nbas-9);

for(i=(nbas+nbas-9);i<(nbas+nbas-1);i++)
{
	double xij = xstart + deltax * 0.5 * (double)i ; // point of expansion
//	printf("\n i = %2i  xij = %+15.10e\n",i,xij);

	for(k=0;k<9;k++) 	// fill in data for Bjorck-Pereyra fit
	{
		alpha[k] = x[k+j] - xij;
		fi[k] = pot[k+j];
//		printf(" k = %2i  alpha[k] = %15.10e  fi[k] = %15.10e\n",k,alpha[k]+xij,fi[k]);
	}

	dvand(9,alpha,fi); // Do Bjorck-Pereyra fit

	ifi=1.0; 												// first integral x^0 = sij * const
	fii[i] = fi[0];
	
	for(k=2;k<9;k+=2)	// Integrate
	{
		ifi *= 0.5 / dtemp2 * (double)(k-1);
		fii[i] +=fi[k]*ifi;
	} 
}

/*
printf("\n Fii vector:\n");
for(i=0;i<nbas+nbas-1;i++)
	printf(" %2i   %+10.15e\n",i,fii[i]);
*/

// start filling
for(i=0;i<nbas;++i)
{
	for(j=0;j<nbas;++j)
	{
		dtemp1 = x[i] - x[j];
		dtemp4 = a * x[i] + a * x[j];
		
/*    Compact code for overlap matrix S
 * S(I,J)=sqrt(sqrt(A(I)*A(J))*TWO/(A(I)+A(J)))*exp(-A(I)*A(J)*(X(I)-X(J))**2/(A(I)+A(J))) */
		s[i][j]=sqrt(2.0 * sqrt(dtemp3) / dtemp2) * exp(- dtemp3 * dtemp1 * dtemp1 / dtemp2);

/*  Original, compact code for kinetic energy
    G(I,J)=-S(I,J)*A(I)*A(J)*(TWO*A(I)*A(J)*(X(I)-X(J))**2-A(I)-A(J))/(A(I)+A(J))**2 */ 
        h[i][j] = - dtemp3 * (2.0*dtemp3*dtemp1*dtemp1 - dtemp2) / (dtemp2 * dtemp2); 

		h[i][j] += fii[i+j];
		h[i][j] *= s[i][j];
	}
}

//END OF FILL

// DIAGONALISE and print

lapack_int info=0;
lapack_int itype=1;
char jobz='V',uplo='U';

info=LAPACKE_dsygv(LAPACK_ROW_MAJOR,itype,jobz,uplo,nbas,&h[0][0],nbas,&s[0][0],nbas,levels);
/* lapack_int LAPACKE_dsygv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, double* a, lapack_int lda,
                          double* b, lapack_int ldb, double* w ); */

if (info != 0) printf("\n!! error in dsygv call, info = %d\n",info);

printf("\n Eigenvalues found:\n");

for(i=0;i<nbas;++i){
	printf("%3d %20.15f\n",i,levels[i]);
}

printf("\n");


printf("\n   Eigenvector printout \n");
printmat1(nbas,h);

// END OF DIAGONALISE

// CREATE the matrix of wavefunction values
// we will save on one matrix (and pay more computation cost), and compute phi(xi) on the fly, in vtemp

for(i=0;i<nbas;++i){
	for(j=0;j<nbas;++j){
		dtemp1=x[i]-x[j];
		vtemp[j]=sqrt(sqrt(2.0 * a / PI)) * exp (- a * dtemp1 *dtemp1);
	}
	for(k=0;k<nbas;++k){
		s[i][k]=0.0;
		for(j=0;j<nbas;++j)
			s[i][k]	+= h[j][k]*vtemp[j];
	}
}		

printf("\n   Wavefunction values at basis centers \n");
printmat1(nbas,h);

// END OF CREATE 
	
// END OF CALCULATE

free(s); free(h);

}

/******************************************************************************/

static inline void dvand (const int n,const double alpha[MAX_POW], double x[MAX_POW] )

/******************************************************************************/
/*
  Purpose:

    DVAND solves a Vandermonde system A' * x = b.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 February 2014

  Author:

    John Burkardt

  Reference:

    Ake Bjorck, Victor Pereyra,
    Solution of Vandermonde Systems of Equations,
    Mathematics of Computation,
    Volume 24, Number 112, October 1970, pages 893-903.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double ALPHA[N], the parameters that define the matrix.
    The values should be distinct.

    Input, double B[N], the right hand side of the linear system.

    Output, double DVAND[N], the solution of the linear system.
*/
{
  int j;
  int k;

//  for (j = 0; j < n; j++ ) x[j] = b[j];

  for ( k = 0; k < n - 1; k++ )
  {
    for ( j = n - 1; k < j; j-- )
    {
      x[j] = ( x[j] - x[j-1] ) / ( alpha[j] - alpha[j-k-1] );
    }
  }

  for ( k = n - 2; 0 <= k; k-- )
  {
    for ( j = k; j < n - 1; j++ )
    {
      x[j] = x[j] - alpha[k] * x[j+1];
    }
  }

}

/******************************************************************************/

static inline double poly_d(const int pow, const double fi[MAX_POW],const double x)

/******************************************************************************/
/*
 * Purpose:
 * 	Returns value of polynom max power pow, coefficients in fi, at point x
*/
{
//"Naive" but working algorhitm
/*double temp1,temp2;
int l;

temp1 = fi[0];
temp2 = 1.0;
for(l=1;l<pow;l++){
	temp2 = temp2 * x;
	temp1 += temp2*fi[l];
}
	
return temp1;
*/
//
// Horner algorhitm
// assumes pow >1 !
//

double temp1;
int l;

temp1 = fi[pow-1];

for (l=pow-2;l>0;l--)
	temp1 = fi[l] + x * temp1;

temp1 = fi[0] + x * temp1;

return temp1; 
}
