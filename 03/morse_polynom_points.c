/*
 * morse_polynom_points.c - test function for interpolating morse potential by various polynoms
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl.h>
//#include <quadmath.h>

#define PI 3.14159265358979323846 	// just pi
#define MAX_POW 32 					// max power of polynomial
static const double C=0.7;     		// basis set broadness

//* Function declaration, print_and_punch.c
void printmat1(const int dim, double mat[][dim]); 
void printmat2(const int dim, double mat[][dim]);

static void printmat2d(const lapack_int dim,const double mat[dim][MAX_POW]); 
static inline double poly_d(const lapack_int pow,const double fi[MAX_POW],const double x);
static void dvand (const int n,const double alpha[MAX_POW], double x[MAX_POW] );

// static void printmat2ld(const lapack_int dim,const __float128 mat[dim][npow]);

/******************************************************************************/

void morse_points(double d, double w, double xstart,double deltax,int nb) 

/******************************************************************************/
{
//double xstart,xend,deltax;     //* basis parameters 
double x[nb],a;             	 //* basis parameters - small, keep them in stack 
double pot[nb]; 				 //* potential values and data for  
double vtemp[nb],levels[nb];     //* temporal array and array with energy levels          
double dtemp1,dtemp2,dtemp3,dtemp4;      //* temporal values 
double ifi;             //* integrals 
// we keep any variable that is connected to array size as lapack_int (mkl_int) - for easy 32/64-bit switch
//int nbas;                      //* amount of basis centers 
//int i,j,k;                     //* loop counters
lapack_int nbas,i,j,k;
int pow=MAX_POW+1;  					 // pow - actual max power of polynomial to be used
double fi[MAX_POW+1],alpha[MAX_POW+1];	 // polynome coefficients
double fii[nb+nb-1];					 // intermediate integrals

// tune sizing
nbas = (lapack_int) nb;
if (nbas<pow) pow=nbas; 

// print banner
printf ("\n");
printf (" !=================================================================================\n");
printf (" !  1D subprogram for Morse potential                   \n");
printf (" !  The potential is interpolated by polynom locally using Bjorck-Pereyra algorhitm    \n");
printf (" !  Integrals with polynom are analytical               \n");
printf ("\n");



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


// FILL in X and A values for equidistant basis

for (i = 0; i < nbas; ++i)
	x[i] = xstart + deltax * i;

a = C * C /(deltax * deltax);

//* Print potential and basis parameters at basis centers

printf("\n   N     X value     A value  Potential(x)\n");
for (i = 0; i < nbas; ++i){
	dtemp1 = d * (exp(-2.0 * w * x[i]) - 2.0* exp(- w * x[i] ));
	printf(" %3d  %10.6f  %10.6f  %12.6f\n",i+1,x[i],a,dtemp1);
	pot[i]=dtemp1;
}

printf("\n Choice of points: \n");

printf("\n  j   Alpha[j]       B[j]\n");

k=nbas/(pow-1);
alpha[0]=x[0];
fi[0]=pot[0];

j=1;
for(i=k-1;j<pow;i+=k)
{
	alpha[j]=x[i];
	fi[j]=pot[i];
	printf(" %2d %+10.5f %+10.5f\n",j,alpha[j],fi[j]);
	j++;
}

dvand(pow,alpha,fi);


printf("\n Power coefficients found: \n");

printf("\n  N              Fit           Analytical \n");

dtemp2 = - w;
dtemp3 = 1.0;
dtemp4 = 1.0;

for(i=0;i<(pow-1);i++) {
    dtemp1 = 2.0 * d * dtemp2 * (dtemp3 - 1.0) * dtemp4;
    dtemp2 *= - w;
    dtemp3 *= 2.0;
    dtemp4 /= (double) i + 2.0;
	printf("F%2d = %+15.10e %+15.10e\n",i+1,fi[i+1],dtemp1);
//	fi[i]=dtemp1;  //use this to try analytic Taylor expansion
}
printf("\n");


// END OF SOLUTION

// PRINT FIT
printf("Fit accuracy : \n");
printf("\n  N            x    V(analyt)       V(fit)         Diff\n");

dtemp1=0.0;

for(i=0;i<nbas;i++){
 dtemp2 = poly_d((pow-1),&fi[1],x[i]);
 dtemp3 = pot[i]-dtemp2+d;
 dtemp1 += dtemp3*dtemp3;
 printf("%3d  %+10.4e  %+10.4e  %+10.4e  %+10.4e\n",i,x[i],pot[i],dtemp2-d,dtemp3);
}

printf("Sigma (sqrt(sum_square_diff)) = %+10.4e\n\n",dtemp1);


printf("Complete points group:\n");
printf("{points\n");
for(i=0;i<nbas;i++)
	printf(" %3i %+20.15e\n",i,pot[i]);
printf("}\n");
//END OF PRINT

//lapack_int info=0;


//DO CALCULATION

// FILL in S,H 

printf("\n Calculating S and H matrices\n");

printf("\n For equidistant basis, sigma = %15e, tenth-wide interval = %15e, while delatx = %15e\n",0.5/sqrt(a),4.29193*0.5/sqrt(a),deltax);

double sigma = 0.0; // variable for sigma - averaged difference value 
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
// print headers
//printf("\n H element comparison in detail:\n");
//printf("\n   i   j        h_analyt         h[i][j]           diff\n"); 

// start filling
for(i=0;i<nbas;++i)
{
	for(j=0;j<nbas;++j)
	{
		dtemp1 = x[i] - x[j];
		dtemp4 = a * x[i] + a * x[j];
		
/*    Compact code for overlap matrix S
 * S(I,J)=sqrt(sqrt(A(I)*A(J))*TWO/(A(I)+A(J)))*exp(-A(I)*A(J)*(X(I)-X(J))**2/(A(I)+A(J))) */
		s[i][j]=sqrt(2.0 * sqrt(dtemp3) / dtemp2) * 
		exp(- dtemp3 * dtemp1 * dtemp1 / dtemp2);

// analytical expression for H element
		double h_analyt = s[i][j]* (- dtemp3 * (2.0*dtemp3*dtemp1*dtemp1 - dtemp2) / (dtemp2 * dtemp2)
		+ d * exp(-2.0 * w * dtemp4 / dtemp2) 
		* ( exp(w * w / dtemp2) - 2.0 * exp(w * (4.0 * dtemp4 + w) / (4.0*dtemp2))));

/*  Original, compact code for kinetic energy
    G(I,J)=-S(I,J)*A(I)*A(J)*(TWO*A(I)*A(J)*(X(I)-X(J))**2-A(I)-A(J))/(A(I)+A(J))**2 */ 
        h[i][j] = - dtemp3 * (2.0*dtemp3*dtemp1*dtemp1 - dtemp2) / (dtemp2 * dtemp2); 

		h[i][j] += fii[i+j];
		h[i][j] *= s[i][j];
// Compare integrals
//		printf(" %3i %3i %+15e %+15e %+15e\n",i,j,h_analyt,h[i][j],(h_analyt-h[i][j]));
		sigma+=(h_analyt-h[i][j])*(h_analyt-h[i][j]);
	}
}

sigma=sqrt(sigma)/(double)(nbas*nbas);
printf("\n Sigma = sqrt(sum(h_analyt-h)^2))/(nbas^2) = %+15.8e\n",sigma);

/*
printf("\n   S matrix printout \n");
printmat1(nbas,s);
		
printf("\n   H matrix printout \n");
printmat1(nbas,h);
*/
//END OF FILL


// DIAGONALISE and print

lapack_int info;
lapack_int itype=1;
char jobz='V',uplo='U';

info=LAPACKE_dsygv(LAPACK_ROW_MAJOR,itype,jobz,uplo,nbas,&h[0][0],nbas,&s[0][0],nbas,levels);
/* lapack_int LAPACKE_dsygv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, double* a, lapack_int lda,
                          double* b, lapack_int ldb, double* w ); */

if (info != 0) printf("\n!! error in dsygv call, info = %d\n",info);

int n_e_max=(int)(sqrt(d * 2.0)/w - 1.0) + 1; // analytical expression for number of levels

printf("\n Eigenvalues found:\n                 Numeric              Analytic\n");

dtemp3 = w * sqrt(d * 2.0);
dtemp4 = w * w;

for(i=0;i<n_e_max;++i){
/*-D+dble(i-0.5)*w*SQRT(D*TWO)-(dble(i-0.5)**2)*(w**2)/TWO*/
    dtemp1 = (double)i + 0.5;
    dtemp2 = -d + dtemp1 * dtemp3 - dtemp1 * dtemp1 * dtemp4 / 2.0;
	printf("%3d %20.15f  %20.15f\n",i,levels[i],dtemp2);
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

void morse_polynom(double d, double w, double xstart,double deltax,int nb) 

/******************************************************************************/
{
//double xstart,xend,deltax;     //* basis parameters 
double x[nb],a[nb];              //* basis arrays - small, keep them in stack 
double pot[nb],fi[MAX_POW];		 //* potential values and polynome coefficients 
double vtemp[nb],levels[nb];     //* temporal array and array with energy levels          
double dtemp1,dtemp2,dtemp3,dtemp4;      //* temporal values 
double ifi[MAX_POW];             //* integrals 
// we keep any variable that is connected to array size as lapack_int (mkl_int) - for easy 32/64-bit switch
//int nbas;                      //* amount of basis centers 
//int i,j,k;                     //* loop counters
lapack_int nbas,i,j,k;
lapack_int info=0,istat[3],pow=MAX_POW;   // int SVD stuff, and pow - actual max power of polynomial to be used
double scale,sva[MAX_POW],stat[7];        // singular values (scale,sva) and other SVD-related stuff
double v[MAX_POW][MAX_POW];               // array for right singular vectors
//__float128 qai,qbi,qfi[MAX_POW],qifi[MAX_POW];  //increased precision variables
//char buf[256];


// tune sizing
nbas = (lapack_int) nb;
if (nbas<MAX_POW) pow=nbas; 

// print banner
printf ("\n");
printf (" !======================================================\n");
printf (" !  1D subprogram for Morse potential                   \n");
printf (" !  The potential is interpolated by polynom globally    \n");
printf (" !  Integrals with polynom are analytical               \n");
printf ("\n");



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


// FILL in X and A values 

for (i = 0; i < nbas; ++i)
	x[i] = xstart + deltax * i;

dtemp1 = x[2] - x[1];
a[0]= C * C / (dtemp1 * dtemp1);
for (i = 1; i < (nbas-1); ++i){
	dtemp1 = x[i+1] - x[i-1];
	a[i]= 4.0 * C * C / (dtemp1 * dtemp1);
}
dtemp1 = x[nbas-1] - x[nbas-2];
a[nbas-1]= C * C / (dtemp1 * dtemp1);


//* Print potential and basis parameters at basis centers

printf("\n   N     X value     A value  Potential(x)\n");
for (i = 0; i < nbas; ++i){
	dtemp1 = d * (exp(-2.0 * w * x[i]) - 2.0* exp(- w * x[i] ));
	printf(" %3d  %10.6f  %10.6f  %12.6f\n",i+1,x[i],a[i],dtemp1);
	pot[i]=dtemp1 + d;
}


// BUILD full Vandermonde martix

for(i=0;i<nbas;++i){
	dvm[i][0]=x[i];
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

printf("LAPACK SVD (dgejsv) completed!\n\n Stat array values = ");
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
for(i=0;i<nbas;i++){
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

/* alternative code - we create single transformation matrix V*Sigma(-1)*U(T) in dvm for future use
// we found it numerically less stable
// dvm(T)= V*Sigma(-1)*U(T)
for(i=0;i<nbas;++i){
	for(j=0;j<pow;++j){
		dvm[i][j]=0.0;
		for(k=0;k<pow;++k){
			dvm[i][j] += v[j][k] * u[i][k] / (sva[k]*scale);
		}
	}
}

//transform pot to coefficients
//fi=dvm(T)*pot
for(i=0;i<pow;++i){
	fi[i]=0.0;
	for(j=0;j<nbas;++j)
		fi[i]+= dvm[j][i]*pot[j];
}
// end of alternative code */

printf("Power coefficients found: \n");

printf("\n  N              Fit           Analytical \n");

dtemp2 = - w;
dtemp3 = 1.0;
dtemp4 = 1.0;

for(i=0;i<pow;i++) {
    dtemp1 = 2.0 * d * dtemp2 * (dtemp3 - 1.0) * dtemp4;
    dtemp2 *= - w;
    dtemp3 *= 2.0;
    dtemp4 /= (double) i + 2.0;
	printf("F%2d = %+15.10e %+15.10e\n",i+1,fi[i],dtemp1);
//	fi[i]=dtemp1;  //use this to try analytic Taylor expansion
}
printf("\n");


// END OF SOLUTION

// PRINT FIT
printf("Fit accuracy : \n");
printf("\n  N            x    V(analyt)       V(fit)         Diff\n");

dtemp1=0.0;

for(i=0;i<nbas;i++){
 dtemp2 = poly_d(pow,fi,x[i]);
 dtemp3 = pot[i]-dtemp2;
 dtemp1 += dtemp3*dtemp3;
 printf("%3d  %+10.4e  %+10.4e  %+10.4e  %+10.4e\n",i,x[i],pot[i]-d,dtemp2-d,dtemp3);
}

printf("Sigma (sqrt(sum_square_diff)) = %+10.4e\n\n",dtemp1);


printf("Complete polynom group:\n");
printf("{polynom\n");
for(i=0;i<pow;i++)
	printf(" %3i %+20.15e\n",i+1,fi[i]);
printf("}\n");
//END OF PRINT

//DO CALCULATION

// FILL in S,H 

for(i=0;i<nbas;++i){
	for(j=0;j<nbas;++j){
		dtemp1 = x[i] - x[j];
		dtemp2 = a[i] + a[j];
		dtemp3 = a[i] * a[j];
		dtemp4 = a[i] * x[i] + a[j] * x[j];
		
/*    Compact code for overlap matrix S
 * S(I,J)=sqrt(sqrt(A(I)*A(J))*TWO/(A(I)+A(J)))*exp(-A(I)*A(J)*(X(I)-X(J))**2/(A(I)+A(J))) */
		s[i][j]=sqrt(2.0 * sqrt(dtemp3) / dtemp2) * 
		exp(- dtemp3 * dtemp1 *dtemp1 / dtemp2);

/* analytical expression for H element
		double h_analyt = s[i][j]* (- dtemp3 * (2.0*dtemp3*dtemp1*dtemp1 - dtemp2) / (dtemp2 * dtemp2)
		+ d * exp(-2.0 * w * dtemp4 / dtemp2) 
		* ( exp(w * w / dtemp2) - 2.0 * exp(w * (4.0 * dtemp4 + w) / (4.0*dtemp2))));
*/
/*  Original, compact code for kinetic energy
    G(I,J)=-S(I,J)*A(I)*A(J)*(TWO*A(I)*A(J)*(X(I)-X(J))**2-A(I)-A(J))/(A(I)+A(J))**2 */ 
//        h[i][j] = - dtemp3 * (2.0*dtemp3*dtemp1*dtemp1 - dtemp2) / (dtemp2 * dtemp2); 

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
		double bi  = 2.0 * dtemp2;		// 2* exponent coefficient
		ifi[0] = xij;                     // integral <phi*x  *phi>
		ifi[1] = 1.0 / bi + ifi[0] * xij; // integral <phi*x^2*phi>

		h[i][j] = (
		-d +  ifi[0] * fi[0]						// potential, const and x part (normally can be ignored)
		+     ifi[1] * fi[1]                     	// potential, x^2 part
		+ 2.0*ifi[1] * a[i] * a[j]               	// kinetic energy , x^2 part
		- 2.0*ifi[0] * a[i] * a[j] * (x[i] + x[j])	// kinetic energy , x   part 
		+ 2.0*         a[i] * a[j] *  x[i] * x[j]);	// kinetic energy , const part

		for(k=2;k<pow;++k){
			ifi[k] = ifi[k-1] * xij + ifi[k-2] / bi * (double) k; 	// integral <phi*x^k*phi>
			h[i][j] += fi[k] * ifi[k];								// potential, x^k part
		} 
		
		h[i][j] *= s[i][j];


/* Use this code to have integration over local Taylor expansion (uncomment below and comment above part)
// analytical local derivatives + analytic integration over taylor series
	h[i][j] = dtemp3 * (dtemp2 - 2.0*dtemp3*dtemp1*dtemp1) / (dtemp2 * dtemp2); //analytical kinetic energy
	double xij = dtemp4 / dtemp2; // point of expansion
	fi[0]=d * (exp(-2.0 * w * xij) - 2.0* exp(- w * xij )); // potential at expansion point
	ifi[0]=1.0; // first integral x^0 = sij
	
	h[i][j] += fi[0]*ifi[0];
	
	dtemp3 = -w;
	dtemp4 = 1.0;
	
	fi[1]=ifi[1]=0.0;
	
	for(k=2;k<9;k++){
		fi[k]=ifi[k]=0.0;
		dtemp3 *= -w / (double) k;
		dtemp4 *= 2.0;
		if (k%2==0){
			fi[k]= 2.0 * dtemp3 * d * (dtemp4 * exp(-2.0 * w * xij) - exp(-w * xij));
			ifi[k] = ifi[k-2] / 2.0 / dtemp2 *(double)(k-1);
			h[i][j] +=fi[k]*ifi[k];
		}
	} 

	h[i][j] *= s[i][j];
*/
	
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

int n_e_max=(int)(sqrt(d * 2.0)/w - 1.0) + 1; // analytical expression for number of levels

printf("\n Eigenvalues found:\n                 Numeric              Analytic\n");

dtemp3 = w * sqrt(d * 2.0);
dtemp4 = w * w;

for(i=0;i<n_e_max;++i){
/*-D+dble(i-0.5)*w*SQRT(D*TWO)-(dble(i-0.5)**2)*(w**2)/TWO*/
    dtemp1 = (double)i + 0.5;
    dtemp2 = -d + dtemp1 * dtemp3 - dtemp1 * dtemp1 * dtemp4 / 2.0;
	printf("%3d %20.15f  %20.15f\n",i,levels[i],dtemp2);
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
		vtemp[j]=sqrt(sqrt(2.0 * a[j] / PI)) * exp (- a[j] * dtemp1 *dtemp1);
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


free (h), free(s);
free (dvm), free(u);

return;
}

/******************************************************************************/

static void printmat2d(const lapack_int dim,const double mat[dim][MAX_POW]) {

/******************************************************************************/
/*
 *  Purpose:
 * 		Prints square matrix, ncols columns per page 
*/
#define NCOLS 16 /* how many columns fit in one page */
int j,k;

	printf("    ");
	for( k=0; k<NCOLS; ++k)
		printf("%10d",k);	
	printf("\n");
	
	for(j=0;j<dim;++j){
		printf(" %3d",j);
		for( k=0; k<NCOLS; ++k)
			printf("%+10.2e",mat[j][k]);	
		printf("\n");
	}	
	printf("\n");

return;
#undef NCOLS
}


/* static void printmat2ld(const lapack_int dim,const __float128 mat[dim][npow]) {
const int ncols=NPOW; 
int j,k;
char buf[16];

	printf("    ");
	for(k=0;k<ncols;++k)
		printf("%10d",k);	
	printf("\n");
	
	for(j=0;j<dim;++j){
		printf(" %3d",j);
		for(k=0;k<ncols;++k){
		    quadmath_snprintf(buf,sizeof buf,"%+10.2QE",mat[j][k]);
			printf("%10s",buf);	
		}
		printf("\n");
	}	
	printf("\n");

return;
} */

/******************************************************************************/

static inline double poly_d(const lapack_int pow, const double fi[MAX_POW], double x)

/******************************************************************************/
/*
 * Purpose:
 * 	Returns value of polynom max power pow, coefficients in fi, at point x
*/
{
double temp1,temp2;
lapack_int l,m;

temp1 = fi[0] * x;
for(l=1;l<pow;l++){
	temp2 = fi[l] * x;
	for(m=0;m<l;m++)	
		temp2 *= x;
	temp1 += temp2;
}
	
return temp1;	
}


/******************************************************************************/

static void dvand (const int n,const double alpha[MAX_POW], double x[MAX_POW] )

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

