// harmonic_and_morse.c - 1D computations on model potentials
//
// Version 0.2, February 2017
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl.h>

// double const w=1.0;     /* frequency of the oscillator */
// double const d=12.0;    /* dissociation limit for morse */
double const c=0.7;     /* basis set broadness */
double const pi=3.14159265358979323846; // just pi

/* Function declaration, print_and_punch.c*/
void printmat1(const int dim,const double mat[][dim]); 
//
//=======================================================================
/******************************************************************************/

void harmonic(double w,double xstart,double deltax,int nb) 

/******************************************************************************/
{

double x[nb],a[nb];              /* basis arrays - small, keep them in stack */
double vtemp[nb],levels[nb];     /* temporal array and array with energy levels */         
double dtemp1,dtemp2,dtemp3;     /* temporal values */
double ai,bi,ifi1,ifi2;          /* integrals */
//int nbas;                      /* amount of basis centers */
//int i,j,k;                     /* loop counters*/
// we keep any variable that is connected to array size as lapack_int (mkl_int) - for easy 32/64-bit switch
lapack_int nbas,i,j,k;

nbas = (lapack_int) nb;

// print banner
printf ("\n");
printf (" !======================================================\n");
printf (" !  Simple 1D subprogram for harmonic potential         \n");
printf (" !  Analytical integration                              \n");
printf ("\n");


/* ALLOCATE memory for large S (overlap) and H (hamiltonian) matrices, in C99 style VLAs, to save space in stack */

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
/*END OF ALLOCATE */


/* FILL in X and A values */

for (i = 0; i < nbas; ++i)
	x[i] = xstart + deltax * i;

dtemp1 = x[2] - x[1];
a[0]= c * c / (dtemp1 * dtemp1);
for (i = 1; i < (nbas-1); ++i){
	dtemp1 = x[i+1] - x[i-1];
	a[i]= 4.0 * c * c / (dtemp1 * dtemp1);
}
dtemp1 = x[nbas-1] - x[nbas-2];
a[nbas-1]= c * c / (dtemp1 * dtemp1);


/* Print potential and basis parameters at basis centers*/

printf("\n   N     X value     A value  Potential(x)\n");
for (i = 0; i < nbas; ++i){
	dtemp1 = w * x[i]*x[i];
	printf(" %3d  %10.6f  %10.6f  %12.6f\n",i+1,x[i],a[i],dtemp1);
}

/* Fill in S */

for(i=0;i<nbas;++i){
	for(j=0;j<nbas;++j){
		dtemp1 = x[i] - x[j];
		dtemp2 = a[i] + a[j];
		dtemp3 = a[i] * a[j];
		
/*    Compact code for overlap matrix S
 * S(I,J)=sqrt(sqrt(A(I)*A(J))*TWO/(A(I)+A(J)))*exp(-A(I)*A(J)*(X(I)-X(J))**2/(A(I)+A(J))) */
		s[i][j]=sqrt(2.0 * sqrt(dtemp3) / dtemp2) * 
		exp(- dtemp3 * dtemp1 *dtemp1 / dtemp2);


/*  Original, compact code for kinetic energy
    G(I,J)=-S(I,J)*A(I)*A(J)*(TWO*A(I)*A(J)*(X(I)-X(J))**2-A(I)-A(J))/(A(I)+A(J))**2 
        h[i][j]=-s[i][j] * dtemp3 * (2.0*dtemp3*dtemp1*dtemp1 - 
        a[i]-a[j]) / (dtemp2 * dtemp2); */

/*  special code for kinetic energy - useful for polynomial expansions - and symmetric
    -(dphi/dx)^2 * 0.5= (2 AI*AJ *x^2 - 2AI*AJ *x *(XI+XJ) + 2(AI*AJ*XI*XJ))* phi^2 
    AI=A(i)*X(i)+A(j)*X(j)
    BI=A(i)+A(j)
    IFI(1)=AI/BI 
    IFI(2)=ONE/(TWO*BI)+IFI(1)*AI/BI - 
*/

		ai = a[i] * x[i] + a[j] * x[j];
		bi = dtemp2;
		ifi1 = ai/bi;                             /*integral <phi*x  *phi>*/
		ifi2 = 1.0 / (2.0 * bi) + ifi1 * ai / bi; /*integral <phi*x^2*phi>*/

		h[i][j] = s[i][j]*
		(      ifi2 * w                           /* harmonic potential */
		+ 2.0* ifi2 * a[i] * a[j]                 /* kinetic energy , x^2 part*/
		- 2.0* ifi1 * a[i] * a[j] * (x[i] + x[j]) /* kinetic energy , x   part */
		+ 2.0*        a[i] * a[j] *  x[i] * x[j]);/* kinetic energy , const part */ 

	}
}

//printf("\n   S matrix printout \n");
//printmat1(nbas,s);
		
//printf("\n   H matrix printout \n");
//printmat1(nbas,h);
/*END OF FILL*/


/* DIAGONALISE and print*/

lapack_int info;
lapack_int itype=1;
char jobz='V',uplo='U';

info=LAPACKE_dsygv(LAPACK_ROW_MAJOR,itype,jobz,uplo,nbas,&h[0][0],nbas,&s[0][0],nbas,levels);
/* lapack_int LAPACKE_dsygv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, double* a, lapack_int lda,
                          double* b, lapack_int ldb, double* w ); */

if (info != 0) printf("\n!! error in dsygv call, info = %d\n",info);

printf("\n Eigenvalues found:\n                 Numeric              Analytic\n");

w=sqrt(2.0*w);

for(i=0;i<nbas;++i)
	printf("%3d %20.15f  %20.15f\n",i,levels[i],w*(double)(i+0.5));

printf("\n");

printf("\n   Eigenvector printout \n");
printmat1(nbas,h);

/*END OF DIAGONALISE*/

// CREATE the matrix of wavefunction values
// we will save on one matrix (and pay more computation cost), and compute phi(xi) on the fly, in vtemp

for(i=0;i<nbas;++i){
	for(j=0;j<nbas;++j){
		dtemp1=x[i]-x[j];
		vtemp[j]=sqrt(sqrt(2.0 * a[j] / pi)) * exp (- a[j] * dtemp1 *dtemp1);
	}
	for(k=0;k<nbas;++k){
		s[i][k]=0.0;
		for(j=0;j<nbas;++j)
			s[i][k]	+= h[j][k]*vtemp[j];
	}
}		

printf("\n   Wavefunction values at basis centers \n");
printmat1(nbas,h);

/* END OF CREATE */
		
free(h), free(s);

return;
}
/******************************************************************************/

void morse(double d,double w,double xstart,double deltax,int nb) {

/******************************************************************************/
double x[nb],a[nb];              /* basis arrays - small, keep them in stack */
double vtemp[nb],levels[nb];     /* temporal array and array with energy levels */         
double dtemp1,dtemp2,dtemp3,dtemp4;      /* temporal values */
double ai,bi,ifi1,ifi2;          /* integrals */
//int nbas;                      /* amount of basis centers */
//int i,j,k;                     /* loop counters*/
// we keep any variable that is connected to array size as lapack_int (mkl_int) - for easy 32/64-bit switch
lapack_int nbas,i,j,k;

nbas = (lapack_int) nb;

// print banner
printf ("\n");
printf (" !======================================================\n");
printf (" !  Simple 1D subprogram for Morse potential            \n");
printf (" !  All integrals are calculated analytically           \n");
printf ("\n");


/* ALLOCATE memory for large S (overlap) and H (hamiltonian) matrices, in C99 style VLAs, to save space in stack */

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
/*END OF ALLOCATE */


/* FILL in X and A values */

for (i = 0; i < nbas; ++i)
	x[i] = xstart + deltax * i;

dtemp1 = x[2] - x[1];
a[0]= c * c / (dtemp1 * dtemp1);
for (i = 1; i < (nbas-1); ++i){
	dtemp1 = x[i+1] - x[i-1];
	a[i]= 4.0 * c * c / (dtemp1 * dtemp1);
}
dtemp1 = x[nbas-1] - x[nbas-2];
a[nbas-1]= c * c / (dtemp1 * dtemp1);


/* Print potential and basis parameters at basis centers*/

printf("\n   N     X value     A value  Potential(x)\n");
for (i = 0; i < nbas; ++i){
	dtemp1 = d * (exp(-2.0 * w * x[i]) - 2.0* exp(- w * x[i] ));
	printf(" %3d  %10.6f  %10.6f  %12.6f\n",i+1,x[i],a[i],dtemp1);
}

/* Fill in S */

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


/*  Original, compact code for kinetic energy
    G(I,J)=-S(I,J)*A(I)*A(J)*(TWO*A(I)*A(J)*(X(I)-X(J))**2-A(I)-A(J))/(A(I)+A(J))**2 */
        h[i][j]=-s[i][j] * dtemp3 * (2.0*dtemp3*dtemp1*dtemp1 - 
        dtemp2) / (dtemp2 * dtemp2); 

/*  original code for potential energy   
 *  rtemp=EXP(-TWO*w*(AI*XI+AJ*XJ)/(AI+AJ))
    fij_analyt=Sij*D*rtemp*(EXP(w**2/(AI+AJ))-TWO*EXP(w*(4.0D0*(AI*XI+AJ*XJ)+w)/(4.0D0*(AI+AJ)))) */
//		dtemp4=
		h[i][j]+=s[i][j]* d * exp(-2.0 * w * dtemp4 / dtemp2) 
		* ( exp(w * w / dtemp2) - 2.0 * exp(w * (4.0 * dtemp4 + w) / (4.0*dtemp2)));

	}
}

//printf("\n   S matrix printout \n");
//printmat1(nbas,s);
		
//printf("\n   H matrix printout \n");
//printmat1(nbas,h);
/*END OF FILL*/


/* DIAGONALISE and print*/

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
    dtemp1 = (double)(i+0.5);
    dtemp2 = -d + dtemp1 * dtemp3 - dtemp1 * dtemp1 * dtemp4 / 2.0;
	printf("%3d %20.15f  %20.15f\n",i,levels[i],dtemp2);
}

printf("\n");

printf("\n   Eigenvector printout \n");
printmat1(nbas,h);

/*END OF DIAGONALISE*/

// CREATE the matrix of wavefunction values
// we will save on one matrix (and pay more computation cost), and compute phi(xi) on the fly, in vtemp

for(i=0;i<nbas;++i){
	for(j=0;j<nbas;++j){
		dtemp1=x[i]-x[j];
		vtemp[j]=sqrt(sqrt(2.0 * a[j] / pi)) * exp (- a[j] * dtemp1 *dtemp1);
	}
	for(k=0;k<nbas;++k){
		s[i][k]=0.0;
		for(j=0;j<nbas;++j)
			s[i][k]	+= h[j][k]*vtemp[j];
	}
}		

printf("\n   Wavefunction values at basis centers \n");
printmat1(nbas,h);

/* END OF CREATE */
		
free(h), free(s);

return;
}





