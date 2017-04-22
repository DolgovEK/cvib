/* gamess.c - subroutines to parse GAMESS(US) output files
 * and perform calculations using these date
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mkl.h>


#define AMU 1822.9					// mass from a.u.(me==1) to a.m.u (1/12 mc == 1)
#define C 0.7						// basis broadness - hardcoded in case of GAMSS(US)
#define PI 3.14159265358979323846 	// just pi
#define CM1 219474.63				// hartree to cm-1 conversion factor

// Forward function declaration
int gamess_1d(const int task_type,const int nmodes,const int npairs,const int nbas);
int gamess_get_nmodes_nbas(int *nmodes,int *nbas);
int parse_gamess_out(const int nmodes, const int npairs, double harmonic[nmodes], double xrange[nmodes], double deltax[nmodes], double fdiag[nmodes][4], double foffdiag[npairs][8]);
void printmat3(const int dim,const double *mat); 


//* Function declaration, print_and_punch.c
void printmat1(const int dim,const double mat[][dim]);

/******************************************************************************/

int gamess_1d(const int task_type,const int nmodes,const int npairs,const int nbas) // task_type currently is not used, but reserved for future use

/******************************************************************************/
{
double dtemp1,levels[nbas],ifi[4];								//temporal values
double a,xrange[nmodes],deltax[nmodes],x[nbas];					//basis parameters
double harmonic[nmodes],fdiag[nmodes][4],foffdiag[npairs][8];	//parameters of potential
// Allocated on stack:
//
// fdiag (manually aligned to 32 byte):
// 0 => Hii
// 1 => Tiii
// 2 => Uiiii
//
// foffdiag (manually aligned to 64 byte):
// 0 => Tiij
// 1 => Tjji
// 2 => Uiiij
// 3 => Ujjji
// 4 => Uiijj
int mode;	// mode counter
int i,j,k,l;// loop counters
// TODO initialize fdiag, foffdiag

// print banner
printf ("\n");
printf (" !=====================================================================================\n");
printf (" !  1D subprogram for polynomial potential                   \n");
printf (" !  The potential is analytic polynom (read as QFF field from GAMESS(US) output file)   \n");
printf (" !  Integrals over polynom are analytical               \n");
printf ("\n");

printf(" Reading output file...\n");

parse_gamess_out(nmodes,npairs,harmonic,xrange,deltax,fdiag,foffdiag);

printf("\n Read results:\n\n");
printf(" nmodes = %3d\n",nmodes);
printf(" npairs = %3d\n",npairs);
printf(" nbas   = %3d\n",nbas);

printf("\n Basis centers:\n");
for (i=0;i<nmodes;i++)
	printf(" Mode = %2d, xrange = %12.6f, deltax = %12.6f\n",i,xrange[i],deltax[i]);


printf("\n Harmonic frequencies (cm-1):\n");
for (i=0;i<nmodes;i++)
	printf(" %2d  %10.1f\n",i,harmonic[i]);

printf("\n Diagonal anharmonic coefficients:\n");
for (i=0;i<nmodes;i++)
{
	printf(" Mode = %2d\n",i);
	printf(" Hii    = %+15.8e\n",fdiag[i][0]);
	printf(" Tiii   = %+15.8e\n",fdiag[i][1]);
	printf(" Uiiii  = %+15.8e\n",fdiag[i][2]);
}

printf("\n Pair-wise anharmonic coefficients:\n");
k=0;
for (i=0;i<(nmodes-1);i++) {
for (j=(i+1);j<nmodes;j++) 
{
	printf(" Mode pair = %2d , %2d\n",i+1,j+1);
	printf(" Tiij  = %+15.8e\n",foffdiag[k][0]);
	printf(" Tjji  = %+15.8e\n",foffdiag[k][1]);
	printf(" Uiiij = %+15.8e\n",foffdiag[k][2]);
	printf(" Ujjji = %+15.8e\n",foffdiag[k][3]);
	printf(" Uiijj = %+15.8e\n",foffdiag[k][4]);
	k++;
}}

// convert to polynom coefficients

for (i=0;i<nmodes;i++)
{
	fdiag[i][0]/=AMU*2.0;
	fdiag[i][1]/=AMU*6.0*sqrt(AMU);
	fdiag[i][2]/=AMU*AMU*24.0;
}

for (i=0;i<npairs;i++)
{
	foffdiag[i][0]/=AMU*6.0*sqrt(AMU);
	foffdiag[i][1]/=AMU*6.0*sqrt(AMU);
	foffdiag[i][2]/=AMU*AMU*24.0;
	foffdiag[i][3]/=AMU*AMU*24.0;
	foffdiag[i][4]/=AMU*AMU*24.0;
}

// ALLOCATE memory for large S (overlap) and H (hamiltonian) matrices, in C99 style VLAs - nbas is known only at runtime 

double (*s)[nbas] = malloc( sizeof(double[nbas][nbas]));
if (s==NULL){
	printf(" Error while allocating memory for s \n");
	exit (1);
}

// todo - add dimension for modes (to keep all WFs for future use)
double (*h)[nmodes][nbas] = malloc( sizeof(double[nmodes][nbas][nbas]));
if (h==NULL){
	printf(" Error while allocating memory for h \n");
	exit (1);
}


// do calculations

for(mode=0;mode<nmodes;mode++)
{
	printf("\n=========================================\n");
	printf(" Starting calculation for mode = %2d\n\n",mode);

// evaluate deltax, xstart
/*
 * QR1=TWO*(AMP+HALF)*PLANC*FREQ(IM)*UNITCONV
            QR2=EIG(I)/AMU
            QRANGE=SQRT(QR1/QR2)
            QMIN   = -QRANGE
            DQ(IM) = (TWO*QRANGE)/(NGRID-1)
 */

// fill in A,X
	for (i = 0; i < nbas; ++i)
		x[i] = -xrange[mode] + deltax[mode] * (double) i;

	a = C * C /(deltax[mode] * deltax[mode]);

// print x and potential
	printf("\n Basis points :\n");
	for (i = 0; i < nbas; ++i)
		printf(" %3d    %12.6f\n",i,x[i]);

// fill H,S

	double a2 = a * a; // a[i] * a[j];
	double bi = 4.0 * a;		// 2* exponent coefficient

	for(i=0;i<nbas;i++){
	for(j=0;j<nbas;j++){

		dtemp1 = x[i] - x[j];
		
//    Compact code for overlap matrix S for equidistant basis
		s[i][j] = exp(- a * dtemp1 * dtemp1 * 0.5);

/*  special code for kinetic energy - useful for polynomial expansions - and symmetric
 *    -(dphi/dx)^2 * 0.5= (2 AI*AJ *x^2 - 2AI*AJ *x *(XI+XJ) + 2(AI*AJ*XI*XJ))* phi^2 */
 
		double xij = 0.5 * (x[i] + x[j]); 		// point of expansion
//        ifi[0] = 1.0;
		ifi[0] = xij;                     			// integral <phi*x  *phi>
		ifi[1] = ifi[0] * xij + 1.0    / bi; 		// integral <phi*x^2*phi>
		ifi[2] = ifi[1] * xij + ifi[0] / bi * 2.0; 	// integral <phi*x^3*phi>
		ifi[3] = ifi[2] * xij + ifi[1] / bi * 3.0; 	// integral <phi*x^4*phi>

		h[mode][i][j] = s[i][j] * (
		      ifi[1] * fdiag[mode][0]            	// potential, x^2 part
		+     ifi[2] * fdiag[mode][1]            	// potential, x^3 part
		+     ifi[3] * fdiag[mode][2]            	// potential, x^4 part
		+ 2.0*ifi[1] * a2               	// kinetic energy , x^2 part
		- 2.0*ifi[0] * a2 * (x[i] + x[j])	// kinetic energy , x   part 
		+ 2.0        * a2 *  x[i] * x[j]);	// kinetic energy , const part
		
//		h[mode][i][j] *= s[i][j];

	}
	}

// print matrices

/*
printf("\n S matrix for mode = %2d \n",mode);
printmat1(nbas,s);

printf("\n H matrix for mode = %2d \n",mode);
printmat3(nbas,&h[mode][0][0]);
*/

// diagonalise

lapack_int info;
lapack_int itype=1;
char jobz='V',uplo='U';

info=LAPACKE_dsygv(LAPACK_ROW_MAJOR,itype,jobz,uplo,nbas,&h[mode][0][0],nbas,&s[0][0],nbas,levels);
/* lapack_int LAPACKE_dsygv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, double* a, lapack_int lda,
                          double* b, lapack_int ldb, double* w ); */

if (info != 0) printf("\n!! error in dsygv call, info = %d\n",info);

// print results	

printf("\n Eigenvalues (hartree) and transition energies (cm-1) for mode = %2d\n\n",mode);

for(i=0;i<nbas;++i){
	printf("%4d %20.15f %20.5f\n",i,levels[i],(levels[i]-levels[0])*CM1);
}

printf("\n");

printf("\n Eigenvectors for mode = %2d \n",mode);
printmat3(nbas,&h[mode][0][0]);

} // end of for(mode=0;mode<nmodes;mode++)


free(h); free(s);
return 0;	
}

/******************************************************************************/

int gamess_get_nmodes_nbas(int *nmodes,int *nbas)

/******************************************************************************/

{
char line[82],*pos;
FILE *fp;
	
		printf("in gamess_get_nmodes_nbas\n");
if ((fp = fopen("output","r"))!=NULL)
	{
	while (fgets(line,80,fp)!=NULL)					// list all file lines till the end
		if ((pos=strstr(line,"NZVAR ="))!=NULL)		// find NZVAR group (can be many)
			{*nmodes=atoi(pos+8);}
		else if ((pos=strstr(line,"NGRID="))!=NULL)// find NGRID group
			{*nbas=atoi(pos+8);}

	} else {
		printf("Can not open \'output\' GAMESS output file, exiting \n");
		fclose(fp);
		exit(1);
	}

fclose(fp);
			
return 0;
}

/******************************************************************************/

int parse_gamess_out(const int nmodes, const int npairs, double harmonic[nmodes], double xrange[nmodes], double deltax[nmodes],double fdiag[nmodes][4], double foffdiag[npairs][8])

/******************************************************************************/
{
char line[82],*pos;
FILE *fp;
int freq_line_count,freq_line_max;
int harmonic_count=0,i,j;

freq_line_count=freq_line_max=(nmodes+6)/5+1; //number of harmonic freqency lines
	
if ((fp = fopen("output","r"))!=NULL)
	{
	while (fgets(line,80,fp)!=NULL)				// list all file lines till the end
	{
		if ((pos=strstr(line,"FREQUENCY:"))!=NULL)	// find FREQUENCY: line (should be many)
		{
			if(freq_line_count==freq_line_max)
			{
//				printf("Found first  HARMONIC line, line count= %d\n",freq_line_count);
				freq_line_count--; // just skip first line
			} 
			else if (freq_line_count==(freq_line_max-1)) // second line
			{
//				printf("Found second HARMONIC line, line count= %d\n",freq_line_count);
				for (i=1;i<5;i++)
					if (harmonic_count<nmodes) 
					{
						harmonic[harmonic_count]=atof(&line[18+12*i]);
						harmonic_count++;
					}

				freq_line_count--; // next line
			} 
			else if (freq_line_count>0)		// third line and till the end
			{
				printf("Found next HARMONIC line, line count= %d\n",freq_line_count);
				for (i=0;i<5;i++)
					if (harmonic_count<nmodes) 
					{
						harmonic[harmonic_count]=atof(&line[18+12*i]);
						harmonic_count++;
					}
				freq_line_count--; // next line
			}
		}
		
		if (strstr(line,"MODE # 1   QRANGE=")!=NULL)	// find QFF diagonal frequency group
		{
			xrange[0]=atof(&line[20]);
			deltax[0]=atof(&line[43]);
			for(i=1;i<nmodes;i++)
			{
				if(fgets(line,80,fp)!=NULL)
				{
					xrange[i]=atof(&line[20]);
					deltax[i]=atof(&line[43]);
				}
			}
		}

		if (strstr(line,"QFF>  MODE=")!=NULL)	// find QFF diagonal frequency group
		{
			for(i=0;i<nmodes;i++)
			{
				if(fgets(line,80,fp)!=NULL)
					fdiag[i][0]=atof(&line[15]);			// Hii
				
				if(fgets(line,80,fp)!=NULL)
					fdiag[i][1]=atof(&line[15]);			// Tiii

				if(fgets(line,80,fp)!=NULL)
					fdiag[i][2]=atof(&line[15]);			// Uiiii
				
				for(j=0;j<15;j++) 
					if(fgets(line,80,fp)!=NULL);// skip to the next mode
			}
		}

		if (strstr(line,"QFF>  MODE (I,J)")!=NULL)	// find QFF frequency group
		{
			for(i=0;i<npairs;i++)
			{
				if(fgets(line,80,fp)!=NULL)
					foffdiag[i][0]=atof(&line[15]);			// Tiij
				
				if(fgets(line,80,fp)!=NULL)
					foffdiag[i][1]=atof(&line[15]);			// Tjji
				
				if(fgets(line,80,fp)!=NULL)
					foffdiag[i][2]=atof(&line[15]);			// Uiiij
				
				if(fgets(line,80,fp)!=NULL)
					foffdiag[i][3]=atof(&line[15]);			// Ujjji
				
				if(fgets(line,80,fp)!=NULL)
					foffdiag[i][4]=atof(&line[15]);			// Uiijj
				
				if(fgets(line,80,fp)!=NULL)
					if(fgets(line,80,fp)!=NULL);			// just skip 2 lines to the next mode pair
			}
		}

		
	}

	} else {
		printf("Can not open \'output\' GAMESS output file, exiting \n");
		fclose(fp);
		exit(1);
	}

fclose(fp);
	
return 0;
}

/******************************************************************************/

void printmat3(const int dim,const double *mat) 

/******************************************************************************/
{
const int ncols=16; /* how many columns fit in one page */
int i,j,k;
int nblocks,nlast;

nblocks = dim / ncols; /* how many blocks of 15 columns to print */
nlast   = dim % ncols; /* how many columns remains in last block*/

for(i=0;i<nblocks;++i)
{
	printf("    ");
	for(k=0;k<ncols;++k)
		printf("%10d",k+i*ncols);	
	printf("\n");
	
	for(j=0;j<dim;++j)
	{
		printf(" %3d",j);
		for(k=0;k<ncols;++k)
			printf(" %+12.5e",mat[j*dim + k + i * ncols]);	
		printf("\n");
	}	
	printf("\n");
}		

if (nlast>0) 
{
	printf("    ");
	for(k=0;k<nlast;++k)
		printf("%10d",k+nblocks*ncols);	
	
	printf("\n");
	
	for(j=0;j<dim;++j)
	{
		printf(" %3d",j);
		for(k=0;k<nlast;++k)
			printf(" %+12.5e",mat[j*dim + k + nblocks * ncols]);	
		
		printf("\n");
	} 	
}
return;
}
