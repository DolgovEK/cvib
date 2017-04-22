// print_and_punch.c -auxillary routines for printing and punching
//
// Version 0.2, February 2017
//
#include <stdio.h>
#include <stdlib.h>

/******************************************************************************/

void print_banner()

/******************************************************************************/
{
printf ("\n");
printf (" !======================================================!\n");
printf (" !           Welcome to program CVIB                    !\n");
printf (" !             Program version 0.2                      !\n");
printf (" !======================================================!\n");

}

/******************************************************************************/

void printmat1(const int dim,const double mat[][dim]) 

/******************************************************************************/
{
const int ncols=16; /* how many columns fit in one page */
int i,j,k;
int nblocks,nlast;

nblocks = dim / ncols; /* how many blocks of 15 columns to print */
nlast   = dim % ncols; /* how many columns remains in last block*/

for(i=0;i<nblocks;++i){
	
	printf("    ");
	for(k=0;k<ncols;++k)
		printf("%10d",k+i*ncols);	
	printf("\n");
	
	for(j=0;j<dim;++j){
		printf(" %3d",j);
		for(k=0;k<ncols;++k)
			printf("%10f",mat[j][k+i*ncols]);	
		printf("\n");
	}	
	printf("\n");
}		

if (nlast>0) {
	printf("    ");
	for(k=0;k<nlast;++k)
		printf("%10d",k+nblocks*ncols);	
	
	printf("\n");
	
	for(j=0;j<dim;++j){
		printf(" %3d",j);
		for(k=0;k<nlast;++k)
			printf("%10f",mat[j][k+nblocks*ncols]);	
		
		printf("\n");
	} 	
}
return;
}

/******************************************************************************/

void printmat2(const int dim,const double mat[][dim]) 

/******************************************************************************/
{
const int ncols=16; /* how many columns fit in one page */
int i,j,k;
int nblocks,nlast;

nblocks = dim / ncols; /* how many blocks of 15 columns to print */
nlast   = dim % ncols; /* how many columns remains in last block*/

for(i=0;i<nblocks;++i){
	
	printf("    ");
	for(k=0;k<ncols;++k)
		printf("%10d",k+i*ncols);	
	printf("\n");
	
	for(j=0;j<dim;++j){
		printf(" %3d",j);
		for(k=0;k<ncols;++k)
			printf("%+10.2e",mat[j][k+i*ncols]);	
		printf("\n");
	}	
	printf("\n");
}		

if (nlast>0) {
	printf("    ");
	for(k=0;k<nlast;++k)
		printf("%10d",k+nblocks*ncols);	
	
	printf("\n");
	
	for(j=0;j<dim;++j){
		printf(" %3d",j);
		for(k=0;k<nlast;++k)
			printf("%+10.2e",mat[j][k+nblocks*ncols]);	
		
		printf("\n");
	} 	
}
return;
}

