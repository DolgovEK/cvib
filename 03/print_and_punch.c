/* print_and_punch.c -auxillary routines for printing and punching
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

/******************************************************************************/

void print_banner()

/******************************************************************************/
{
printf ("\n");
printf (" !======================================================!\n");
printf (" !           Welcome to program CVIB                    !\n");
printf (" !             Program version 0.3                      !\n");
printf (" !======================================================!\n");

}

/******************************************************************************/

void printmat1(const int dim, double mat[][dim]) 

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

void printmat2(const int dim, double mat[][dim]) 

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

/******************************************************************************/

void printmat3(const int dim, double *mat) 

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
		printf("%13d",k+i*ncols);	
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
		printf("%13d",k+nblocks*ncols);	
	
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
