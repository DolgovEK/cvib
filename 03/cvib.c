/*
 * cvib.c - main file for cvib program
 *
 * Version 0.3, July 2017
 * 
 * Copyright 2017 Evgeny Dolgov <dolgov.evgeny@gmail.com>
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

// Version history
//
// Version 0.1 do calculation for 1D harmonic oscillator problem
// it reads xstart and deltax as command line arguments
// invocation: ~#./cvib.x task_type xstart deltax nbas > cvib.out
//
// Version 0.2 adds global polynom fit via SVD, 
// local polynom fit using Bjorck-Pereyra,
// 1D general solver with polynom and pointwise potentials,
// 1D morse oscillator,
// input via cvib.inp,
// GAMESS(US) output parser via friend=gamess
// invocation: ~#./cvib.x cvib.inp > cvib.out
//
// Version 0.3 adds vibrational CAS-CI,
// and vibrational SCF followed by MP2, all based on QFF potential
// invocation: ~#./cvib.x cvib.inp > cvib.out
// GAMESS(US) output file should be called just "output"
//


/* Function declarations, print_and_punch.c*/
void printmat1(const int dim,const double mat[][dim]);
void print_banner();

/* Function declarations, harmonic_and_morse.c*/
void harmonic(double w,double xstart,double xend,int nb);
void morse(double d,double w,double xstart,double deltax,int nb);

/* Function declarations, morse_polynom_points.c*/
void morse_polynom(double d,double w,double xstart,double deltax,int nb);
void morse_points(double d,double w,double xstart,double deltax,int nb);

/* Function declarations, 1d_polynom_points.c*/
void one_d_points(const char *filename, const double xstart, const double deltax, const int nb); 
void one_d_polynom(const char *filename, const double xstart, const double deltax, const int nb, const int npow); 
void one_d(const char *filename, const double xstart, const double deltax, const int nb, const int npow); 

/* Function declaration, parser.c*/
void parse_control(const char *filename,int *task_type,int *pot_type,int *friend);
void parse_harmonic(const char *filename,double *w,int *nbas,double *xstart,double *xend);
void parse_morse(const char *filename,double *d,double *w,int *nbas,double *xstart,double *xend);
void parse_1d(const char *filename,int *npow,int *nbas,double *xstart,double *xend); 
void parse_size(const char *filename, int *nmodes, int *nbas);

/* Function declaration, gamess.c*/
int gamess_get_nmodes_nbas(int *nmodes,int *nbas);

/* Function declaration, vscf_mp2.c*/
int vscf_qff(const char *filename, const int nmodes,const int npairs,const int nbas);

/* Function declaration, ci.c*/
int ci_qff(const char *filename, const int nmodes,const int npairs,const int nbas);

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// main function of the program
//
//
/******************************************************************************/

int main (int argc, char *argv[]) 

/******************************************************************************/
{
int task_type=0; /* =0 -> harmonic (default)
                    =1 -> morse 
                    =2 -> 1D      
                    =3 -> VCI */
int pot_type=0;  /* =0 -> analytic (default)
                    =1 -> polynom 
                    =2 -> points */
int friend=0;    /* =0 -> none (default)
                    =1 -> gamess  */
                    
double xstart=-15.5,xend=15.5,deltax=0.2;     	// basis parameters 
double d=12.0,w=0.5;						// potential parameters for harmonic and morse 1D tasks 
int nbas=32,npow=16,nmodes=3;

/* Print banner */

print_banner();

if (argc>0)

parse_control(argv[1],&task_type,&pot_type,&friend);

//	printf("friend = %2d\n",friend);


if (friend==1)
{
//	printf("in friend\n");
	gamess_get_nmodes_nbas(&nmodes,&nbas);

	parse_size(argv[1],&nmodes,&nbas); 

	switch (task_type) 
	{
	case 0: 
	case 1: 
	case 2:
//TODO: restore simple 1D calculation for GAMESS(US)	 
//	parse_1d(argv[1],&npow,&nbas,&xstart,&xend);
//	break;  
	case 3: // VCI
		ci_qff(argv[1],nmodes,nmodes*(nmodes-1)/2,nbas);	
	break;
	case 4: // VSCF
		vscf_qff(argv[1],nmodes,nmodes*(nmodes-1)/2,nbas);	
	break;
	}
	
} else {  // if no friend

parse_size(argv[1],&nmodes,&nbas); 

switch (task_type) 
{
case 0: 
	parse_harmonic(argv[1],&w,&nbas,&xstart,&xend);
	break;  
case 1: 
	parse_morse(argv[1],&d,&w,&nbas,&xstart,&xend);
	break;  
case 2: 
	parse_1d(argv[1],&npow,&nbas,&xstart,&xend);
	break;  
}

deltax = (xend - xstart)/ ((double) (nbas - 1));

printf ("\n Subprogram choice:\n");	
printf (" task   = ");

switch (task_type) 
{
case 0: 
	printf("harmonic\n w      = %10g\n",w);
	break;
case 1: 
	printf("morse\n d      = %10g\n w      = %10g\n",d,w);
	break;
case 2: 
	printf("1D\n");
	break;
case 3: 
	printf("VCI\n");
	break;
case 4: 
	printf("VSCF\n");
	break;
}

printf (" potential   = ");

switch (pot_type) 
{
case 0: 
	printf("analytic\n");
	break;
case 1: 
	printf("polynom\n");
	break;
case 2: 
	printf("points\n");
	break;
}
 
if (task_type < 3) printf ("\n Starting parameters: xstart = %10g\n xend   = %10g\n deltax = %10g\n nbas   = %10d\n",xstart,xend,deltax,nbas);	


switch (task_type) 
{
case 0: //harmonic 
    harmonic(w,xstart,deltax,nbas); // call to this function is irrelevant to pot_type value
	break;
case 1: //morse
	switch (pot_type) {
	case 0:
		morse(d,w,xstart,deltax,nbas);           // analytic evaluation
		break;
	case 1:
		morse_polynom(d,w,xstart,deltax,nbas);   // global polynomial interpolation
		break;
	case 2:
		morse_points(d,w,xstart,deltax,nbas);    // local interpolation
		break; 
	}
	break;
case 2: //1d
	switch (pot_type) {
	case 0:
		one_d(argv[1],xstart,deltax,nbas,npow); // analytic integrals, polynom is red from file
		break;
	case 1:
		one_d_polynom(argv[1],xstart,deltax,nbas,npow);// analytic integrals, global polynomial interpolation 
		break;
	case 2:
		one_d_points(argv[1],xstart,deltax,nbas);		// local interpolation 
		break;
	}
case 3: //VCI 
    ci_qff(argv[1],nmodes,nmodes*(nmodes-1)/2,nbas); // call to this function is irrelevant to pot_type value
	break;
case 4: // VSCF
	vscf_qff(argv[1],nmodes,nmodes*(nmodes-1)/2,nbas);	
	break;
}

} // end of if (friend)

printf("\n All done!\n\n Arevoir!\n\n");

return 0;
}

