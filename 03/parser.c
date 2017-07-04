/* parser.c - functions to parse cvib.inp
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

/******************************************************************************/

void parse_control(const char *filename,int *task_type,int *pot_type,int *friend) 

/******************************************************************************/
/*
 * Purpose
 *  to parse {control} group
 *  in the form {control: task harmonic, potential polynom, friend none}
 * 
 * */
 {
char line[82];
FILE *fp;
	
if ((fp = fopen(filename,"r"))!=NULL){
	
	while (fgets(line,80,fp)!=NULL){ 				// list all file lines till the end
	if (strstr(line,"{control")!=NULL){				// find start of {control group (can be many)
		do {										// list all file lines till the end (inside the group)

			if (strstr(line,"task")!=NULL){			// parse task
				if (strstr(line,"harmonic")!=NULL)
					*task_type=0;
				if (strstr(line,"morse")!=NULL)
					*task_type=1;
				if (strstr(line,"1D")!=NULL)
					*task_type=2;
				if (strstr(line,"VCI")!=NULL)
					*task_type=3;			
				if (strstr(line,"VSCF")!=NULL)
					*task_type=4;			}
				
			if (strstr(line,"potential")!=NULL){	// parse potential
				if (strstr(line,"analytic")!=NULL)
					*pot_type=0;
				if (strstr(line,"polynom")!=NULL)
					*pot_type=1;
				if (strstr(line,"points")!=NULL)
					*pot_type=2;
			}

			if (strstr(line,"friend")!=NULL){		//parse friend
				if (strstr(line,"none")!=NULL)
					*friend=0;
				if (strstr(line,"gamess")!=NULL)
					*friend=1;
			}
			if (strstr(line,"}")!=NULL)				// if end of the group - break 
				break; 

		} while (fgets(line,80,fp)!=NULL);
	}
	}	
	fclose(fp);
} else 
	// can not open input file
	// message about fileopen failure
	printf("\n Cannot open input file %s, using default parameters for {control}\n",filename);
	

// PRINT entire {control} group

printf("\n \n {control: task");

switch (*task_type) {
	case 0:
		printf(" harmonic, potential");
		break;
	case 1:
		printf(" morse, potential");
		break;
	case 2:
		printf(" 1D, potential");
		break;
	case 3:
		printf(" VCI, potential");
		break;
	case 4:
		printf(" VSCF, potential");
		break;
}
switch (*pot_type) {
	case 0:
		printf(" analytic, friend");
		break;
	case 1:
		printf(" polynom, friend");
		break;
	case 2:
		printf(" points, friend");
		break;
}
switch (*friend) {
	case 0:
		printf(" none}\n");
		break;
	case 1:
		printf(" gamess}\n");
		break;
}

return;
}

/******************************************************************************/

void parse_size(const char *filename, int *nmodes, int *nbas)

/******************************************************************************/

/*
 * Purpose
 *  to parse {sizing} group
 *  in the form {sizing nmodes = 16, nbas = 64}
 * 
*/
{
char line[82],*pos;
FILE *fp;
	
if ((fp = fopen(filename,"r"))!=NULL)
{
	while (fgets(line,80,fp)!=NULL)
	{ 				// list all file lines till the end
	if (strstr(line,"{sizing")!=NULL)
	{			// find start of {sizing group (can be many)
		do 
		{										// list all file lines till the end (inside the group)
			pos=strtok(line,"=");
			while ( pos!=NULL )
			{
				if (strstr(pos,"nmodes")!=NULL)
				{												// parse nmodes
					pos=strtok(NULL,"=");
					*nmodes=atoi(pos);
				} else if (strstr(pos,"nbas")!=NULL){			// parse nbas
					pos=strtok(NULL,"=");
					*nbas=atoi(pos);
				} else 
				  pos=strtok(NULL,"=");
			}
			if (strstr(line,"}")!=NULL)				// if end of the group - break 
				break;

		} while (fgets(line,80,fp)!=NULL);
	}
	}	
	fclose(fp);
} else 
	// can not open input file
	// message about fileopen failure
	printf("\n Cannot open input file %s, using default parameters for {harmonic}\n",filename);
	

// PRINT will be done elsewhere

return;
}



/******************************************************************************/

int parse_power(const char *filename, const int npow, double fi[]) 

/******************************************************************************/
/*
 * Purpose
 *  to parse {polynom} group
 *  in the form 
 * 	{polynom
 *    0   0.0
 * 	  1   0.0
 * 	  2   1.0
 * 	...
 * 	}
 */
{
char line[82],*pos;
FILE *fp;
int i=0;
	
if ((fp = fopen(filename,"r"))!=NULL)
{
	while (fgets(line,80,fp)!=NULL)
	{ 				// list all file lines till the end
	if (strstr(line,"{polynom")!=NULL)
	{			// find start of {harmonic group (can be many)
		for (i=0;i<npow;i++) 
		{										// list all file lines till the end (inside the group)
			if (fgets(line,80,fp)==NULL) break;
			fi[i]=atof(&line[4]);
			if (strstr(line,"}")!=NULL)	break; // if end of the group - break
		} 
	}
	}	
	fclose(fp);
} else {
	// can not open input file
	// message about fileopen failure
	printf("\n Cannot open input file %s, using default parameters for {polynom}\n",filename);
	return 1;
}
// The group will be printed in the calculation program

return 0;
}

/******************************************************************************/

int parse_points(const char *filename, const int nbas, double pot[]) 

/******************************************************************************/
/*
 * Purpose
 *  to parse {points} group
 *  in the form 
 * 	{points
 * 	  0   0.0
 * 	  1   1.0
 * 	...
 * 	}

*/
{
char line[82],*pos;
FILE *fp;
int i=0;
	
if ((fp = fopen(filename,"r"))!=NULL)
{
	while (fgets(line,80,fp)!=NULL)
	{ 				// list all file lines till the end
	if (strstr(line,"{points")!=NULL)
	{			// find start of {harmonic group (can be many)
		for (i=0;i<nbas;i++) 
		{										// list all file lines till the end (inside the group)
			if (fgets(line,80,fp)==NULL) break;
			pot[i]=atof(&line[4]);
//			printf("line = %80s, i=%4i pot[i]=%20.15f\n",line,i,pot[i]);
			if (strstr(line,"}")!=NULL)	break; // if end of the group - break
		} 
	}
	}	
	fclose(fp);
} else {
	// can not open input file
	// message about fileopen failure
	printf("\n Cannot open input file %s, using default parameters for {points}\n",filename);
	return 1;
}
// The group will be printed in the calculation program

return 0;
}

/******************************************************************************/

void parse_1d(const char *filename,int *npow,int *nbas,double *xstart,double *xend) 

/******************************************************************************/
/*
 * Purpose
 *  to parse {1D} group
 *  in the form {1D: npow = 16, nbas = 64, xstart = -10.0, xend = 10.0}
 * 
*/
{
char line[82],*pos;
FILE *fp;
	
if ((fp = fopen(filename,"r"))!=NULL)
{
	while (fgets(line,80,fp)!=NULL)
	{ 				// list all file lines till the end
	if (strstr(line,"{1D")!=NULL)
	{			// find start of {harmonic group (can be many)
		do 
		{										// list all file lines till the end (inside the group)
			pos=strtok(line,"=");
			while ( pos!=NULL )
			{
				if (strstr(pos,"npow")!=NULL)
				{												// parse npow
					pos=strtok(NULL,"=");
					*npow=atoi(pos);
				} else if (strstr(pos,"nbas")!=NULL){			// parse nbas
					pos=strtok(NULL,"=");
					*nbas=atoi(pos);
				} else if (strstr(pos,"xstart")!=NULL){			// parse xstart
					pos=strtok(NULL,"=");
					*xstart=atof(pos);
				} else if (strstr(pos,"xend")!=NULL){			// parse xend
					pos=strtok(NULL,"=");
					*xend=atof(pos);
				} else 
				  pos=strtok(NULL,"=");
			}
			if (strstr(line,"}")!=NULL)				// if end of the group - break 
				break;

		} while (fgets(line,80,fp)!=NULL);
	}
	}	
	fclose(fp);
} else 
	// can not open input file
	// message about fileopen failure
	printf("\n Cannot open input file %s, using default parameters for {harmonic}\n",filename);
	

// PRINT entire {1D} group

printf("\n Entire {1D} line: \n {1D: ");
printf("npow = %5i, nbas = %5i, xstart = %15.10f, xend = %15.10f}\n",*npow,*nbas,*xstart,*xend);

return;
}


/******************************************************************************/

void parse_harmonic(const char *filename,double *w,int *nbas,double *xstart,double *xend) 

/******************************************************************************/
/*
 * Purpose
 *  to parse {harmonic} group
 *  in the form {harmonic: W = 1.0, nbas = 64, xstart = -10.0, xend = 10.0}
 * 
*/
{
char line[82],*pos;
FILE *fp;
	
if ((fp = fopen(filename,"r"))!=NULL){
	
	while (fgets(line,80,fp)!=NULL){ 				// list all file lines till the end
	if (strstr(line,"{harmonic")!=NULL){			// find start of {harmonic group (can be many)
		do {										// list all file lines till the end (inside the group)
			pos=strtok(line,"=");
			while ( pos!=NULL ){
				if (strstr(pos,"W")!=NULL){						// parse w
					pos=strtok(NULL,"=");
					*w=atof(pos);
				} else if (strstr(pos,"nbas")!=NULL){			// parse nbas
					pos=strtok(NULL,"=");
					*nbas=atoi(pos);
				} else if (strstr(pos,"xstart")!=NULL){			// parse xstart
					pos=strtok(NULL,"=");
					*xstart=atof(pos);
				} else if (strstr(pos,"xend")!=NULL){			// parse xend
					pos=strtok(NULL,"=");
					*xend=atof(pos);
				} else 
				  pos=strtok(NULL,"=");
			}
			if (strstr(line,"}")!=NULL)				// if end of the group - break 
				break;

		} while (fgets(line,80,fp)!=NULL);
	}
	}	
	fclose(fp);
} else 
	// can not open input file
	// message about fileopen failure
	printf("\n Cannot open input file %s, using default parameters for {harmonic}\n",filename);
	

// PRINT entire {harmonic} group

printf("\n Entire {harmonic} line: \n {harmonic: ");
printf("W = %15.10f, nbas = %5i, xstart = %15.10f, xend = %15.10f}\n",*w,*nbas,*xstart,*xend);

return;
}

/******************************************************************************/

void parse_morse(const char *filename,double *d,double *w,int *nbas,double *xstart,double *xend) 

/******************************************************************************/
/*
 * Purpose
 *  to parse {morse} group
 *  in the form {morse: D = 12.0, W = 1.0, nbas = 64, xstart = -4.0, xend = 15.0}
 * 
*/
{
char line[82],*pos;
FILE *fp;
	
if ((fp = fopen(filename,"r"))!=NULL){
	
	while (fgets(line,80,fp)!=NULL){ 				// list all file lines till the end
	if (strstr(line,"{morse")!=NULL){				// find start of {morse group (can be many)
		do {										// list all file lines till the end (inside the group)
			pos=strtok(line,"=");
			while ( pos!=NULL ){
				if (strstr(pos,"D")!=NULL){						// parse d
					pos=strtok(NULL,"=");
					*d=atof(pos);
				} else if (strstr(pos,"W")!=NULL){				// parse w
					pos=strtok(NULL,"=");
					*w=atof(pos);
				} else if (strstr(pos,"nbas")!=NULL){			// parse nbas
					pos=strtok(NULL,"=");
					*nbas=atoi(pos);
				} else if (strstr(pos,"xstart")!=NULL){			// parse xstart
					pos=strtok(NULL,"=");
					*xstart=atof(pos);
				} else if (strstr(pos,"xend")!=NULL){			// parse xend
					printf("found xend!\n");
					pos=strtok(NULL,"=");
					*xend=atof(pos);
				} else {
				  pos=strtok(NULL,"=");
				}
			}
			if (strstr(line,"}")!=NULL)				// if end of the group - break 
				break;

		} while (fgets(line,80,fp)!=NULL);
	}
	}	
	fclose(fp);
} else 
	// can not open input file
	// message about fileopen failure
	printf("\n Cannot open input file %s, using default parameters for {morse}\n",filename);
	

// PRINT entire {morse} group

printf("\n Entire {morse} line: \n {morse: ");
printf("D = %15.10f, W = %15.10f, nbas = %5i, xstart = %15.10f, xend = %15.10f}\n",*d,*w,*nbas,*xstart,*xend);

return;
}
