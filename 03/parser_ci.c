/*
 * parser_ci.c - functions to parse more groups in cvib.c
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
#include <string.h>

/******************************************************************************/

int parse_excit(const char *filename, const int nmodes, int excit[]) 

/******************************************************************************/
/*
 * Purpose
 *  to parse {excit} group
 *  in the form 
 * 	{excit mode max_excit
 *            0         4
 * 	          1         4
 * 	          2         4
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
	if (strstr(line,"{excit")!=NULL)
	{			// find start of {excit group (can be many)
		for (i=0;i<nmodes;i++) 
		{										// list all file lines till the end (inside the group)
			if (fgets(line,80,fp)==NULL) break;
			excit[i]=atoi(&line[12]);
			if (strstr(line,"}")!=NULL)	break; // if end of the group - break
		} 
	}
	}	
	fclose(fp);
} else {
	// can not open input file
	// message about fileopen failure
	printf("\n Cannot open input file %s, using default parameters for {excit}\n",filename);
	return 1;
}
// The group will be printed in the calling function

return 0;
}

/******************************************************************************/

int parse_ci(const char *filename, int *blocksize, int *maxit, double *convtol)

/******************************************************************************/
/*
 * Purpose
 *  to parse {ci} group
 *  in the form {ci blocksize =  16, maxit =  30, convtol = 1.00e-06}
 * 
*/
{
char line[82],*pos;
FILE *fp;
	
if ((fp = fopen(filename,"r"))!=NULL){
	
	while (fgets(line,80,fp)!=NULL){ 				// list all file lines till the end
	if (strstr(line,"{ci")!=NULL){					// find start of {ci group (can be many)
		do {										// list all file lines till the end (inside the group)
			pos=strtok(line,"=");
			while ( pos!=NULL ){
				if (strstr(pos,"blocksize")!=NULL){					// parse blocksize
					pos=strtok(NULL,"=");
					*blocksize=atoi(pos);
				} else if (strstr(pos,"maxit")!=NULL){				// parse maxit
					pos=strtok(NULL,"=");
					*maxit=atoi(pos);
				} else if (strstr(pos,"convtol")!=NULL){			// parse convtol
					pos=strtok(NULL,"=");
					*convtol=atof(pos);
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
	printf("\n Cannot open input file %s, using default parameters for {ci}\n",filename);
	

// PRINT will be done in calling function

return 0;
}

/******************************************************************************/

int parse_basis(const char * filename,const int nmodes, double * xrange, double * deltax)

/******************************************************************************/
/*
 * Purpose
 *  to parse {basis} group
 *  in the form 
 * 	{basis  mode        xrange       deltax
 *             0 +33.148356000 +4.419781000
 *             1 +33.523792000 +4.469839000
 *             2 +42.998096000 +5.733079000 
 *  	...
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
	if (strstr(line,"{basis")!=NULL)
	{			// find start of {basis group (can be many)
		for (i=0;i<nmodes;i++) 
		{										// list all file lines till the end (inside the group)
			if (fgets(line,80,fp)==NULL) break;
			xrange[i]=atof(&line[13]);
			deltax[i]=atof(&line[27]);
			if (strstr(line,"}")!=NULL)	break; // if end of the group - break
		} 
	}
	}	
	fclose(fp);
} else {
	// can not open input file
	// message about fileopen failure
	printf("\n Cannot open input file %s, using default parameters for {basis}\n",filename);
	return 1;
}
// The group will be printed in the calling function

return 0;
}

/******************************************************************************/

int parse_qff(const char * filename,const int nmodes,double fdiag[][4],double fpair[][6])

/******************************************************************************/
/*
 * Purpose
 *  to parse {qff} group
 *  in the form 
 * {qff 
 * Mode =  0
 * Hii    = +3.39701450e-01
 * Tiii   = +1.35559800e-06
 * Uiiii  = +1.73603670e+00
 * ...
 * Mode pair =  0 ,  1; pair #   0
 * Tiij  = -8.03824730e-01
 * Tjji  = +6.45822200e-07
 * Uiiij = +9.33730410e-06
 * Ujjji = +8.56116890e-06
 * Uiijj = +1.59735330e+00
 * ...
 * } 
 */
{
char line[82],*pos;
FILE *fp;
int i=0;
int npair = nmodes * (nmodes-1);
	
if ((fp = fopen(filename,"r"))!=NULL)
{
	while (fgets(line,80,fp)!=NULL)
	{ 				// list all file lines till the end
	if (strstr(line,"{qff")!=NULL)
	{			// find start of {harmonic group (can be many)
		for (i=0;i<nmodes;i++) 
		{										// list all file lines till the end (inside the group)
			if (fgets(line,80,fp)==NULL) break; // skip first line
			if (fgets(line,80,fp)==NULL) break;
			fdiag[i][0]=atof(&line[8]); // Hii
			if (fgets(line,80,fp)==NULL) break;
			fdiag[i][1]=atof(&line[8]); // Tiii
			if (fgets(line,80,fp)==NULL) break;
			fdiag[i][2]=atof(&line[8]); // Uiiii
			if (strstr(line,"}")!=NULL)	break; // if end of the group - break
		} 
		for (i=0;i<npair;i++) 
		{										// list all file lines till the end (inside the group)
			if (fgets(line,80,fp)==NULL) break; // skip first line
			if (fgets(line,80,fp)==NULL) break;
			fpair[i][0]=atof(&line[8]); // Tiij
			if (fgets(line,80,fp)==NULL) break;
			fpair[i][1]=atof(&line[8]); // Tjji
			if (fgets(line,80,fp)==NULL) break;
			fpair[i][2]=atof(&line[8]); // Uiiij
			if (fgets(line,80,fp)==NULL) break;
			fpair[i][3]=atof(&line[8]); // Ujjji
			if (fgets(line,80,fp)==NULL) break;
			fpair[i][4]=atof(&line[8]); // Uiijj
			if (strstr(line,"}")!=NULL)	break; // if end of the group - break
		} 
	}
	}	
	fclose(fp);
} else {
	// can not open input file
	// message about fileopen failure
	printf("\n Cannot open input file %s, using default parameters for {qff}\n",filename);
	return 1;
}
// The group will be printed in the calling function

return 0;
}
