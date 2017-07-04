/* gamess.c - subroutines to parse GAMESS(US) output files
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
#include <string.h>
#include <math.h>

/******************************************************************************/

int gamess_get_nmodes_nbas(int *nmodes,int *nbas)

/******************************************************************************/

{
char line[82],*pos;
FILE *fp;
	
//		printf("in gamess_get_nmodes_nbas\n");
if ((fp = fopen("output","r"))!=NULL)
	{
	while (fgets(line,80,fp)!=NULL)					// list all file lines till the end
		if ((pos=strstr(line,"NZVAR ="))!=NULL)		// find NZVAR group (can be many)
			{*nmodes=atoi(pos+8);}
		else if ((pos=strstr(line,"NGRID="))!=NULL)// find NGRID group
			{*nbas=atoi(pos+8);}

	} else {
		printf("Can not open \'output\' GAMESS output file, exiting \n");
//		fclose(fp);
		exit(1);
	}

fclose(fp);
			
return 0;
}



/******************************************************************************/

int parse_gamess_out(const int nmodes, const int npairs, double harmonic[nmodes], double xrange[nmodes], double deltax[nmodes],double fdiag[nmodes][4], double foffdiag[npairs][6])

/******************************************************************************/
{
char line[82],*pos;
FILE *fp;
int ex;
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
//				printf("Found next HARMONIC line, line count= %d\n",freq_line_count);
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
				{
					fdiag[i][0]=atof(&line[15]);			// Hii
					ex=atoi(&line[28]);
					fdiag[i][0]*=pow(10,ex);
				}
				if(fgets(line,80,fp)!=NULL)
				{	
					fdiag[i][1]=atof(&line[15]);			// Tiii
					ex=atoi(&line[28]);
					fdiag[i][1]*=pow(10,ex);
				}

				if(fgets(line,80,fp)!=NULL)
				{
					fdiag[i][2]=atof(&line[15]);			// Uiiii
					ex=atoi(&line[28]);
					fdiag[i][2]*=pow(10,ex);
				}
				
				for(j=0;j<15;j++) 
					if(fgets(line,80,fp)!=NULL);// skip to the next mode
			}
		}

		if (strstr(line,"QFF>  MODE (I,J)")!=NULL)	// find QFF frequency group
		{
			for(i=0;i<npairs;i++)
			{
				if(fgets(line,80,fp)!=NULL)
				{
					foffdiag[i][0]=atof(&line[15]);			// Tiij
					ex=atoi(&line[28]);
					foffdiag[i][0]*=pow(10,ex);
				}
				
				if(fgets(line,80,fp)!=NULL)
				{
					foffdiag[i][1]=atof(&line[15]);			// Tjji
					ex=atoi(&line[28]);
					foffdiag[i][1]*=pow(10,ex);
				}
				
				if(fgets(line,80,fp)!=NULL)
				{
					foffdiag[i][2]=atof(&line[15]);			// Uiiij
					ex=atoi(&line[28]);
					foffdiag[i][2]*=pow(10,ex);
				}
				
				if(fgets(line,80,fp)!=NULL)
				{
					foffdiag[i][3]=atof(&line[15]);			// Ujjji
					ex=atoi(&line[28]);
					foffdiag[i][3]*=pow(10,ex);
				}
				
				if(fgets(line,80,fp)!=NULL)
				{
					foffdiag[i][4]=atof(&line[15]);			// Uiijj
					ex=atoi(&line[28]);
					foffdiag[i][4]*=pow(10,ex);
				}
				
				if(fgets(line,80,fp)!=NULL)
					if(fgets(line,80,fp)!=NULL);			// just skip 2 lines to the next mode pair
			}
		}

		
	}

	fclose(fp);
	
	} else {
		printf("Can not open \'output\' GAMESS output file, data should be provided in the input file \n");
//		fclose(fp);
//		exit(1);
	}

	
return 0;
}
