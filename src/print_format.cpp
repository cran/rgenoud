/*

  RGENOUD

  Walter R. Mebane, Jr.
  Cornell University
  http://macht.arts.cornell.edu/wrm1
  <wrm1@macht.arts.cornell.edu>

  Jasjeet Singh Sekhon 
  UC Berkeley
  http://sekhon.polisci.berkeley.edu
  <sekhon@berkeley.edu>

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/print_format.cpp,v 2.15 2005/10/29 06:14:44 jsekhon Exp jsekhon $

*/


#include "genoud.h"


/********************************************************************************/
/*  ReadPopulation():                                                           */
/*                                                                              */
/*  This function reads in an old population file and initializes the           */
/*  newpopulation with it.                                                      */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

long ReadPopulation(double **Data, long NewPopSize, long NewVars, FILE *output, FILE *fp)
{
  char ctmp[MAXPATH];
  long generation, PopSize, nvars, UsePopSize;
  long i, j, ltmp;
  double **OldData;
  short trip=0;

  /* This reads the "Generations:" name */

  while(!(feof(fp))) {

    /* pos = ftell(fp); */

    fscanf(fp, "%s", ctmp);
    fscanf(fp, " %d", &generation);
    fprintf(output, "Old PopFile:\n");
    fprintf(output, "Generation: %d\n", generation);
    /* This reads the "Population" name */
    fscanf(fp, "%s", ctmp);
    /* This reads the "Size:" name */
    fscanf(fp, "%s", ctmp);
    fscanf(fp, " %d", &PopSize);
    fprintf(output, "Population Size: %d\n", PopSize); 
    /* This reads the "Variables:" name */
    fscanf(fp, "%s", ctmp);
    fscanf(fp, " %d", &nvars);
    fprintf(output, "Nvars: %d\n", nvars); 

    if (trip==0) {
      OldData = JaMatrixAllocate(PopSize+2, nvars+2);

      if (nvars!=NewVars) return(0);
      trip++;
    }
    
    /* loop over the main data part */
    for (i=1; i<=PopSize; i++) {
      fscanf(fp,"%d",&ltmp);
      for (j=0; j<=nvars; j++) {
	fscanf(fp,"%lf", &OldData[i][j]);
      }
    }

  }

  /* Map OldData to Data */
  if (NewPopSize < PopSize) UsePopSize = NewPopSize;
  else UsePopSize=PopSize;
  for (i=1; i<=UsePopSize; i++) {
    Data[i][nvars+1] = 0.0;
    for (j=0; j<=nvars; j++) {
      Data[i][j] = OldData[i][j];
    }
  }

  /* let's print the population file */
  fprintf(output, "\nRead in Population: UsePopSize: %d\n", UsePopSize); 
  for (i=1; i<=UsePopSize; i++) {
    fprintf(output, "%d \t", i); 
    for (j=0; j<=nvars; j++) {
      fprintf(output, "%e \t", Data[i][j]);
    }
    fprintf(output, "\n");
  }
  fflush(output);

  JaMatrixFree(OldData, PopSize);
  return(PopSize);
} /* end Read Population */



/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   print_domains()                              */
/*                                                                              */
/*           SYNOPSIS          :   void print_domains(equal,t_equ)              */
/*                                                                              */
/*           DESCRIPTION       :   This function prints the matrix passed, on to*/
/*                                  the standard output, in the format of       */
/*                                  domains.                                    */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   main()                                       */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/




void print_domains(MATRIX equal, int t_equ, short DataType, FILE *output)
     /*
       MATRIX equal;   the domains matrix, with the upper and lower limits
       int t_equ;      *the total number of domains
     */
{
  int i,j;

  fprintf(output,"Domains:\n");
  //Integer
  if (DataType==1) 
  {
      for(i=1; i<=t_equ; i++)
      {
	  for(j=1; j<=3; j++)
	  {
	      if(j == 2)
		  fprintf(output,"  <=  X%-2d  <=   ",(int)equal[i][j]);
	      else
		  fprintf(output," %d ",(int) equal[i][j]);
	  }
	  fprintf(output,"\n");
      }
  } else {
      for(i=1; i<=t_equ; i++)
      {
	  for(j=1; j<=3; j++)
	  {
	      if(j == 2)
		  fprintf(output,"  <=  X%-2d  <=   ",(int)equal[i][j]);
	      else
		  fprintf(output," %e ",equal[i][j]);
	  }
	  fprintf(output,"\n");
      }
  }
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   print_matrix()                               */
/*                                                                              */
/*           SYNOPSIS          :   void print_matrix(lr,ur,lc,uc,mat)           */
/*                                                                              */
/*           DESCRIPTION       :   This function prints a given double matrix,   */
/*                                  on to the standard output                   */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   main(),                                      */
/*                                 optimization().                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/


void print_matrix(int lr, int ur, int lc, int uc, MATRIX mat, FILE *output)
     /*
       int lr,ur,lc,uc;
       MATRIX mat;
     */
{
  int i,j;

  for(i=lr; i<=ur; i++)
    {
      fprintf(output,"\n");
      for(j=lc; j<=uc; j++)
        fprintf(output,"%5.2f\t",mat[i][j]);
    }
  fprintf(output,"\n\n");
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   print_population()                           */
/*                                                                              */
/********************************************************************************/

void print_population(long popsize, long nvars, long generation, long lexical, double **foo, FILE *out)
{
  long i,j;

  if (lexical < 2)
    {
      fprintf(out,"Generation: %d \t Population Size: %d \t Fit Values: 1 \t Variables: %d\n\n", 
	      generation, popsize, nvars);
      for(i = 1; i <= popsize; i++)
	{
	  fprintf(out,"%d \t %e \t",i, foo[i][0]);
	  for (j = 1; j <= nvars; j++)
	    {
	      fprintf(out,"%e \t ",foo[i][j]);
	    }
	  fprintf(out,"\n");
	}
      fprintf(out,"\n\n");
    }
  else 
    {
      long lexical_end = lexical-1+nvars+2;

      fprintf(out,"Generation: %d \t Population Size: %d \t Fit Values: %d \t Variables: %d\n\n", 
	      generation, popsize, lexical, nvars);
      for(i = 1; i <= popsize; i++)
	{
	  fprintf(out,"%d \t ", i);

	  /* print lexical fit values */
	  fprintf(out,"%e \t ",foo[i][0]);
	  for(j=(nvars+2);j<lexical_end;j++)
	    {
	      fprintf(out,"%e \t ",foo[i][j]);
	    }		      
	  
	  /* print variables */
	  for (j = 1; j <= nvars; j++)
	    {
	      fprintf(out,"%e \t ",foo[i][j]);
	    }
	  fprintf(out,"\n");
	}
      fprintf(out,"\n\n");
    }
} /* end */


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   print_vector()                               */
/*                                                                              */
/*           SYNOPSIS          :   void print_vector(arr,l,u)                   */
/*                                                                              */
/*           DESCRIPTION       :   This function prints a given double vector,   */
/*                                  on to the standard output                   */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   main()                                       */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void print_vector(VECTOR arr, int l, int u, FILE *output)
     /*
       VECTOR arr;
       int l,u;
     */
{
  int i;

  for(i=l; i<=u; i++)
    fprintf(output,"%5.2f\t",arr[i]);
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   print_ivector()                              */
/*                                                                              */
/*           SYNOPSIS          :   void print_ivector(arr,l,u)                  */
/*                                                                              */
/*           DESCRIPTION       :   This function prints a given integer vector, */
/*                                  on to the standard output                   */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :                                                */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/


void print_ivector(IVECTOR arr, int l, int u, FILE *output)
     /*
       int l,u;
       IVECTOR arr;
     */
{
  int i;

  for(i=l; i<=u; i++)
    fprintf(output,"%d\t",arr[i]);
  fprintf(output,"\n\n");
}

