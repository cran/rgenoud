/*

  RGENOUD (limited version)

  Walter R. Mebane, Jr.
  Cornell University
  http://macht.arts.cornell.edu/wrm1
  wrm1@macht.arts.cornell.edu

  Jasjeet Singh Sekhon 
  Harvard University and Lamarck, Inc.
  http://jsekhon.fas.harvard.edu/
  jsekhon@fas.harvard.edu

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/genoud.cpp,v 1.18 2002/10/17 03:45:20 jsekhon Exp $

*/

#include "genoud.h"
#include "unif.h"

/* unif.h integer definition */
long NewUnifSeed[MAXTHREADS];
long RandIntSeed[MAXTHREADS];
long ThreadNumber;

extern double func4g(double *X);

double genoud(struct GND_IOstructure *Structure)
{

  extern long NewUnifSeed[MAXTHREADS];
  extern long RandIntSeed[MAXTHREADS];
  extern long ThreadNumber;

  MATRIX 
    domains,      /*Matrix for Domains*/
    final_mat;    /*The final Domain*/

  VECTOR 
    LowerBounds,
    UpperBounds,
    X;            /*Initial values for the variables*/

  INDEX 
    fin;          /*Size of final matrix*/

  int 
    i,            /*Counter variable*/
    nvars;        /*Remaining variables after p-variables were eliminated*/

  time_t start_time,
         stop_time;
  double delta_time;
  long   hours, minutes, seconds;
  char   time_str[27];

  double peak_val;

  static long BaseNewUnifSeed=81282,
    BaseRandIntSeed=53058;
  static short firsttime=1;

  FILE *output;

  /* Lamarck Child Test Variables */
  // char LVMchar[1000];
  long LVMreturn;

/********************************************************************************/

  LVMreturn=0;
  time(&start_time);
  strcpy(time_str, ctime(&start_time));

  /* Fault Tolerant MinMax */
  if (Structure->MinMax > 0) Structure->MinMax=1;
  else Structure->MinMax=0;

  if (Structure->OutputType==0) {
    output=stdout;
  }
  else if (Structure->OutputType==1) {
    if((output = fopen(Structure->OutputPath, "w")) == NULL) {
      fprintf(output,"%s", Structure->OutputPath);

      if (Structure->MinMax==1)
	return(ERROR_CODE);
      else
	return(-1*ERROR_CODE);
    }
  }
  else if (Structure->OutputType==2) {
    if((output = fopen(Structure->OutputPath, "a")) == NULL) {
      fprintf(output,"%s", Structure->OutputPath);
      if (Structure->MinMax==1)
	return(ERROR_CODE);
      else
	return(-1*ERROR_CODE);
    }
  }
  else {
    if (Structure->MinMax==1)
      return(ERROR_CODE);
    else
      return(-1*ERROR_CODE);
  }

  fprintf(output,"\n\n%s",time_str);

  ThreadNumber=Structure->ThreadNumber;
  if (ThreadNumber > MAXTHREADS) {
    fprintf(output,"\nERROR: NO MORE THAN %d THREADS ALLOWED\n\n", MAXTHREADS);
    if (Structure->MinMax==1)
      return(ERROR_CODE);
    else
      return(-1*ERROR_CODE);
  }
  if (Structure->ProvideSeeds == 1) {
    /*
      Only toggle the instance number if we have threads! */
    NewUnifSeed[ThreadNumber] = Structure->UnifSeed;
    RandIntSeed[ThreadNumber] = Structure->IntSeed;
  }
  else {
    /* If a Seed is NOT provided, use the base random number and run from that base!
       In other words, might as well the ThreadNumber equal to 0 
    */
    if (firsttime==1) {
      NewUnifSeed[0] = BaseNewUnifSeed;
      RandIntSeed[0] = BaseRandIntSeed;	
      firsttime=0;
    }
    ThreadNumber = 0;
  }

  fin.r =   Structure->nvars;            /*total number of inequalities + domains*/
  fin.c =   Structure->nvars+2;          /*x2 variables + lower limits + upper limits*/

  nvars = Structure->nvars;
  if (nvars > MAX_VAR)
    fprintf(output, "Too many variables - Increase MAX_VAR in header file\n");

  /* Boundary Enforcement: 1: no trespassing; 0: trespassing but no camping; -1 anything goes! */

  /*Allocating memory for all the vectors and matrices*/
  final_mat = matrix(1,fin.r,1,fin.c);
  domains = matrix(1,nvars,1,3);
  LowerBounds = Gvector(1, nvars);
  UpperBounds = Gvector(1, nvars);
  X = Gvector(1,nvars);

  /* SETUP DOMAINS */

  /* alter the domains when we are using integers because of the "open
     set" problem.  We only extend the UPPER domain bound */
  if (Structure->DataType==1) {
    for(i=0; i<Structure->nvars; i++)
      Structure->Domains[i][1] = Structure->Domains[i][1] + 0.99999;
  }
  
  for(i=1; i<=Structure->nvars; i++)
    {
      domains[i][1] = Structure->Domains[i-1][0];
      domains[i][2] =  (double) i;
      domains[i][3] = Structure->Domains[i-1][1];
    }
  
  for (i=1; i<=nvars; i++)
    {
      LowerBounds[i] = domains[i][1];
      UpperBounds[i] = domains[i][3];
    }

  /*Initialization*/
  print_domains(domains,nvars,Structure->DataType, output);


#ifdef SQL_DEFINE
  // Initialize the SQL database for network evaluation
  if (Structure->Network==1)
    {
      short sval;

      sval = (short) InitializeMySQL(NULL, Structure->DBname);
      if (sval < 0)
	{
	  Structure->Status = -1;
	  return(sval);
	}
      
      // Let's setup our databases
      sval = setup_database(Structure->DBname, Structure->AgentName);
      if (sval < 0)
	{
	  Structure->Status = -1;
	  return(sval);
	}
    }
#endif

  //  initialize(final_mat,fin);
  // find_final_mat1(LowerBounds,UpperBounds,final_mat,nvars,fin.c);

  /* --- Lamarck Child Test ---- */
	//  strcpy(LVMchar,"Genoud.StartUp");
	// LVM_TextToGenome(LVMchar, Structure->pvm->GenomeBuffer, TRUE);
	// LVMreturn = LVM_SubEval(Structure->pvm, "Genoud Child Test", 0);
#ifdef MS_WINDOWS
  //LVMreturn = LVM_EvalSubScript(Structure->pvm, "GENOUD Startup Logo", "Genoud.startup end", 0);
#endif
	

  /*This procedure initializes the initial population with the values generated*/
  /*and applies the genetic operators and reproduces a new generation; evaluates*/
  /*each agent of the new generation, and again goes through the cycle of*/
  /*reproduction and evaluation, for the number of times, user has specified*/
  /*and prints out a final population with best agents*/
  if (Structure->DataType==1) {
    peak_val = JaIntegerOptimization(Structure, X, domains, output);
  }
  else {
    peak_val = optimization(Structure, X, domains, output);
  }

  /* print out the hessian matrix */
  /*
    dohessians(nvars, runcases, 9, X+1, func4g, func4g);
  */

  /* free memory */
  free_matrix(final_mat,1,fin.r,1);
  free_matrix(domains, 1, nvars, 1);
  free_vector(LowerBounds,1);
  free_vector(UpperBounds,1);
  free_vector(X,1);
  
  /* print final numbers and the time this has taken */
  fprintf(output, "\n");
  fprintf(output, "Solution Found Generation %d\n", Structure->oPeakGeneration);
  fprintf(output, "Number of Generations Run %d\n", Structure->oGenerations);

  time(&stop_time);

  strcpy(time_str, ctime(&stop_time));
  fprintf(output,"\n\n%s",time_str);

  delta_time = difftime(stop_time, start_time);
  hours   = (int) (delta_time/3600);
  minutes = (int) (delta_time - (hours*3600))/60;
  seconds = (int) delta_time - (hours*3600) - (minutes*60);
  
  fprintf(output,"Total run time : %d hours %d minutes and %d seconds\n", 
	  hours, minutes, seconds);

  fflush(output);

  if (Structure->OutputType==1) fclose(output);

  if (peak_val==ERROR_CODE)
    {
    if (Structure->MinMax==1)
      return(ERROR_CODE);
    else
      return(-1*ERROR_CODE);
    }
  else
    return(peak_val);
}

