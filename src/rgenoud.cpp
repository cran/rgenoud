/*

  RGENOUD

  Walter R. Mebane, Jr.
  Cornell University
  http://macht.arts.cornell.edu/wrm1
  wrm1@macht.arts.cornell.edu

  Jasjeet Singh Sekhon 
  Harvard University
  http://jsekhon.fas.harvard.edu/
  jsekhon@fas.harvard.edu

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/rgenoud.cpp,v 1.31 2005/03/01 06:36:36 jsekhon Exp $

*/

#include <R.h>
#include "genoud.h"

void *global_fn;
void *global_rho;
void *global_fn_optim;
int global_nvars;
int global_max;

double genoud(struct GND_IOstructure *Structure);
double EvalFitness(double *X, long *Status);
double setupGenoud(struct GND_IOstructure *MainStructure, int nvars);

extern "C" 
{

  // <Rdefines.h> must appear inside the {extern "C"} declaration
  // because this header file, unlike <R.h>, does not have an {#ifdef
  // __cplusplus} statment included.
#include <Rdefines.h>

  // mkanswer
  SEXP mkans(double fit, double *oResults, double *oGradients, long *oP, long oGenerations,
	     long oPeakGeneration, long oPopSize, long nvars)
  { 
    SEXP ans;
    long length, i, indx, operators;
    
    operators=9;
    length=(nvars*2) + 4 + operators;

    PROTECT(ans=allocVector(REALSXP,length));
    REAL(ans)[0] = fit;
    REAL(ans)[1] = (double) oGenerations;
    REAL(ans)[2] = (double) oPeakGeneration;
    REAL(ans)[3] = (double) oPopSize;
    // include results
    indx = 3;
    for (i=0; i<nvars; i++) {
      indx++;
      REAL(ans)[indx] = oResults[i];
    }
    // include gradients
    for (i=0; i<nvars; i++) {
      indx++;
      REAL(ans)[indx] = oGradients[i];
    }
    // include the actual operator count
    for (i=0; i<operators; i++) {
      indx++;
      REAL(ans)[indx] = oP[i];
    }
    UNPROTECT(1);

    return(ans);
  } // end of mkans


  double feval(SEXP fn, SEXP rho, double *X, int parameters)
  {
    SEXP R_fcall, x;
    double ret;
    int i;
    
    PROTECT(x = allocVector(REALSXP, parameters));

    for (i=0; i<parameters; i++)
      {
	REAL(x)[i] = X[i];
      }

    PROTECT(R_fcall = lang2(fn, R_NilValue));
    SETCADR(R_fcall, x);
    ret =  REAL(eval(R_fcall, rho))[0];
    UNPROTECT(2);

    return(ret);
  } // end of feval

  double genoud_optim(double *X, int parameters)
  {
    extern void *global_rho;
    extern void *global_fn_optim;

    SEXP fn_optim, rho;
    SEXP ans, R_fcall, x;
    double fit;
    int i;

    fn_optim= (SEXP) global_fn_optim;
    rho= (SEXP) global_rho;
    
    PROTECT(x = allocVector(REALSXP, parameters));

    for (i=0; i<parameters; i++)
      {
	REAL(x)[i] = X[i];
      }

    PROTECT(R_fcall = lang2(fn_optim, R_NilValue));
    SETCADR(R_fcall, x);

    ans = eval(R_fcall, rho);
    fit = REAL(ans)[0];

    for(i=0; i<parameters; i++)
      {
	X[i] = REAL(ans)[i+1];
      }

    UNPROTECT(2);
    return(fit);
  } // end of genoud_optim()



  // See "Writing R Extensions" page 31.
  SEXP rgenoud(SEXP fn, SEXP rho,
	       SEXP nvars, SEXP pop_size, SEXP max_generations, SEXP wait_generations,
	       SEXP n_starting_values, SEXP starting_values,
	       SEXP P, SEXP Domains, 
	       SEXP max, SEXP gradient_check, SEXP boundary_enforcement,
	       SEXP solution_tolerance, SEXP BFGS, SEXP data_type_int,
	       SEXP provide_seeds, SEXP unif_seed, SEXP int_seed,
	       SEXP print_level, SEXP share_type, SEXP instance_number,
	       SEXP MemoryMatrix, SEXP Debug,
	       SEXP output_path, SEXP output_type, SEXP project_path,
	       SEXP hard_generation_limit,
	       SEXP fn_optim, SEXP optim) 
{

  // Let's do the R checking right here

  extern void *global_fn;
  extern void *global_rho;
  extern void *global_fn_optim;
  extern int global_nvars;
  extern int global_max;

  SEXP ret;
  double solution;
  long parameters, i, j;

  double *Results, *Gradients;

  if(!isEnvironment(rho)) 
    error ("`rho' should be an environment");

  parameters = asInteger(nvars);
  global_fn       = fn;
  global_rho      = rho;
  global_fn_optim = fn_optim;
  global_nvars    = parameters;
  global_max      = asInteger(max);


  // setup GENOUD
  struct GND_IOstructure *MainStructure;
  MainStructure = (struct GND_IOstructure *) malloc(sizeof(struct GND_IOstructure));

  double **domains;
  domains = (double **) malloc(parameters*sizeof(double));
  for (i=0; i<parameters; i++) {
    domains[i] = (double *) malloc(2*sizeof(double));
  }

  for (j=0; j<2; j++) {
    for (i=0; i<parameters; i++) {
      domains[i][j] = REAL(Domains)[i + j*parameters];
    }
  }

  // starting values
  double **StartingValues;
  int nStartingValues;
  nStartingValues = asInteger(n_starting_values);
  if (nStartingValues > 0) {
    StartingValues = (double **) malloc(nStartingValues*sizeof(double));
    for (i=0; i<nStartingValues; i++) {
      StartingValues[i] = (double *) malloc(parameters*sizeof(double));
    }

    for(i=0; i<parameters; i++) {
      StartingValues[0][i] = REAL(starting_values)[i];
    }
  }

  MainStructure->Optim=asInteger(optim);
  MainStructure->nvars=parameters;
  MainStructure->PopSize=asInteger(pop_size);
  MainStructure->MaxGenerations=asInteger(max_generations);
  MainStructure->WaitGenerations=asInteger(wait_generations);
  MainStructure->HardGenerationLimit=asInteger(hard_generation_limit);
  MainStructure->nStartingValues=nStartingValues;
  MainStructure->StartingValues=StartingValues;
  MainStructure->P[0]=REAL(P)[0];
  MainStructure->P[1]=REAL(P)[1];
  MainStructure->P[2]=REAL(P)[2];
  MainStructure->P[3]=REAL(P)[3];
  MainStructure->P[4]=REAL(P)[4];
  MainStructure->P[5]=REAL(P)[5];
  MainStructure->P[6]=REAL(P)[6];
  MainStructure->P[7]=REAL(P)[7];
  MainStructure->P[8]=REAL(P)[8];
  MainStructure->Domains=domains;
  MainStructure->MinMax=asInteger(max);
  MainStructure->GradientCheck=asInteger(gradient_check);
  MainStructure->BoundaryEnforcement=asInteger(boundary_enforcement);
  MainStructure->SolutionTolerance=asReal(solution_tolerance);
  MainStructure->UseBFGS=asInteger(BFGS);

  MainStructure->MemoryUsage=asInteger(MemoryMatrix);
  MainStructure->Debug=asInteger(Debug);

  MainStructure->InstanceNumber=asInteger(instance_number);

  MainStructure->ProvideSeeds=asInteger(provide_seeds);
  MainStructure->UnifSeed=asInteger(unif_seed);
  MainStructure->IntSeed=asInteger(int_seed);
  MainStructure->PrintLevel=asInteger(print_level);
  MainStructure->DataType=asInteger(data_type_int);

  /* 
     Share Type:
     (0) no reading of the existing project file and no looking at the public population file
     (1) reading of any existing project file, but no examining of public population file
     (2) NO reading of any existing project file but examination of public population file
     (3) BOTH reading of any existing project file AND examination of public population file
  */
  MainStructure->ShareType=asInteger(share_type);

  //Paths
  char OutputPath[1000], ProjectPath[1000];
  strcpy(OutputPath,STRING_VALUE(output_path));
  strcpy(ProjectPath,STRING_VALUE(project_path));
  MainStructure->OutputPath=OutputPath;
  MainStructure->ProjectPath=ProjectPath;
  MainStructure->OutputType=asInteger(output_type);

  //Fixed Stuff
  MainStructure->Network=0;
    

  /* output data structures */
  Results = (double *)  malloc(parameters*sizeof(double));
  Gradients = (double *)  malloc(parameters*sizeof(double));
  
  MainStructure->oResults=Results;
  MainStructure->oGradients=Gradients;

  solution = setupGenoud(MainStructure, parameters);

  ret = mkans(solution, 
	      MainStructure->oResults, MainStructure->oGradients, MainStructure->oP,
	      MainStructure->oGenerations, MainStructure->oPeakGeneration,
	      MainStructure->oPopSize, MainStructure->nvars);

  // Free memory
  free(MainStructure);
  for (i=0; i<parameters; i++) 
    free(domains[i]);
  free(domains);
  free(Results);
  free(Gradients);
    

  if (nStartingValues > 0) {
    free(StartingValues[0]);
  free(StartingValues);
  }
  
  // return the solution
  return(ret);
  
  /*
    SEXP R_fcall, ans, tmp;
    double dtmp;
    
    if(!isFunction(fn)) error("The first argument, `fn' must be a function");
    if(!isEnvironment(rho)) error ("`rho' should be an environment");
    
    PROTECT(R_fcall = lang2(fn, R_NilValue));
    PROTECT(ans = allocVector(VECSXP, 1));
    VECTOR(ans)[0] = eval(R_fcall, rho);
    dtmp = REAL(ans)[0];
    printf("ans[0]: %lf\n", dtmp);
    UNPROTECT(2);
    return(ans);
  */
} // end of rgenoud()
  
} // end of extern "C"


double EvalFitness(double *X, long *Status)
{
  extern void *global_fn;
  extern void *global_rho;
  extern int  global_nvars;
  extern int  global_max;

  double fit;
  int isFinite=0;

  fit = feval( (SEXP) global_fn, (SEXP) global_rho, X, global_nvars);  

  *Status = 2;

  isFinite = R_finite(fit);
  if (!isFinite)
    {
      if (global_max)
	{
	  return(-1*BIGNUMBER);
	}
      else
	{
	  return(BIGNUMBER);
	}
    }

  return(fit);
} // end of EvalFitness


double setupGenoud(struct GND_IOstructure *MainStructure, int nvars)
{
  int i, 
    pop_size, maxGenerations, waitGenerations;
  double P0, P1, P2, P3, P4, P5, P6, P7, P8;
  double **Domains;
  double *Results, *Gradients, SolutionTolerance;
  double solution;
  double (*p)(double *X, long *Status);
  long   (*pM)(double **Population, long npop, long nvars);
  char public_population_path[1000];
  short int output_type;
  int Status=0;
  
  short int MinMax;
  short int BoundaryEnforcement;
  short int GradientCheck;
  short int UseBFGS;
  short int DataType, PrintLevel, ShareLevel;

  SolutionTolerance=0.0001;
  MinMax = 1;
  
  output_type = 0; /* to stndout */

  /* the function to optimize */
  /* p = test_abs; */
  p = EvalFitness;
  pM = NULL;
  
  /*setup memory */
  Results = (double *)  malloc(nvars*sizeof(double));
  Gradients = (double *)  malloc(nvars*sizeof(double));
  Domains = (double **) malloc(nvars*sizeof(double));
  for (i=0; i<nvars; i++) {
    Domains[i] = (double *) malloc(2*sizeof(double));
  }

  for (i=0; i<nvars; i++) {
    Domains[i][0] = -2;
    Domains[i][1] =  2;
  }
    
  BoundaryEnforcement = 2; /* 0=anything goes, 1: regular; 2: no trespassing! */
  GradientCheck = 1;
  UseBFGS = 1;
  PrintLevel=2;
  /* 
     ShareLevel, not yet used.
     (n) Dump (or read) this many people to/from the public population file
  */
  ShareLevel = 1;

  DataType=0;
  /* the population Size */
  pop_size = 15;
  /* the number of generation */
  maxGenerations=100;
  waitGenerations=10;
  /* operators */
  P0=10;
  P1=10;
  P2=10;
  P3=10;
  P4=10;
  P5=10;
  P6=10;
  P7=10;
  P8=0;

  MainStructure->AgentFit=(double (*)(double *, long int *)) *p;
  MainStructure->AgentFitMatrix=(long int (*)(double **, long int, long int))*pM;
  // MainStructure->nvars=nvars;
  // MainStructure->PopSize=pop_size;
  // MainStructure->MaxGenerations=maxGenerations;
  // MainStructure->WaitGenerations=waitGenerations;
  // MainStructure->P[0]=P0;
  //MainStructure->P[1]=P1;
  //MainStructure->P[2]=P2;
  //MainStructure->P[3]=P3;
  //MainStructure->P[4]=P4;
  //MainStructure->P[5]=P5;
  //MainStructure->P[6]=P6;
  //MainStructure->P[7]=P7;
  //MainStructure->P[8]=P8;
  //MainStructure->Domains=Domains;
  //MainStructure->OutputType=output_type;
  //MainStructure->MinMax=MinMax;
  //MainStructure->GradientCheck=GradientCheck;
  //MainStructure->BoundaryEnforcement=BoundaryEnforcement;
  //MainStructure->SolutionTolerance=SolutionTolerance;
  MainStructure->Status=Status;
  //MainStructure->UseBFGS=UseBFGS;

  /* output data structures */
  // MainStructure->oResults=Results;
  // MainStructure->oGradients=Gradients;
  // MainStructure->oP[0]=oP[0]; */
  MainStructure->oGenerations=0;
  MainStructure->oPeakGeneration=0;
  MainStructure->oPopSize=0;

  MainStructure->ThreadNumber=0;
  // MainStructure->ProvideSeeds=0;

  // MainStructure->UnifSeed=81282;
  // MainStructure->IntSeed=53058;

  /* Paths */
  MainStructure->PublicPopulationPath=public_population_path;

  /* printlevel: 0=minimal, 1=normal, 2=detailed */
  // MainStructure->PrintLevel=PrintLevel;

  MainStructure->AllowDynamicUpdating=0;
  // MainStructure->DataType=DataType;
  MainStructure->ShareLevel=ShareLevel;

  solution = genoud(MainStructure);

  // free memory
  free(Results);
  free(Gradients);
  for (i=0; i<nvars; i++)
    free(Domains[i]);
  free(Domains);

  return(solution);
} // end of setupGenoud
