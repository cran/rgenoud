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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/evaluate.cpp,v 1.23 2004/02/02 08:01:29 jsekhon Exp $

*/

#include "genoud.h"
#include "gradient.h"

#ifndef OPTIM
extern double genoud_optim(double *X, int nvars);
#endif

long Gnvars[MAXINSTANCES];
struct GND_IOstructure *ExternStructure;

int JaIntegerCMP(double **a, double **b) 
{
  extern long Gnvars[MAXINSTANCES];
  extern struct GND_IOstructure *ExternStructure;


  long i = 0;
  long nvars;

  nvars=Gnvars[ExternStructure->InstanceNumber];

  for (i=1; i<=nvars; i++) {
    if ( (int) a[0][i] != (int) b[0][i])
      break;
  }

  if ( (int) a[0][i] >  (int) b[0][i]) i = 1;
  else if ( (int) a[0][i] <  (int) b[0][i]) i = -1;

  return i;
} /* end of JaIntegerCMP */


int JaDoubleCMP(double **a, double **b) 
{
  extern long Gnvars[MAXINSTANCES];
  extern struct GND_IOstructure *ExternStructure;

  long i = 0;
  long nvars;

  nvars=Gnvars[ExternStructure->InstanceNumber];

  for (i=1; i<=nvars; i++) {
    if ( a[0][i] != b[0][i])
      break;
  }

  if ( a[0][i] > b[0][i]) i = 1;
  else if ( a[0][i] < b[0][i]) i = -1;

  return i;
} /* end of JaCMP */


/* Cummulative probability on crossover */
/* Random probability on mutation       */
/* NO multiple hits per agent possible  */


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   optimization()                               */
/*                                                                              */
/*           SYNOPSIS          :   double optimization(X,x1,x2,fin_mat,rc,tot_eq) */
/*                                                                              */
/*           DESCRIPTION       :   This procedure initializes the population    */
/*                                  with the values X passed from main, and     */
/*                                  evaluates them.  After assigning weight     */
/*                                  for each member of the populaiton, a group  */
/*                                  of them are chosen to reproduce and a group */
/*                                  is chosen to die.  Genetic operators are    */
/*                                  applied and a new generation is produced    */
/*                                  to replace the members that died.  This     */
/*                                  cycle continues for the number of times,    */
/*                                  user specifies in the input file            */
/*                                                                              */
/*           FUNCTIONS CALLED  :   assign_probab(),                             */
/*                                 evaluate(),                                  */
/*                                 find_cum_probab(),                           */
/*                                 find_live_die(),                             */
/*                                 find_parent(),                               */
/*                                 ivector(),                                   */
/*                                 matrix(),                                    */
/*                                 oper1(),                                     */
/*                                 oper2(),                                     */
/*                                 oper3(),                                     */
/*                                 oper4(),                                     */
/*                                 oper5(),                                     */
/*                                 oper6(),                                     */
/*                                 print_population(),                          */
/*                                 sort(),                                      */
/*                                 Gvector() was vector().                                    */
/*                                                                              */
/*           CALLING FUNCITONS :   main()                                       */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

/*
		    int nvars, int pop_size, int MaxGenerations, int WaitGenerations,
		    int P1, int P2, int P3, int P4, int P5, int P6, int P7, int P8, 
		    FILE *output, short int MinMax, short int GradientCheck, short int BoundaryEnforcement, 
		    double SolutionTolerance, double *Results, double *Gradients,
		    int *Status
*/

double optimization(struct GND_IOstructure *Structure, VECTOR X, 
		    MATRIX domains, FILE *output)
{

  extern struct GND_IOstructure *ExternStructure;

  MATRIX 
    new_genera,   /*Temporary storage for the new generation*/
    population,   /*Population of x2 variables*/
    temp;

  VECTOR probab,       /*Probability of agents to die or live*/
         cum_probab,   /*Cumilative probability of agents*/
         t_vec;

  IVECTOR live;


  long count_gener= 1; /*Counter to keep track of the number of generations*/
  unsigned long peak_cnt;

  int                     /*Total number of agents chosen to reproduce*/
    j1,
    j2,
    j3,
    j4,
    j5,
    j6,
    j7,
    j8,
    oper,
    ocnt,
    B,                     /*Parameter for the 3rd operator - nonuniform mutation*/
    STEP,                  /*Parameter for the 5th operator - simple arithmetical crossover*/
    first_live=0,            /*Index of the two parents for crossover parents*/
    second_live=0,
    first_die,             /*Index of the two parents for crossover death*/
    second_die,
    die_now,               /*index of agent to replace in current operation*/
    i,
    j,
    s;


  double Q,                   /*Probability of the best agent*/
         Teval=0,               /*Evaluation of the best agent*/
         peak_val;

  FLAG  same;
  double **Jnew; 

  int evaliter;
  double bfgsfit, evalgtol;

  double *grad, *evalX, *finalhessin, *bfgsoutX;

  int nochange_gen=0;

  double oldfitvalue=0;

  int IncreaseGenerations;
  short int GradientTrigger=0;
  short int BoundaryTrigger;
  long InstanceNumber;

  /* Strucutre fixup! */
  long nvars, MaxGenerations, WaitGenerations;
  long *Status;
  long pop_size, P, P0, P1, P2, P3, P4, P5, P6, P7, P8;
  short int MinMax, GradientCheck, BoundaryEnforcement, UseBFGS;
  double SolutionTolerance, *Results, *Gradients;
  short PrintLevel, HardGenerationLimit;

  /* Old variables which may change when SetRunTimeParameters is run during a run! */
  long pop_size_old;
  double **population_old;

  /* Summary Statistics (mean, variance etc) */
  double popmean, popvar, popwrk, popstat;

  /* Population Print population*/
  FILE *popout;
  long *tobs;
  double *mean, *var, *skew, *kur;

  /* Stuff for the Unique Stuff (how's that for an informative comment! */
  /* A big Matrix which remembers all of our past evaluations. It's
     maximum memory is set in genoud.h */
  extern long Gnvars[MAXINSTANCES];
  double **Memory;
  long MemorySize=0, UniqueCount, OldUniqueCount=0;

  /* LVM calls */
  long LVMreturn;

  LVMreturn = 0;
  /* fine two unique parents count */
  long SameCount, UniquePairs;

  // NetworkEvaluate() Stuff
  long NetworkNumber; // number of individuals to be evaluated

  ExternStructure=Structure;
  Status=&(Structure->Status);

  Results=Structure->oResults;
  Gradients=Structure->oGradients;

  /* Structure Done */
  SetRunTimeParameters(Structure, 1,
		       &pop_size, &nvars, &MaxGenerations, &WaitGenerations,
		       &MinMax, &GradientCheck, &BoundaryEnforcement, &UseBFGS, &SolutionTolerance,
		       &InstanceNumber, &P, &P0, &P1, &P2, &P3, &P4, &P5, &P6, &P7, &P8, 
		       &PrintLevel, &HardGenerationLimit, output);

  /*Space allocation for all the vectors and matrices involved*/
  population    = JaMatrixAllocate(pop_size+2, nvars+2);
  new_genera    = JaMatrixAllocate(pop_size+2, nvars+2);

  /* new_genera = JaMatrixAllocate(pop_size+2, nvars+2); */
  temp       = matrix(1,2,0,nvars);
  probab     = Gvector(1,pop_size);
  t_vec      = Gvector(1,nvars);
  cum_probab = Gvector(1,pop_size);
  live       = ivector(1,pop_size);

  Gnvars[Structure->InstanceNumber]=nvars;

  if (Structure->MemoryUsage==1)
    {
      if (HardGenerationLimit==0)
	MemorySize=(MaxGenerations+1)*pop_size+1+pop_size;
      else
	MemorySize=(MaxGenerations+1)*pop_size+1+pop_size;
      
      Memory = JaMatrixAllocate(MemorySize, nvars+2);
    }

  grad = (double *) malloc((nvars)*sizeof(double));
  evalX = (double *) malloc((nvars)*sizeof(double));
  finalhessin = (double *) malloc(((nvars*nvars)+(nvars))*sizeof(double));
  bfgsoutX = (double *) malloc((nvars+1)*sizeof(double));

  /* populationstats variables */
  mean = (double *) malloc((nvars+1)*sizeof(double));
  var = (double *) malloc((nvars+1)*sizeof(double));
  skew = (double *) malloc((nvars+1)*sizeof(double));
  kur = (double *) malloc((nvars+1)*sizeof(double));
  tobs = (long *) malloc((nvars+1)*sizeof(long));

  Q=0.5;
  B=6;
  STEP=10;

  fprintf(output,"\n\n");

  switch(MinMax) {
  case 0:
    fprintf(output,"Minimization Problem.\n\n");  
    break;
  case 1:
    fprintf(output,"Maximization Problem.\n\n");  
    break;
  }

  fprintf(output,"Parameter B (hardcoded): %d\n", B);
  fprintf(output,"Parameter Q (hardcoded): %f\n", Q);
  fprintf(output,"\n");

  fflush(output);

  peak_val = 0;
  peak_cnt = 0;

  pop_size_old=0;
  if (Structure->ShareType == 1 || Structure->ShareType == 3) {
    fprintf(output, "Using old population file to initialize new population\n");
    if((popout = fopen(Structure->ProjectPath, "r")) == NULL) {
      fprintf(output,"WARNING: Unable to open the old project file: %s\n", 
	      Structure->ProjectPath);
      fprintf(output,"         Generating new population\n");
    }
    else {
      pop_size_old=ReadPopulation(population, pop_size, nvars, output, popout);
      fclose(popout);
      if (pop_size_old<2) {
	fprintf(output,
		"WARNING: The old population file appears to be from the run of a different model!\n");
	pop_size_old=0;
      }
    }
    if (PrintLevel==2) {
      if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
	fprintf(output,"Unable to open the project file: %s", 
		Structure->ProjectPath);
	
	/* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);
	
	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	
        /* free numeric.c allocations */
	if (Structure->MemoryUsage==1)
	  JaMatrixFree(Memory, MemorySize);

	JaMatrixFree(population, pop_size+2);
	JaMatrixFree(new_genera, pop_size+2);
	
	free_matrix(temp, 1, 2, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);

	return(ERROR_CODE);
      }
      fclose(popout);
    }
  } /* end of ShareType 0 */
  else {
    if (PrintLevel==2) {
      if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
	fprintf(output,"Unable to open the project file: %s", 
		Structure->ProjectPath);

	/* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);
	
	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	
        /* free numeric.c allocations */
	if (Structure->MemoryUsage==1)
	  JaMatrixFree(Memory, MemorySize);

	JaMatrixFree(population, pop_size+2);
	JaMatrixFree(new_genera, pop_size+2);
	
	free_matrix(temp, 1, 2, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);

	return(ERROR_CODE);
      }
      fclose(popout);
    }
  }

  /* The new initial value matrix: setting a new initial value for every individual */
  if (ExternStructure->nStartingValues > 0) 
    {
      fprintf(output,"\nSTARTING VALUES\n\n");
      // seed the starting values until we run out of population or starting values!
      j = pop_size_old;
      for(s=0; s<ExternStructure->nStartingValues; s++) {
	j++;
	for(i=1; i<=nvars; i++) {
	  population[j][i] = ExternStructure->StartingValues[s][i-1];
	  population[j][nvars+1] = -1.0;
	}
      } // end of for loop
      pop_size_old = j;

      // randomly add on people if we still have population left over!
      for(j=pop_size_old+1; j<=pop_size; j++) {
	for(i=1; i<=nvars; i++) { 
	  population[j][i] = frange_ran(domains[i][1], domains[i][3]); 
	  population[j][nvars+1] = -1.0;
	}
      }
    } // end of we have starting values!
  else 
    {
      for(j=pop_size_old+1; j<=pop_size; j++) {
	for(i=1; i<=nvars; i++) { 
	  population[j][i] = frange_ran(domains[i][1], domains[i][3]); 
	  population[j][nvars+1] = -1.0;
	}
      }
    } // end of else

  if (Structure->MemoryUsage==1)
    {
      /* BINARY SEARCH.  (Knuth 3:409). */
      /* We have sorted Memory by the "key" */
      OldUniqueCount=UniqueCount=0;
      
      /* BINARY SEARCH.  (Knuth 3:409). */
      /* We have sorted Memory by the "key" */
      /* JaIntegerSort(population, pop_size, nvars+2);  */
      JaDoubleSort(population, pop_size, nvars+2); 
      
      OldUniqueCount=UniqueCount;

      JaDoubleMemoryMatrix_Gen0(Structure,
				Memory, population, X,
				&UniqueCount, OldUniqueCount, pop_size, nvars, 
				output, Status);

      if (*Status < 0)
	{
	  fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);

	  /* free memory */
	  if (Structure->MemoryUsage==1)
	    JaMatrixFree(Memory, MemorySize);
	      
	      /* free populationstats stuff */
	  free(mean);
	  free(var);
	  free(skew);
	  free(kur);
	  free(tobs);
	      
	  free(bfgsoutX);
	  free(finalhessin);
	  free(evalX);
	  free(grad);
	      
	  /* free numeric.c allocations */
	  JaMatrixFree(population, pop_size+2);
	  JaMatrixFree(new_genera, pop_size+2);
	      
	  free_matrix(temp, 1, 2, 0);
	  free_vector(probab, 1);
	  free_vector(t_vec, 1);
	  free_vector(cum_probab, 1);
	  free_ivector(live, 1);
	      
	  return(ERROR_CODE);
	}
      if ( (UniqueCount+pop_size) >= MemorySize )
	{
	  Structure->MemoryUsage=0;
	  fprintf(output,"\nWARNING: Turning Off MemoryMatrix because memory usage is too great.\n\n");
	} /* end of if */
    } // end of Memory based evaluation
  else
    {
      NetworkNumber=0;
      for (i=1; i<=pop_size; i++) 
	{
	  if (Structure->DynamicPopulation==2)
	    {
	      JaDynamicPopulationCheck(Structure, population, i, pop_size, nvars, output, Status);
	      if (*Status < 0) {
		fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
		
		/* free memory */
		
		/* free populationstats stuff */
		free(mean);
		free(var);
		free(skew);
		free(kur);
		free(tobs);
		
		free(bfgsoutX);
		free(finalhessin);
		free(evalX);
		free(grad);
		
		/* free numeric.c allocations */
		JaMatrixFree(population, pop_size+2);
		JaMatrixFree(new_genera, pop_size+2);
		
		free_matrix(temp, 1, 2, 0);
		free_vector(probab, 1);
		free_vector(t_vec, 1);
		free_vector(cum_probab, 1);
		free_ivector(live, 1);
		
		return(ERROR_CODE);		
	      } // end of Status < 0
	    } // end of DynamicPopulation==2

	  if (population[i][nvars+1]==-1.0 || population[i][nvars+1]==11.0)
	    {
	      
	      if (Structure->Network==1)
		{
		  NetworkNumber++;
		  population[i][0] = EVALUATE;
		}
	      else
		{
		  for(j=1; j<=nvars; j++)
		    X[j] = population[i][j];
		  
		  population[i][0] = evaluate(Structure->AgentFit, X, nvars, Status);
		  
		  if (*Status < 0) {
		    fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
		    
		    /* free memory */
		    
		    /* free populationstats stuff */
		    free(mean);
		    free(var);
		    free(skew);
		    free(kur);
		    free(tobs);
		    
		    free(bfgsoutX);
		    free(finalhessin);
		    free(evalX);
		    free(grad);
		    
		    /* free numeric.c allocations */
		    JaMatrixFree(population, pop_size+2);
		    JaMatrixFree(new_genera, pop_size+2);
		    
		    free_matrix(temp, 1, 2, 0);
		    free_vector(probab, 1);
		    free_vector(t_vec, 1);
		    free_vector(cum_probab, 1);
		    free_ivector(live, 1);
		    
		    return(ERROR_CODE);
		  }
		} // else
	    }
	} //end of i loop
      if (Structure->Network==1)
	{
	  NetworkEvaluate(Structure->AgentFit, Structure->DBname, Structure->AgentName, 
			  population, pop_size, nvars, NetworkNumber, Status, 10);
	  
	  if (*Status < 0) {
	    fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
	    
	    /* free memory */
	    
	    /* free populationstats stuff */
	    free(mean);
	    free(var);
	    free(skew);
	    free(kur);
	    free(tobs);
	    
	    free(bfgsoutX);
	    free(finalhessin);
	    free(evalX);
	    free(grad);
	    
	    /* free numeric.c allocations */
	    JaMatrixFree(population, pop_size+2);
	    JaMatrixFree(new_genera, pop_size+2);
	    
	    free_matrix(temp, 1, 2, 0);
	    free_vector(probab, 1);
	    free_vector(t_vec, 1);
	    free_vector(cum_probab, 1);
	    free_ivector(live, 1);
	    
	    return(ERROR_CODE);
	  }	      
	} 
    } // end of default evaluation


  /*Sort the initial inidivduals based on their evaluation function*/
  sort(MinMax,population,pop_size,0);

  switch(MinMax) {
  case 0:
    Teval = population[1][0];
    peak_cnt = count_gener;
    peak_val = population[1][0];
    break;
  case 1:
    Teval = population[1][0];
    peak_cnt = count_gener;
    peak_val = population[1][0];
    break;
  }

  fprintf(output,"\nThe 2 best initial individuals are\n");
  for(i=1; i<3; i++) {
    print_vector(population[i],1,nvars,output);
    fprintf(output,"\nfitness = %e", population[i][0]);
    fprintf(output,"\n\n");
  }

  fprintf(output,"\nThe worst fit of the population is: %e\n", 
      population[pop_size][0]);
  fprintf(output,"\n\n");

  fprintf(output,"\n\nGeneration#\t    Solution Value\n");
  fprintf(output,"\n%7d \t%e\n", 0, population[1][0]);

  /* compute and print mean and variance of population */
  if (PrintLevel==2) {
      fprintf(output,"GENERATION: 0 (initializing the population)\n");
      populationstats(population, pop_size, nvars, mean, var, skew, kur, tobs);
      for (i=0; i<=nvars; i++) {
	  if (i==0) {
	      fprintf(output, "Fitness Value... %e\n", population[1][i]);
	      fprintf(output, "mean............ %e\n", mean[i]);
	      fprintf(output, "var............. %e\n", var[i]);
	      fprintf(output, "skewness........ %e\n", skew[i]);
	      fprintf(output, "kurtosis........ %e\n", kur[i]);
	      fprintf(output, "#null........... %d\n", pop_size-tobs[i]);
	      if(Structure->MemoryUsage==1)
		fprintf(output, "#unique......... %d, #Total UniqueCount: %d\n", 
			UniqueCount-OldUniqueCount, UniqueCount);
	      fprintf(output, "tobs............ %d\n", tobs[i]);
	  }
	  else {
	      fprintf(output, "var %d:\n", i);
	      fprintf(output, "best............ %e\n", population[1][i]);
	      fprintf(output, "mean............ %e\n", mean[i]);
	      fprintf(output, "var............. %e\n", var[i]);
	      fprintf(output, "skewness........ %e\n", skew[i]);
	      fprintf(output, "kurtosis........ %e\n", kur[i]);
	      fprintf(output, "#null........... %d\n", pop_size-tobs[i]);
	      fprintf(output, "tobs............ %d\n", tobs[i]);
	  }
      }
  } /* end of printlevel if */

  if (PrintLevel==1) {
    popmean = popvar = 0.0 ;
    popwrk = 1.0 / pop_size ;
    for(i=1; i<=pop_size; i++) {
      popmean += population[i][0] ;
    }
    popmean *= popwrk ; 
    for(i=1; i<=pop_size; i++) {
      popstat =  population[i][0] - popmean ;
      popvar += (popstat*popstat) ;
    }
    popvar *= popwrk ;
    fprintf(output, "   mean = %e, variance = %e\n\n", popmean, popvar);
  }

  fflush(output);

  /* Print the population file */
  if ( PrintLevel == 1 ) {
    if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
      fprintf(output,"Unable to open the project file: %s", 
	      Structure->ProjectPath);

      /* free populationstats stuff */
      free(mean);
      free(var);
      free(skew);
      free(kur);
      free(tobs);
      
      free(bfgsoutX);
      free(finalhessin);
      free(evalX);
      free(grad);
      
      /* free numeric.c allocations */
      if (Structure->MemoryUsage==1)
	JaMatrixFree(Memory, MemorySize);

      JaMatrixFree(population, pop_size+2);
      JaMatrixFree(new_genera, pop_size+2);
      
      free_matrix(temp, 1, 2, 0);
      free_vector(probab, 1);
      free_vector(t_vec, 1);
      free_vector(cum_probab, 1);
      free_ivector(live, 1);
			  
      return(ERROR_CODE);
    }
    print_population(pop_size, nvars, 0, population, popout);
    fclose(popout);
  } /* end of PrintLevel if */
  if ( PrintLevel == 2 ) {
    if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
      fprintf(output,"Unable to open the project file: %s", 
	      Structure->ProjectPath);

      /* free populationstats stuff */
      free(mean);
      free(var);
      free(skew);
      free(kur);
      free(tobs);
      
      free(bfgsoutX);
      free(finalhessin);
      free(evalX);
      free(grad);
      
      /* free numeric.c allocations */
      if (Structure->MemoryUsage==1)
	JaMatrixFree(Memory, MemorySize);

      JaMatrixFree(population, pop_size+2);
      JaMatrixFree(new_genera, pop_size+2);
      
      free_matrix(temp, 1, 2, 0);
      free_vector(probab, 1);
      free_vector(t_vec, 1);
      free_vector(cum_probab, 1);
      free_ivector(live, 1);

      return(ERROR_CODE);
    }
    print_population(pop_size, nvars, 0, population, popout);
    fflush(popout);
    fclose(popout);
  }

  /*Assigning probability of survival for each of the agent, with the*/
  /*probability provided by the user for the best agent*/
  assign_probab(probab,pop_size,Q); 

  /*Finding the cumulative probability of the agents*/
  find_cum_probab(cum_probab,probab,pop_size);

  /*Reproducing and evaluating for the total number of generations times*/
  do
    {

      /*Initializing the live vector*/
      for(j=1; j<=pop_size; j++)
        {
          live[j] = 0;
          for(i=0; i<=nvars; i++)
            new_genera[j][i] = population[j][i];
        }

      /*Finding the agents that will die and the agents that will reproduce*/
      find_live(cum_probab,live,pop_size,P);
      /* set die_now counter to start replacements with the worst agent.
         use of die_now is okay if the entire population (except the previous
         best) is to be replaced in each generation */
      die_now = pop_size;

      j1=j2=j3=j4=j5=j6=j7=j8=0;

      UniquePairs= UniqueCount-OldUniqueCount;
      UniquePairs= (int) (0.5*(UniquePairs*UniquePairs-UniquePairs));
      if ( MAX_OPER_UNIQUE_TRY < UniquePairs)
	UniquePairs = MAX_OPER_UNIQUE_TRY;

      /* main operator loop */
      while(j1+j2+j3+j4+j4+j5+j5+j6+j7+j7+j8 < P)
        {
          oper = irange_ran(1,8);
          switch (oper)
            {
              case 1:
		/* JS Description: Uniform Mutation */
                     /*Applying the first operator, uniform mutation*/
                    if (j1 < P1)
                      {
			/*Find one parent for mutation operator 1*/
			first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  fprintf(output, "No agents to be replaced\n");
			  exit(1);
			}

			new_genera[die_now][nvars+1] = 1.0;
			for(i=1; i<=nvars; i++)
			  t_vec[i] = population[first_live][i];
			for (ocnt = irange_ran(1,nvars); ocnt>0; ocnt--)
			  oper1(t_vec,domains,nvars);
			for(i=1; i<=nvars; i++)
			  new_genera[die_now][i] = t_vec[i];
			die_now--;
			j1++;
		      }
                    break;
              case 2:
		/* JS Description: Boundary Mutation */
                    /*Applying the second operator, boundary mutation*/
                    if (j2 < P2)
                      {
                        /*Find one parent for mutation operator 2*/
                        first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  fprintf(output, "No agents to be replaced\n");
			  exit(1);
			}

                        new_genera[die_now][nvars+1] = 2.0;
                        for(i=1; i<=nvars; i++)
                          t_vec[i] = population[first_live][i];
			oper2(t_vec,domains,nvars);
                        for(i=1; i<=nvars; i++)
                          new_genera[die_now][i] = t_vec[i];
			die_now--;
                        j2++;
                      }
                    break;
              case 3:
		/* JS Description: Non-uniform Mutation */
                    /*Applying the third operator, non-uniform mutation*/
                    if (j3 < P3)
                      {
                        /*Find one parent for mutation operator 3*/
                        first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  fprintf(output, "No agents to be replaced\n");
			  exit(1);
			}

                        new_genera[die_now][nvars+1] = 3.0;
                        for(i=1; i<=nvars; i++)
                          t_vec[i] = population[first_live][i];
			for (ocnt = irange_ran(1,nvars); ocnt>0; ocnt--)
			  oper3(t_vec,domains,nvars,MaxGenerations,count_gener,B);
                        for(i=1; i<=nvars; i++)
                          new_genera[die_now][i] = t_vec[i];
			die_now--;
                        j3++;
                      }
                    break;

              case 4:
		/* JS Description: Polytope Crossover */
                    /*Applying the fourth operator, whole arithmetical crossover*/
                    if (j4 < (int) P4/2)
                      {
                        /*Find two distinct parents for crossover operator 4*/
                        same = TRUE;
			SameCount=0;
			while (same==TRUE) {
			  SameCount++;

			  first_live  = find_parent(live,pop_size);
			  second_live = find_parent(live,pop_size);

			  if (SameCount >= (UniquePairs) ) 
			    break;
			  
			  for(i=1; i<=nvars; i++)
			    if (population[first_live][i] != population[second_live][i])
			      same = FALSE;
			} /* end of while same==TRUE loop */
			/* check that agents to replace are in range */
			if (die_now < 3) {
			  fprintf(output,"Not enough agents to be replaced\n");
			  exit(1);
			}
			live[first_live]--;
			live[second_live]--;
			first_die   = die_now-- ;
			second_die  = die_now-- ;
			new_genera[first_die][nvars+1]  = 4.0;
			new_genera[second_die][nvars+1] = 4.0;
                        if (!same)
                          {
                            for(i=1; i<=nvars; i++)
                              {
                                temp[1][i] = population[first_live][i];
                                temp[2][i] = population[second_live][i];
                              }
                            oper4(temp[1],temp[2],nvars);
                            for(i=1; i<=nvars; i++)
                              {
                                new_genera[first_die][i]  = temp[1][i];
                                new_genera[second_die][i] = temp[2][i];
			      }
                          }
			else {
			  /* copy agent chosen twice into two new indivs */
			  for(i=1; i<=nvars; i++) {
			    new_genera[first_die][i]  = 
			      population[first_live][i];
			    new_genera[second_die][i] = 
			      population[second_live][i];
			  }
			}
                        j4++;
                      }
                    break;
              case 5:
		/* JS Description: Multiple Point Simple Crossover
		   Applying the fifth operator, simple arithmetical crossover*/
                    if (j5 < (int) P5/2)
                      {
                        /*Find two distinct parents for crossover operator 5*/
                        same = TRUE;
			SameCount=0;
			while (same==TRUE) {
			  SameCount++;
			  
			  first_live  = find_parent(live,pop_size);
			  second_live = find_parent(live,pop_size);

			  if (SameCount >= (UniquePairs) ) 
			    break;

			  for(i=1; i<=nvars; i++)
			    if (population[first_live][i] != population[second_live][i])
			      same = FALSE;
			} /* end of while same==TRUE loop */
			/* check that agents to replace are in range */
			if (die_now < 3) {
			  fprintf(output,"Not enough agents to be replaced\n");
			  exit(1);
			}
			live[first_live]--;
			live[second_live]--;
			first_die   = die_now-- ;
			second_die  = die_now-- ;
			new_genera[first_die][nvars+1]  = 5.0;
			new_genera[second_die][nvars+1] = 5.0;
                        if (!same)
                          {
                            for(i=1; i<=nvars; i++)
                              {
                                temp[1][i] = population[first_live][i];
                                temp[2][i] = population[second_live][i];
                              }
                            oper5(temp[1],temp[2],STEP,domains,nvars);
                            for(i=1; i<=nvars; i++)
                              {
                                new_genera[first_die][i]  = temp[1][i];
                                new_genera[second_die][i] = temp[2][i];
                              }
                          }
			else {
			  /* copy agent chosen twice into two new indivs */
			  for(i=1; i<=nvars; i++) {
			    new_genera[first_die][i]  = 
			      population[first_live][i];
			    new_genera[second_die][i] = 
			      population[second_live][i];
			  }
			}
                        j5++;
                      }
                    break;
              case 6:
		/* JS Description: Whole Non-uniform Mutation */
                    /*Applying the sixth operator, whole non-uniform mutation*/
                    if (j6 < P6)
                      {
                        /*Find one parent for mutation operator 6*/
                        first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  fprintf(output, "No agents to be replaced\n");
			  exit(1);
			}

                        new_genera[die_now][nvars+1] = 6.0;
                        for(i=1; i<=nvars; i++)
                          t_vec[i] = population[first_live][i];
                        oper6(t_vec,domains,nvars,MaxGenerations,count_gener,B);
                        for(i=1; i<=nvars; i++)
                          new_genera[die_now][i] = t_vec[i];
			die_now--;
                        j6++;
                      }
                    break;
              case 7:
		/* JS Description: Heuristic Crossover */
                    /*Applying the seventh operator*/
                    if (j7 < (int) P7/2)
                      {
                        /*Find two distinct parents for operator 7*/
                        same = TRUE;
			SameCount=0;
			while (same==TRUE) {
			  SameCount++;

			  first_live  = find_parent(live,pop_size);
			  second_live = find_parent(live,pop_size);
			  
			  if (SameCount >= (UniquePairs) ) 
			    break;

			  for(i=1; i<=nvars; i++)
			    if (population[first_live][i] != population[second_live][i])
			      same = FALSE;
			} /* end of while same==TRUE loop */
			/* check that agents to replace are in range */
			if (die_now < 3) {
			  fprintf(output,"Not enough agents to be replaced\n");
			  exit(1);
			}
			live[first_live]--;
			live[second_live]--;
			first_die   = die_now-- ;
			second_die  = die_now-- ;
			new_genera[first_die][nvars+1]  = 7.0;
			new_genera[second_die][nvars+1] = 7.0;
                        if (!same) {
			  if (first_live < second_live)
			    /* first agent is better agent */
			    for(i=1; i<=nvars; i++) {
			      temp[2][i] = population[first_live][i];
			      temp[1][i] = population[second_live][i];
			    }
			  else
			    /* second agent is better agent */
			    for(i=1; i<=nvars; i++) {
			      temp[2][i] = population[second_live][i];
			      temp[1][i] = population[first_live][i];
			    }
			  oper7(temp[1],temp[2],domains,nvars);
			  for(i=1; i<=nvars; i++)
			    new_genera[first_die][i]  = temp[1][i];
			  if (first_live < second_live)
			    /* first agent is better agent */
			    for(i=1; i<=nvars; i++) {
			      temp[2][i] = population[first_live][i];
			      temp[1][i] = population[second_live][i];
			    }
			  else
			    /* second agent is better agent */
			    for(i=1; i<=nvars; i++) {
			      temp[2][i] = population[second_live][i];
			      temp[1][i] = population[first_live][i];
			    }
			  oper7(temp[1],temp[2],domains,nvars);
			  for(i=1; i<=nvars; i++)
			    new_genera[second_die][i]  = temp[1][i];
			}
			else {
			  /* copy agent chosen twice into two new indivs */
			  for(i=1; i<=nvars; i++) {
			    new_genera[first_die][i]  = 
			      population[first_live][i];
			    new_genera[second_die][i] = 
			      population[second_live][i];
			  }
			}
                        j7++;
                      }
              case 8:
		/* JS Description: Local-Minimum Crossover */
                     /*Applying the eighth operator, homotopy (BFGS) */
                    if (j8 < P8)
                      {
                        /*Find one parent for BFGS operator 1*/
                        first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  fprintf(output, "No agents to be replaced\n");
			  exit(1);
			}

                        new_genera[die_now][nvars+1] = 8.0;
                        for(i=1; i<=nvars; i++)
                          t_vec[i] = population[first_live][i];
                        oper8(Structure->AgentFit, t_vec, domains, SolutionTolerance, 
			      Structure->Optim, nvars, 
			      BoundaryEnforcement, MinMax, InstanceNumber, output, 
			      Status, PrintLevel);
			if (*Status < 0) {
			  fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
			  
			  /* free populationstats stuff */
			  free(mean);
			  free(var);
			  free(skew);
			  free(kur);
			  free(tobs);

			  free(bfgsoutX);
			  free(finalhessin);
			  free(evalX);
			  free(grad);
			  
				/* free numeric.c allocations */
			  if (Structure->MemoryUsage==1)
			    JaMatrixFree(Memory, MemorySize);

			  JaMatrixFree(population, pop_size+2);
			  JaMatrixFree(new_genera, pop_size+2);

			  free_matrix(temp, 1, 2, 0);
			  free_vector(probab, 1);
			  free_vector(t_vec, 1);
			  free_vector(cum_probab, 1);
			  free_ivector(live, 1);
			  
			  return(ERROR_CODE);
			}
			for(i=1; i<=nvars; i++)
			  new_genera[die_now][i] = t_vec[i];
			die_now--;
			j8++;
                      }
                    break;
            }
        }
      
      /*Replace the population with the new generation */
      Jnew = new_genera;
      new_genera = population;
      population = Jnew;

      if (Structure->DynamicPopulation==1 || Structure->DynamicPopulation==2)
	{
	  FILE *DynamicInput;
	  double tmp;

	  if((DynamicInput = fopen(Structure->DynamicPopulationPath, "r")) == NULL) {
	    fprintf(output,"WARNING: Unable to open the DynamicPopulationPath: %s\n", 
		    Structure->DynamicPopulationPath);
	    fprintf(output,"         Continuing the process.\n");

	    Structure->DynamicPopulation=0;
	    break;
	  }
	  
	  fprintf(output,"\nDynamically Reading in Individuals from: %s\n", Structure->DynamicPopulationPath);

	  j = pop_size;
	  i=0;
	  while (fscanf(DynamicInput,"%lf", &tmp)==1)
	    {
	      i++;
	      if (i>nvars)
		{
		  j--;
		  if (j<2)
		    {
		      fprintf(output,
			      "\nWARNING: Dynamic Population input includes more individuals than the population size.\n");
		      fprintf(output,
			      "         Keeping the current best individual and discarding extra individuals.\n");
		      fprintf(output,
			      "         See file: %s\n\n", Structure->DynamicPopulationPath);
		      break;
		    }
		  i = 1;
		} // end of if

	      population[j][i] = tmp;
	      population[j][nvars+1] = 10.0;

	      if (Structure->Debug==1)
		{
		  printf("\nDEBUG: WE READ IN THE FOLLOWING INDIVIDUALS:\n");
		  printf("population[%d][%d]: %e\n", j, i, population[j][i]);
		}

	    } // end of while loop
	  fclose(DynamicInput);

	  fprintf(output,"   Read in %d Individuals\n\n", pop_size-j+1);
	  
	  Structure->DynamicPopulation=0;
	} // end of DynamicPopulation


      if (Structure->MemoryUsage==1)
	{
	  /* BINARY SEARCH.  (Knuth 3:409). */
	  /* We have sorted Memory by the "key" */
	  /* JaIntegerSort(population, pop_size, nvars+2);  */
	  JaDoubleSort(population, pop_size, nvars+2); 
	  JaDoubleSort(Memory, UniqueCount, nvars+2);
	  
	  OldUniqueCount=UniqueCount;

	  JaDoubleMemoryMatrix(Structure,
			       Memory, population, X,
			       &UniqueCount, OldUniqueCount,
			       pop_size, nvars, output, Status);

	  if (*Status < 0) {
	    fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
		
	    /* free memory */
	    if (Structure->MemoryUsage==1)
	      JaMatrixFree(Memory, MemorySize);
		
	    /* free populationstats stuff */
	    free(mean);
	    free(var);
	    free(skew);
	    free(kur);
	    free(tobs);
		
	    free(bfgsoutX);
	    free(finalhessin);
	    free(evalX);
	    free(grad);
		
	    /* free numeric.c allocations */
	    JaMatrixFree(population, pop_size+2);
	    JaMatrixFree(new_genera, pop_size+2);
		
	    free_matrix(temp, 1, 2, 0);
	    free_vector(probab, 1);
	    free_vector(t_vec, 1);
	    free_vector(cum_probab, 1);
	    free_ivector(live, 1);
		
	    return(ERROR_CODE);
	  }	  
	  if ( (UniqueCount+pop_size) >= MemorySize )
	    {
	      Structure->MemoryUsage=0;
	      fprintf(output,"\nWARNING: Turning Off MemoryMatrix because memory usage is too great.\n\n");
	    } /* end of if */
	} // end of MemoryUsage==1
      else
	{
	  NetworkNumber=0;
	  for (i=1; i<=pop_size; i++) 
	    {
	      if (i > 1 && Structure->DynamicPopulation==2)
		{
		  JaDynamicPopulationCheck(Structure, population, i, pop_size, nvars, output, Status);
		  if (*Status < 0) {
		    fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
		
		    /* free memory */
		
		    /* free populationstats stuff */
		    free(mean);
		    free(var);
		    free(skew);
		    free(kur);
		    free(tobs);
		
		    free(bfgsoutX);
		    free(finalhessin);
		    free(evalX);
		    free(grad);
		
		    /* free numeric.c allocations */
		    JaMatrixFree(population, pop_size+2);
		    JaMatrixFree(new_genera, pop_size+2);
		
		    free_matrix(temp, 1, 2, 0);
		    free_vector(probab, 1);
		    free_vector(t_vec, 1);
		    free_vector(cum_probab, 1);
		    free_ivector(live, 1);
		
		    return(ERROR_CODE);		
		  } // end of Status < 0
		} // end of DynamicPopulation==2
	      
	      if (population[i][nvars+1]!=0)
		{
		  if (Structure->Network==1)
		    {
		      NetworkNumber++;
		      population[i][0] = EVALUATE;
		    }
		  else
		    {
		      for(j=1; j<=nvars; j++)
			X[j] = population[i][j];
		      
		      population[i][0] = evaluate(Structure->AgentFit, X, nvars, Status);
		      
		      if (*Status < 0) {
			fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
			
			/* free memory */
			
			/* free populationstats stuff */
			free(mean);
			free(var);
			free(skew);
			free(kur);
			free(tobs);
			
			free(bfgsoutX);
			free(finalhessin);
			free(evalX);
			free(grad);
			
			/* free numeric.c allocations */
			JaMatrixFree(population, pop_size+2);
			JaMatrixFree(new_genera, pop_size+2);
			
			free_matrix(temp, 1, 2, 0);
			free_vector(probab, 1);
			free_vector(t_vec, 1);
			free_vector(cum_probab, 1);
			free_ivector(live, 1);
			
			return(ERROR_CODE);
		      }
		    }
		}
	    } //end of i loop	  
	  if (Structure->Network==1)
	    {
	      NetworkEvaluate(Structure->AgentFit, Structure->DBname, Structure->AgentName, 
			      population, pop_size, nvars, NetworkNumber, Status, 10);
	      
	      if (*Status < 0) {
		fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
		
		/* free memory */
		
		/* free populationstats stuff */
		free(mean);
		free(var);
		free(skew);
		free(kur);
		free(tobs);
		
		free(bfgsoutX);
		free(finalhessin);
		free(evalX);
		free(grad);
		
		/* free numeric.c allocations */
		JaMatrixFree(population, pop_size+2);
		JaMatrixFree(new_genera, pop_size+2);
		
		free_matrix(temp, 1, 2, 0);
		free_vector(probab, 1);
		free_vector(t_vec, 1);
		free_vector(cum_probab, 1);
		free_ivector(live, 1);
		
		return(ERROR_CODE);
	      }	    
	    }  
	} //end of default evaluation scheme

      /*Sort the new population based on their evaluation function*/
      sort(MinMax,population,pop_size,0);

      /* apply the bfgs to the best individual */
      if (UseBFGS != 0) {
	for (i=1; i<=nvars; i++)
	  {
	    bfgsoutX[i-1]=population[1][i];
	  }

	if (Structure->Optim==0)
	  {
	    evalgtol=SolutionTolerance;

	    dfgsmin(Structure->AgentFit, bfgsoutX, nvars, evalgtol, &evaliter, &bfgsfit, 
		    finalhessin, MinMax, BoundaryEnforcement, InstanceNumber, domains, 
		    Status, PrintLevel, output);

	    if (*Status < 0) {
	      /* Free Memory */
	    
	      /* free populationstats stuff */
	      free(mean);
	      free(var);
	      free(skew);
	      free(kur);
	      free(tobs);
	    
	      free(bfgsoutX);
	      free(finalhessin);
	      free(evalX);
	      free(grad);
	    
	      /* free numeric.c allocations */
	      if (Structure->MemoryUsage==1)
		JaMatrixFree(Memory, MemorySize);

	      JaMatrixFree(population, pop_size+2);
	      JaMatrixFree(new_genera, pop_size+2);
	    
	      free_matrix(temp, 1, 2, 0);
	      free_vector(probab, 1);
	      free_vector(t_vec, 1);
	      free_vector(cum_probab, 1);
	      free_ivector(live, 1);
	    
	      return(ERROR_CODE);
	    }
	
	    if (MinMax==1) bfgsfit=-1*bfgsfit;
	  } // Optim==0
	else
	  {
	    bfgsfit = genoud_optim(bfgsoutX, nvars);
	  }
	
	switch(MinMax) {
	case 0:
	  if (population[1][0] > bfgsfit) /* minimize */
	    {
	      /* is the BFGS individual in the bounds? */
	      BoundaryTrigger=0; /* outside of bounds ? */
	      for (i=0; i<nvars; i++) {
		j = i+1;
		if (bfgsoutX[i] < domains[j][1]) {
		  BoundaryTrigger=1;
		  fprintf(output,
			  "\nWARNING: BFGS hit on best individual produced Out of Boundary individual.\n");
		  fprintf(output,"WARNING: Generation: %d \t Parameter: %d \t Value: %e\n\n", 
			  count_gener, i+1, bfgsoutX[i]);
		  fprintf(output,"WARNING: Fit: %e\n\n", bfgsfit);
		}
		if (bfgsoutX[i] > domains[j][3]) {
		  BoundaryTrigger=1;
		  fprintf(output,
			  "\nWARNING: BFGS hit on best individual produced Out of Boundary individual.\n");
		  fprintf(output,"WARNING: Generation: %d \t Parameter: %d \t Value: %e\n", 
			  count_gener, i+1, bfgsoutX[i]);
		  fprintf(output,"WARNING: Fit: %e\n\n", bfgsfit);
		}
	      } /* end for loop */
	      
	      /* if we we use out of bounds individuals then proceed */
	      /* 0=anything goes, 1: regular; 2: no trespassing! */
	      if (BoundaryEnforcement==0) {
		for(i=1;i<=nvars;i++) 
		  {
		    population[1][i]=bfgsoutX[i-1];
		  }
		population[1][0]=bfgsfit;
	      }
	      else if (BoundaryTrigger==0) {
		for(i=1;i<=nvars;i++) 
		  {
		    population[1][i]=bfgsoutX[i-1];
		  }
		population[1][0]=bfgsfit;
	      }
	    } /* end if (population[1][0] > bfgs) */
	case 1:
	  if (population[1][0] < bfgsfit) /* maximize */
	    {
	      /* is the BFGS individual in the bounds? */
	      BoundaryTrigger=0; /* outside of bounds ? */
	      for (i=0; i<nvars; i++) {
		j = i+1;
		if (bfgsoutX[i] < domains[j][1]) {
		  BoundaryTrigger=1;
		  fprintf(output,
			  "\nWARNING: BFGS hit on best individual produced Out of Boundary individual.\n");
		  fprintf(output,"WARNING: Generation: %d \t Parameter: %d \t Value: %e\n\n", 
			  count_gener, i+1, bfgsoutX[i]);
		}
		if (bfgsoutX[i] > domains[j][3]) {
		  BoundaryTrigger=1;
		  fprintf(output,
			  "\nWARNING: BFGS hit on best individual produced Out of Boundary individual.\n");
		  fprintf(output,"WARNING: Generation: %d \t Parameter: %d \t Value: %e\n\n", 
			  count_gener, i+1, bfgsoutX[i]);
		}
	      } /* end for loop */
	      
	      /* if we we use out of bounds individuals then proceed */
	      /* 0=anything goes, 1: regular; 2: no trespassing! */
	      if (BoundaryEnforcement==0) {
		for(i=1;i<=nvars;i++) 
		  {
		    population[1][i]=bfgsoutX[i-1];
		  }
		population[1][0]=bfgsfit;
	      }
	      else if (BoundaryTrigger==0) {
		for(i=1;i<=nvars;i++) 
		  {
		    population[1][i]=bfgsoutX[i-1];
		  }
		population[1][0]=bfgsfit;
	      }
	    } /* end if (population[1][0] < bfgsfit) */
	} /* end switch */
	/* end of bfgs stuff */
      } /* end of UseBFGS */  

      switch(MinMax)
        {
	case 0:
	  if(Teval > population[1][0])
	    {
	      Teval = population[1][0];
	      fprintf(output,"%7lu \t%e\n",
		      count_gener,population[1][0]); 
	      fflush(output);
	      peak_cnt = count_gener;
	      peak_val = population[1][0];
	    }
	  break;
	case 1:
	  if(Teval < population[1][0])
	    {
	      Teval = population[1][0];
	      fprintf(output,"%7lu \t%e\n",
		      count_gener,population[1][0]);
	      fflush(output);
	      peak_cnt = count_gener;
	      peak_val = population[1][0];
	    }
	  break;
        }
      
      /* compute and print mean and variance of population */
      if (PrintLevel==2) {
	fprintf(output,"\nGENERATION: %d\n", count_gener);
	populationstats(population, pop_size, nvars, mean, var, skew, kur, tobs);
	for (i=0; i<=nvars; i++) {
	  if (i==0) {
	      fprintf(output, "Fitness Value... %e\n", population[1][i]);
	      fprintf(output, "mean............ %e\n", mean[i]);
	      fprintf(output, "var............. %e\n", var[i]);
	      fprintf(output, "skewness........ %e\n", skew[i]);
	      fprintf(output, "kurtosis........ %e\n", kur[i]);
	      fprintf(output, "#null........... %d\n", pop_size-tobs[i]);
	      if(Structure->MemoryUsage==1)
		fprintf(output, "#unique......... %d, #Total UniqueCount: %d\n", 
			UniqueCount-OldUniqueCount, UniqueCount);
	      fprintf(output, "tobs............ %d\n", tobs[i]);
	  }
	  else {
	      fprintf(output, "var %d:\n", i);
	      fprintf(output, "best............ %e\n", population[1][i]);
	      fprintf(output, "mean............ %e\n", mean[i]);
	      fprintf(output, "var............. %e\n", var[i]);
	      fprintf(output, "skewness........ %e\n", skew[i]);
	      fprintf(output, "kurtosis........ %e\n", kur[i]);
	      fprintf(output, "#null........... %d\n", pop_size-tobs[i]);
	      fprintf(output, "tobs............ %d\n", tobs[i]);
	  }
	}
      } /* end of printlevel if */

      if (PrintLevel==1) {
	popmean = popvar = 0.0 ;
	popwrk = 1.0 / pop_size ;
	for(i=1; i<=pop_size; i++) {
	  popmean += population[i][0] ;
	}
	popmean *= popwrk ;
	for(i=1; i<=pop_size; i++) {
	  popstat =  population[i][0] - popmean ;
	  popvar += (popstat*popstat) ;
	}
	popvar *= popwrk ;
	fprintf(output, "   mean = %e, variance = %e\n\n", popmean, popvar);
      }

      fflush(output);

      /*
#ifdef MS_WINDOWS
      sort(MinMax,population,pop_size,1);
      // Let's graph the current population! 
      for (i=1; i<=pop_size; i++) 
	{
	  // x1 
	  Structure->pvm->vmArgs[2*(i-1)].value = population[i][1];
	  Structure->pvm->vmArgs[2*(i-1)].type  = 0;

	  // fit 
	  Structure->pvm->vmArgs[2*(i-1)+1].value = population[i][0];
	  Structure->pvm->vmArgs[2*(i-1)+1].type  = 0;
	} // end of i loop 

      LVMreturn = LVM_EvalSubScript(Structure->pvm, "GenoudChart", "AllArgs Genoud.Graph1 end", pop_size*2);
	  fprintf(output,"Done LVM_EvalSub\n");
	  fprintf(output,"LVMreturn: %d\n", LVMreturn);	
      sort(MinMax,population,pop_size,0);
#endif
      */
	
      /* Print the population file */
      if ( PrintLevel == 1 ) {
	if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
	  fprintf(output,"Unable to open the project file: %s", 
		  Structure->ProjectPath);

	  /* free populationstats stuff */
	  free(mean);
	  free(var);
	  free(skew);
	  free(kur);
	  free(tobs);
	  
	  free(bfgsoutX);
	  free(finalhessin);
	  free(evalX);
	  free(grad);
	  
	  /* free numeric.c allocations */
	  if (Structure->MemoryUsage==1)
	      JaMatrixFree(Memory, MemorySize);

	  JaMatrixFree(population, pop_size+2);
	  JaMatrixFree(new_genera, pop_size+2);
	  
	  free_matrix(temp, 1, 2, 0);
	  free_vector(probab, 1);
	  free_vector(t_vec, 1);
	  free_vector(cum_probab, 1);
	  free_ivector(live, 1);

	  return(ERROR_CODE);
	}
	print_population(pop_size, nvars, count_gener, population, popout);
	fclose(popout);
      } /* end of PrintLevel if */
      if ( PrintLevel == 2 ) {
	if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
	  fprintf(output,"Unable to open the project file: %s", 
		  Structure->ProjectPath);

	  /* free populationstats stuff */
	  free(mean);
	  free(var);
	  free(skew);
	  free(kur);
	  free(tobs);
	  
	  free(bfgsoutX);
	  free(finalhessin);
	  free(evalX);
	  free(grad);
	  
	  /* free numeric.c allocations */
	  if (Structure->MemoryUsage==1)
	    JaMatrixFree(Memory, MemorySize);

	  JaMatrixFree(population, pop_size+2);
	  JaMatrixFree(new_genera, pop_size+2);
	  
	  free_matrix(temp, 1, 2, 0);
	  free_vector(probab, 1);
	  free_vector(t_vec, 1);
	  free_vector(cum_probab, 1);
	  free_ivector(live, 1);

	  return(ERROR_CODE);
	}
	print_population(pop_size, nvars, count_gener, population, popout);
	fflush(popout);
	fclose(popout);
      }
      
      if (count_gener == 1) {
	oldfitvalue=population[1][0];
      }
      
      switch(MinMax)
	{
	case 0:
	  if (oldfitvalue - SolutionTolerance > population[1][0]) {
	    nochange_gen=0;
	    oldfitvalue=population[1][0];
	  }
	  else nochange_gen++;
	  break;
	case 1:
	  if (oldfitvalue + SolutionTolerance < population[1][0]) {
	    nochange_gen=0;
	    oldfitvalue=population[1][0];
	  }
	  else nochange_gen++;	      
	  break;
	}
      
      if (nochange_gen > (WaitGenerations)) {
	/* increase the number of WaitGenerations if the gradients are NOT zero! */	  
	if (GradientCheck==0) {
	  fprintf(output,"\nSoft Generation Wait Limit Hit.\n");
	  fprintf(output,"No Improvement in %d Generations\n", nochange_gen-1);
	  fflush(output);
	  MaxGenerations = 0;
	  nochange_gen=0;
	}
	else  {
	  for (i=1; i<=nvars; i++)
	    {
		  bfgsoutX[i-1]=population[1][i];
	    }
	  gradient(Structure->AgentFit, bfgsoutX, grad, nvars, MinMax, BoundaryEnforcement, InstanceNumber, 
		   domains, Status);
	  if (*Status < 0) {
	    /* Free Memory */
	      
	      /* free populationstats stuff */
	      free(mean);
	      free(var);
	      free(skew);
	      free(kur);
	      free(tobs);
	      
	      free(bfgsoutX);
	      free(finalhessin);
	      free(evalX);
	      free(grad);
	      
	      /* free numeric.c allocations */
	      if (Structure->MemoryUsage==1)
		JaMatrixFree(Memory, MemorySize);

	      JaMatrixFree(population, pop_size+2);
	      JaMatrixFree(new_genera, pop_size+2);
	      
	      free_matrix(temp, 1, 2, 0);
	      free_vector(probab, 1);
	      free_vector(t_vec, 1);
	      free_vector(cum_probab, 1);
	      free_ivector(live, 1);
	      
	      return(ERROR_CODE);
	  }
	  GradientTrigger = 0;
	  for (i=0; i<nvars; i++) {
	    if (fabs(grad[i]) > SolutionTolerance) {
	      GradientTrigger = 1;
	      break;
		}
	  } /* end for loop */
	  if (GradientTrigger==1) {
	    IncreaseGenerations = WaitGenerations;
	    WaitGenerations += IncreaseGenerations;
		fprintf(output,
			"\nDoubling Soft Maximum Wait Generation Limit to %d (from %d).\n", 
			WaitGenerations, IncreaseGenerations);
		fprintf(output,"I'm doing this because at least one gradient is too large.\n");
		fprintf(output,"G[%d]: %e\t Solution Tolerance: %e\n\n", 
			i+1, grad[i], SolutionTolerance);
	  }
	  else {
	    fprintf(output,"\nSoft Generation Wait Limit Hit.\n");
	    fprintf(output,"No Improvement in %d Generations\n", nochange_gen-1);
		fflush(output);
		MaxGenerations = 0;
		nochange_gen=0;
	  }
	}/* end if loop */
      }
      
      if ( (count_gener == MaxGenerations) && (GradientTrigger==1) ) 
	{
	  if (HardGenerationLimit==0)
	    {
	      IncreaseGenerations = MaxGenerations;
	      MaxGenerations += IncreaseGenerations;
	      fprintf(output,
		      "\nIncreasing Soft Maximum Generation Limit by %d (MaxGenerations) to %d.\n", 
		      IncreaseGenerations, MaxGenerations);
	      fprintf(output,"I'm doing this because at least one gradient is too large.\n\n");
	    } // if (Structure->HardGenerationLimit==0)
	  else
	    {
	      fprintf(output,"\nSTOPPING: HARD MAXIMUM GENERATION LIMIT HIT\n");
	      fprintf(output,"          At least one gradient is still too large\n");
	    } // else
	} // if ( (count_gener == MaxGenerations) && (GradientTrigger==1) ) 


      /* increase the number of generations if fitness has been improving */
      if ( (count_gener == MaxGenerations) &&  (nochange_gen < WaitGenerations) ) {
	if (HardGenerationLimit==0)
	  {
	    if (WaitGenerations > MaxGenerations) {
	      IncreaseGenerations = WaitGenerations;
	      MaxGenerations += (int) (IncreaseGenerations);
	      fprintf(output,
		      "\nIncreasing Soft Maximum Generation Limit by %d (WaitGenerations) to %d\n", 
		      IncreaseGenerations, MaxGenerations);
	      fprintf(output,"I'm doing this because the fitness is still impoving.\n\n");
	    }
	    else {
	      IncreaseGenerations = MaxGenerations;
	      MaxGenerations += (int) (IncreaseGenerations);
	      fprintf(output,
		      "\nIncreasing Soft Maximum Generation Limit by %d (MaxGenerations) to %d.\n", 
		      IncreaseGenerations, MaxGenerations);
	      fprintf(output,"I'm doing this because the fitness is still improving.\n\n");
	    }
	  } // if (Structure->HardGenerationLimit==0)
	else
	  {
	    fprintf(output,"\nSTOPPING: HARD MAXIMUM GENERATION LIMIT HIT\n");
	    fprintf(output,"          But fitness is still improving\n");
	  }
      } // if ( (count_gener == MaxGenerations) &&  (nochange_gen < WaitGenerations) )
      
      fflush(output);

      /* Should we recheck the main data structure for changes to the operator set? If so, let's
       do it now */
      if (Structure->AllowDynamicUpdating==1) {
	pop_size_old = pop_size;
	fprintf(output,"\nUpdating Main Data Structure:\n");
	SetRunTimeParameters(Structure, 0,
			     &pop_size, &nvars, &MaxGenerations, &WaitGenerations,
			     &MinMax, &GradientCheck, &BoundaryEnforcement, &UseBFGS, &SolutionTolerance,
			     &InstanceNumber, &P, &P0, &P1, &P2, &P3, &P4, &P5, &P6, &P7, &P8, 
			     &PrintLevel, &HardGenerationLimit, output);
	if (pop_size > pop_size_old) {
	  population_old = JaMatrixAllocate(pop_size_old+2, nvars+2);
	  
	  for (i=1; i<=pop_size_old;i++) {
	    for (j=0; j<=nvars+1; j++) {
	      population_old[i][j] = population[i][j];
	    }
	  }
	  JaMatrixFree(population, pop_size_old+2);
	  population    = JaMatrixAllocate(pop_size+2, nvars+2);

	  for (i=1; i<=pop_size_old;i++) {
	    for (j=0; j<=nvars+1; j++) {
	      population[i][j] = population_old[i][j];
	    }
	  }	  
	  JaMatrixFree(population_old, pop_size_old+2);

	  /* we need to add individuals to population! */
	  for (j=(pop_size_old+1); j<=pop_size; j++) {
	    for (i=1; i<=nvars; i++) {
	      population[j][i] = frange_ran(domains[i][1], domains[i][3]); 
	      population[j][nvars+1] = 100.0; /* 100.0=new! */
	    }
	  }
	} /* end of if popsize > pop_size_old */

	JaMatrixFree(new_genera, pop_size_old+2);	
	free_vector(probab, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);

	new_genera    = JaMatrixAllocate(pop_size+2, nvars+2);
	temp       = matrix(1,2,0,nvars);
	probab     = Gvector(1,pop_size);
	t_vec      = Gvector(1,nvars);
	cum_probab = Gvector(1,pop_size);
	live       = ivector(1,pop_size);

	Structure->AllowDynamicUpdating=0;
      } // end of if AllowDynamicUpdating==1
      
    } /* end of do loop */
  /*Increment iteration count and test whether all generations are done*/
  while (++count_gener <= MaxGenerations);
  
  fprintf(output,"\nBest Fit Found at Generation %lu\nFit Value = %e\n",peak_cnt,peak_val);
  fprintf(output,"\n\nParameters at the Solution (value, gradient):\n\n");

  /* output data structure */
  Structure->oPeakGeneration=peak_cnt;
  Structure->oGenerations=count_gener-1;

  /* obtain gradients */
  if (GradientCheck==0 && UseBFGS==0) {
    fprintf(output,"\nNot Obtaining Gradient Information\n");
    for (i=0; i< nvars; i++) {
      grad[i]=-1.0;
    }
  }
  else {
    for (i=1; i<=nvars; i++)
      {
	bfgsoutX[i-1]=population[1][i];
      }
    gradient(Structure->AgentFit, bfgsoutX, grad, nvars, MinMax, BoundaryEnforcement, 
	     InstanceNumber, domains, Status);
    if (*Status < 0) {
	/* Free Memory */
	
	/* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);
	
	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	
	/* free numeric.c allocations */
	free_matrix(population, 1, pop_size,0);
	free_matrix(new_genera, 1, pop_size, 0);
	free_matrix(temp, 1, 2, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);
	
	return(ERROR_CODE);
    }
  }
  
  /* print best solution */
  for(j = 1; j <= nvars; j++) {
    i = j-1;
    fprintf(output," X[%2d] :\t%e\tG[%2d] :\t%e\n",j,population[1][j],j,grad[i]);
    Results[i] = population[1][j];
    Gradients[i] = grad[i];
  }
  
  /* free memory */

  /* free populationstats stuff */
  free(mean);
  free(var);
  free(skew);
  free(kur);
  free(tobs);
  
  free(bfgsoutX);
  free(finalhessin);
  free(evalX);
  free(grad);

  /* free numeric.c allocations */
  if (Structure->MemoryUsage==1)
    JaMatrixFree(Memory, MemorySize);

  JaMatrixFree(population, pop_size+2);
  JaMatrixFree(new_genera, pop_size+2);

  free_matrix(temp, 1, 2, 0);
  free_vector(probab, 1);
  free_vector(t_vec, 1);
  free_vector(cum_probab, 1);
  free_ivector(live, 1);

  return(peak_val);

}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   sort()                                       */
/*                                                                              */
/*           SYNOPSIS          :   void sort(MinMax, population,pop_size,       */
/*                                 variable)                                    */
/*                                                                              */
/*           DESCRIPTION       :   This function sorts the population, in the   */
/*                                  ascending or the descending order of the    */
/*                                  evaluation function, depending on whether   */
/*                                  it is a maximization or a minimization      */
/*                                  function, respectively.                     */
/*                                                                              */
/*                                  As an alternative, the sortq function below */
/*                                  can be used, That sorting function uses     */
/*                                  the quicksort algorithm.                    */
/*                                                                              */
/*                                                                              */
/*           FUNCTIONS CALLED  :   swap()                                       */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/********************************************************************************/
void sort(short int MinMax, MATRIX  population, int pop_size,
	  long nvar)
     /*
       short int MinMax;      Tells whether it is a maximizaiton or a minimization function
       int pop_size;          Population size
       MATRIX population;     Array of population
     */
{
  int i,j;


  /*If MinMax is 0 sorts in the descending order, and*/
  /*if it is 1 sorts in the ascending order*/
  /*Sorted in ascending or descending order, based on*/
  /*the evaluated values of each of the agents*/
  switch(MinMax)
    {
    case 0 :
      for(i=1; i<=pop_size; i++)
        for(j=i+1; j<=pop_size; j++)
          if(population[i][nvar] > population[j][nvar])
            swap(&population[i],&population[j]);
      break;

    case 1 :
      for(i=1; i<=pop_size; i++)
        for(j=i+1; j<=pop_size; j++)
          if(population[i][nvar] < population[j][nvar])
            swap(&population[i],&population[j]);
      break;
    }
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   swap()                                       */
/*                                                                              */
/*           SYNOPSIS          :   void swap(x,y)                               */
/*                                                                              */
/*           DESCRIPTION       :   This function interchanges the values of     */
/*                                  x and y.                                    */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   sort()                                       */
/*                                                                              */
/*                                                                              */
/********************************************************************************/



void swap(double **x, double **y)
     /* double **x,**y; */
{
  double *temp;

  temp = *x;
  *x = *y;
  *y = temp;
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_parent()                                */
/*                                                                              */
/*           SYNOPSIS          :   int find_parent(live,pop_size)               */
/*                                                                              */
/*           DESCRIPTION       :   This function returns the index of the       */
/*                                  agent in the population, which is to be     */
/*                                  chosen for reproduction.                    */
/*                                                                              */
/*           FUNCTIONS CALLED  :   irange_ran()                                 */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

int find_parent(IVECTOR live, int pop_size)
     /*
       int pop_size;    Population size
       IVECTOR live;    Vector containing the number of times each agent
                        is going to reproduce
     */
{
  int i,temp,t1=0,tot=0;

  /*Finding the total number of parents to reproduce*/
  for(i=1; i<=pop_size; i++)
    tot = tot + live[i];
  if(tot==0)
    {
      // printf("No agents to select");
      exit(1);
    }

  /*Choosing one of them randomly*/
  temp = irange_ran(1,tot);

  tot = 0;
  i = 1;
  do{
    if(live[i]!=0)
      t1 = i;
    tot = tot + live[i++];
  }while(tot<temp);

  /*Decrementing the number of times the parent chosen is going to reproduce*/
  // live[t1]--;
  return(t1);
}
/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   assign_probab()                              */
/*                                                                              */
/*           SYNOPSIS          :   void assign_probab(probab,pop_size,Q)        */
/*                                                                              */
/*           DESCRIPTION       :   This function assigns probability of survival*/
/*                                  to each of the agents determined by the     */
/*                                  value provided by the user for the          */
/*                                  probability of the best agnet.              */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/********************************************************************************/


void assign_probab(VECTOR probab, int pop_size, double Q)
     /*
       int pop_size;     Population size
       double Q;         The probability of survival of the best agent
       VECTOR probab;    Array to contain the probability of survival
                         of each of the agents
     */
{
  int i;

  /* Q, Q(1-Q)^1, Q(1-Q)^2 ... Q(1-Q)^n */
  for(i=1; i<=pop_size; i++)
    probab[i] = Q * x_pow_y(1-Q,i-1);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   x_pow_y()                                    */
/*                                                                              */
/*           SYNOPSIS          :   double x_pow_y(x,y)                           */
/*                                                                              */
/*           DESCRIPTION       :   This function returns the value of x to the  */
/*                                  power of y.                                 */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   assign_probab()                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/********************************************************************************/



double x_pow_y(double x, int y)
     /*
       double x;
       int y;
     */
{
  int i;
  double tot = 1.0;

  for(i=0; i < y; i++)
    tot = tot * x;
  return(tot);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_cum__probab()                           */
/*                                                                              */
/*           SYNOPSIS          :   void find_cum__probab(cum_probab,probab,     */
/*                                                                     pop_size)*/
/*                                                                              */
/*           DESCRIPTION       :   This function finds the cumulative           */
/*                                  probability of each of the agents, from the */
/*                                  individual probability found earlier.       */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/********************************************************************************/



void find_cum_probab(VECTOR cum_probab, VECTOR probab, int pop_size)
     /*
       int pop_size;     Population size
       VECTOR probab,      Individual probability of survival of each of the agent
       cum_probab;         Cumulative probability of survival of each of the agent
     */
{
  int i;

  cum_probab[1] = probab[1];

  for(i=2; i<=pop_size; i++)
    cum_probab[i] = cum_probab[i-1] + probab[i];

}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_live()                                  */
/*                                                                              */
/*           SYNOPSIS          :   void find_live(cum_probab,live,pop_size,P4+P5*/
/*                                                                              */
/*           DESCRIPTION       :   This function finds the agents from the      */
/*                                  population, who are going to live - those   */
/*                                  who are going to reproduce, which is done   */
/*                                  based on the cumulative probability of      */
/*                                  survival of each of the agents.             */
/*                                                                              */
/*           FUNCTIONS CALLED  :   frange_ran()                                 */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/********************************************************************************/


void find_live(VECTOR cum_probab, IVECTOR live, int pop_size, int P)
     /*
       VECTOR cum_probab;  Cumulative probability
       IVECTOR live;       Agents that are going to reproduce
       int pop_size,       Population size
       P;                  Total number of parents needed to reproduce
     */
{
  double random;
  int count=0,/*Count of the number of agents chosen to live*/
      i;

  do
    {
      /*Choosing a random cumulative probability*/
      random = frange_ran(0.0,1.0);
      i=0;
      /*Finding the agent with the chosen cumulative probability*/
      do{
        i++;
        }while((random > cum_probab[i]) && (i< pop_size));

      /*Chosing the parent with that probability to reproduce*/
      if(count < P)
        {
          live[i]++;
          count++;
        }
    }while(count < P);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_die()                                   */
/*                                                                              */
/*           SYNOPSIS          :   void find_die(cum_probab,die,pop_size,P4+P5) */
/*                                                                              */
/*           DESCRIPTION       :   This function finds the agents from the      */
/*                                  population, who are going to die.           */
/*                                                                              */
/*           FUNCTIONS CALLED  :   frange_ran()                                 */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/********************************************************************************/



int find_die(VECTOR cum_probab, IVECTOR die, int pop_size)
     /*
       VECTOR cum_probab; Cumulative probability
       IVECTOR die;       Agents that are going to die
       int pop_size;      Population size
     */
{
  double random;
  int i;
  int done = FALSE;

  do
    {
      /*Choosing a random cumulative probability*/
      random = frange_ran(0.0,1.0);
      i=0;
      /*Finding the agent with the chosen cumulative probability*/
      do{
        i++;
        }
      while((random > cum_probab[i]) && (i< pop_size));

      /*Chosing the agent to die*/
      if ((die[pop_size+1-i] == 0) && (i < pop_size))
        done = TRUE;
    }
  while(!done);
  return(pop_size+1-i);
}



void SetRunTimeParameters(struct GND_IOstructure *Structure, 
			  short FirstTime,
			  long *PopSize, long *nvars, long *MaxGenerations, long *WaitGenerations,
			  short *MinMax, short *GradientCheck, short *BoundaryEnforcement, short *UseBFGS,
			  double *SolutionTolerance,
			  long *InstanceNumber, long *P, long *P0, long *P1, long *P2, long *P3, long *P4, long *P5, 
			  long *P6, long *P7, long *P8, short *PrintLevel, 
			  short *HardGenerationLimit, FILE *output)
{
  double tP;
  int i;

  *PopSize=Structure->PopSize;
  *nvars=Structure->nvars;
  *MaxGenerations=Structure->MaxGenerations;
  *WaitGenerations=Structure->WaitGenerations;

  if (FirstTime==1)
    *HardGenerationLimit=Structure->HardGenerationLimit;

  *MinMax=Structure->MinMax;
  *GradientCheck=Structure->GradientCheck;
  *BoundaryEnforcement=Structure->BoundaryEnforcement;
  *UseBFGS=Structure->UseBFGS;
  *InstanceNumber=Structure->InstanceNumber;

  *SolutionTolerance=Structure->SolutionTolerance;
  *PrintLevel=Structure->PrintLevel;

  /* Check to make sure that all operators are positve numbers! */
  if (Structure->P[0] < 0 ) {
    fprintf(output,"\n\nWARNING: Operator 0 (Cloning) was Assigned an Illegal Value: %d\n", 
	    Structure->P[1]);
    fprintf(output,"WARNING: I'm Assigning Zero Operators of Type 0.\n");
    Structure->P[0]=0.0;
  }
  if (Structure->P[1] < 0 ) {
    fprintf(output,"\n\nWARNING: Operator 1 (Uniform Mutation) was Assigned an Illegal Value: %d\n", 
	    Structure->P[1]);
    fprintf(output,"WARNING: I'm Assigning Zero Operators of Type 1.\n");
    Structure->P[1]=0.0;
  }
  if (Structure->P[2] < 0 ) {
    fprintf(output,"\n\nWARNING: Operator 2 (Boundary Mutation) was Assigned an Illegal Value: %d\n", 
	    Structure->P[2]);
    fprintf(output,"WARNING: I'm Assigning Zero Operators of Type 2.\n");
    Structure->P[2]=0;
  }
  if (Structure->P[3] < 0 ) {
    fprintf(output,"\n\nWARNING: Operator 3 (Non-Uniform Mutation) was Assigned an Illegal Value: %d\n", 
	    Structure->P[3]);
    fprintf(output,"WARNING: I'm Assigning Zero Operators of Type 3.\n");
    Structure->P[3]=0;
  }
  if (Structure->P[4] < 0 ) {
    fprintf(output,"\n\nWARNING: Operator 4 (Polytope Crossover) was Assigned an Illegal Value: %d\n", 
	    Structure->P[4]);
    fprintf(output,"WARNING: I'm Assigning Zero Operators of Type 4.\n");
    Structure->P[4]=0;
  }
  if (Structure->P[5] < 0 ) {
    fprintf(output,"\n\nWARNING: Operator 5 (Multiple Point Simple Crossover) was Assigned an Illegal Value: %d\n", 
	    Structure->P[5]);
    fprintf(output,"WARNING: I'm Assigning Zero Operators of Type 5.\n");
    Structure->P[5]=0;
  }
  if (Structure->P[6] < 0 ) {
    fprintf(output,"\n\nWARNING: Operator 6 (Whole Non-Uniform Mutation) was Assigned an Illegal Value: %d\n", 
	    Structure->P[6]);
    fprintf(output,"WARNING: I'm Assigning Zero Operators of Type 6.\n");
    Structure->P[6]=0;
  }
  if (Structure->P[7] < 0 ) {
    fprintf(output,"\n\nWARNING: Operator 7 (Heuristic Crossover) was Assigned an Illegal Value: %d\n", 
	    Structure->P[7]);
    fprintf(output,"WARNING: I'm Assigning Zero Operators of Type 7.\n");
    Structure->P[7]=0;
  }

  /* if DataType==1 (i.e., integer) we are not giong to use any gradient information etc. */
  if (Structure->DataType==1) {
    *UseBFGS=0;
    *GradientCheck=0;

    if (Structure->P[8] > 0) {
      fprintf(output,"\n\nWARNING: Operator 8 (Local-Minimum Crossover) was Assigned an Illegal Value: %d\nThis is an illegal value because we are working with integer data", 
	      Structure->P[8]);
      fprintf(output,"WARNING: I'm Assigning Zero Operators of Type 8.\n");
      Structure->P[8]=0;
    } /* end of if */
  }
  else {
    if (Structure->P[8] < 0 ) {
      fprintf(output,"\n\nWARNING: Operator 8 (Local-Minimum Crossover) was Assigned an Illegal Value: %d\n", 
	      Structure->P[8]);
      fprintf(output,"WARNING: I'm Assigning Zero Operators of Type 8.\n");
      Structure->P[8]=0;
    }
  } /* end of else */

  /* Let's figure out the number of operators we need.  Move stuff to absolute space */
  tP = 0;
  for (i=0; i<9; i++) {
    tP = tP + Structure->P[i] ;
  }
  
  if (tP > 0) {
    *P0 = Iround(  (Structure->P[0] /  tP) * (double) (*PopSize-2) );
    *P1 = Iround(  (Structure->P[1] /  tP) * (*PopSize-2) );
    *P2 = Iround(  (Structure->P[2] /  tP) * (*PopSize-2) );
    *P3 = Iround(  (Structure->P[3] /  tP) * (*PopSize-2) );
    *P4 = Iround(  (Structure->P[4] /  tP) * (*PopSize-2) );
    *P5 = Iround(  (Structure->P[5] /  tP) * (*PopSize-2) );
    *P6 = Iround(  (Structure->P[6] /  tP) * (*PopSize-2) );
    *P7 = Iround(  (Structure->P[7] /  tP) * (*PopSize-2) );
    *P8 = Iround(  (Structure->P[8] /  tP) * (*PopSize-2) );
  }
  else {
    *P0 = 0;
    *P1 = 0;
    *P2 = 0;
    *P3 = 0;
    *P4 = 0;
    *P5 = 0;
    *P6 = 0;
    *P7 = 0;
    *P8 = 0;
  }

  /* Check to make sure that all operators (i.e., 4, 5, 7) which have to be even numbers are */
  if (fmod(*P4,2) > 0.0) {
    fprintf(output,"\nWARNING: Operator 5 (Polytope Crossover) may only be started\n");
    fprintf(output,"WARNING: an even number of times.  I am increasing this operator by one.\n");
    *P4=*P4+1;
  }
  if (fmod(*P5,2) > 0.0) {
    fprintf(output,"\nWARNING: Operator 6 (Multiple Point Simple Crossover) may only be started\n");
    fprintf(output,"WARNING: an even number of times.  I am increasing this operator by one.\n");
    *P5=*P5+1;
  }
  if (fmod(*P7,2) > 0.0) {
    fprintf(output,"\nWARNING: Operator 8 (Heuristic Crossover) may only be started\n");
    fprintf(output,"WARNING: an even number of times.  I am increasing this operator by one.\n");
    *P7=*P7+1;
  }

  /*P is the total number of parents needed for applying all the operators*/
  *P = *P1 + *P2 + *P3 + *P4 + *P5 + *P6 + *P7 + *P8;
  if(*P > *PopSize)
    {
      fprintf(output,"\nWARNING: The total number of operators greater than population size\n");

      if (fmod(*P+1,2) > 0.0) {
	*PopSize = *P+2;
      fprintf(output,"WARNING: I'm increasing the population size to %d (operators+2).\n", *PopSize);
      }
      else {
	*PopSize = *P+1;
	fprintf(output,"WARNING: I'm increasing the population size to %d (operators+1).\n", *PopSize);
      }
    }
  else if ( *P== *PopSize) {
      fprintf(output,"\nWARNING: The total number of operators equal to the population size\n");

      if (fmod(*P+1,2) > 0.0) {
	*PopSize = *P+2;
      fprintf(output,"WARNING: I'm increasing the population size to %d (operators+2).\n", *PopSize);
      }
      else {
	*PopSize = *P+1;
	fprintf(output,"WARNING: I'm increasing the population size to %d (operators+1).\n", *PopSize);
      }
  }

  if (fmod(*PopSize,2) > 0.0) {
    fprintf(output,"WARNING: population size is not an even number.\n");
    fprintf(output,"WARNING: Increasing population size by 1\n");
    *PopSize=*PopSize+1;
  }

  /* return PopSize and P values to the return data structure */
  *P0 = *PopSize-*P-1;
  Structure->oP[0]=*P0;
  Structure->oP[1]=*P1;
  Structure->oP[2]=*P2;
  Structure->oP[3]=*P3;
  Structure->oP[4]=*P4;
  Structure->oP[5]=*P5;
  Structure->oP[6]=*P6;
  Structure->oP[7]=*P7;
  Structure->oP[8]=*P8;
  Structure->oPopSize=*PopSize;

  fprintf(output, "\n");
  if (Structure->DataType==1) fprintf(output, "Data Type: Integer\n");
  else fprintf(output, "Data Type: Floating Point\n");
  fprintf(output,"Operators (code number, name, population) \n");
  fprintf(output,"\t(1) Cloning........................... \t%d\n", *P0);
  fprintf(output,"\t(2) Uniform Mutation.................. \t%d\n", *P1);
  fprintf(output,"\t(3) Boundary Mutation................. \t%d\n", *P2);
  fprintf(output,"\t(4) Non-Uniform Mutation.............. \t%d\n", *P3);
  fprintf(output,"\t(5) Polytope Crossover................ \t%d\n", *P4);
  fprintf(output,"\t(6) Multiple Point Simple Crossover... \t%d\n", *P5);
  fprintf(output,"\t(7) Whole Non-Uniform Mutation........ \t%d\n", *P6);
  fprintf(output,"\t(8) Heuristic Crossover............... \t%d\n", *P7);
  fprintf(output,"\t(9) Local-Minimum Crossover........... \t%d\n\n", *P8);
  if (*HardGenerationLimit==0)
    fprintf(output,"SOFT Maximum Number of Generations: %lu\n", *MaxGenerations);
  else
    fprintf(output,"HARD Maximum Number of Generations: %lu\n", *MaxGenerations);
  fprintf(output,"Maximum Nonchanging Generations: %lu\n", *WaitGenerations);
  fprintf(output,"Population size       : %d\n", *PopSize);
  fprintf(output,"Convergence Tolerance: %e\n", *SolutionTolerance);

  fprintf(output, "\n");
  if (*UseBFGS !=0) {
    fprintf(output,
	    "Using the BFGS Derivative Based Optimizer on the Best Individual Each Generation.\n");
  }
  else {
    fprintf(output,
	    "Not Using the BFGS Derivative Based Optimizer on the Best Individual Each Generation.\n");
  }
  if (*GradientCheck==0)
    fprintf(output,"Not Checking Gradients before Stoping\n");
  else 
    fprintf(output,"Checking Gradients before Stoping\n");

  if (*BoundaryEnforcement==0) 
    fprintf(output,"Using Out of Bounds Individuals.\n\n");
  else if (*BoundaryEnforcement==1) 
    fprintf(output,"Not Using Out of Bounds Individuals But Allowing Trespassing.\n\n");
  else if (*BoundaryEnforcement==2) 
    fprintf(output,"Not Using Out of Bounds Individuals and Not Allowing Trespassing.\n\n");

  // some more consistency checks
  if (Structure->Network==1 && Structure->MemoryUsage==1)
    {
      fprintf(output,"Turning Memory Usage off because Network GENOUD has been requested\n"); 
      Structure->MemoryUsage=0;
    }
    
  fflush(output);

} /* End SetOperators */


/********************************************************************************/
/*  JaIntegerOptimization:                                                      */
/*                                                                              */
/*  This function assumes that the X variables are integers.                    */
/*                                                                              */
/*  The cross over operators are different!                                     */
/*                                                                              */
/********************************************************************************/

double JaIntegerOptimization(struct GND_IOstructure *Structure, VECTOR X, 
			     MATRIX domains, FILE *output)
{
  extern struct GND_IOstructure *ExternStructure;
  
  MATRIX new_genera,   /*Temporary storage for the new generation*/
         population,   /*Population of x2 variables*/
         temp;

  VECTOR probab,       /*Probability of agents to die or live*/
         cum_probab,   /*Cumilative probability of agents*/
         t_vec;

  IVECTOR live;


  long count_gener= 1; /*Counter to keep track of the number of generations*/
  unsigned long peak_cnt;

  int                     /*Total number of agents chosen to reproduce*/
    j1,
    j2,
    j3,
    j4,
    j5,
    j6,
    j7,
    j8,
    oper,
    ocnt,
    B,                     /*Parameter for the 3rd operator - nonuniform mutation*/
    STEP,                  /*Parameter for the 5th operator - simple arithmetical crossover*/
    first_live=0,          /*Index of the two parents for crossover parents*/
    second_live=0,
    first_die,             /*Index of the two parents for crossover death*/
    second_die,
    die_now,               /*index of agent to replace in current operation*/
    i,
    j,
    s;


  double Q,                   /*Probability of the best agent*/
         Teval=0,             /*Evaluation of the best agent*/
         peak_val;

  FLAG  same;
  double **Jnew;

  double *grad, *evalX, *finalhessin, *bfgsoutX;

  int nochange_gen=0;

  double oldfitvalue=0;

  int IncreaseGenerations;
  short int GradientTrigger=0;
  long InstanceNumber;

  /* Strucutre fixup! */
  long nvars, MaxGenerations, WaitGenerations;
  long *Status;
  long pop_size, P, P0, P1, P2, P3, P4, P5, P6, P7, P8;
  short int MinMax, GradientCheck, BoundaryEnforcement, UseBFGS;
  double SolutionTolerance, *Results, *Gradients;
  short PrintLevel, HardGenerationLimit;

  /* Old variables which may change when SetRunTimeParameters is run during a run! */
  long pop_size_old;
  double **population_old, **PopulationTmp;

  /* Summary Statistics (mean, variance etc) */
  double popmean, popvar, popwrk, popstat;

  /* Population Print population*/
  FILE *popout;
  long *tobs;
  double *mean, *var, *skew, *kur;

  /* Stuff for the Unique Stuff (how's that for an informative comment! */
  /* A big Matrix which remembers all of our past evaluations. It's
     maximum memory is set in genoud.h */
  extern long Gnvars[MAXINSTANCES];
  double **Memory;
  long UniqueCount, OldUniqueCount=0, MemorySize=0;
  // FLAG UniqueFlag, Redundant;
  // long upper, lower, midpoint;

  /* LVM calls */
  // long LVMreturn;

  /* fine two unique parents count */
  long SameCount, UniquePairs;

  // NetworkEvaluate() Stuff
  long NetworkNumber; // number of individuals to be evaluated

  ExternStructure=Structure;
  Status=&(Structure->Status);

  Results=Structure->oResults;
  Gradients=Structure->oGradients;

  /* Structure Done */
  SetRunTimeParameters(Structure, 1,
		       &pop_size, &nvars, &MaxGenerations, &WaitGenerations,
		       &MinMax, &GradientCheck, &BoundaryEnforcement, &UseBFGS, &SolutionTolerance,
		       &InstanceNumber, &P, &P0, &P1, &P2, &P3, &P4, &P5, &P6, &P7, &P8, 
		       &PrintLevel, &HardGenerationLimit, output);

  /* The population matrix is used from 0 to popsize (inclusive) and o to nvars+1 (inclusive) */
  /* population           = matrix(0,pop_size+2,0,nvars+1); */
  population    = JaMatrixAllocate(pop_size+2, nvars+2);
  PopulationTmp = JaMatrixAllocate(pop_size+1, nvars+2);
  new_genera    = JaMatrixAllocate(pop_size+2, nvars+2);

  temp       = matrix(1,2,0,nvars);
  probab     = Gvector(1,pop_size);
  t_vec      = Gvector(1,nvars);
  cum_probab = Gvector(1,pop_size);
  live       = ivector(1,pop_size);

  Gnvars[Structure->InstanceNumber]=nvars;

  if (Structure->MemoryUsage==1)
    {
      if (HardGenerationLimit==0)
	MemorySize=(MaxGenerations+1)*pop_size+1+pop_size;
      else
	MemorySize=(MaxGenerations+1)*pop_size+1+pop_size;
      
      Memory = JaMatrixAllocate(MemorySize, nvars+2);
    }

  grad = (double *) malloc((nvars)*sizeof(double));
  evalX = (double *) malloc((nvars)*sizeof(double));
  finalhessin = (double *) malloc(((nvars*nvars)+(nvars))*sizeof(double));
  bfgsoutX = (double *) malloc((nvars+1)*sizeof(double));

  /* populationstats variables */
  mean = (double *) malloc((nvars+1)*sizeof(double));
  var = (double *) malloc((nvars+1)*sizeof(double));
  skew = (double *) malloc((nvars+1)*sizeof(double));
  kur = (double *) malloc((nvars+1)*sizeof(double));
  tobs = (long *) malloc((nvars+1)*sizeof(long));

  Q=0.2;
  B=6;
  STEP=10;

  /* Read the Population File if we are told to.  For the moment we
     shall assume that we are supposed to read the file */

  fprintf(output,"\n\n");

  switch(MinMax) {
  case 0:
    fprintf(output,"Minimization Problem.\n\n");  
    break;
  case 1:
    fprintf(output,"Maximization Problem.\n\n");  
    break;
  }

  fprintf(output,"Parameter B (hardcoded): %d\n", B);
  fprintf(output,"Parameter Q (hardcoded): %f\n", Q);
  fprintf(output,"\n");

  fflush(output);

  peak_val = 0;
  peak_cnt = 0;

  pop_size_old=0;
  if (Structure->ShareType == 1 || Structure->ShareType == 3) {
    fprintf(output, "Using old population file to initialize new population\n");
    if((popout = fopen(Structure->ProjectPath, "r")) == NULL) {
      fprintf(output,"WARNING: Unable to open the old project file: %s\n", 
	      Structure->ProjectPath);
      fprintf(output,"         Generating new population\n");
    }
    else {
      pop_size_old=ReadPopulation(population, pop_size, nvars, output, popout);
      fclose(popout);

      for (i=1; i<=pop_size_old; i++) {
		  for (j=1; j<=nvars; j++) {
			population[i][j] = (int) population[i][j];
		  }
	  }

      if (pop_size_old<2) {
	fprintf(output,
		"WARNING: The old population file appears to be from the run of a different model!\n");
	pop_size_old=0;
      }
    }
    if (PrintLevel==2) {
      if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
	fprintf(output,"Unable to open the project file: %s", 
		Structure->ProjectPath);

	/* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);
	
	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	
	/* free numeric.c allocations */
	if (Structure->MemoryUsage==1)
	  JaMatrixFree(Memory, MemorySize);
	JaMatrixFree(population, pop_size+2);
	JaMatrixFree(new_genera, pop_size+2);
	
	free_matrix(temp, 1, 2, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);

	return(ERROR_CODE);
      }
      fclose(popout);
    }
  } /* end of ShareType 0 */
  else {
    if (PrintLevel==2) {
      if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
	fprintf(output,"Unable to open the project file: %s", 
		Structure->ProjectPath);
	
	/* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);
	
	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	
	/* free numeric.c allocations */
	if (Structure->MemoryUsage==1)
	  JaMatrixFree(Memory, MemorySize);
	JaMatrixFree(population, pop_size+2);
	JaMatrixFree(new_genera, pop_size+2);
	
	free_matrix(temp, 1, 2, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);
	
	return(ERROR_CODE);
      }
      fclose(popout);
    }
  }

  /* The new initial value matrix: setting a new initial value for every individual */
  if (ExternStructure->nStartingValues > 0) 
    {
      // seed the starting values until we run out of population or starting values!
      j = pop_size_old;
      for(s=0; s<ExternStructure->nStartingValues; s++) {
	j++;
	for(i=1; i<=nvars; i++) {
	  population[j][i] = (int) ExternStructure->StartingValues[s][i-1];
	  population[j][nvars+1] = -1.0;
	}
      } // end of for loop
      pop_size_old = j;

      // randomly add on people if we still have population left over!
      for(j=pop_size_old+1; j<=pop_size; j++) {
	for(i=1; i<=nvars; i++) { 
	  population[j][i] = (int) frange_ran(domains[i][1], domains[i][3]); 
	  population[j][nvars+1] = -1.0;
	}
      }
    } // end of we have starting values!
  else 
    {
      for(j=pop_size_old+1; j<=pop_size; j++) {
	for(i=1; i<=nvars; i++) { 
	  population[j][i] = (int) frange_ran(domains[i][1], domains[i][3]); 
	  population[j][nvars+1] = -1.0;
	}
      }
    } // end of else


  if (Structure->MemoryUsage==1)
    {
      /* BINARY SEARCH.  (Knuth 3:409). */
      OldUniqueCount=UniqueCount=0;
      
      /* BINARY SEARCH.  (Knuth 3:409). */
      /* We have sorted Memory by the "key" */
      /* JaIntegerSort(population, pop_size, nvars+2);  */
      JaIntegerSort(population, pop_size, nvars+2); 
      
      OldUniqueCount=UniqueCount;

      JaIntMemoryMatrix_Gen0(Structure,
				Memory, population, X,
				&UniqueCount, OldUniqueCount, pop_size, nvars, output, Status);
      if (*Status < 0) {
	fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
	  
	  /* free memory */
	if (Structure->MemoryUsage==1)
	  JaMatrixFree(Memory, MemorySize);

	  /* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);

	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	  
	/* free numeric.c allocations */
	JaMatrixFree(population, pop_size+2);
	JaMatrixFree(new_genera, pop_size+2);
	  
	free_matrix(temp, 1, 2, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);
	  
	return(ERROR_CODE);
      }      
      if ( (UniqueCount+pop_size) >= MemorySize )
	{
	  Structure->MemoryUsage=0;
	  fprintf(output,"\nWARNING: Turning Off MemoryMatrix because memory usage is too great.\n\n");
	} /* end of if */
    } // end of Memory based evaluation
  else
    {
      NetworkNumber=0;
      for (i=1; i<=pop_size; i++) 
	{

	  if (Structure->DynamicPopulation==2)
	    {
	      JaDynamicPopulationCheck(Structure, population, i, pop_size, nvars, output, Status);
	      if (*Status < 0) {
		fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
		
		/* free memory */
		
		/* free populationstats stuff */
		free(mean);
		free(var);
		free(skew);
		free(kur);
		free(tobs);
		
		free(bfgsoutX);
		free(finalhessin);
		free(evalX);
		free(grad);
		
		/* free numeric.c allocations */
		JaMatrixFree(population, pop_size+2);
		JaMatrixFree(new_genera, pop_size+2);
		
		free_matrix(temp, 1, 2, 0);
		free_vector(probab, 1);
		free_vector(t_vec, 1);
		free_vector(cum_probab, 1);
		free_ivector(live, 1);
		
		return(ERROR_CODE);
	      } // end of Status < 0
	    } // end of DynamicPopulation==2
	  
	  if (population[i][nvars+1]==-1.0 || population[i][nvars+1]==11.0)
	    {
	      if (Structure->Network==1)
		{
		  NetworkNumber++;
		  population[i][0] = EVALUATE;
		}
	      else
		{
		  for(j=1; j<=nvars; j++)
		    X[j] = population[i][j];
		  
		  population[i][0] = evaluate(Structure->AgentFit, X, nvars, Status);

		  if (*Status < 0) {
		    fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
		    
		    /* free memory */
		    
		    /* free populationstats stuff */
		    free(mean);
		    free(var);
		    free(skew);
		    free(kur);
		    free(tobs);
		    
		    free(bfgsoutX);
		    free(finalhessin);
		    free(evalX);
		    free(grad);
		    
		    /* free numeric.c allocations */
		    JaMatrixFree(population, pop_size+2);
		    JaMatrixFree(new_genera, pop_size+2);
		    
		    free_matrix(temp, 1, 2, 0);
		    free_vector(probab, 1);
		    free_vector(t_vec, 1);
		    free_vector(cum_probab, 1);
		    free_ivector(live, 1);
		    
		    return(ERROR_CODE);
		  }
		}
	    }
	} //end of i loop
      if (Structure->Network==1)
	{
	  NetworkEvaluate(Structure->AgentFit, Structure->DBname, Structure->AgentName, 
			  population, pop_size, nvars, NetworkNumber, Status, 10);
	  
	  if (*Status < 0) {
	    fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
	    
	    /* free memory */
	    
	    /* free populationstats stuff */
	    free(mean);
	    free(var);
	    free(skew);
	    free(kur);
	    free(tobs);
	    
	    free(bfgsoutX);
	    free(finalhessin);
	    free(evalX);
	    free(grad);
	    
	    /* free numeric.c allocations */
	    JaMatrixFree(population, pop_size+2);
	    JaMatrixFree(new_genera, pop_size+2);
	    
	    free_matrix(temp, 1, 2, 0);
	    free_vector(probab, 1);
	    free_vector(t_vec, 1);
	    free_vector(cum_probab, 1);
	    free_ivector(live, 1);
	    
	    return(ERROR_CODE);
	  }
	}
    } // end of default evaluation
  
  /*Sort the initial inidivduals based on their evaluation function*/
  sort(MinMax,population,pop_size,0);

  switch(MinMax) {
  case 0:
    Teval = population[1][0];
    peak_cnt = count_gener;
    peak_val = population[1][0];
    break;
  case 1:
    Teval = population[1][0];
    peak_cnt = count_gener;
    peak_val = population[1][0];
    break;
  }

  fprintf(output,"\nThe 2 best initial individuals are\n");
  for(i=1; i<3; i++) {
    print_vector(population[i],1,nvars,output);
    fprintf(output,"\nfitness = %e", population[i][0]);
    fprintf(output,"\n\n");
  }

  fprintf(output,"\nThe worst fit of the population is: %e\n", 
      population[pop_size][0]);
  fprintf(output,"\n\n");

  fprintf(output,"\n\nGeneration#\t    Solution Value\n");
  fprintf(output,"\n%7d \t%e\n", 0, population[1][0]);

  /* compute and print mean and variance of population */
  if (PrintLevel==2) {
      fprintf(output,"GENERATION: 0 (initializing the population)\n");
      populationstats(population, pop_size, nvars, mean, var, skew, kur, tobs);
      for (i=0; i<=nvars; i++) {
	  if (i==0) {
	      fprintf(output, "Fitness Value... %e\n", population[1][i]);
	      fprintf(output, "mean............ %e\n", mean[i]);
	      fprintf(output, "var............. %e\n", var[i]);
	      fprintf(output, "skewness........ %e\n", skew[i]);
	      fprintf(output, "kurtosis........ %e\n", kur[i]);
	      fprintf(output, "#null........... %d\n", pop_size-tobs[i]);
	      if(Structure->MemoryUsage==1)
		fprintf(output, "#unique......... %d, #Total UniqueCount: %d\n", 
			UniqueCount-OldUniqueCount, UniqueCount);
	      fprintf(output, "tobs............ %d\n", tobs[i]);
	  }
	  else {
	      fprintf(output, "var %d:\n", i);
	      fprintf(output, "best............ %e\n", population[1][i]);
	      fprintf(output, "mean............ %e\n", mean[i]);
	      fprintf(output, "var............. %e\n", var[i]);
	      fprintf(output, "skewness........ %e\n", skew[i]);
	      fprintf(output, "kurtosis........ %e\n", kur[i]);
	      fprintf(output, "#null........... %d\n", pop_size-tobs[i]);
	      fprintf(output, "tobs............ %d\n", tobs[i]);
	  }
      }
  } /* end of printlevel if */

  if (PrintLevel==1) {
    popmean = popvar = 0.0 ;
    popwrk = 1.0 / pop_size ;
    for(i=1; i<=pop_size; i++) {
      popmean += population[i][0] ;
    }
    popmean *= popwrk ; 
    for(i=1; i<=pop_size; i++) {
      popstat =  population[i][0] - popmean ;
      popvar += (popstat*popstat) ;
    }
    popvar *= popwrk ;
    fprintf(output, "   mean = %e, variance = %e\n\n", popmean, popvar);
  }

  fflush(output);

  /* Print the population file */
  if ( PrintLevel == 1 ) {
    if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
      fprintf(output,"Unable to open the project file: %s", 
	      Structure->ProjectPath);
      
      /* free populationstats stuff */
      free(mean);
      free(var);
      free(skew);
      free(kur);
      free(tobs);
      
      free(bfgsoutX);
      free(finalhessin);
      free(evalX);
      free(grad);
      
      /* free numeric.c allocations */
      if (Structure->MemoryUsage==1)
	JaMatrixFree(Memory, MemorySize);
      JaMatrixFree(population, pop_size+2);
      JaMatrixFree(new_genera, pop_size+2);
      
      free_matrix(temp, 1, 2, 0);
      free_vector(probab, 1);
      free_vector(t_vec, 1);
      free_vector(cum_probab, 1);
      free_ivector(live, 1);
      
      return(ERROR_CODE);
    }
    print_population(pop_size, nvars, 0, population, popout);
    fclose(popout);
  } /* end of PrintLevel if */
  if ( PrintLevel == 2 ) {
    if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
      fprintf(output,"Unable to open the project file: %s", 
	      Structure->ProjectPath);
      
      /* free populationstats stuff */
      free(mean);
      free(var);
      free(skew);
      free(kur);
      free(tobs);
      
      free(bfgsoutX);
      free(finalhessin);
      free(evalX);
      free(grad);
      
      /* free numeric.c allocations */
      if (Structure->MemoryUsage==1)
	JaMatrixFree(Memory, MemorySize);
      JaMatrixFree(population, pop_size+2);
      JaMatrixFree(new_genera, pop_size+2);
      
      free_matrix(temp, 1, 2, 0);
      free_vector(probab, 1);
      free_vector(t_vec, 1);
      free_vector(cum_probab, 1);
      free_ivector(live, 1);
      
      return(ERROR_CODE);
    }
    print_population(pop_size, nvars, 0, population, popout);
    fflush(popout);
    fclose(popout);
  }

  /*Assigning probability of survival for each of the agent, with the*/
  /*probability provided by the user for the best agent*/
  assign_probab(probab,pop_size,Q); 

  /*Finding the cumulative probability of the agents*/
  find_cum_probab(cum_probab,probab,pop_size);

  /*Reproducing and evaluating for the total number of generations times*/
  do
    {

      /*Initializing the live vector*/
      for(j=1; j<=pop_size; j++)
        {
          live[j] = 0;
          for(i=0; i<=nvars; i++)
            new_genera[j][i] = population[j][i];
        }

      /*Finding the agents that will die and the agents that will reproduce*/
      find_live(cum_probab,live,pop_size,P);
      /* set die_now counter to start replacements with the worst agent.
         use of die_now is okay if the entire population (except the previous
         best) is to be replaced in each generation */
      die_now = pop_size;

      j1=j2=j3=j4=j5=j6=j7=j8=0;

      UniquePairs= UniqueCount-OldUniqueCount;
      UniquePairs= (int) (0.5*(UniquePairs*UniquePairs-UniquePairs));
      if ( MAX_OPER_UNIQUE_TRY < UniquePairs)
	UniquePairs = MAX_OPER_UNIQUE_TRY;

      /* main operator loop */
      while(j1+j2+j3+j4+j4+j5+j5+j6+j7+j7 < P)
        {
          oper = irange_ran(1,7);
          switch (oper)
            {
	    case 1:
	      /* JS Description: Uniform Mutation */
	      /*Applying the first operator, uniform mutation*/
	      if (j1 < P1)
		{
		  /*Find one parent for mutation operator 1*/
		  first_live  = find_parent(live,pop_size);
		  live[first_live]--;
		  /* check that agent to replace is in range */
		  if (die_now < 2) {
		    fprintf(output, "No agents to be replaced\n");
		    exit(1);
		  }
		  
		  new_genera[die_now][nvars+1] = 1.0;
		  for(i=1; i<=nvars; i++)
		    t_vec[i] = population[first_live][i];
		  for (ocnt = irange_ran(1,nvars); ocnt>0; ocnt--)
		    JaIntegerOper1(t_vec,domains,nvars);
		  for(i=1; i<=nvars; i++)
		    new_genera[die_now][i] = t_vec[i];
		  die_now--;
		  j1++;
		}
	      break;
	    case 2:
	      /* JS Description: Boundary Mutation */
	      /*Applying the second operator, boundary mutation*/
	      if (j2 < P2)
		{
		  /*Find one parent for mutation operator 2*/
		  first_live  = find_parent(live,pop_size);
		  live[first_live]--;
		  /* check that agent to replace is in range */
		  if (die_now < 2) {
		    fprintf(output, "No agents to be replaced\n");
		    exit(1);
		  }
		  
		  new_genera[die_now][nvars+1] = 2.0;
		  for(i=1; i<=nvars; i++)
		    t_vec[i] = population[first_live][i];
		  JaIntegerOper2(t_vec,domains,nvars);
		  for(i=1; i<=nvars; i++)
		    new_genera[die_now][i] = t_vec[i];
		  die_now--;
		  j2++;
		}
	      break;
	    case 3:
	      /* JS Description: Non-uniform Mutation */
	      /*Applying the third operator, non-uniform mutation*/
	      if (j3 < P3)
		{
		  /*Find one parent for mutation operator 3*/
		  first_live  = find_parent(live,pop_size);
		  live[first_live]--;
		  /* check that agent to replace is in range */
		  if (die_now < 2) {
		    fprintf(output, "No agents to be replaced\n");
			  exit(1);
		  }
		  
		  new_genera[die_now][nvars+1] = 3.0;
		  for(i=1; i<=nvars; i++)
		    t_vec[i] = population[first_live][i];
		  for (ocnt = irange_ran(1,nvars); ocnt>0; ocnt--)
		    JaIntegerOper3(t_vec,domains,nvars,MaxGenerations,count_gener,B);
		  for(i=1; i<=nvars; i++)
		    new_genera[die_now][i] = t_vec[i];
		  die_now--;
		  j3++;
		}
	      break;
	      
	    case 4:
	      /* JS Description: Polytope Crossover */
	      /*Applying the fourth operator, whole arithmetical crossover*/
	      if (j4 < (int) P4/2)
		{
		  /*Find two distinct parents for crossover operator 4*/
		  same = TRUE;
		  SameCount=0;
		  while (same==TRUE) {
		    SameCount++;
		    
		    first_live  = find_parent(live,pop_size);
		    second_live = find_parent(live,pop_size);

		    if (SameCount >= (UniquePairs) ) 
		      break;

		    for(i=1; i<=nvars; i++)
		      if ( (int) population[first_live][i] != (int) population[second_live][i])
			same = FALSE;
		  } /* end of while same==TRUE loop */
		  /* check that agents to replace are in range */
		  if (die_now < 3) {
		    fprintf(output,"Not enough agents to be replaced\n");
		    exit(1);
		  }
		  live[first_live]--;
		  live[second_live]--;
		  first_die   = die_now-- ;
		  second_die  = die_now-- ;
		  new_genera[first_die][nvars+1]  = 4.0;
		  new_genera[second_die][nvars+1] = 4.0;
		  if (!same)
		    {
		      for(i=1; i<=nvars; i++)
			{
			  temp[1][i] = population[first_live][i];
			  temp[2][i] = population[second_live][i];
			}
		      JaIntegerOper4(temp[1],temp[2],nvars);
		      for(i=1; i<=nvars; i++)
			{
			  new_genera[first_die][i]  = temp[1][i];
			  new_genera[second_die][i] = temp[2][i];
			}
		    }
		  else {
		    /* copy agent chosen twice into two new indivs */
		    for(i=1; i<=nvars; i++) {
		      new_genera[first_die][i]  = 
			population[first_live][i];
		      new_genera[second_die][i] = 
			population[second_live][i];
		    }
		  }
		  j4++;
		}
	      break;
	    case 5:
	      /* JS Description: Multiple Point Simple Crossover
		 Applying the fifth operator, simple arithmetical crossover*/
	      if (j5 < (int) P5/2)
		{
		  /*Find two distinct parents for crossover operator 5*/
		  same = TRUE;
		  SameCount=0;
		  while (same==TRUE) {
		    SameCount++;
		    
		    first_live  = find_parent(live,pop_size);
		    second_live = find_parent(live,pop_size);

		    if (SameCount >= (UniquePairs) ) 
		      break;

		    for(i=1; i<=nvars; i++)
		      if ((int) population[first_live][i] != (int) population[second_live][i])
			same = FALSE;
		  } /* end of while same==TRUE loop */
		  /* check that agents to replace are in range */
		  if (die_now < 3) {
		    fprintf(output,"Not enough agents to be replaced\n");
		    exit(1);
		  }
		  live[first_live]--;
		  live[second_live]--;
		  first_die   = die_now-- ;
		  second_die  = die_now-- ;
		  new_genera[first_die][nvars+1]  = 5.0;
		  new_genera[second_die][nvars+1] = 5.0;
		  if (!same)
		    {
		      for(i=1; i<=nvars; i++)
			{
			  temp[1][i] = population[first_live][i];
			  temp[2][i] = population[second_live][i];
			}
		      JaIntegerOper5(temp[1],temp[2],STEP,domains,nvars);
		      for(i=1; i<=nvars; i++)
			{
			  new_genera[first_die][i]  = temp[1][i];
			  new_genera[second_die][i] = temp[2][i];
			}
		    }
		  else {
		    /* copy agent chosen twice into two new indivs */
		    for(i=1; i<=nvars; i++) {
		      new_genera[first_die][i]  = 
			      population[first_live][i];
		      new_genera[second_die][i] = 
			population[second_live][i];
		    }
		  }
		  j5++;
		}
	      break;
	    case 6:
	      /* JS Description: Whole Non-uniform Mutation */
	      /*Applying the sixth operator, whole non-uniform mutation*/
                    if (j6 < P6)
                      {
                        /*Find one parent for mutation operator 6*/
                        first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  fprintf(output, "No agents to be replaced\n");
			  exit(1);
			}
			
                        new_genera[die_now][nvars+1] = 6.0;
                        for(i=1; i<=nvars; i++)
                          t_vec[i] = population[first_live][i];
                        JaIntegerOper6(t_vec,domains,nvars,MaxGenerations,count_gener,B);
                        for(i=1; i<=nvars; i++)
                          new_genera[die_now][i] = t_vec[i];
			die_now--;
                        j6++;
                      }
                    break;
	    case 7:
	      /* JS Description: Heuristic Crossover */
	      /*Applying the seventh operator*/
	      if (j7 < (int) P7/2)
		{
		  /*Find two distinct parents for operator 7*/
		  same = TRUE;
		  SameCount=0;
		  while (same==TRUE) {
		    SameCount++;
		    first_live  = find_parent(live,pop_size);
		    second_live = find_parent(live,pop_size);

		    if (SameCount >= (UniquePairs) ) 
		      break;

		    for(i=1; i<=nvars; i++)
		      if (population[first_live][i] != population[second_live][i])
			same = FALSE;
		  } /* end of while same==TRUE loop */
		  /* check that agents to replace are in range */
		  if (die_now < 3) {
		    fprintf(output,"Not enough agents to be replaced\n");
		    exit(1);
		  }
		  live[first_live]--;
		  live[second_live]--;
		  first_die   = die_now-- ;
		  second_die  = die_now-- ;
		  new_genera[first_die][nvars+1]  = 7.0;
		  new_genera[second_die][nvars+1] = 7.0;
		  if (!same) {
		    if (first_live < second_live)
		      /* first agent is better agent */
		      for(i=1; i<=nvars; i++) {
			temp[2][i] = population[first_live][i];
			temp[1][i] = population[second_live][i];
			    }
		    else
		      /* second agent is better agent */
		      for(i=1; i<=nvars; i++) {
			temp[2][i] = population[second_live][i];
			temp[1][i] = population[first_live][i];
		      }
		    JaIntegerOper7(temp[1],temp[2],domains,nvars);
		    for(i=1; i<=nvars; i++)
		      new_genera[first_die][i]  = temp[1][i];
		    if (first_live < second_live)
		      /* first agent is better agent */
		      for(i=1; i<=nvars; i++) {
			temp[2][i] = population[first_live][i];
			temp[1][i] = population[second_live][i];
			    }
		    else
		      /* second agent is better agent */
		      for(i=1; i<=nvars; i++) {
			temp[2][i] = population[second_live][i];
			temp[1][i] = population[first_live][i];
		      }
		    JaIntegerOper7(temp[1],temp[2],domains,nvars);
		    for(i=1; i<=nvars; i++)
		      new_genera[second_die][i]  = temp[1][i];
		  }
		  else {
		    /* copy agent chosen twice into two new indivs */
		    for(i=1; i<=nvars; i++) {
			    new_genera[first_die][i]  = 
			      population[first_live][i];
			    new_genera[second_die][i] = 
			      population[second_live][i];
		    }
		  }
		  j7++;
		}
	    }
	}
      
      /*Replace the population with the new generation */
      Jnew = new_genera;
      new_genera = population;
      population = Jnew;

      if (Structure->DynamicPopulation==1 || Structure->DynamicPopulation==2)
	{
	  FILE *DynamicInput;
	  double tmp;

	  if((DynamicInput = fopen(Structure->DynamicPopulationPath, "r")) == NULL) {
	    fprintf(output,"WARNING: Unable to open the DynamicPopulationPath: %s\n", 
		    Structure->DynamicPopulationPath);
	    fprintf(output,"         Continuing the process.\n");

	    Structure->DynamicPopulation=0;
	    break;
	  }

	  fprintf(output,"\nDynamically Reading in Individuals from: %s\n", Structure->DynamicPopulationPath);
	  
	  j = pop_size;
	  i=0;
	  while (fscanf(DynamicInput,"%lf", &tmp)==1)
	    {
	      i++;
	      if (i>nvars)
		{
		  j--;
		  if (j<2)
		    {
		      fprintf(output,
			      "\nWARNING: Dynamic Population input includes more individuals than the population size.\n");
		      fprintf(output,
			      "         Keeping the current best individual and discarding extra individuals.\n");
		      fprintf(output,
			      "         See file: %s\n\n", Structure->DynamicPopulationPath);
		      break;
		    }
		  i = 1;
		} // end of if

	      population[j][i] = (int) tmp;
	      population[j][nvars+1] = 10.0;

	      if (Structure->Debug==1)
		{
		  fprintf(output,"\nDEBUG: WE READ IN THE FOLLOWING INDIVIDUALS:\n");
		  fprintf(output,"population[%d][%d]: %e\n", j, i, population[j][i]);
		}
	    } // end of while loop
	  
	  fclose(DynamicInput);

	  fprintf(output,"   Read in %d Individuals\n\n", pop_size-j+1);

	  Structure->DynamicPopulation=0;
	} // end of DynamicPopulation

      if (Structure->MemoryUsage==1)
	{
	  /* BINARY SEARCH.  (Knuth 3:409). */
	  /* We have sorted Memory by the "key" */
	  /* JaIntegerSort(population, pop_size, nvars+2);  */
	  JaIntegerSort(population, pop_size, nvars+2); 
	  JaIntegerSort(Memory, UniqueCount, nvars+2);

	  OldUniqueCount=UniqueCount;

	  JaIntMemoryMatrix(Structure,
			    Memory, population, X,
			    &UniqueCount, OldUniqueCount, pop_size, nvars, output, Status);

	  if (*Status < 0) {
	    fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
	      
	    /* free memory */
	    if (Structure->MemoryUsage==1)
	      JaMatrixFree(Memory, MemorySize);
	    
	    /* free populationstats stuff */
	    free(mean);
	    free(var);
	    free(skew);
	    free(kur);
	    free(tobs);
	      
	    free(bfgsoutX);
	    free(finalhessin);
	    free(evalX);
	    free(grad);
	      
	    /* free numeric.c allocations */
	    JaMatrixFree(population, pop_size+2);
	    JaMatrixFree(new_genera, pop_size+2);
	      
	    free_matrix(temp, 1, 2, 0);
	    free_vector(probab, 1);
	    free_vector(t_vec, 1);
	    free_vector(cum_probab, 1);
	    free_ivector(live, 1);
	      
	    return(ERROR_CODE);
	  }
	  if ( (UniqueCount+pop_size) >= MemorySize )
	    {
	      Structure->MemoryUsage=0;
	      fprintf(output,"\nWARNING: Turning Off MemoryMatrix because memory usage is too great.\n\n");
	    } /* end of if */
	} // end of Memory based evaluation
        else
	  {
	    NetworkNumber=0;
	    for (i=1; i<=pop_size; i++) 
	      {

		if (i > 1 && Structure->DynamicPopulation==2)
		{
		  JaDynamicPopulationCheck(Structure, population, i, pop_size, nvars, output, Status);
		  if (*Status < 0) {
		    fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
		
		      /* free memory */
		
		      /* free populationstats stuff */
		    free(mean);
		    free(var);
		    free(skew);
		    free(kur);
		    free(tobs);
		
		    free(bfgsoutX);
		    free(finalhessin);
		    free(evalX);
		    free(grad);
		
		    /* free numeric.c allocations */
		    JaMatrixFree(population, pop_size+2);
		    JaMatrixFree(new_genera, pop_size+2);
		
		    free_matrix(temp, 1, 2, 0);
		    free_vector(probab, 1);
		    free_vector(t_vec, 1);
		    free_vector(cum_probab, 1);
		    free_ivector(live, 1);
		
		    return(ERROR_CODE);
		  } // end of Status < 0
		} // end of DynamicPopulation==2
	  
		if (population[i][nvars+1]!=0)
		  {
		  if (Structure->Network==1)
		    {
		      NetworkNumber++;
		      population[i][0] = EVALUATE;
		    }
		  else
		    {
		      for(j=1; j<=nvars; j++)
			X[j] = population[i][j];
		      
		      population[i][0] = evaluate(Structure->AgentFit, X, nvars, Status);

		      if (*Status < 0) {
			fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
			
			/* free memory */
			
			/* free populationstats stuff */
			free(mean);
			free(var);
			free(skew);
			free(kur);
			free(tobs);
			
			free(bfgsoutX);
			free(finalhessin);
			free(evalX);
			free(grad);
			
			/* free numeric.c allocations */
			JaMatrixFree(population, pop_size+2);
			JaMatrixFree(new_genera, pop_size+2);
			
			free_matrix(temp, 1, 2, 0);
			free_vector(probab, 1);
			free_vector(t_vec, 1);
			free_vector(cum_probab, 1);
			free_ivector(live, 1);
			
			return(ERROR_CODE);
		      }	      

		    }
		  }
	      } //end of i loop
	    
	    if (Structure->Network==1)
	      {
		NetworkEvaluate(Structure->AgentFit, Structure->DBname, Structure->AgentName, 
				population, pop_size, nvars, NetworkNumber, Status, 10);
		
		if (*Status < 0) {
		  fprintf(output,"EVALUATE.C: Elegantly exiting GENOUD because of exit (status) code: %d\n", *Status);
		  
		  /* free memory */
		  
		  /* free populationstats stuff */
		  free(mean);
		  free(var);
		  free(skew);
		  free(kur);
		  free(tobs);
		  
		  free(bfgsoutX);
		  free(finalhessin);
		  free(evalX);
		  free(grad);
		  
		  /* free numeric.c allocations */
		  JaMatrixFree(population, pop_size+2);
		  JaMatrixFree(new_genera, pop_size+2);
		  
		  free_matrix(temp, 1, 2, 0);
		  free_vector(probab, 1);
		  free_vector(t_vec, 1);
		  free_vector(cum_probab, 1);
		  free_ivector(live, 1);
		  
		  return(ERROR_CODE);
		}	      
	      }
	  } // end of default evaluation

      /*Sort the new population based on their evaluation function*/
      sort(MinMax,population,pop_size,0);

      switch(MinMax)
        {
	case 0:
	  if(Teval > population[1][0])
	    {
	      Teval = population[1][0];
	      fprintf(output,"%7lu \t%e\n",
		      count_gener,population[1][0]); 
	      fflush(output);
	      peak_cnt = count_gener;
	      peak_val = population[1][0];
	    }
	  break;
          case 1:
            if(Teval < population[1][0])
              {
                Teval = population[1][0];
                fprintf(output,"%7lu \t%e\n",
			count_gener,population[1][0]);
		fflush(output);
                peak_cnt = count_gener;
                peak_val = population[1][0];
              }
            break;
        }

      /* compute and print mean and variance of population */
      if (PrintLevel==2) {
	  fprintf(output,"\nGENERATION: %d\n", count_gener);
	  populationstats(population, pop_size, nvars, mean, var, skew, kur, tobs);
	  for (i=0; i<=nvars; i++) {
	      if (i==0) {
		  fprintf(output, "Fitness Value... %e\n", population[1][i]);
		  fprintf(output, "mean............ %e\n", mean[i]);
		  fprintf(output, "var............. %e\n", var[i]);
		  fprintf(output, "skewness........ %e\n", skew[i]);
		  fprintf(output, "kurtosis........ %e\n", kur[i]);
		  fprintf(output, "#null........... %d\n", pop_size-tobs[i]);
		  if(Structure->MemoryUsage==1)
		    fprintf(output, "#unique......... %d, #Total UniqueCount: %d\n", 
			    UniqueCount-OldUniqueCount, UniqueCount);
		  fprintf(output, "tobs............ %d\n", tobs[i]);
	      }
	      else {
		  fprintf(output, "var %d:\n", i);
		  fprintf(output, "best............ %e\n", population[1][i]);
		  fprintf(output, "mean............ %e\n", mean[i]);
		  fprintf(output, "var............. %e\n", var[i]);
		  fprintf(output, "skewness........ %e\n", skew[i]);
		  fprintf(output, "kurtosis........ %e\n", kur[i]);
		  fprintf(output, "#null........... %d\n", pop_size-tobs[i]);
		  fprintf(output, "tobs............ %d\n", tobs[i]);
	      }
	  }
      } /* end of printlevel if */
      
      if (PrintLevel==1) {
	popmean = popvar = 0.0 ;
	popwrk = 1.0 / pop_size ;
	for(i=1; i<=pop_size; i++) {
	  popmean += population[i][0] ;
	}
	popmean *= popwrk ;
	for(i=1; i<=pop_size; i++) {
	  popstat =  population[i][0] - popmean ;
	  popvar += (popstat*popstat) ;
	}
	popvar *= popwrk ;
	fprintf(output, "   mean = %e, variance = %e\n\n", popmean, popvar);
      }

      fflush(output);
	
      /* Print the population file */
      if ( PrintLevel == 1 ) {
	if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
	  fprintf(output,"Unable to open the project file: %s", 
		  Structure->ProjectPath);

	  /* free populationstats stuff */
	  free(mean);
	  free(var);
	  free(skew);
	  free(kur);
	  free(tobs);
	  
	  free(bfgsoutX);
	  free(finalhessin);
	  free(evalX);
	  free(grad);
	  
	  /* free numeric.c allocations */
	  if (Structure->MemoryUsage==1)
	    JaMatrixFree(Memory, MemorySize);
	  JaMatrixFree(population, pop_size+2);
	  JaMatrixFree(new_genera, pop_size+2);
	  
	  free_matrix(temp, 1, 2, 0);
	  free_vector(probab, 1);
	  free_vector(t_vec, 1);
	  free_vector(cum_probab, 1);
	  free_ivector(live, 1);

	  return(ERROR_CODE);
	}
	print_population(pop_size, nvars, count_gener, population, popout);
	fclose(popout);
      } /* end of PrintLevel if */
      if ( PrintLevel == 2 ) {
	if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
	  fprintf(output,"Unable to open the project file: %s", 
		  Structure->ProjectPath);
	  
	  /* free populationstats stuff */
	  free(mean);
	  free(var);
	  free(skew);
	  free(kur);
	  free(tobs);
	  
	  free(bfgsoutX);
	  free(finalhessin);
	  free(evalX);
	  free(grad);
	  
	  /* free numeric.c allocations */
	  if (Structure->MemoryUsage==1)
	    JaMatrixFree(Memory, MemorySize);
	  JaMatrixFree(population, pop_size+2);
	  JaMatrixFree(new_genera, pop_size+2);
	  
	  free_matrix(temp, 1, 2, 0);
	  free_vector(probab, 1);
	  free_vector(t_vec, 1);
	  free_vector(cum_probab, 1);
	  free_ivector(live, 1);
	  
	  return(ERROR_CODE);
	}
	print_population(pop_size, nvars, count_gener, population, popout);
	fflush(popout);
	fclose(popout);
      }
      
      if (count_gener == 1) {
	oldfitvalue=population[1][0];
      }
      
      /*
	if (oldfitvalue == population[1][0]) {
	nochange_gen++;
	}
      */
      
      switch(MinMax)
	{
	case 0:
	  if (oldfitvalue - SolutionTolerance > population[1][0]) {
	    nochange_gen=0;
	    oldfitvalue=population[1][0];
	  }
	  else nochange_gen++;
	  break;
	case 1:
	  if (oldfitvalue + SolutionTolerance < population[1][0]) {
	    nochange_gen=0;
	    oldfitvalue=population[1][0];
	  }
	  else nochange_gen++;	      
	  break;
	}
      
      if (nochange_gen > (WaitGenerations)) {
	/* increase the number of WaitGenerations if the gradients are NOT zero! */	  
	if (GradientCheck==0) {
	  fprintf(output,"\nSoft Generation Wait Limit Hit.\n");
	  fprintf(output,"No Improvement in %d Generations\n", nochange_gen-1);
	  fflush(output);
	  MaxGenerations = 0;
	  nochange_gen=0;
	}
	else  {
	  for (i=1; i<=nvars; i++)
	    {
		  bfgsoutX[i-1]=population[1][i];
	    }
	  gradient(Structure->AgentFit, bfgsoutX, grad, nvars, MinMax, BoundaryEnforcement, InstanceNumber, 
		   domains, Status);
	  if (*Status < 0) {
	    /* Free Memory */
	    if (Structure->MemoryUsage==1)
	      JaMatrixFree(Memory, MemorySize);

	    /* free populationstats stuff */
	    free(mean);
	    free(var);
	    free(skew);
	    free(kur);
	    free(tobs);

	    free(bfgsoutX);
	    free(finalhessin);
	    free(evalX);
	    free(grad);
	    
	    /* free numeric.c allocations */
	    JaMatrixFree(population, pop_size+2); 
	    JaMatrixFree(new_genera,pop_size+2); 

	    free_matrix(temp, 1, 2, 0);
	    free_vector(probab, 1);
	    free_vector(t_vec, 1);
	    free_vector(cum_probab, 1);
	    free_ivector(live, 1);
	    
	    return(ERROR_CODE);
	  }
	  GradientTrigger = 0;
	  for (i=0; i<nvars; i++) {
	    if (fabs(grad[i]) > SolutionTolerance) {
	      GradientTrigger = 1;
	      break;
		}
	  } /* end for loop */
	  if (GradientTrigger==1) {
	    IncreaseGenerations = WaitGenerations;
	    WaitGenerations += IncreaseGenerations;
		fprintf(output,
			"\nDoubling Soft Maximum Wait Generation Limit to %d (from %d).\n", 
			WaitGenerations, IncreaseGenerations);
		fprintf(output,"I'm doing this because at least one gradient is too large.\n");
		fprintf(output,"G[%d]: %e\t Solution Tolerance: %e\n\n", 
			i+1, grad[i], SolutionTolerance);
	  }
	  else {
	    fprintf(output,"\nSoft Generation Wait Limit Hit.\n");
	    fprintf(output,"No Improvement in %d Generations\n", nochange_gen-1);
		fflush(output);
		MaxGenerations = 0;
		nochange_gen=0;
	  }
	}/* end if loop */
      }
      
      if ( (count_gener == MaxGenerations) && (GradientTrigger==1) ) 
	{
	  if (HardGenerationLimit==0)
	    {
	      IncreaseGenerations = MaxGenerations;
	      MaxGenerations += IncreaseGenerations;
	      fprintf(output,
		      "\nIncreasing Soft Maximum Generation Limit by %d (MaxGenerations) to %d.\n", 
		      IncreaseGenerations, MaxGenerations);
	      fprintf(output,"I'm doing this because at least one gradient is too large.\n\n");
	    } // if (Structure->HardGenerationLimit==0)
	  else
	    {
	      fprintf(output,"\nSTOPPING: HARD MAXIMUM GENERATION LIMIT HIT\n");
	      fprintf(output,"          At least one gradient is still too large\n");
	    } // else
	} // if ( (count_gener == MaxGenerations) && (GradientTrigger==1) ) 


      /* increase the number of generations if fitness has been improving */
      if ( (count_gener == MaxGenerations) &&  (nochange_gen < WaitGenerations) ) {
	if (HardGenerationLimit==0)
	  {
	    if (WaitGenerations > MaxGenerations) {
	      IncreaseGenerations = WaitGenerations;
	      MaxGenerations += (int) (IncreaseGenerations);
	      fprintf(output,
		      "\nIncreasing Soft Maximum Generation Limit by %d (WaitGenerations) to %d\n", 
		      IncreaseGenerations, MaxGenerations);
	      fprintf(output,"I'm doing this because the fitness is still impoving.\n\n");
	    }
	    else {
	      IncreaseGenerations = MaxGenerations;
	      MaxGenerations += (int) (IncreaseGenerations);
	      fprintf(output,
		      "\nIncreasing Soft Maximum Generation Limit by %d (MaxGenerations) to %d.\n", 
		      IncreaseGenerations, MaxGenerations);
	      fprintf(output,"I'm doing this because the fitness is still improving.\n\n");
	    }
	  } // if (Structure->HardGenerationLimit==0)
	else
	  {
	    fprintf(output,"\nSTOPPING: HARD MAXIMUM GENERATION LIMIT HIT\n");
	    fprintf(output,"          But fitness is still improving\n");
	  }
      } // if ( (count_gener == MaxGenerations) &&  (nochange_gen < WaitGenerations) )
      
      fflush(output);

      /* Should we recheck the main data structure for changes to the operator set? If so, let's
       do it now */
      if (Structure->AllowDynamicUpdating==1) {
	pop_size_old = pop_size;
	fprintf(output,"\nUpdating Main Data Structure:\n");
	SetRunTimeParameters(Structure, 0,
			     &pop_size, &nvars, &MaxGenerations, &WaitGenerations,
			     &MinMax, &GradientCheck, &BoundaryEnforcement, &UseBFGS, &SolutionTolerance,
			     &InstanceNumber, &P, &P0, &P1, &P2, &P3, &P4, &P5, &P6, &P7, &P8, 
			     &PrintLevel, &HardGenerationLimit, output);
	if (pop_size > pop_size_old) {
	  population_old = JaMatrixAllocate(pop_size_old+2, nvars+2);
	  
	  for (i=1; i<=pop_size_old;i++) {
	    for (j=0; j<=nvars; j++) {
	      population_old[i][j] = population[i][j];
	    }
	  }
	  JaMatrixFree(population, pop_size_old+2);
	  population    = JaMatrixAllocate(pop_size+2, nvars+2);

	  for (i=1; i<=pop_size_old;i++) {
	    for (j=0; j<=nvars; j++) {
	      population[i][j] = population_old[i][j];
	    }
	  }	  
	  JaMatrixFree(population_old, pop_size_old+2);

	  /* we need to add individuals to population! */
	  for (j=(pop_size_old+1); j<=pop_size; j++) {
	    for (i=1; i<=nvars; i++) {
	      population[j][i] = frange_ran(domains[i][1], domains[i][3]); 
	      population[j][nvars+1] = 1.0;
	    }
	  }
	} /* end of if popsize > pop_size_old */
	
	JaMatrixFree(new_genera, pop_size_old+2);
	free_vector(probab, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);

	new_genera    = JaMatrixAllocate(pop_size+2, nvars+2);

	temp       = matrix(1,2,0,nvars);
	probab     = Gvector(1,pop_size);
	t_vec      = Gvector(1,nvars);
	cum_probab = Gvector(1,pop_size);
	live       = ivector(1,pop_size);

	Structure->AllowDynamicUpdating=0;
      } // end of if AllowDynamicUpdating==1
      
    } /* end of do loop */
  /*Increment iteration count and test whether all generations are done*/
  while (++count_gener <= MaxGenerations);
  
  fprintf(output,"\nBest Fit Found at Generation %lu\nFit Value = %e\n",peak_cnt,peak_val);
  fprintf(output,"\n\nParameters at the Solution (value, gradient):\n\n");

  /* output data structure */
  Structure->oPeakGeneration=peak_cnt;
  Structure->oGenerations=count_gener-1;

  /* obtain gradients */
  if (GradientCheck==0 && UseBFGS==0) {
    fprintf(output,"\nNot Obtaining Gradient Information\n");
    for (i=0; i< nvars; i++) {
      grad[i]=-1.0;
    }
  }
  else {
    for (i=1; i<=nvars; i++)
      {
	bfgsoutX[i-1]=population[1][i];
      }
    gradient(Structure->AgentFit, bfgsoutX, grad, nvars, MinMax, BoundaryEnforcement, 
	     InstanceNumber, domains, Status);
    if (*Status < 0) {
	/* Free Memory */
      if (Structure->MemoryUsage==1)
	JaMatrixFree(Memory, MemorySize);
	
	/* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);
	
	
	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	
	/* free numeric.c allocations */
	JaMatrixFree(population, pop_size+2); 
	JaMatrixFree(new_genera,pop_size+2); 
	
	free_matrix(temp, 1, 2, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);
	
	return(ERROR_CODE);
    }
  }
  
  /* print best solution */
  for(j = 1; j <= nvars; j++) {
    i = j-1;
    fprintf(output," X[%2d] :\t%e\tG[%2d] :\t%e\n",j,population[1][j],j,grad[i]);
    Results[i] = population[1][j];
    Gradients[i] = grad[i];
  }
  
  /* free memory */
  if (Structure->MemoryUsage==1)
    JaMatrixFree(Memory, MemorySize); 

  /* free populationstats stuff */
  free(mean);
  free(var);
  free(skew);
  free(kur);
  free(tobs);
  
  
  free(bfgsoutX);
  free(finalhessin);
  free(evalX);
  free(grad);

  /* free numeric.c allocations */
  /* free_matrix(population, 0, pop_size+1,0); */
  JaMatrixFree(population, pop_size+2); 
  JaMatrixFree(new_genera,pop_size+2); 

  free_matrix(temp, 1, 2, 0);
  free_vector(probab, 1);
  free_vector(t_vec, 1);
  free_vector(cum_probab, 1);
  free_ivector(live, 1);

  return(peak_val);

} /* end of JaIntegerOptimization */


/********************************************************************************/
/*  JaIntegerSort():                                                            */
/*                                                                              */
/*  This function sorts a double** on an integer basis.                         */
/*  The function also assumes that the double** is indexed from 1 in its rows   */
/*  and from zero in its columns.                                               */
/*                                                                              */
/********************************************************************************/

void JaIntegerSort(double **InMatrix, long n, long k)
{
  /* extern int JaIntegerCMP(); */
  long i, j;
  double **Tmp;
  extern long Gnvars[MAXINSTANCES];
  long nvars;

  Tmp = JaMatrixAllocate(n, k);

  nvars=Gnvars[ExternStructure->InstanceNumber];

  for (i=1; i<=n; i++) {
    for (j=0; j<k; j++) {
      Tmp[i-1][j] = InMatrix[i][j];
    }
  }

#ifdef MS_WINDOWS 
     qsort(Tmp, n, sizeof(double *), 
	   (int (__cdecl *)(const void *,const void *)) JaIntegerCMP);
#else
     qsort(Tmp, n, sizeof(double *), 
	   (int (*)(const void *, const void *)) JaIntegerCMP);
#endif

  for (i=1; i<=n; i++) {
    for (j=0; j<k; j++) {
      InMatrix[i][j] = Tmp[i-1][j];
    }
  }

  JaMatrixFree(Tmp, n);
} /* end of JaIntegerSort */


/********************************************************************************/
/*  JaDoubleSort():                                                             */
/*                                                                              */
/*  This function sorts a double** on an double  basis.                         */
/*  The function also assumes that the double** is indexed from 1 in its rows   */
/*  and from zero in its columns.                                               */
/*                                                                              */
/********************************************************************************/

void JaDoubleSort(double **InMatrix, long n, long k)
{
  /* extern int JaDoubleCMP(); */
  long i, j;
  double **Tmp;
  extern long Gnvars[MAXINSTANCES];
  long nvars;

  Tmp = JaMatrixAllocate(n, k);

  nvars=Gnvars[ExternStructure->InstanceNumber];

  for (i=1; i<=n; i++) {
    for (j=0; j<k; j++) {
      Tmp[i-1][j] = InMatrix[i][j];
    }
  }

#ifdef MS_WINDOWS 
     qsort(Tmp, n, sizeof(double *), 
	   (int (__cdecl *)(const void *,const void *)) JaDoubleCMP);
#else
     qsort(Tmp, n, sizeof(double *), 
	   (int (*)(const void *, const void *)) JaDoubleCMP);
#endif
  for (i=1; i<=n; i++) {
    for (j=0; j<k; j++) {
      InMatrix[i][j] = Tmp[i-1][j];
    }
  }

  JaMatrixFree(Tmp, n);
} /* end of JaDoubleSort */

void JaDoubleMemoryMatrix_Gen0(struct GND_IOstructure *Structure, 
			       double **Memory, double **population, double *X,
			       long *UniqueCount, long OldUniqueCount,
			       int pop_size, int nvars, 
			       FILE *output, long *Status)
{
  int i, j;
  FLAG UniqueFlag;

  /* search over i in population */
  for(i=1; i<=pop_size; i++) {

    if (Structure->DynamicPopulation==2)
      {
	JaDynamicPopulationCheck(Structure, population, i, pop_size, nvars, output, Status);
	if (*Status < 0) {
	  return;
	} // end of Status < 0
      } // end of DynamicPopulation==2
    
    UniqueFlag=FALSE;
    if (i>1) {
      /* 
	 Take only unique people out of the population[][] matrix 
      */
      for (j=1; j<=nvars; j++) {
	/* Unique in this population? */
	if (population[i][j] != population[i-1][j]) {
	  UniqueFlag=TRUE;
	  break;
	}
      } /* end of j loop */
    } /* end of if */
    else UniqueFlag=TRUE;
    
    if (UniqueFlag) {
      ++*UniqueCount; /* *UniqueCount counts from 1 */
	  
      for (j=0; j<=nvars; j++) 
	Memory[*UniqueCount][j] = population[i][j];

      if (population[i][nvars+1]==-1.0 || population[i][nvars+1]==11.0) {
	for(j=1; j<=nvars; j++)
	  X[j] = Memory[*UniqueCount][j];
	  
	Memory[*UniqueCount][0] = evaluate(Structure->AgentFit, X, nvars, Status);
	if (*Status < 0) {
	  return;
	}
	Memory[*UniqueCount][nvars+1] = 0.0;
	    
	population[i][0] = Memory[*UniqueCount][0];
      } /* end of if */
      else { 
	Memory[*UniqueCount][0] = population[i][0];
	Memory[*UniqueCount][nvars+1] = 0.0;
      }
    } /* end of main UniqueFlag */
    else {
      for (j=0; j<=nvars+1; j++) 
	population[i][j] = population[i-1][j];
    } /* end of else */
  } /* end of i loop */  

} //end of JaDoubleMemoryMatrix_Gen0

void JaDoubleMemoryMatrix(struct GND_IOstructure *Structure, 
			  double **Memory, double **population, double *X,
			  long *UniqueCount, long OldUniqueCount,
			  int pop_size, int nvars, FILE *output, long *Status)
{
  int i, j;
  FLAG UniqueFlag, Redundant;
  long upper, lower, midpoint;

  /* search over i in population */
  for(i=1; i<=pop_size; i++) {

    if (i > 1 && Structure->DynamicPopulation==2)
      {
	JaDynamicPopulationCheck(Structure, population, i, pop_size, nvars, output, Status);
	if (*Status < 0) {
	  return;
	} // end of Status < 0
      } // end of DynamicPopulation==2

    UniqueFlag=FALSE;
    Redundant=FALSE;
    if (population[i][nvars+1]==0) {
      Redundant=TRUE;
    }
    else if (i>1) {
      /* 
	 Take only unique people out of the population[][] matrix 
      */
      for (j=1; j<=nvars; j++) {
	/* the integer cast here is required for the integer version! */
	/* Unique in this population? */
	if (population[i][j] != population[i-1][j]) {
	  UniqueFlag=TRUE;
	  break;
	}
      } /* end of j loop */
    } /* end of if */
    else UniqueFlag=TRUE;
	    
    if (UniqueFlag) {
      /* B1 initialize the upper and lower bounds */
      upper=OldUniqueCount;
      lower=1;
	      
      while(!(upper < lower)) {
	UniqueFlag=FALSE;
	/* B2 obtain the approximate midpoint between upper and lower */
	midpoint = upper+lower;
	midpoint= (int) midpoint / 2;
		
	/* B3 is k < k_{i} or k > k_{i} or is k==k_{i}? */
	for (j=1; j<=nvars; j++) {
	  /* B4 */
	  if (population[i][j] < Memory[midpoint][j]) {
	    upper=midpoint-1;
	    UniqueFlag=TRUE;
	  } /* end of if */
	  /* B4 */
	  else if (population[i][j] > Memory[midpoint][j]) {
	    lower=midpoint+1;
	    UniqueFlag=TRUE;
	  } /* end of if */
	  if (UniqueFlag) break;
	} /* end of j loop */
	if (!UniqueFlag) {
	  population[i][0] = Memory[midpoint][0];
	  population[i][nvars+1] = 0.0;
	  upper=0;
	}
      } /* end of midpoint while */
	      
      /* we have a unique individual */
      if (UniqueFlag) {
	++*UniqueCount; /* *UniqueCount counts from 1 */
		
	for (j=0; j<=nvars; j++) 
	  Memory[*UniqueCount][j] = population[i][j];
		
	for(j=1; j<=nvars; j++)
	  X[j] = Memory[*UniqueCount][j];
		
	Memory[*UniqueCount][0] = evaluate(Structure->AgentFit, X, nvars, Status);
		
	if (*Status < 0) {
	  return;
	}
		
	Memory[*UniqueCount][nvars+1] = 0.0;
		
	for (j=0; j<=nvars+1; j++) 
	  population[i][j] = Memory[*UniqueCount][j];
		
	/* This matrix indexes the position in the population matrix which is equal to the 
	   unique individual in the population matrix */
	/* fprintf(output,"!match: i: %d, *UniqueCount: %d\n", i, *UniqueCount); */	  
      } /* end of if */
    } /* end of main UniqueFlag */
    else if (!Redundant) {
      for (j=0; j<=nvars+1; j++) 
	population[i][j] = population[i-1][j];
    } /* end of else */
  } /* end of i loop */
      
} // end of JaDoubleMemoryMatrix


void JaIntMemoryMatrix_Gen0(struct GND_IOstructure *Structure, 
			       double **Memory, double **population, double *X,
			       long *UniqueCount, long OldUniqueCount,
			       int pop_size, int nvars, FILE *output, long *Status)
{
  int i, j;
  FLAG UniqueFlag;

  /* search over i in population */
  for(i=1; i<=pop_size; i++) {

    if (Structure->DynamicPopulation==2)
      {
	JaDynamicPopulationCheck(Structure, population, i, pop_size, nvars, output, Status);
	if (*Status < 0) {
	  return;
	} // end of Status < 0
      } // end of DynamicPopulation==2

    UniqueFlag=FALSE;
    if (i>1) {
      /* 
	 Take only unique people out of the population[][] matrix 
      */
      for (j=1; j<=nvars; j++) {
	/* the integer cast here is required for the integer version! */
	/* Unique in this population? */
	if ((int) population[i][j] != (int) population[i-1][j]) {
	  UniqueFlag=TRUE;
	  break;
	}
      } /* end of j loop */
    } /* end of if */
    else UniqueFlag=TRUE;
	
    if (UniqueFlag) {
      ++*UniqueCount; /* *UniqueCount counts from 1 */
      
      for (j=0; j<=nvars; j++) 
	Memory[*UniqueCount][j] = population[i][j];

      if (population[i][nvars+1]==-1.0 || population[i][nvars+1]==11.0) {
	for(j=1; j<=nvars; j++)
	  X[j] = Memory[*UniqueCount][j];
      
	Memory[*UniqueCount][0] = evaluate(Structure->AgentFit, X, nvars, Status);
	if (*Status < 0) {
	  return;
	}
	Memory[*UniqueCount][nvars+1] = 0.0;
	
	population[i][0] = Memory[*UniqueCount][0];
      } /* end of if */
      else { 
	Memory[*UniqueCount][0] = population[i][0];
	Memory[*UniqueCount][nvars+1] = 0.0;
      }
    } /* end of main UniqueFlag */
    else {
      for (j=0; j<=nvars+1; j++) 
	population[i][j] = population[i-1][j];
    } /* end of else */
  } /* end of i loop */

} //end of JaIntMemoryMatrix_Gen0


void JaIntMemoryMatrix(struct GND_IOstructure *Structure, 
			  double **Memory, double **population, double *X,
			  long *UniqueCount, long OldUniqueCount,
			  int pop_size, int nvars, FILE *output, long *Status)
{
  int i, j;
  FLAG UniqueFlag, Redundant;
  long upper, lower, midpoint;

  /* search over i in population */
  for(i=1; i<=pop_size; i++) {

    if (i > 1 && Structure->DynamicPopulation==2)
      {
	JaDynamicPopulationCheck(Structure, population, i, pop_size, nvars, output, Status);
	if (*Status < 0) {
	  return;
	} // end of Status < 0
      } // end of DynamicPopulation==2

    UniqueFlag=FALSE;
    Redundant=FALSE;
    if (population[i][nvars+1]==0) {
      Redundant=TRUE;
    }
    else if (i>1) {
      /* 
	 Take only unique people out of the population[][] matrix 
      */
      for (j=1; j<=nvars; j++) {
	/* the integer cast here is required for the integer version! */
	/* Unique in this population? */
	if ((int) population[i][j] != (int) population[i-1][j]) {
	  UniqueFlag=TRUE;
	  break;
	}
      } /* end of j loop */
    } /* end of if */
    else UniqueFlag=TRUE;
	
    if (UniqueFlag) {
      /* B1 initialize the upper and lower bounds */
      upper=OldUniqueCount;
      lower=1;
	  
      while(!(upper < lower)) {
	UniqueFlag=FALSE;
	/* B2 obtain the approximate midpoint between upper and lower */
	midpoint = upper+lower;
	midpoint= (int) midpoint / 2;
	    
	/* B3 is k < k_{i} or k > k_{i} or is k==k_{i}? */
	for (j=1; j<=nvars; j++) {
	  /* B4 */
	  if ((int) population[i][j] < (int) Memory[midpoint][j]) {
	    upper=midpoint-1;
	    UniqueFlag=TRUE;
	  } /* end of if */
	  /* B4 */
	  else if ((int) population[i][j] > (int) Memory[midpoint][j]) {
	    lower=midpoint+1;
	    UniqueFlag=TRUE;
	  } /* end of if */
	  if (UniqueFlag) break;
	} /* end of j loop */
	if (!UniqueFlag) {
	  population[i][0] = Memory[midpoint][0];
	  population[i][nvars+1] = 0.0;
	  upper=0;
	}
      } /* end of midpoint while */
	  
      /* we have a unique individual */
      if (UniqueFlag) {
	++*UniqueCount; /* *UniqueCount counts from 1 */
	    
	for (j=0; j<=nvars; j++) 
	  Memory[*UniqueCount][j] = (int) population[i][j];
	    
	for(j=1; j<=nvars; j++)
	  X[j] = Memory[*UniqueCount][j];
	    
	Memory[*UniqueCount][0] = evaluate(Structure->AgentFit, X, nvars, Status);
	    
	if (*Status < 0) {
	  return;
	}
	    
	Memory[*UniqueCount][nvars+1] = 0.0;
	    
	for (j=0; j<=nvars+1; j++) 
	  population[i][j] = Memory[*UniqueCount][j];
	    
	/* This matrix indexes the position in the population matrix which is equal to the 
	   unique individual in the population matrix */
	/* fprintf(output,"!match: i: %d, *UniqueCount: %d\n", i, *UniqueCount); */	  
      } /* end of if */
    } /* end of main UniqueFlag */
    else if (!Redundant) {
      for (j=0; j<=nvars+1; j++) 
	population[i][j] = population[i-1][j];
    } /* end of else */
  } /* end of i loop */
      
} //end of JaIntMemoryMatrix

void JaDynamicPopulationCheck(struct GND_IOstructure *Structure,
				    double **population, int location, int pop_size, int nvars, 
				    FILE *output, long *Status)
{

  FILE *DynamicInput;
  double tmp;
  long i, j, indx, maximum;
  int Dobs;

  maximum = pop_size-location+1;


  if((DynamicInput = fopen(Structure->DynamicPopulationPath, "r")) == NULL) {
    fprintf(output,"WARNING: Unable to open the DynamicPopulationPath: %s\n", 
	    Structure->DynamicPopulationPath);
    fprintf(output,"         Continuing the process.\n");
    
    Structure->DynamicPopulation=0;
    return;
  }

  fprintf(output,"\nDynamically Reading in Individuals from: %s\n", Structure->DynamicPopulationPath);

  j=0;
  i=0;
  while (fscanf(DynamicInput,"%lf", &tmp)==1)
    {
      i++;
      if (i>nvars)
	{
	  j++;
	  if (j>maximum)
	    {
	      fprintf(output,
		      "\nWARNING: Dynamic Population input includes more individuals than can be added at this time.\n");
	      fprintf(output,
		      "         If you want more freedom in adding individuals into the population please use do not use the MemoryMatrix feature.\n");
	      break;
	    }
	  i = 1;
	} // end of if

      indx = location+j;
      population[indx][i] = tmp;
      population[indx][nvars+1] = 11.0;
      
      if (Structure->Debug==1)
	{
	  fprintf(output, "\nDEBUG: WE READ IN THE FOLLOWING INDIVIDUALS:\n");
	  fprintf(output, "population[%d][%d]: %e\n", indx, i, population[indx][i]);
	}
		    
    } // end of while loop
  fclose(DynamicInput);
		
  Dobs = j+1;
		
  fprintf(output,"   Read in %d Individuals\n\n", Dobs);		
		
  Structure->DynamicPopulation=0; 
} // end of DynamicPopulationCheck

