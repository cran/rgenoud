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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/operators.cpp,v 1.20 2002/11/06 02:12:02 jsekhon Exp $

*/

#include "genoud.h"

#ifndef OPTIM
extern double genoud_optim(double *X, int nvars);
#endif

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper1()                                      */
/*                                 Uniform Mutation                             */
/*                                                                              */
/*           SYNOPSIS          :   void oper1(parent,domains,nvars)             */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, uniform mutation.             */
/*                                                                              */
/*           FUNCTIONS CALLED  :   find_range(),                                */
/*                                 frange_ran(),                                */
/*                                 irange_ran(),                                */
/*                                 Gvector().                                   */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void oper1(VECTOR parent, double **domains, int nvars)
     /* VECTOR parent;  The parent vector*/
{
  int comp;
  double llim,ulim;/*Lower and Upper limits of the value to be mutated*/

  FLAG SAME;
  double tmp;
  long count;

  count=0;
  SAME=TRUE;
  while (SAME==TRUE) 
    {
      count++;

      comp = irange_ran(1,nvars);
      
      /*Finding the lower and upper limits between which the values are to be mutated*/
      find_range(&llim,&ulim,comp,domains,nvars,parent);
      
      /*Find a random value between the lower and the upper limits, to substitute*/
      /*for the old value*/
      tmp = frange_ran(llim,ulim);

      if ( parent[comp] != tmp)
	SAME=FALSE;
      else if (count >= MAX_OPER_UNIQUE_TRY)
	SAME=FALSE;
    } /* end of while */

  parent[comp] = tmp;
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper2()                                      */
/*                                 Boundary Mutatation                          */
/*                                 No Uniqueness checking here                  */
/*                                 Don't use this oper often!                   */
/*                                                                              */
/*           SYNOPSIS          :   void oper2(parent,fin_mat,nvars)             */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, boundary mutation.            */
/*                                                                              */
/*           FUNCTIONS CALLED  :   find_range(),                                */
/*                                 flip(),                                      */
/*                                 irange_ran(),                                */
/*                                 Gvector().                                   */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/


void oper2(VECTOR parent, double **domains, int nvars)
     /* VECTOR parent;  The parent vector*/
     /* MATRIX fin_mat; The final matrix*/
{
  int comp;
  double llim,ulim;/*Lower and Upper limits of the value to be mutated*/

  FLAG SAME;
  double tmp;
  long count;

  count=0;
  SAME=TRUE;
  while (SAME==TRUE) 
    {
      count++;
      
      comp = irange_ran(1,nvars);
      
      /*Finding the lower and upper limits between which the values are to be mutated*/
      find_range(&llim,&ulim,comp,domains,nvars,parent);
      
      /*Replace either the lower limit or the upper limit at random,*/
      /*for the old value*/
      tmp = (flip() == TAIL) ? llim : ulim;

      if ( tmp != parent[comp])
	SAME=FALSE;
      else if (count >= MAX_OPER_UNIQUE_TRY)
	SAME=FALSE;
    } /* end of while */
  
  parent[comp] = tmp;
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper3()                                      */
/*                                                                              */
/*           SYNOPSIS          :   void oper3(parent,fin_mat,r,c,T,t,B)         */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, non-uniform mutation.         */
/*                                                                              */
/*           FUNCTIONS CALLED  :   find_range(),                                */
/*                                 flip(),                                      */
/*                                 get_F(),                                     */
/*                                 irange_ran(),                                */
/*                                 Gvector().                                   */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void oper3(VECTOR parent, double **domains, int nvars, int T, int t, int B)
  /* VECTOR parent; */
  /* unsigned long T;   Total number of generations*/
  /* unsigned long t;   Current generation number*/
  /* int B; */
{
  int comp;
  double llim,ulim;

  FLAG SAME;
  double tmp;
  long count;

  count=0;

  SAME=TRUE;
  while (SAME==TRUE) {
    count++;

    comp = irange_ran(1,nvars);
    find_range(&llim,&ulim,comp,domains,nvars,parent);
    
    /*From the current value of the component to be mutated, chooose at random*/
    /*whether to mutate with a lesser value or a greater value*/
    /*Then find a value lesser or greater than the original value from the*/
    /*function get_f()*/
    tmp = (flip() == TAIL) ? parent[comp]-get_F(T,t,parent[comp]-llim,B) :
      parent[comp]+get_F(T,t,ulim-parent[comp],B);

    if ( parent[comp] != tmp)
      SAME=FALSE;
    else if (count >= MAX_OPER_UNIQUE_TRY)
      SAME=FALSE;
  } /* end of while */

  parent[comp] = tmp;
}



/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper4()                                      */
/*                                 Polytope Crossover                           */
/*                                                                              */
/*           SYNOPSIS          :   void oper4(p1,p2,x2_vari)                    */
/*                                                                              */
/*           DESCRIPTION       :   This function returns two new vectors        */
/*                                  generated after whole arithmetical          */
/*                                  crossover, from the two parent vectors.     */
/*                                                                              */
/*           FUNCTIONS CALLED  :   matrix()                                     */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void oper4(VECTOR p1, VECTOR p2, int nvars)
     /* VECTOR p1,p2;  The two parents chosen for crossover*/
     /* int nvars;   Length of the vector*/
{
  double **child;
  long    i;
  double  A;

  FLAG SAME;
  long count, tcount;

  child = JaMatrixAllocate(3, nvars+1); 

  count=0;

  SAME=TRUE;
  while (SAME==TRUE) {
    count++;

    do
      A = frange_ran(0.0,1.0);
    while (A==0);                   /* insure A is above 0 */
    
    for(i=1; i<=nvars; i++)
      {
	child[1][i] = p1[i] * A + p2[i] * (1.0-A);
	child[2][i] = p2[i] * A + p1[i] * (1.0-A);
      }

    if (count >= MAX_OPER_UNIQUE_TRY)
      {
	SAME=FALSE;
	break;
      }
    
    /* Are the two new individuals unique? */    
    tcount=0;
    for (i=1; i<=nvars; i++) {
      if ( (int) child[1][i] != (int) p1[i] )
	tcount++;
      
      if ( (int) child[2][i] != (int) p2[i] )
	tcount++;
    } /* end of i loop */
    
    if (tcount==(nvars*2)) SAME=FALSE;
  } /* end of while */

  for(i=1; i<=nvars; i++)
    {
      p1[i] = child[1][i];
      p2[i] = child[2][i];
    }
  
  JaMatrixFree(child, 3);  
} /* end of oper4() */


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper5()                                      */
/*                                 Multiple Point Simple Crossover              */
/*                                                                              */
/*           SYNOPSIS          :   void oper5(p1,p2,STEP,nvars,fin_mat,X,x2)    */
/*                                                                              */
/*           DESCRIPTION       :   This function returns two new vectors        */
/*                                  generated after simple arithmetical         */
/*                                  crossover, from the two parent vectors.     */
/*                                                                              */
/*           FUNCTIONS CALLED  :   irange_ran()                                 */
/*                                 matrix(),                                    */
/*                                 satis_con()                                  */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void oper5(VECTOR p1, VECTOR p2, int STEP, double **domains, int nvars)
     /* VECTOR p1,p2;   *The two parents for crossing over*/
     /* int    STEP;    *Parameter for the crossover*/
{
  MATRIX child;
  FLAG BFLAG1 = FALSE,/*Check to see if the newly created vectors satisfies the*/
       BFLAG2 = FALSE;/*set of constraints*/
  int i,n=1,cut;

  /* unique check variables */
  FLAG SAME;
  long count, tcount;

  child = matrix(1,2,1,nvars);

  count=0;
  SAME=TRUE;
  while (SAME==TRUE)
    {
      count++;

      /*Get a random spot on the vector for crossover*/
      cut = irange_ran(1,nvars);
      /*Copy the parent vectors on to the child vectors*/
      for(i=1; i<=cut; i++)
	{
	  child[1][i] = p1[i];
	  child[2][i] = p2[i];
	}
      do
	{
	  /*Cross the two vectors*/
	  for(i=cut + 1; i<=nvars; i++)
	    {
	      child[1][i] = p1[i] * (double)n/(double)STEP + p2[i] * (1.0-(double)n/(double)STEP);
	      child[2][i] = p2[i] * (double)n/(double)STEP + p1[i] * (1.0-(double)n/(double)STEP);
	    }
	  
	  /*Check to see if they satisfy the constraints*/
	  BFLAG1 = InBounds(child[1],domains,nvars);
	  BFLAG2 = InBounds(child[2],domains,nvars);
	  n++;
	  /*If the constraints not satisfied, then generate another*/
	  /*set of crossed over values*/
	}while((n<=STEP) && ((BFLAG1 == FALSE) || (BFLAG2 == FALSE)));
      
      /* Are the two new individuals unique? */
      if (count >= MAX_OPER_UNIQUE_TRY)
	{
	  SAME=FALSE;
	  break;
	}

      tcount=0;
      for (i=1; i<=nvars; i++) {
	if ( child[1][i] != p1[i] )
	  tcount++;
	
	if ( child[2][i] != p2[i] )
	  tcount++;
      } /* end of i loop */

      if (tcount==(nvars*2)) SAME=FALSE;

    } /* end of while (SAME==TRUE) */

  if (BFLAG1==TRUE && BFLAG2==TRUE)
    {
      for(i=1; i<=nvars; i++)
	{
	  p1[i] = child[1][i];
	  p2[i] = child[2][i];
	}
    }
  
  free_matrix(child,1,2,1);
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper6()                                      */
/*                                 Whole Non-Uniform Mutation                   */
/*                                                                              */
/*           SYNOPSIS          :   void oper6(parent,fin_mat,nvars,T,t,B)       */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, whole non-uniform mutation.   */
/*                                                                              */
/*           FUNCTIONS CALLED  :   find_range(),                                */
/*                                 flip(),                                      */
/*                                 get_F(),                                     */
/*                                 irange_ran(),                                */
/*                                 Gvector().                                   */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/********************************************************************************/


void oper6(VECTOR parent, double **domains, int nvars, int T, int t, int B)
     /* VECTOR parent;
	MATRIX fin_mat;
	unsigned long T;    Total number of generations
	unsigned long t;    Current generation number
	int B; */
{
  int  comp=0,
       i,
      *next;
  double llim,ulim;

  /* unique check variables */
  FLAG SAME;
  long count;
  double tmp=0;

  count=0;
  SAME=TRUE;
  while (SAME==TRUE)
    {
      count++;
      
      next = ivector(1, nvars);
      for(i=1; i<=nvars; i++)
	next[i] = 0;
      
      for (i=1; i<=nvars; i++)
	{
	  do
	    comp = irange_ran(1, nvars);
	  while (next[comp] == 1);
	  next[comp] = 1;
	  
	  find_range(&llim,&ulim,comp,domains,nvars,parent);
	  
	  /*From the current value of the component to be mutated, chooose at random*/
	  /*whether to mutate with a lesser value or a greater value*/
	  /*Then find a value lesser or greater than the original value from the*/
	  /*function get_f()*/
	  tmp = (flip() == TAIL) ? parent[comp]-get_F(T,t,parent[comp]-llim,B) :
	    parent[comp]+get_F(T,t,ulim-parent[comp],B);
	}

      if ( parent[comp] != tmp)
	SAME=FALSE;
      else if (count >= MAX_OPER_UNIQUE_TRY)
	SAME=FALSE;
    } /* end of while loop */

  free_ivector(next,1);
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper7()                                      */
/*                                 Heuristic Crossover                          */
/*                                                                              */
/*           SYNOPSIS          :   void oper7(p1,p2,nvars,fin_mat,X,x2)         */
/*                                                                              */
/*           DESCRIPTION       :   This function returns one new vector         */
/*                                                                              */
/*           FUNCTIONS CALLED  :   frange_ran()                                 */
/*                                 Gvector(),                                   */
/*                                 satis_con()                                  */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void oper7(VECTOR p1, VECTOR p2, double **domains, int nvars)
{
  VECTOR child;
  FLAG BFLAG = FALSE;/*Check to see if the newly created vector satisfies the*/
                      /*set of constraints*/
  int i,n=2,tries=MAX_OPER_UNIQUE_TRY;
  double A;

  /* unique check variables */
  FLAG SAME;
  long count;

  child = Gvector(1,nvars);

  count=0;
  SAME=TRUE;
  while (SAME==TRUE)
    {
      count++;
      
      do
	{
	  A = frange_ran(0.0,1.0);
	  for(i=1; i<=nvars; i++)
	    child[i] = (  A * (p2[i] - p1[i]) + p2[i] );
	  
	  /*Check to see if it satisfies the constraints*/
	  BFLAG = InBounds(child,domains,nvars);
	  n++;
	  /*If the constraints not satisfied, then try again */
	}
      while((n<=tries) && (BFLAG == FALSE));

      /* Is the new individual unique? */
      if (count >= MAX_OPER_UNIQUE_TRY)
	{
	  SAME=FALSE;
	  break;
	}
      
      for (i=1; i<=nvars; i++) {
	if ( child[i] != p1[i] ) {
	  SAME=FALSE;
	  break;
	}
      }

    } /* end of while SAME loop */

  if (BFLAG==TRUE)
    {
      for(i=1; i<=nvars; i++)
	p1[i] = child[i];
    }


  free_vector(child,1);
} /* end of oper7() */


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_range()                                 */
/*                                                                              */
/*           SYNOPSIS          :   void find_range(llim,ulim,domains,nvars,     */
/*                                                 parent                       */
/*                                                                              */
/*           DESCRIPTION       :   This function finds the upper and lower      */
/*                                  limits, within which the mutation should    */
/*                                  occur.                                      */
/*                                                                              */
/*           FUNCTIONS CALLED  :   none()                                       */
/*                                                                              */
/*           CALLING FUNCITONS :   oper1()                                      */
/*                                 oper2()                                      */
/*                                 oper3()                                      */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void find_range(double *llim, double *ulim, int comp, double **domains, int nvars, VECTOR parent)
     /* double *llim,*ulim; Upper and lower limits*/
     /* int comp;           Component of the vector to be mutated*/
     /* VECTOR parent;      The vector with the values of the variables*/
{
  double A, B;

  A = frange_ran(0.0,1.0);
  B = 1.0 - A;

  *llim = (A*domains[comp][1]) + B* parent[comp];

  A = frange_ran(0.0,1.0);
  B = 1.0 - A;

  *ulim = B * parent[comp] + (A*domains[comp][3]);
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   irange_ran()                                 */
/*                                                                              */
/*           SYNOPSIS          :   int irange_ran(llim,ulim)                    */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a random integer       */
/*                                  between the llim and ulim.                  */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   find_parent(),                               */
/*                                 oper1().                                     */
/*                                 oper2().                                     */
/*                                 oper3().                                     */
/*                                 oper5().                                     */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/
int irange_ran(int llim, int ulim)
{
  int num;

  do
    num =  llim + ((int) ((newrand()*(long)(ulim-llim+1))/(long) 65535));
  while ((num < llim) || (num > ulim));
  return(num);
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   get_F()                                      */
/*                                                                              */
/*           SYNOPSIS          :   double get_F(T,t,y,B)                         */
/*                                                                              */
/*           DESCRIPTION       :   This function returns the double value which  */
/*                                  is the evaluated value of the function,     */
/*                                  needed for the operators 3 and 6            */
/*                                                                              */
/*           FUNCTIONS CALLED  :   none()                                       */
/*                                                                              */
/*           CALLING FUNCITONS :   oper3()                                      */
/*           CALLING FUNCITONS :   oper6()                                      */
/*                                                                              */
/********************************************************************************/

double get_F(int T, int t, double y, int B)
     /*
       unsigned long t,T;
       int           B;
       double         y;
     */
{
  double factor;

  factor =  (double) pow(1.0 - (double)t/(double)T,(double)B);
  factor = factor * frange_ran(0.0,1.0);
  if (factor < 0.00001)
    factor = 0.00001;
  return(y * factor);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper8()                                      */
/*                                 Local-Minimum Crossover: bfgs                */
/*                                                                              */
/*           SYNOPSIS          :   void oper8(parent,fin_mat,rc)                */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, BFGS.                         */
/*                                                                              */
/*           FUNCTIONS CALLED  :   dfgsmin(),                                   */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void oper8(double (*VMfunction)(double *LX, long *LStatus),
	   VECTOR parent, MATRIX domains, 
	   double SolutionTolerance, short Optim, int nvars, 
	   short int MinMax, short BoundaryEnforcement, 
	   long InstanceNumber, FILE *output, long *LVMstatus,
	   short PrintLevel)
{

  double *parm, *hessin, *work;
  int i, j, evaliter, btest;
  double bfgsfit, evalgtol;
  double A, B;

  parm  = (double *) malloc((nvars)*sizeof(double)); 
  work  = (double *) malloc((nvars+1)*sizeof(double));
  hessin  = (double *) malloc((nvars*nvars+nvars)*sizeof(double));    

  A = frange_ran(0.0,1.0);
  B = 1.0 - A;

  for (i=0; i<nvars; i++) {
    parm[i] = parent[i+1];
  }

  if(Optim==0)
    {
      evalgtol=SolutionTolerance;
      dfgsmin(VMfunction, parm, nvars, evalgtol, &evaliter, &bfgsfit, hessin, MinMax,
	      BoundaryEnforcement, InstanceNumber, domains, LVMstatus, PrintLevel, output);
      if (*LVMstatus < 0) {
	free(hessin);
	free(work);
	free(parm);
    
	return ;
      }
    }
  else
    {
      bfgsfit = genoud_optim(parm, nvars);
    }
  

  if (BoundaryEnforcement<0) {
    for(i=1; i<=nvars; i++) {
      parent[i] = A * parm[i-1] + B * parent[i];
    }
  }
  else {
    for (j=0; j<20; j++) 
      {
	btest = 0;
	for (i=1; i<=nvars; i++) 
	  {
	    work[i] = A * parm[i-1] + B * parent[i]; 

	    btest = (domains[i][1] > work[i]) || (work[i] > domains[i][3]) ;
	    /* shrink point until all parameters are in bounds */
	    if (btest) 
	      {
		fprintf(output, "WARNING: killing boundary trigger in bfgs oper(8). fit:%10.8lf\n",bfgsfit);
		fprintf(output, "WARNING: oper(8) Parameter: %d \t Value: %e\n\n", i, work[i]);
	      }
	  }
	if (btest==0) break;
	A *= 0.5 ;
	B = 1.0 - A;
      }
    if (j<20) 
      { 
	/* leave parent unchanged if not in boundary after 20 halvings */
	for (i=1; i<=nvars; i++) 
	  {
	    parent[i] = work[i];
	  }
      }
  }

  free(hessin);
  free(work);
  free(parm);

  return ;
}

/********************************************************************************/
/* Integer Operators!                                                           */
/*                                                                              */
/********************************************************************************/

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegerOper1()                             */
/*                                 Uniform Mutation                             */
/*                                                                              */
/*           SYNOPSIS          :   void oper1(parent,fin_mat,rc)                */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, uniform mutation.             */
/*                                                                              */
/*           FUNCTIONS CALLED  :   find_range(),                                */
/*                                 frange_ran(),                                */
/*                                 irange_ran(),                                */
/*                                 Gvector().                                   */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void JaIntegerOper1(VECTOR parent, double **domains, int nvars)
     /* VECTOR parent;  The parent vector*/
     /* MATRIX fin_mat; The final matrix*/
     /* INDEX rc;       Row and column of the final matrix*/
{
  long comp;
  double llim,ulim;/*Lower and Upper limits of the value to be mutated*/
  FLAG SAME;
  double tmp;
  long count;

  count=0;

  SAME=TRUE;
  while (SAME==TRUE) 
    {
      count++;
      comp = irange_ran(1,nvars);
        
      /*Finding the lower and upper limits between which the values are to be mutated*/
      find_range(&llim,&ulim,comp,domains,nvars,parent);
      
      /*Find a random value between the lower and the upper limits, to substitute*/
      /*for the old value*/
      tmp =  (int) frange_ran(llim,ulim);

      if ( (int) parent[comp] != (int) tmp)
	SAME=FALSE;
      else if (count >= MAX_OPER_UNIQUE_TRY)
	SAME=FALSE;
    } /* end of while */

  parent[comp] = (int) tmp;
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegerOper2()                             */
/*                                 Boundary Mutatation                          */
/*                                 No Uniqueness checking here                  */
/*                                 Don't use this oper often!                   */
/*                                                                              */
/*           SYNOPSIS          :   void oper2(parent,fin_mat,rc)                */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, boundary mutation.            */
/*                                                                              */
/*           FUNCTIONS CALLED  :   find_range(),                                */
/*                                 flip(),                                      */
/*                                 irange_ran(),                                */
/*                                 Gvector().                                   */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void JaIntegerOper2(VECTOR parent, double **domains, int nvars)
     /* VECTOR parent;  The parent vector*/
     /* MATRIX fin_mat; The final matrix*/
     /* INDEX rc;       Row and column of the final matrix*/
{
  int comp;
  double llim,ulim;/*Lower and Upper limits of the value to be mutated*/

  FLAG SAME;
  double tmp;
  long count;

  count=0;
  SAME=TRUE;
  while (SAME==TRUE) 
    {
      count++;

      comp = irange_ran(1,nvars);
      
      /*Finding the lower and upper limits between which the values are to be mutated*/
      find_range(&llim,&ulim,comp,domains,nvars,parent);
      
      /*Replace either the lower limit or the upper limit at random,*/
      /*for the old value*/
      tmp =  (int) (flip() == TAIL) ? llim : ulim;

      if ( (int) tmp != (int) parent[comp])
	SAME=FALSE;
      else if (count >= MAX_OPER_UNIQUE_TRY)
	SAME=FALSE;
    } /* end of while */

  parent[comp] = (int) tmp;
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegeroper3()                             */
/*                                 Non-Uniform Mutation                         */
/*                                                                              */
/*           SYNOPSIS          :   VECTOR oper3(parent,fin_mat,r,c,T,t,B)       */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, non-uniform mutation.         */
/*                                                                              */
/*           FUNCTIONS CALLED  :   find_range(),                                */
/*                                 flip(),                                      */
/*                                 get_F(),                                     */
/*                                 irange_ran(),                                */
/*                                 Gvector().                                   */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void JaIntegerOper3(VECTOR parent, double **domains, int nvars, int T, int t, int B)
     /* VECTOR parent;
	MATRIX fin_mat;
	INDEX rc; */
     /* unsigned long T;   Total number of generations*/
     /* unsigned long t;   Current generation number*/
     /* int B; */
{
  int comp;
  double llim,ulim;

  FLAG SAME;
  double tmp;
  long count;

  count=0;

  SAME=TRUE;
  while (SAME==TRUE) {
    count++;
    comp = irange_ran(1,nvars);
 
    find_range(&llim,&ulim,comp,domains,nvars,parent);

    /*From the current value of the component to be mutated, chooose at random*/
    /*whether to mutate with a lesser value or a greater value*/
    /*Then find a value lesser or greater than the original value from the*/
    /*function get_f()*/
    tmp = 
      (int) (flip() == TAIL) ? parent[comp]-get_F(T,t,parent[comp]-llim,B) :
      parent[comp]+get_F(T,t,ulim-parent[comp],B);

    if ( (int) parent[comp] != (int) tmp)
      SAME=FALSE;
    else if (count >= MAX_OPER_UNIQUE_TRY)
      SAME=FALSE;
  } /* end of while */
  parent[comp] = (int) tmp;
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegeroper4()                             */
/*                                 Polytope Crossover                           */
/*                                                                              */
/*           SYNOPSIS          :   void oper3(parent,fin_mat,r,c,T,t,B)         */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, non-uniform mutation.         */
/*                                                                              */
/*           FUNCTIONS CALLED  :   find_range(),                                */
/*                                 flip(),                                      */
/*                                 get_F(),                                     */
/*                                 irange_ran(),                                */
/*                                 Gvector().                                   */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/


void JaIntegerOper4(VECTOR p1, VECTOR p2, int nvars)
     /* VECTOR p1,p2;  The two parents chosen for crossover*/
     /* int nvars;   Length of the vector*/
{
  double **child;
  long    i;
  double  A;

  FLAG SAME;
  long count, tcount;

  child = JaMatrixAllocate(3, nvars+1); 

  count=0;

  SAME=TRUE;
  while (SAME==TRUE) {
    count++;

    do
      A = frange_ran(0.0,1.0);
    while (A==0);                   /* insure A is above 0 */
    
    for(i=1; i<=nvars; i++)
      {
	child[1][i] = ( p1[i] * A + p2[i] * (1.0-A) );
	child[2][i] = ( p2[i] * A + p1[i] * (1.0-A) );
      }

    if (count >= MAX_OPER_UNIQUE_TRY)
      {
	SAME=FALSE;
	break;
      }

    /* Are the two new individuals unique? */    
    tcount=0;
    for (i=1; i<=nvars; i++) {
      if ( (int) child[1][i] != (int) p1[i] )
	tcount++;

      if ( (int) child[2][i] != (int) p2[i] )
	tcount++;
    } /* end of i loop */
	 
    if (tcount==(nvars*2)) SAME=FALSE;
  } /* end of SAME while */

    for(i=1; i<=nvars; i++)
      {
	p1[i] = (int) child[1][i];
	p2[i] = (int) child[2][i];
      }


  JaMatrixFree(child, 3);  
} /* end of JaIntegerOper4 */


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegerOper5()                             */
/*                                 Multiple Point Simple Crossover              */
/*                                                                              */
/*           SYNOPSIS          :   void oper5(p1,p2,STEP,rc,fin_mat,X,x2)       */
/*                                                                              */
/*           DESCRIPTION       :   This function returns two new vectors        */
/*                                  generated after simple arithmetical         */
/*                                  crossover, from the two parent vectors.     */
/*                                                                              */
/*           FUNCTIONS CALLED  :   irange_ran()                                 */
/*                                 matrix(),                                    */
/*                                 satis_con()                                  */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void JaIntegerOper5(VECTOR p1, VECTOR p2, int STEP, double **domains, int nvars)
     /* VECTOR p1,p2;   *The two parents for crossing over*/
     /* INDEX rc;       *Row and column of the final matrix*/
     /* MATRIX fin_mat; *The final matrix*/
     /* int    STEP;    *Parameter for the crossover*/
{
  MATRIX child;
  FLAG BFLAG1 = FALSE,/*Check to see if the newly created vectors satisfies the*/
       BFLAG2 = FALSE;/*set of constraints*/
  int i,n=1,cut;

  /* unique check variables */
  FLAG SAME;
  long count, tcount;


  child = matrix(1,2,1,nvars);

  count=0;
  SAME=TRUE;
  while (SAME==TRUE)
    {
      count++;
      /*Get a random spot on the vector for crossover*/
      cut = irange_ran(1,nvars);
      /*Copy the parent vectors on to the child vectors*/
      for(i=1; i<=cut; i++)
	{
	  child[1][i] = p1[i];
	  child[2][i] = p2[i];
	}
      do
	{
	  /*Cross the two vectors*/
	  for(i=cut + 1; i<=nvars; i++)
	    {
	      child[1][i] = p1[i] * (double)n/(double)STEP + p2[i] * (1.0-(double)n/(double)STEP);
	      child[2][i] = p2[i] * (double)n/(double)STEP + p1[i] * (1.0-(double)n/(double)STEP);
	    }
	  
	  /*Check to see if they satisfy the constraints*/
	  BFLAG1 = InBounds(child[1],domains,nvars);
	  BFLAG2 = InBounds(child[2],domains,nvars);
	  n++;
	  /*If the constraints not satisfied, then generate another*/
	  /*set of crossed over values*/
	}while((n<=STEP) && ((BFLAG1 == FALSE) || (BFLAG2 == FALSE)));

      /* Are the two new individuals unique? */
      if (count >= MAX_OPER_UNIQUE_TRY)
	{
	  SAME=FALSE;
	  break;
	}

      tcount=0;
      for (i=1; i<=nvars; i++) {
	if ( (int) child[1][i] != (int) p1[i] )
	  tcount++;
	
	if ( (int) child[2][i] != (int) p2[i] )
	  tcount++;
      } /* end of i loop */
      
      if (tcount==(nvars*2)) SAME=FALSE;
      
    } /* end of while (SAME==TRUE); */

  if (BFLAG1==TRUE && BFLAG2==TRUE)
    {
      for(i=1; i<=nvars; i++)
	{
	  p1[i] = (int) child[1][i];
	  p2[i] = (int) child[2][i];
	}
    }

  free_matrix(child,1,2,1);
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegerOper6()                             */
/*                                 Whole Non-Uniform Mutation                   */
/*                                                                              */
/*           SYNOPSIS          :   void oper6(parent,fin_mat,nvars,T,t,B)       */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, whole non-uniform mutation.   */
/*                                                                              */
/*           FUNCTIONS CALLED  :   find_range(),                                */
/*                                 flip(),                                      */
/*                                 get_F(),                                     */
/*                                 irange_ran(),                                */
/*                                 Gvector().                                   */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/********************************************************************************/


void JaIntegerOper6(VECTOR parent, double **domains, int nvars, int T, int t, int B)
     /* VECTOR parent;
	unsigned long T;    Total number of generations
	unsigned long t;    Current generation number
	int B; */
{
  int  comp=0,
       i,
      *next;
  double llim,ulim;

  /* unique check variables */
  FLAG SAME;
  long count;
  double tmp=0;

  count=0;
  SAME=TRUE;
  while (SAME==TRUE)
    {
      count++;

      next = ivector(1, nvars);
      for(i=1; i<=nvars; i++)
	next[i] = 0;
      
      for (i=1; i<=nvars; i++)
	{
	  do
	    comp = irange_ran(1, nvars);
	  while (next[comp] == 1);
	  next[comp] = 1;
	  
	  find_range(&llim,&ulim,comp,domains,nvars,parent);
	  
	  /*From the current value of the component to be mutated, chooose at random*/
	  /*whether to mutate with a lesser value or a greater value*/
	  /*Then find a value lesser or greater than the original value from the*/
	  /*function get_f()*/
	  tmp = (int) (flip() == TAIL) ? parent[comp]-get_F(T,t,parent[comp]-llim,B) :
	    parent[comp]+get_F(T,t,ulim-parent[comp],B);
	}

      if ( (int) parent[comp] != (int) tmp)
	SAME=FALSE;
      else if (count >= MAX_OPER_UNIQUE_TRY)
	SAME=FALSE;
    } /* end of SAME loop */

  parent[comp] = (int) tmp;

  free_ivector(next,1);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegerOper7()                             */
/*                                 Heuristic Crossover                          */
/*                                                                              */
/*           DESCRIPTION       :   This function returns one new vector         */
/*                                                                              */
/*           FUNCTIONS CALLED  :   frange_ran()                                 */
/*                                 Gvector(),                                   */
/*                                 satis_con()                                  */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void JaIntegerOper7(VECTOR p1, VECTOR p2, double **domains, int nvars)
{
  VECTOR child;
  FLAG BFLAG = FALSE;/*Check to see if the newly created vector satisfies the*/
                      /*set of constraints*/
  int i,n=2,tries=MAX_OPER_UNIQUE_TRY;
  double A;

  /* unique check variables */
  FLAG SAME;
  long count;

  child = Gvector(1,nvars);

  count=0;
  SAME=TRUE;
  while (SAME==TRUE)
    {
      count++;
      
      do
	{
	  A = frange_ran(0.0,1.0);
	  for(i=1; i<=nvars; i++)
	    child[i] = (int) (  A * (p2[i] - p1[i]) + p2[i] );
	  
	  /*Check to see if it satisfies the constraints*/
	  BFLAG = InBounds(child,domains,nvars);
	  n++;
	  /*If the constraints not satisfied, then try again */
	}
      while((n<=tries) && (BFLAG == FALSE));

      /* Is the new individual unique? */
      if (count >= MAX_OPER_UNIQUE_TRY)
	{
	  SAME=FALSE;
	  break;
	}
      
      for (i=1; i<=nvars; i++) {
	if ( (int) child[i] != (int) p1[i] ) {
	  SAME=FALSE;
	  break;
	}
      }

    } /* end of while SAME loop */

  if (BFLAG==TRUE) 
    {
      for(i=1; i<=nvars; i++)
	p1[i] = (int) child[i];
    }


  free_vector(child,1);
} /* end of JaIntegerOper7() */


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   InBounds(child, domains, nvars)              */
/*                                                                              */
/*           DESCRIPTION       :   This function returns TRUE or FALSE depending*/
/*                                  on whether the vector passed satisfies all  */
/*                                  the constraints or not.                     */
/*                                                                              */
/*           FUNCTIONS CALLED  :   none()                                       */
/*                                                                              */
/*           CALLING FUNCITONS :   oper5()                                      */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/
FLAG InBounds(VECTOR child, double **domains, int nvars)
     /* VECTOR child;   The vector to be checked for violation of constriants*/
{
  int i;
  

  for(i=1; i<=nvars; i++)
    {
      if( (child[i] < domains[i][1]) || (child[i] > domains[i][3]) )
	return(FALSE);
    }

  return(TRUE);
}


