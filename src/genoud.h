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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/genoud.h,v 1.20 2002/11/06 02:11:35 jsekhon Exp $

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include <stdarg.h>
#include "LmkGenoud.h"

#ifdef OPTIM
extern "C"
{
  extern double genoud_optim(double *X, int parameters);
}
#endif

extern "C"
{

void MyRprintf(FILE *foo, char *out, ...);

inline void MyRprintf(FILE *foo, char *out, ...)
{
  extern void Rprintf(char*, ...);
  extern void Rvprintf(const char *format, va_list arg);
  
  va_list(ap); /*will point to each unnamed argument in turn*/
  va_start(ap, out);
  if (foo == stdout) {
    Rvprintf(out,ap); 
  }
  else if (foo == stderr) {
    Rvprintf(out,ap); 
  }
  else {
    vfprintf(foo, out, ap); 
  }
  va_end(ap); 
} // end of MyRprintf

#define fprintf MyRprintf  
}

#define M(ROW,COL,NCOLS) (((ROW)*(NCOLS))+(COL))
#define EVALUATE -645271937
#define MAXPATH 1000
#define MAXTHREADS 20
#define MAXINSTANCES 20
#define ERROR_CODE -99999
#define BIGNUMBER 1e300
#define MAX_OPER_UNIQUE_TRY 200

#define TRUE 1
#define FALSE 0

#define DOS_SYS  FALSE           /* set to true for dos, false for unix */
#define UNIX_SYS TRUE            /* set to true for unix, false for dos */

#define flip()  ((int) ((newrand()*(long)2)/(long) 65535))
#define MIN -32768
#define MAX 32768
#define BIGINT 100000000
#define HEAD 1
#define TAIL 0
#define TRIES 1000
#define MULT 25173
#define INCR 13849
#define MOD ((long int) 65536)
#define SHUFFLE 256   /* size of random number shuffle array */
#define EST_OFFSET 0  /* offset in hours from Eastern Standard Time zone)  */

#define NOTDONE_ADD 0.25
/* number of generations with no changes to be treated as convergence at
   generation limit */
#define NOTDONE_LIM 50

typedef double **MATRIX;
typedef double *VECTOR;
typedef int **IMATRIX;
typedef int *IVECTOR;
typedef int FLAG;
typedef int TOSS;
typedef struct {int r; int c;}INDEX;



struct GND_IOstructure
{
  /* --- Basic Input Parameters ---- */
  double	(*AgentFit)(double *LX, long *LStatus);
  long      (*AgentFitMatrix)(double **Population, long population, long nvars);
  long		nvars;
  long		PopSize;
  long		MaxGenerations;
  long		WaitGenerations; 
  double	P[9];		  /* Operators */
  double	**Domains;
  short		MinMax;
  short		GradientCheck;
  short		BoundaryEnforcement;  /* 0=anything goes, 1: regular; 2: no trespassing! */
  double	SolutionTolerance;
  long		Status;
  long		ThreadNumber;	/* indexed from zero */
  long          InstanceNumber; /* indexed from zero, the number of parallel runs */
  short		UseBFGS;        /* Use BFGS on the Best Individual 1= True, 0=False */
  short         AllowDynamicUpdating; /* T, F */
  short         DataType;       /* 1== integer, everything else equals float */
  short         DynamicPopulation; /* T, F */
  short         MemoryUsage;    /* 0=minimize, 1=normal */
  short         Debug;          /* T, F */
  short         HardGenerationLimit; // T, F
  short         Optim;               // T, F
  short         Network;             // T, F

  /* Starting Values (if we want to provide any) */
  double        **StartingValues; /* a matrix of starting values (each set consists of a row) */
  long          nStartingValues;  /* number of starting values */

  /* Random Number Stuff */
  short         ProvideSeeds; /* 0: no, 1: yes */
  long          UnifSeed;
  long          IntSeed;

  /* --- Ouput Diagnostics --- */
  double	*oResults; 
  double	*oGradients;
  long		oP[9];				/* operators used */
  long		oGenerations;
  long		oPeakGeneration;
  long		oPopSize;

  /* Output Files */
  char*		OutputPath;
  char*         ProjectPath;
  char*         PublicPopulationPath;

  /* I/O types */
  short		OutputType;
  short         PrintLevel;

  /* Parallel Processing Stuff */
  short         ShareLevel;
  short         ShareType;
  char*		DynamicPopulationPath;
  char*         DBname;         // The Name of the Agent to optimize
  char*         AgentName;      // The Name of the Agent to optimize


  T_VMRECORD	*pvm; /* read only.  Required for rentrant call to LamarckVM */

};

/* bfgs.c */
void dfgsmin(double (*VMfunction)(double *LX, long *LStatus),
	     double *p, int n, double gtol, int *iter, double *fret, double *hessian,
	     short int MinMax, short int BoundaryEnforcement, long InstanceNumber,
	     double **Domains, long *LVMstatus, short PrintLevel, FILE *output);

/* change_order.c file */
void get_var_order(IVECTOR tot, IVECTOR cart, IMATRIX var_order);
void find_x1_x2(int tot, IMATRIX var_order, IVECTOR x1, IVECTOR x2);
void find_ac1_ac2(int t1, int t2, int t3, IVECTOR x1, IVECTOR x2, MATRIX mat, MATRIX ac1, MATRIX ac2);
void find_lu1_lu2(IVECTOR tot, IVECTOR x1, IVECTOR x2, VECTOR dom, VECTOR dom1, VECTOR dom2);
void find_limits(int tot, MATRIX domains, VECTOR llim, VECTOR ulim);
void find_new_in_eq(VECTOR a1b, MATRIX a1a2, VECTOR ll, VECTOR ul, INDEX rc, MATRIX newin);
void find_org_in_eq(VECTOR a1_b, MATRIX a1_a2, VECTOR vec_d, MATRIX c1, MATRIX c2, int c1row,
		    INDEX a1a2, MATRIX org_ineq);
void initialize(MATRIX mat, INDEX rc);
void find_final_mat1(VECTOR l2, VECTOR u2, MATRIX finmat, int row, int col);
void find_final_mat2(MATRIX newin, int r, int c, int finr, MATRIX finmat);
void find_final_mat3(MATRIX orgin, int r, int c, int finr, MATRIX finmat);

/* eval.c */
double evaluate(double (*VMfunction)(double *LX, long *LStatus),
		VECTOR X, int nvars, long *Status);
void EvaluateMatrix(long (*VMfunctionMatrix)(double **Population, long population, long nvars),
		    MATRIX Population, long npopulation, long nvars, long *Status);

#ifdef SQL_DEFINE
long NetworkEvaluate(double (*VMfunction)(double *LX, long *LStatus), 
		     char *DBname, char *AgentName, 
		     MATRIX population, long pop_size, long nvars, long NetworkNumber,
		     long *Status, long PackageSize);
short setup_database(char *DBname, char *AgentName);
#else
inline long NetworkEvaluate(double (*VMfunction)(double *LX, long *LStatus), 
		     char *DBname, char *AgentName, 
		     MATRIX population, long pop_size, long nvars, long NetworkNumber,
		     long *Status, long PackageSize)
{
  long a;
  a=1;
  return(a);
  
}
#endif


/* evaluate.c */
double optimization(struct GND_IOstructure *Structure, VECTOR X, 
		    MATRIX domains, FILE *output);
void sort(short int MinMax, MATRIX  population, int pop_size,
	  long nvar);
void swap(double **x, double **y);
int find_parent(IVECTOR live, int pop_size);
void assign_probab(VECTOR probab, int pop_size, double Q);
double x_pow_y(double x, int y);
void find_cum_probab(VECTOR cum_probab, VECTOR probab, int pop_size);
void find_live(VECTOR cum_probab, IVECTOR live, int pop_size, int P);
int find_die(VECTOR cum_probab, IVECTOR die, int pop_size);
void SetRunTimeParameters(struct GND_IOstructure *Structure, 
			  short FirstTime,
			  long *PopSize, long *nvars, long *MaxGenerations, long *WaitGenerations,
			  short *MinMax, short *GradientCheck, short *BoundaryEnforcement, short *UseBFGS,
			  double *SolutionTolerance,
			  long *InstanceNumber, long *P, long *P0, long *P1, long *P2, long *P3, long *P4, long *P5, 
			  long *P6, long *P7, long *P8, short *PrintLevel, 
			  short *HardGenerationLimit, FILE *output);
double JaIntegerOptimization(struct GND_IOstructure *Structure, VECTOR X, 
			     MATRIX domains, FILE *output);
void JaIntegerSort(double **InMatrix, long n, long k);
int JaIntegerCMP(double **a, double **b) ;
void JaDoubleSort(double **InMatrix, long n, long k);
int JaDoubleCMP(double **a, double **b) ;
void JaDoubleMemoryMatrix_Gen0(struct GND_IOstructure *Structure, 
			       double **Memory, double **population, double *X,
			       long *UniqueCount, long OldUniqueCount,
			       int pop_size, int nvars, 
			       FILE *output, long *Status);
void JaDoubleMemoryMatrix(struct GND_IOstructure *Structure, 
			  double **Memory, double **population, double *X,
			  long *UniqueCount, long OldUniqueCount,
			  int pop_size, int nvars, FILE *output, long *Status);
void JaIntMemoryMatrix_Gen0(struct GND_IOstructure *Structure, 
			       double **Memory, double **population, double *X,
			       long *UniqueCount, long OldUniqueCount,
			       int pop_size, int nvars, FILE *output, long *Status);
void JaIntMemoryMatrix(struct GND_IOstructure *Structure, 
			  double **Memory, double **population, double *X,
			  long *UniqueCount, long OldUniqueCount,
			  int pop_size, int nvars, FILE *output, long *Status);
void JaDynamicPopulationCheck(struct GND_IOstructure *Structure,
			      double **population, int location, int pop_size, int nvars, 
			      FILE *output, long *Status);

/* frange_ran.c */
double newunif(void);
double frange_ran(double llim, double ulim);
unsigned int randint(void);
unsigned int newrand(void);

/* math.c */
/* Not needed here.  In here for completeness! */
void add(double *in1, double *in2, double *out, int row, int col);
void copy(double *in, double *target, int row, int col);
void multi(double *in1, double *in2, double *out,
	   int row1, int col1, int row2, int col2, int outrowcol[2],
	   FILE *output);
void scalarmulti(double scalar, double *in1, double *out, int row, int col) ;
void scalarmultioffdiag(double scalar, double *in1, double *out, int row, int col) ;
void subtract(double *in1, double *in2, double *out, int row, int col);
double trace(double *a, int n);
void transpose(double *orig_matrix, double *t_matrix, int orig_rows, int orig_columns);
void copy_matrix(MATRIX mat1, MATRIX mat2, int lr, int ur, int lc, int uc);
int Iround(double in);
void samplestats(double **obsdata, int numobsv, int novarsv, int weightflag, 
		 double *weightdata, FILE *output);
void populationstats(double **popdata, int numobsv, int novarsv, 
		     double *mean, double *var, double *skew, double *kur,
		     long *tobs);

/* multiply.c */
void mmprod(int m, int nm, int n, MATRIX mul_cm, MATRIX mul_am, MATRIX mul_bm);
void mvprod(int m, int nm, VECTOR cm, MATRIX am, VECTOR bm);

/* numerics.c */
double **JaMatrixAllocate(long n, long k);
void JaMatrixFree(double **M, long k);
short **JaShortMatrixAllocate(long nobs, long nvars);
void JaShortMatrixFree(double **M, long nobs);
MATRIX matrix(int nrl, int nrh, int ncl, int nch);
void nrerror(char error_text[]);
double *Gvector(int nl, int nh);
int **imatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);
void free_vector(double *v, int nl);
void free_ivector(int *v, int nl);
void free_matrix(double **m, int nrl, int nrh, int ncl);
void free_imatrix(int **m, int nrl, int nrh, int ncl);


/* operators.c file */
void oper1(VECTOR parent, double **domains, int nvars);
void oper2(VECTOR parent, double **domains, int nvars);
void oper3(VECTOR parent, double **domains, int nvars, int T, int t, int B);
void oper4(VECTOR p1, VECTOR p2, int x2_vari);
void oper5(VECTOR p1, VECTOR p2, int STEP, double **domains, int nvars);
void oper6(VECTOR parent, double **domains, int nvars, int T, int t, int B);
void oper7(VECTOR p1, VECTOR p2, double **domains, int nvars);
void oper8(double (*VMfunction)(double *LX, long *LStatus),
	   VECTOR parent, MATRIX domains, 
	   double SolutionTolerance, short Optim, int nvars, 
	   short int MinMax, short BoundaryEnforcement, 
	   long InstanceNumber, FILE *output, long *LVMstatus,
	   short PrintLevel);
void find_range(double *llim, double *ulim, int comp, double **domains, int nvars, VECTOR parent);
int irange_ran(int llim, int ulim);
double get_F(int T, int t, double y, int B);
void JaIntegerOper1(VECTOR parent, double **domains, int nvars);
void JaIntegerOper2(VECTOR parent, double **domains, int nvars);
void JaIntegerOper3(VECTOR parent, double **domains, int nvars, int T, int t, int B);
void JaIntegerOper4(VECTOR p1, VECTOR p2, int nvars);
void JaIntegerOper5(VECTOR p1, VECTOR p2, int STEP, double **domains, int nvars);
void JaIntegerOper6(VECTOR parent, double **domains, int nvars, int T, int t, int B);
void JaIntegerOper7(VECTOR p1, VECTOR p2, double **domains, int nvars);
FLAG InBounds(VECTOR child, double **domains, int nvars);

/*print_format.c */
long ReadPopulation(double **Data, long NewPopSize, long NewVars, FILE *output, FILE *fp);
void print_domains(MATRIX equal, int t_equ, short DataType, FILE *output);
void print_matrix(int lr, int ur, int lc, int uc, MATRIX mat, FILE *output);
void print_population(int popsize, int nvars, int generation, double **foo, FILE *out);
void print_vector(VECTOR arr, int l, int u, FILE *output);
void print_ivector(IVECTOR arr, int l, int u, FILE *output);


// mysql.cpp
#ifdef SQL_DEFINE

int InitializeMySQL(char *CnfFileName, char *dbname);
int gDestroy_MySQL(void);
long mysql_insert_pop( char *DBname, double *vec_parms, long nvars );
long mysql_check_completion(char *DBname);
double mysql_extract_value(char *DBname, long rownumb);
short mysql_noreturn_command(char *command);

#endif

// if we don't have any Numerical Recipes code

#ifdef NONR
inline void dfgsmin(double (*VMfunction)(double *LX, long *LStatus),
	     double *p, int n, double gtol, int *iter, double *fret, double *hessian,
	     short int MinMax, short int BoundaryEnforcement, long InstanceNumber,
	     double **Domains, long *LVMstatus, short PrintLevel, FILE *output)
{
  
}
#endif
