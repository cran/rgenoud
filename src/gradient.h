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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/gradient.h,v 1.19 2002/10/19 08:29:54 jsekhon Exp $

*/

#define EPS 0.000000000000001   /* machine precision */

void gradient(double (*VMfunction)(double *LX, long *LStatus),
	      double *p, double *g, long nparms, short MinMax, short BoundaryEnforcement,
	      long InstanceNumber, double **Domains, long *LVMstatus);

double func4g(double (*VMfunction)(double *LX, long *LStatus),
	      double *X, long nvars, short MinMax, short BoundaryEnforcement, 
	      double **Domains, long *Status);

void numgrad(double (*VMfunction)(double *LX, long *LStatus),
	     double *epsacc, double *optint,
	     int nparms, double *invals, double *grads, double *wrk,
	      double (*func)(double (*VMfunction)(double *LX, long *LStatus),
			     double *X, int nvars, short int MinMax, long *Status), 
	     short MinMax, long *LVMstatus);

void numgradc(double (*VMfunction)(double *LX, long *LStatus),
	      double *epsacc, double *optint,	      
	      int nparms, double *invals, double *grads, double *wrk,
	      double (*func)(double (*VMfunction)(double *LX, long *LStatus),
			     double *X, long nvars, short MinMax, 
			     short BoundaryEnforcement, double **Domains, long *Status), 
	      short MinMax, short BoundaryEnforcement, double **Domains, long *LVMstatus);

extern double *numopgc(double *epsacc, double *optint,
		int nparms, int nobs, double *invals, double *opg, double *wrk,
		int (*func)(double *X, double *outvec));

double **eaccuracy(double (*VMfunction)(double *LX, long *LStatus),
		   int nparms, int ndiffs, double h, double *invals,
		   double *wrk, 
		   double (*func)(double (*VMfunction)(double *LX, long *LStatus),
				  double *X, long nvars, short MinMax, 
				  short BoundaryEnforcement, double **Domains, long *Status), 
		   short MinMax, short BoundaryEnforcement, double **Domains, long *LVMstatus);

struct estints {
  int nparms;
  int *errors;  /* 0 == OK, >=1 == error */
  double
    *hf,   /* interval */
    *phi,  /* f' (first derivative) */
    *phic,  /* f' (first derivative, central-difference) */
    *phi2, /* f'' (second derivative) */
    *ef,   /* error bound */
    *hessian;  /* hessian matrix (lower triangle) */
};

struct estints *algfd(double (*VMfunction)(double *LX, long *LStatus),
		      int nparms, double *eps, double *invals, double *wrk,
		      double (*func)(double (*VMfunction)(double *LX, long *LStatus),
				     double *X, long nvars, short MinMax, 
				     short BoundaryEnforcement, double **Domains, long *Status), 
		      short MinMax, short BoundaryEnforcement, double **Domains, long *LVMstatus);

void fdestimates(double (*VMfunction)(double *LX, long *LStatus),
		 int parm, double fvalue, double *invals, double *wrk,
		 double eps, double h,
		 double *fplus, double *fminus,
		 double *phif, double *phib, double *phic, double *phi2,
		 double *cf, double *cb, double *c2,
		 double (*func)(double (*VMfunction)(double *LX, long *LStatus),
				double *X, long nvars, short MinMax, 
				short BoundaryEnforcement, double **Domains, long *Status), 
		 int nparms, short MinMax, short BoundaryEnforcement, double **Domains, 
		 long *LVMstatus);

struct estints *numhessian(struct estints *instruc, double *invals, double *wrk,
			   double (*func)(double *));

struct estints *numhessianc(double (*VMfunction)(double *LX, long *LStatus),
			    struct estints *instruc, double *invals, double *wrk,
			    double (*func)(double (*VMfunction)(double *LX, long *LStatus),
					   double *X, long nvars, short MinMax, 
					   short BoundaryEnforcement, double **Domains, long *Status), 
			    short MinMax, short BoundaryEnforcement, double **Domains, long *LVMstatus);

void estoptint(double (*VMfunction)(double *LX, long *LStatus),
	       double *epsacc, double *optint,
	       int nparms, int ndiffs, int pflag, double *invals,
	       double (*func)(double (*VMfunction)(double *LX, long *LStatus),
			      double *X, long nvars, short MinMax, short BoundaryEnforcement,
			      double **Domains, long *Status), 
	       short MinMax, short BoundaryEnforcement, double **Domains, long *LVMstatus);

void dohessians(double (*VMfunction)(double *LX, long *LStatus),
		double *epsacc, 
		int nparms, int nobs, int ndiffs, double *invals,
		double (*func)(double (*VMfunction)(double *LX, long *LStatus),
			       double *X, long nvars, short MinMax, 
			       short BoundaryEnforcment, double **Domains, long *Status), 
		double (*funco)(double *, double *),
		short MinMax, short BoundaryEnforcement, double **Domains, long *LVMstatus);
