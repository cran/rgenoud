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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/frange_ran.cpp,v 1.20 2002/11/06 02:11:25 jsekhon Exp $

*/

#include "genoud.h"

/* replacements for random number functions in GENOCOP's frange_ran.c */

#include "unif.h"
#include "urans.h"

double newunif(void)
{
  extern long NewUnifSeed[MAXTHREADS];
  extern long ThreadNumber;
  /*  static integer aux[TLPAUXSIZE], iseed = NEWUNIFSEED ;*/
  static integer aux[TLPAUXSIZE];
  double wrkout;
  double wrk;

  /* get a random uniform double on (0,1] */
  ruxorv (&NewUnifSeed[ThreadNumber], 1, &wrk, aux);
  wrkout = (double) wrk;
  return(wrkout);
}

double frange_ran(double llim, double ulim)
     /*
       double ulim;
       double llim;
     */
{

  /*llim, ulim:  The upper and lower limits between which the random
      number is to be generated*/

  double num1, diff = ulim - llim;

  if (diff == 0)
    return(llim);
  else if(diff < 0.0001)
    return((flip() == TAIL) ? llim : ulim);
  do
    num1 = llim +  newunif()*(ulim-llim) ;
  while((num1<llim)||(num1>ulim));
  return(num1);
}

unsigned int randint (void)
{
  extern long RandIntSeed[MAXTHREADS];
  extern long ThreadNumber;
  /*  static integer aux[TLPAUXSIZE], iseed = RANDINTSEED ; */
  static integer aux[TLPAUXSIZE];
  integer wrk;
  int num;

  /* get a signed 32-bit number from the TLP generator */
  tlpseq (&RandIntSeed[ThreadNumber], 1, &wrk, aux);
  /* truncate to 16 bits */
  num = wrk%65535;
  return (num);
}

unsigned int newrand (void)
{
  return (randint());
}

