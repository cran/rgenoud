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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/unif.h,v 1.20 2002/11/06 02:12:21 jsekhon Exp $

*/


/* Tausworthe-Lewis-Payne generator */

#ifndef DEFTLPGEN
typedef long int integer;

#define TLPAUXSIZE 1300
#define TLPDBSIZ 2000
#define ZEROI ((integer)0)
#define TLPMOD 2147483647

void ruxorv (integer *, int, double *, integer *);
void tlpseq (integer *, int, integer *, integer *);
void tauint (integer *);
void tlpcor (integer, int, integer *, integer *);

#define DEFTLPGEN
#endif
