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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/eval.cpp,v 1.18 2002/10/17 03:45:19 jsekhon Exp $

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "genoud.h"

#ifdef SQL_DEFINE
#include "mysql.h"
#endif

#ifdef UNIX
  #include <sys/time.h>
  #include <sys/types.h>
  #include <unistd.h>   
#endif

#ifdef UNIX
	#include <unistd.h>   
	#ifndef os_getpid
		#define os_getpid getpid
	#endif
#else
        #include <winsock.h>
	#include <process.h>
	#ifndef os_getpid
		#define os_getpid _getpid
	#endif
#endif

double evaluate(double (*VMfunction)(double *LX, long *LStatus),
		VECTOR X, int nvars, long int *Status)
{
  double fit;

  fit = VMfunction(X+1, Status);

  return(fit);
}


void EvaluateMatrix(long (*VMfunctionMatrix)(double **Population, long population, long nvars),
		    MATRIX Population, long npopulation, long nvars, long *Status)
{

  *Status = VMfunctionMatrix(Population, npopulation, nvars);

}


#ifdef SQL_DEFINE
/*
  This function stacks up evaluation requests for parallel processing.
  But this does not allow for Memory Matrix.  That must be turned off!

  For the time being let's evaluate everything in the population matrix.
*/
long NetworkEvaluate(double (*VMfunction)(double *LX, long *LStatus), 
		     char *DBname, char *AgentName, 
		     MATRIX population, long pop_size, long nvars, long NetworkNumber,
		     long *Status, long PackageSize)
{
  double *X, tmp;
  long *VRowNumber, *VPnumber, RowNumber=0;
  int i, p, retval;
  struct timeval tv;                        // Timeout value
  fd_set rd_set;

  //printf("Network Evaluation\n");

  // Memory Allocation
  X = (double *) malloc((nvars)*sizeof(double));
  VRowNumber = (long *) malloc((pop_size)*sizeof(long));
  VPnumber = (long *) malloc((pop_size)*sizeof(long));

  // dispatch the individuals needed to be evaluated
  for (p=1; p<=pop_size; p++)
    {
      if (population[p][0] == EVALUATE)
	{
	  for (i=0; i<nvars; i++)
	    {
	      X[i] = population[p][i+1];
	    }
	  VPnumber[RowNumber]=p;

	  VRowNumber[RowNumber] = mysql_insert_pop(DBname, X, nvars);
	  if (VRowNumber[RowNumber] < 0)
	    {
	      *Status = -1;
	      retval = VRowNumber[RowNumber];
	      //printf("VRowNumber[RowNumber] < 0: %d\n", VRowNumber[RowNumber]);
	      goto CleanUP;
	    }
	  RowNumber++;
	} // end of if (population[p][0] == EVALUATE)
    } // end of for

  // use select for OS independent microsecond sleep
  // sample timeout of 0 seconds, needs to be reset each time!
  FD_ZERO(&rd_set);          // initialize
  FD_SET(893,&rd_set);	

  i=0;
  for (;;)
    {
      tmp = VMfunction(X, Status);
      if (*Status < 0)
	{
	  retval = -1;
	  goto CleanUP;
	}

      tv.tv_sec =  1;       //seconds
      tv.tv_usec = 0;     //microseconds (1/1,000,000)

      retval = select(1, &rd_set, NULL, NULL, &tv);

      retval = mysql_check_completion(DBname);
      //printf("remaining agents: %d, %d \n", retval, ++i);
      if (retval==0)
	break;
    }

  // retrieve results
  for (i=0; i<RowNumber; i++)
    {
      p=VPnumber[i];
      // printf("retrieving row: %d\n", VRowNumber[i]);
      population[p][0] = mysql_extract_value(DBname, VRowNumber[i]);
      // printf("got it\n");
      if (population[p][0]==EVALUATE)
	{
	  *Status = -1;
	  retval = -1;
	  goto CleanUP;
	}
    } // end of i loop;
  
  retval = 0;

 CleanUP:

  free(VPnumber);
  free(VRowNumber);
  free(X);
  return(retval);
} // end of EvaluateStacked


short setup_database(char *DBname, char *AgentName)
{

  const int  MAX_bbSQL_INSERT = 10240;
  char *command;
  short retval;
  int   p       = 0;
  
  /* create SQL command string  */
  if ((command = (char *) malloc(MAX_bbSQL_INSERT)) == NULL) {
    printf("Failed to allocated command in setup_database()\n");
    return -1;
  }

  // printf("p\n");
  sprintf(command,"DROP DATABASE %s", DBname);
  retval = mysql_noreturn_command(command);
  // printf("p\n");

  // printf("p\n");
  sprintf(command,"CREATE DATABASE %s", DBname);
  retval = mysql_noreturn_command(command);
  // printf("p\n");
  if (retval < 0)
    return retval;

  p = 0;
  p += sprintf(command+p,"CREATE TABLE %s.status ", DBname);
  p += sprintf(command+p,"( ");
  p += sprintf(command+p,"jobname       varchar(255) NOT NULL, ");
  p += sprintf(command+p,"tabname       varchar(20)  NOT NULL, ");
  p += sprintf(command+p,"row           int(11)      NOT NULL, ");
  p += sprintf(command+p,"computername  varchar(255) , ");
  p += sprintf(command+p,"ipaddress     varchar(16)  , ");
  p += sprintf(command+p,"port          int(11)      , ");
  p += sprintf(command+p,"pid           int(11)      , ");
  p += sprintf(command+p,"time1         datetime     , ");
  p += sprintf(command+p,"time2         datetime     , ");
  p += sprintf(command+p,"jobtime       double       , ");
  p += sprintf(command+p,"status        int(11)      , ");
  p += sprintf(command+p,"attempts      int(11)      , ");
  p += sprintf(command+p,"PRIMARY KEY (jobname,tabname,row), ");
  p += sprintf(command+p,"INDEX tidx (computername), ");
  p += sprintf(command+p,"INDEX tidx2 (status) ");
  p += sprintf(command+p,")");
  retval = mysql_noreturn_command(command);
  if (retval < 0)
    return retval;

  p = 0;
  p += sprintf(command+p,"CREATE TABLE %s.list ", DBname);
  p += sprintf(command+p,"( ");
  p += sprintf(command+p,"arg1 double, ");
  p += sprintf(command+p,"ret1 double, ");
  p += sprintf(command+p,"row integer AUTO_INCREMENT, ");
  p += sprintf(command+p,"PRIMARY KEY (row), ");
  p += sprintf(command+p,"INDEX tidx (row)  ");
  p += sprintf(command+p,")");
  retval = mysql_noreturn_command(command);
  if (retval < 0)
    return retval;

  p = 0;
  p += sprintf(command+p,"CREATE TABLE %s.job ",DBname);
  p += sprintf(command+p,"( ");
  p += sprintf(command+p,"jobname varchar(255) NOT NULL, ");
  p += sprintf(command+p,"program text, ");
  p += sprintf(command+p,"argc	  integer, ");
  p += sprintf(command+p,"retc	  integer, ");
  p += sprintf(command+p,"INDEX  tidx (jobname), ");
  p += sprintf(command+p,"PRIMARY KEY (jobname) ");
  p += sprintf(command+p,")");
  retval = mysql_noreturn_command(command);
  // printf("retval: %d\n", retval);
  if (retval < 0)
    return retval;

  p = 0;
  p += sprintf(command+p,"INSERT INTO pgenoudTest.job ",DBname);
  p += sprintf(command+p," VALUES ('pgenoud','arg1 %s',1,1)",AgentName);
  // printf("sending:%s*\n", command);
  retval = mysql_noreturn_command(command);
  if (retval < 0)
    return retval;

} // endof setup_database
#endif
