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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/LmkGenoud.h,v 1.22 2003/07/12 04:57:03 jsekhon Exp $

*/

#ifdef MS_WINDOWS 
     #include <limits.h>
#else
     #include <inttypes.h>
#endif

#ifdef MS_WINDOWS
     #include <basetsd.h>			/*  __int64  LARGE_INTEGER */
#endif 

/* ----- Lamarck Variable types (also T_FLOAT above) ----- */
#ifdef SINGLE_PRECISION
	typedef	float T_FLOAT;

	#define  L_OVERFLOW		FLT_MAX	    /* max value = number that is returned by VM when overflow errors occur. */
	#define  L_UNDERFLOW	FLT_MIN 	/* min value = number that is returned by VM when underflow errors occur */
	#define  L_EPSILON		FLT_EPSILON /* Smallest number such that 1.0 + L_EPSILON != 1.0 */
	#define  L_DIG			FLT_DIG		/* number of decimal digits of precision */
#else
	typedef	double T_FLOAT;

	#define  L_OVERFLOW		DBL_MAX	    /* max value = number that is returned by VM when overflow errors occur. */
	#define  L_UNDERFLOW	DBL_MIN 	/* min value = number that is returned by VM when underflow errors occur */
	#define  L_EPSILON		DBL_EPSILON /* Smallest number such that 1.0 + L_EPSILON != 1.0 */
	#define  L_DIG			DBL_DIG		/* number of decimal digits of precision */
#endif

typedef char			T_BOOLEAN;		/* Lamarck Boolean is FALSE = 0, TRUE <> 0, usually TRUE = 1 */
typedef char			T_CHAR;			/* storage class for (ASCII) character set */
typedef char			T_BYTE;
typedef short			T_ERRNO;		/* Lamarck Error Number */
typedef short			T_DEPTH;		/* depth or complexity of a logic structure */
typedef	short			T_SHORT;
typedef long			T_CODON;		/* storage type for "codons", individual instruction numbers */
typedef long			T_GENOME_INDEX;	/* storage type for indices of Lamarck Agent Genomes (refer to L_MAX_GENOME_LENGTH) */
typedef long			T_DATA_INDEX;	/* storage type for indices to data stacks and memory (refer to L_MAX_DATA_STACK) */
typedef long			T_LONG;
typedef	unsigned long	T_COUNTER;		/* Lamarck counter: steps, errors, etc. */
typedef unsigned long	T_MEM_ADDRESS;
typedef float			T_SINGLE;
typedef double          T_DOUBLE;

#ifdef MS_WINDOWS
  typedef __int64		T_TIME;
#else
  typedef int64_t               T_TIME;
#endif

typedef struct {
	T_FLOAT			value;
	T_CODON			type;
	T_GENOME_INDEX	origin;
	T_GENOME_INDEX  dependency;
} T_ITEM;

/* -- Virtual Machine Record Structure -- */
/* === Control Flow structures ==== */
typedef struct {
    T_GENOME_INDEX	ElseBit;
	T_GENOME_INDEX	EndBit;
    T_BOOLEAN		TestResult;
} T_IFENDIF;


typedef struct {
    T_FLOAT			value;
	T_FLOAT			BegVal;
	T_FLOAT			EndVal;
	T_FLOAT			stepSize;
    T_COUNTER		LoopCount;
	T_GENOME_INDEX	indexReturn;
	T_GENOME_INDEX  indexNext;
	T_DEPTH			entryDepth;
	T_DEPTH			entryDepthDL;
	T_DEPTH			entryDepthIF;
} T_FORNEXT;


typedef struct {
	T_COUNTER		LoopCount;
	T_GENOME_INDEX	indexReturn;
	T_DEPTH			entryDepth;
	T_DEPTH			entryDepthFN;
	T_DEPTH			entryDepthIF;
} T_DOLOOP;


typedef struct {
	T_CODON			codon;
	T_DATA_INDEX	ArgC;
	T_ITEM			*ArgV;
} T_CODON_PARMS;


/* --------------------------------------------------------------------------------------------------
	Any changes to T_VMRECORD need to be echoed in the vmFields information structure
	and in the Registration of vmRecord.___ functions
----------------------------------------------------------------------------------------------------- */
typedef struct {
	T_BOOLEAN		flagDataMode;		/* treat codons as data.  Put instruction numbers (codon#s) on stack */
	T_BOOLEAN		flagDebug;			/* stop execution after each codon */
	T_BOOLEAN		flagBreakOnError;	/* stop execution after each error */
	T_BOOLEAN		flagShareMemory;	/* share memory with sub-agent */
	T_BOOLEAN		flagShareStack;		/* share stack with sub-agent */
	T_BOOLEAN		flagAgentCodon;		/* current codon "behaves like" an agent (e.g. C_CALLME) */
	T_BOOLEAN		flagRedirection;	/* allow redirection of codons */
	T_BOOLEAN		flagSwitchToSubAgent;	/* transfer control to sub-agent */
	T_BOOLEAN		flagGrammarAction;	/* grammar repair in process */
	T_BOOLEAN		flagErrorAction;	/* error repair in process */
	T_BOOLEAN		flagSuccessAction;
	T_BOOLEAN		flagAbort;			/* abort evaluation */
	T_BOOLEAN		flagSubExpression;	/* SubExpression in process (vs. SubAgent) */
	T_BOOLEAN		flagSilentMode;
	T_BOOLEAN		flagSaveStack;	
	T_BOOLEAN		flagCircularLogic;	/* Check for circular logic before invoking agents */ 
	T_ERRNO			errorNumber;	/* Error Number */
	T_DEPTH         depth;			/* VM depth current */
	T_DEPTH			depthEntry;		/* VM record number of "Process" entry point */
	T_DEPTH			depthParent;	/* VM record number of Parent Agent (records can be interlaced) */
	T_DEPTH			depthChild;		/* VM record number of Child Agent */
    T_DEPTH			depthFN;		/* Nesting of For..Next */
    T_DEPTH			depthDL;		/* Nesting of Do..Loop */
    T_DEPTH			depthIF;		/* Nesting of If..EndepthIF */
	T_DEPTH			maxDepthMe;		/* Maximum depth for this agent */
	T_DEPTH			maxDepthSubs;	/* Maximum depth of sub agents */
	T_GENOME_INDEX  genome_length; 
	T_GENOME_INDEX	indexCurrent;		/* (last) index into genome */
	T_GENOME_INDEX	indexNext;
	T_GENOME_INDEX	indexResume;      /* If an error occurs, jump back to indexResume if L_RESUME issued */
	T_GENOME_INDEX	indexDependency;
	T_CODON			codon;			/* (last) genetic code bit */
	T_CODON			codonLast;		
	T_CODON			codonNext;
	T_CODON			codonErrorAction;	/* Jump execution to this marker codon, Tag[1] .. Tag[9] if an error is encountered */
	T_CODON			codonGrammarAction;		/* Transfer execution to this codon if insufficient arguments for a function (Grammar Interrupt) */
	T_CODON			codonSuccessAction;
	T_DATA_INDEX	dp;				/* Data stack pointer */
	T_DATA_INDEX	mp;				/* number of items in memory */
	T_DATA_INDEX 	nArgsAgent;		/* number of arguments passed to agent (args copied into A[]. */
	T_DATA_INDEX	nArgsCodon;		/* number of arguments for current codon */
	T_DATA_INDEX	nRets;			/* return values from current codon (or agent at C_END) */
	T_DATA_INDEX 	lastArgC;
	T_DATA_INDEX	maxStack;
	T_DATA_INDEX	maxMem;
	T_DATA_INDEX	maxArgs;
	T_COUNTER		steps;				/* Total number of instructions */
	T_COUNTER		errorCount;			/* Counts the number of times the error handling routines are invoked */
	T_COUNTER		errorCountAtStart;	/* ErrorCount at Initialization */
    T_COUNTER		maxSteps;
	T_TIME          maxTimeMe;
	T_TIME          maxTimeAgent;
	T_TIME			timeStart;		/* time Agent started */
	T_TIME			timeEnd;		/* time when system time last called.  Initially, timeEnd = timeStart */
	T_TIME			timeElapsed;	/* cumulative runtime.  Adjusts for time spent sitting in debug window */
	T_FLOAT			ActionThreshold;
	T_ITEM	        Evaluation;		    /* Evaluation */
    /* ------ Simple pointers (vectors non-allocated on per vmRecord basis)  ----- */
	T_CHAR			*AgentName;		/* pointer to Agent name (passed to evaluate()).  Need for "SaveAgent" */
	T_CODON			*GenomeBuffer;
	T_CODON			*ptrActVarType;
	T_FLOAT		    *ptrActVarValue;
	T_FLOAT			*ptrBrkOnChgVar;
	T_ITEM			*D;				/* ptr to start of Agent data stack, located in public D */
    T_ITEM			*A;				/* ptr to start of Agent argument space (Args passed to Agent) */
	T_ITEM			*M;				/* ptr to start of Agent Memory, located in public M */
	T_ITEM			*vmArgs;
	T_ITEM			*vmRets;
	/* ------- Pointers to allocated memory ----------- */
	T_ERRNO			*Intron;		/* Use tor to determine trimmed, deadwood-free, genome */
    T_CODON			*Genome;		/* Agent Program.  String of Ints with termination element, L_END */
	T_ITEM			*Var;			/* Private variables */
	T_FORNEXT		*FN;			/* Data structures for FOR..NEXT */
    T_DOLOOP		*DL;			/* Data structures for WHILE (expr) DO..LOOP */
    T_IFENDIF		*IF;			/* data structure for IF..ELSE ..ENdepthIF */
	T_CODON_PARMS	*AV;			/* data structure for Parameters of Active Variable */
} T_VMRECORD;


/* -------------- Lamarck Virtual Machines functions --------------------- */
/* 
#ifdef MS_WINDOWS
     #define CALL_STYLE __stdcall
#else 
     #define CALL_STYLE __attribute__((stdcall))
#endif

T_DATA_INDEX	CALL_STYLE	LVM_SubEval(T_VMRECORD *pvmParent, T_CHAR *SubAgent, T_DATA_INDEX nArgs);
T_DATA_INDEX	CALL_STYLE	LVM_EvalSubScript(T_VMRECORD *pvmParent, T_CHAR *SubAgentName, T_CHAR *AgentScript, T_DATA_INDEX nArgs);
T_GENOME_INDEX	CALL_STYLE	LVM_TextToGenome(T_CHAR *AgentText, T_CODON *Genome, T_BOOLEAN AutoRegister);
T_GENOME_INDEX  CALL_STYLE	LVM_Max_Genome_Length(void);
*/
