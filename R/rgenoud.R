#
#
#  RGENOUD (limited version)
#
#  Walter R. Mebane, Jr.
#  Cornell University
#  http://macht.arts.cornell.edu/wrm1
#  wrm1@macht.arts.cornell.edu
#
#  Jasjeet Singh Sekhon 
#  Harvard University and Lamarck, Inc.
#  http://jsekhon.fas.harvard.edu/
#  jsekhon@fas.harvard.edu
#
#  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/rgenoud.R,v 1.20 2002/11/06 02:12:13 jsekhon Exp $
#


genoud <- function(fn, nvars, max=FALSE, pop.size=1000, max.generations=100, wait.generations=10,
                  hard.generation.limit=TRUE, starting.values=NULL, MemoryMatrix=NULL, Debug=FALSE, 
                  Domains=NULL, default.domains=10,
                  gradient.check=TRUE, boundary.enforcement=2,
                  solution.tolerance=0.001, BFGS=TRUE, data.type.int=FALSE, hessian=FALSE,
                  unif.seed=812821, int.seed=53058,
                  print.level=2, share.type=0, instance.number=0,
                  output.path="stdout", output.append=FALSE, project.path="genoud.pro",
                  P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0)
{
  # let's load rgenoud.so if it is not already loaded.
#  if (!is.loaded("rgenoud"))
#    {
#      dyn.load("rgenoud.so");
#    }
#  else if (instance.number==0)
#    {
#      dyn.unload("rgenoud.so");
#      dyn.load("rgenoud.so");
#    }
  #we can only use R's optim() with this version of rgenoud().
  roptim <- TRUE;

  #do we have stating values?
  if (is.null(starting.values)) {
    nStartingValues <- 0;
  }
  else {
    nStartingValues <- 1;
  }

  #MemoryMatrix
  if (is.null(MemoryMatrix)) {
    MemoryMatrix <- TRUE;

    if (nvars > 20)
      {
        cat("\nWARNING: Since the number of parameters is greater than 20,\nWARNING: MemoryMatrix has been turned off by default.\nWARNING: You may turn it on using the MemoryMatrix flag.\nWARNING: This option increases speed at the cost of extra memory usage.\n\n")
        MemoryMatrix <- FALSE;
      }
  }

  #set output.type
  if (output.path=="stdout")
    {
      output.type <- 0;
    }
  else
    {
      if (output.append)
        {
          output.type <- 2;
        }
      else
        {
          output.type <- 1;
        }
    }

  # let's create the Domains if none have been passed.
  if (!(is.matrix(Domains)))
    {
      Domains <- matrix(nrow=nvars, ncol=2);
      for (i in 1:nvars)
        {
          Domains[i,1] <- -1*default.domains;
          Domains[i,2] <- default.domains;
        } # end of for loop
    } # end of Domains if
  # create the P vector
  P <- vector(length=9, mode="numeric");
  P[1] <- P1; P[2] <- P2; P[3] <- P3; P[4] <- P4;
  P[5] <- P5; P[6] <- P6; P[7] <- P7; P[8] <- P8;
  P[9] <- P9;

  # has the user provided any seeds?
  if (unif.seed==812821 && int.seed==53058)
    provide.seeds <- FALSE
  else
    provide.seeds <- TRUE;

  if (max==FALSE)
        {
          g.scale <- 1;
        }
  else
    {
      g.scale <- -1;
    }

  #optim st
  genoud.optim.wrapper101 <- function(foo.vals)
    {
      ret <- optim(foo.vals, fn=as.function(fn), method="BFGS",
                  control=list(fnscale=g.scale));
      return(c(ret$value,ret$par));
    } # end of genoud.optim.wrapper101


  gout <- .Call("rgenoud", as.function(fn), new.env(),
               as.integer(nvars), as.integer(pop.size), as.integer(max.generations),
               as.integer(wait.generations),
               as.integer(nStartingValues), as.vector(starting.values),
               as.vector(P), as.matrix(Domains),
               as.integer(max), as.integer(gradient.check), as.integer(boundary.enforcement),
               as.double(solution.tolerance), as.integer(BFGS), as.integer(data.type.int),
               as.integer(provide.seeds), as.integer(unif.seed), as.integer(int.seed),
               as.integer(print.level), as.integer(share.type), as.integer(instance.number),
               as.integer(MemoryMatrix), as.integer(Debug),
               as.character(output.path), as.integer(output.type), as.character(project.path),
               as.integer(hard.generation.limit),
               as.function(genoud.optim.wrapper101), as.integer(roptim));

  if (hessian==TRUE)
    {
      hes <- optim(gout[5:(nvars+4)], fn, method="BFGS", hessian=TRUE,
                  control=list(fnscale=g.scale));
      
      hes <- hes$hessian;

      ret <- list(value=gout[1], generations=gout[2], peakgeneration=gout[3], popsize=gout[4],
                 par=gout[5:(nvars+4)], gradients=gout[(nvars+5):(nvars+nvars+4)],
                 operators=gout[(nvars+nvars+5):(nvars+nvars+9+4)],
                 hessian=hes);
    }
  else
    {
      ret <- list(value=gout[1], generations=gout[2], peakgeneration=gout[3], popsize=gout[4],
                 par=gout[5:(nvars+4)], gradients=gout[(nvars+5):(nvars+nvars+4)],
                 operators=gout[(nvars+nvars+5):(nvars+nvars+9+4)]);
    }

  return(ret);
} #end of genoud()

