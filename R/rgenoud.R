#
#
#  RGENOUD
#
#  Walter R. Mebane, Jr.
#  Cornell University
#  http://macht.arts.cornell.edu/wrm1
#  <wrm1@macht.arts.cornell.edu>
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.polisci.berkeley.edu
#  <sekhon@berkeley.edu>
#
#  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/rgenoud.R,v 2.0 2005/09/19 03:58:47 jsekhon Exp jsekhon $
#


genoud <- function(fn, nvars, max=FALSE, pop.size=1000, max.generations=100, wait.generations=10,
                   hard.generation.limit=TRUE, starting.values=NULL, MemoryMatrix=TRUE, 
                   Domains=NULL, default.domains=10, solution.tolerance=0.001,
                   gr=NULL, boundary.enforcement=0, lexical=FALSE, gradient.check=TRUE, BFGS=TRUE, 
                   data.type.int=FALSE, hessian=FALSE, unif.seed=812821, int.seed=53058,
                   print.level=2, share.type=0, instance.number=0,
                   output.path="stdout", output.append=FALSE, project.path=NULL,
                   P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
                   cluster=FALSE, balance=FALSE, debug=FALSE, ...)
{
  fn1 <- function(par) fn(par, ...)
  gr1 <- if (!is.null(gr)) {
    function(par) gr(par, ...)
  } else {
    gr1 <- NULL
  }

  #setpath to tempdir
  if(is.null(project.path))
    {
      project.path=paste(tempdir(),"/genoud.pro",sep="")
    }
  
  #do we have stating values?
  if (is.null(starting.values)) {
    nStartingValues <- 0;
  }
  else {
    nStartingValues <- 1;
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
      ret <- optim(foo.vals, fn=fn1, gr=gr1, method="BFGS",
                   control=list(fnscale=g.scale));
      return(c(ret$value,ret$par));
    } # end of genoud.optim.wrapper101

  #if lexical==TRUE we need to know how many items will be returned
  if (lexical < 0)
    {
      warning("lexical < 0.  Resetting to FALSE\n")
      lexical <- 0
    }
  if (lexical==1)
    {
      foo <- fn1(Domains[,1])
      lexical = length(as.vector(foo))
    }
  if (lexical > 0)
    {
      #All derivative stuff is turned off if we are going to do lexical sorting
      BFGS=FALSE
      gradient.check=FALSE
      hessian=FALSE
      P9 = 0
    }
  if (lexical==0)
    lexical <- 1

  # create the P vector
  P <- vector(length=9, mode="numeric");
  P[1] <- P1; P[2] <- P2; P[3] <- P3; P[4] <- P4;
  P[5] <- P5; P[6] <- P6; P[7] <- P7; P[8] <- P8;
  P[9] <- P9;

  clustertrigger=1
  if (is.logical(cluster))
    {
      if (cluster==FALSE)  {
        clustertrigger=0
      } else {
        stop("cluster option must be either FALSE, an object of the 'cluster' class (from the 'snow' package) or a list of machines so 'genoud' can create such an object")
      }
    }

  if(clustertrigger) {
    snow.exists = require("snow")
    if (!snow.exists) {
      stop("The 'cluster' feature cannot be used unless the package 'snow' can be loaded.")
    }
  }

  if(clustertrigger)
    {
      if(!MemoryMatrix)
        {
          MemoryMatrix = TRUE
          warning("MemoryMatrix has been set to 'TRUE' because the cluster option has been requested.")
        }

      GENclusterExport <- function (cl, list) 
        {
          gets <- function(n, v) {
            assign(n, v)
            NULL
          }
          for (name in list) {
            clusterCall(cl, gets, name, get(name))
          }
        }
      
      if (class(cluster)[1]=="SOCKcluster" | class(cluster)[1]=="PVMcluster" | class(cluster)[1]=="spawnedMPIcluster" | class(cluster)[1]=="MPIcluster") {
        clustertrigger=1
        cl <- cluster
        GENclusterExport(cl, "fn")
        GENclusterExport(cl, "fn1")
      } else {
        clustertrigger=2
        cluster <- as.vector(cluster)
        cat("You will now be prompted for passwords so your cluster can be setup.\n")
        cl <- makeSOCKcluster(cluster)
        GENclusterExport(cl, "fn")        
        GENclusterExport(cl, "fn1")
      }      
    }

  fnLexicalSort <- function(mat, parms)
    {
      # parms = (1) MinMax, (2) nvars, (3) lexical_end, (4) type (0, vars, 1, obj)

      decreasing=FALSE
      if(parms[1]==1)
        decreasing=TRUE

      if(parms[4]==0)
         {
           #on nvars not on obj function
           foo = "indx <- order("
           for(i in (2:(parms[2]+1)) )
             {
               foo = paste(foo,"mat[,",i,"], ",sep="")
             }
           foo = paste(foo,"mat[,",parms[2]+2,"], ",sep="")
           foo = paste(foo,"decreasing=FALSE)",sep="")
           eval(parse(text=foo))
           mat = mat[indx,]           
         } else {
           #lexical on obj function
           foo = "indx <- order(mat[,1], "
           for(i in (parms[2]+3):parms[3] )
             {
               foo = paste(foo,"mat[,",i,"], ",sep="")
             }
           foo = paste(foo,"decreasing=",decreasing,")",sep="")
           eval(parse(text=foo))
           mat = mat[indx,]
         }
      return(mat)
    }

  fnMemoryMatrixEvaluate <- function(Memory, population, parms)
    {
      EVALUATE = -93813381

      MinMax = parms[1]
      nvars = parms[2]
      lexical = parms[3]
      
      lexical.end = ncol(population)
      pop.size = nrow(population)

      vars.indx <- 2:(nvars+1)
      lexical.indx <- c(1,(nvars+3):lexical.end)

      FIRSTTIME = TRUE
      if (nrow(Memory) > 1)
        FIRSTTIME = FALSE

      nevaluate = 0

      mfunc <- function(pop,memory)
        {
          return(row(memory)[match(pop,memory)])
        }

      population.mat = matrix(population[,vars.indx], ncol=nvars)
      if (!FIRSTTIME)
        {
          Memory.mat = matrix(Memory[,vars.indx], ncol=nvars)
          match.matrix.indx = matrix(mfunc(population.mat,Memory.mat),ncol=nvars)

          for (i in 1:pop.size)
            {
              match.indx = match.matrix.indx[i,]
              if (is.finite(sum(match.indx)))
                {
                  if (sum(match.indx==rep(match.indx[1],nvars))==nvars)
                    {
                      #found match
                      population[i,] <- Memory[match.indx[1],]
                    }
                } else {
                  if (population[i,nvars+2] != 0) {
                    population[i,nvars+2] = EVALUATE
                    nevaluate = nevaluate+1
                  }
                }
            }
        } else {
          for (i in 1:pop.size)
            {
              population[i,nvars+2] = EVALUATE
              nevaluate = nevaluate+1
            }
        }

      #evaluation loop
      if (nevaluate > 0)
        {
          eval.indx <- population[,nvars+2]==EVALUATE
          ret=0
          if (clustertrigger==0)
            {
              in.mat = matrix(population.mat[eval.indx,], ncol=nvars)
              ret <- matrix(t(apply(in.mat, 1, fn1)), ncol=lexical)
            } else {
              if (balance==TRUE) {
                in.mat = t(matrix(population.mat[eval.indx,], ncol=nvars))
                cl.in <- as.list(as.data.frame(in.mat))
                cl.out <- clusterApplyLB(cl, cl.in, fn1, ...)
                try(ret <- matrix(t(data.frame(cl.out)), ncol=lexical), TRUE)
                if (!is.matrix(ret)) {
                  if (!debug) {
                    stop("Cluster returned an object which could not be turned into a data.frame.  Cluster may be in bad state.  Please consider restarting it. To see what the cluster returned please turn on the 'debug' option.")
                  } else {
                    cat("Cluster returned an object which could not be turned into a data.frame:\n")
                    print(cl.out)
                    stop("Cluster may be in bad state.  Please consider restarting it.")
                  }
                }
              } else {
                in.mat = matrix(population.mat[eval.indx,], ncol=nvars)
                cl.out = parRapply(cl, in.mat, fn1, ...)
                try(ret <- matrix(cl.out, byrow=TRUE, ncol=lexical), TRUE)
                if (!is.matrix(ret)) {
                  if (!debug) {
                    stop("Cluster returned an object which could not be turned into a data.frame.  Cluster may be in bad state.  Please consider restarting it. To see what the cluster returned please turn on the 'debug' option.")
                  } else {
                    cat("Cluster returned an object which could not be turned into a data.frame:\n")
                    print(cl.out)
                    stop("Cluster may be in bad state.  Please consider restarting it.")
                  }
                } 
              }
            } # else clustertrigger==0

          if (lexical < 2)
            {
              population[eval.indx,1] = ret[,1]
            } else {
              population[eval.indx,lexical.indx] = ret
            }
          population[eval.indx,nvars+2] = 0
          
          if(!FIRSTTIME)
            {
              Memory = rbind(Memory,population[eval.indx,])
            } else {
              Memory = matrix(population[eval.indx,], ncol=lexical.end)
              FIRSTTIME = FALSE
            }
        }#end of nevaluate
      
      if (lexical > 1)
        {
          population <- fnLexicalSort(population, c(MinMax,nvars,lexical.end,1))
        } else {
          if (MinMax==0)
            {
              population <- population[order(population[,1]),]
            } else {
              population <- population[order(population[,1], decreasing=TRUE),]
            }
        }

      return(as.vector(c(nrow(Memory), Memory, population)))
    } #end of fnMemoryMatrixEvaluate

  if (!is.null(gr))
    {
      UserGradient = TRUE
      gr1func <- gr1
    } else {
      UserGradient = FALSE      
      gr1func <- function() {}
    }

  if(data.type.int)
    {
      BFGS = FALSE
      gradient.check=FALSE
    }

  gout <- .Call("rgenoud", as.function(fn1), new.env(),
                as.integer(nvars), as.integer(pop.size), as.integer(max.generations),
                as.integer(wait.generations),
                as.integer(nStartingValues), as.vector(starting.values),
                as.vector(P), as.matrix(Domains),
                as.integer(max), as.integer(gradient.check), as.integer(boundary.enforcement),
                as.double(solution.tolerance), as.integer(BFGS), as.integer(data.type.int),
                as.integer(provide.seeds), as.integer(unif.seed), as.integer(int.seed),
                as.integer(print.level), as.integer(share.type), as.integer(instance.number),
                as.integer(MemoryMatrix), as.integer(debug),
                as.character(output.path), as.integer(output.type), as.character(project.path),
                as.integer(hard.generation.limit),
                as.function(genoud.optim.wrapper101), 
                as.integer(lexical), as.function(fnLexicalSort), as.function(fnMemoryMatrixEvaluate),
                as.integer(UserGradient), as.function(gr1func), 
                PACKAGE="rgenoud");

  indx1 <- 4;
  indx2 <- (indx1+lexical-1);
  value=gout[indx1:indx2];

  indx1 <- indx2+1
  indx2 <- indx1+nvars-1
  par = gout[indx1:indx2]

  indx1 <- indx2+1
  indx2 <- indx1+nvars-1
  if(!gradient.check & !BFGS )
    {
      gradients= gout[indx1:indx2]
      gradients = rep(NA, length(gradients))
    } else {
      gradients = gout[indx1:indx2]
    }

  indx1 <- indx2+1
  indx2 <- indx1+8
  operators=gout[indx1:indx2]

  if (hessian==TRUE)
    {
      hes <- optim(gout[5:(nvars+4)], fn=fn1, gr=gr1, method="BFGS", hessian=TRUE,
                   control=list(fnscale=g.scale));
      
      hes <- hes$hessian;
      
      ret <- list(value=value, par=par, gradients=gradients,
                  generations=gout[1], peakgeneration=gout[2], popsize=gout[3],
                  operators=operators,
                  hessian=hes);
    }
  else
    {
      ret <- list(value=value, par=par, gradients=gradients,
                  generations=gout[1], peakgeneration=gout[2], popsize=gout[3],
                  operators=operators);
    }

  if (clustertrigger==2)
    stopCluster(cl)
  
  return(ret)  
} #end of genoud()

