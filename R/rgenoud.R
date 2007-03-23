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
#


genoud <- function(fn, nvars, max=FALSE, pop.size=1000, max.generations=100, wait.generations=10,
                   hard.generation.limit=TRUE, starting.values=NULL, MemoryMatrix=TRUE, 
                   Domains=NULL, default.domains=10, solution.tolerance=0.001,
                   gr=NULL, boundary.enforcement=0, lexical=FALSE, gradient.check=TRUE, BFGS=TRUE, 
                   data.type.int=FALSE, hessian=FALSE, unif.seed=812821, int.seed=53058,
                   print.level=2, share.type=0, instance.number=0,
                   output.path="stdout", output.append=FALSE, project.path=NULL,
                   P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
                   P9mix=NULL, BFGSfn=NULL, BFGShelp = NULL,
                   cluster=FALSE, balance=FALSE, debug=FALSE, ...)
{
  if(!is.null(BFGShelp) && !is.function(BFGShelp)) stop("'BFGShelp' must be NULL or a function")

  if(!is.null(P9mix) && !is.real(P9mix))  {
    stop("'P9mix' must be NULL or a number between 0 and 1")
  } else {
    if(is.null(P9mix)) {
      P9mix <- -1
    } else {
      if(! ( (1 >= P9mix) && (P9mix > 0) ))
        stop("'P9mix' must be NULL or a number between 0 and 1 (it may be equal to 1)")
    }
  }

  if (max==FALSE)
    {
      g.scale <- 1;
      FiniteBadFitValue <- .Machine$double.xmax
    } else  {
      g.scale <- -1;
      FiniteBadFitValue <- -.Machine$double.xmax
    }

  if(!lexical & !is.null(BFGSfn))
    {
      stop("'BFGSfn' can only be provided with lexical optimization")
    }
  if (!is.null(BFGSfn) & BFGS==FALSE)
    {
      if (!is.function(BFGSfn))
        stop("IF 'BFGSfn' is not a function, it must be NULL")
      warning("setting BFGS==TRUE because 'BFGSfn' is not null")
      BFGS <- TRUE
    }

  fn1 <- function(par) {
    fit <- fn(par, ...)

    if(is.null(fit))
      fit <- FiniteBadFitValue

    if(length(fit)==1)
      if(!is.finite(fit))
        fit <- FiniteBadFitValue

    return(fit)
  }#end of fn1
  
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

  #if lexical==TRUE we need to know how many items will be returned
  if (lexical < 0)
    {
      warning("lexical < 0.  Resetting to FALSE\n")
      lexical <- 0
    }
  if (lexical>=1) # BG: check regardless in case someone inputs the wrong value for lexical
    {
      if(nStartingValues)
        {
          foo <- fn1(starting.values)          
        } else {
          rfoo <- runif(nrow(Domains), Domains[,1], Domains[,2])
          if(data.type.int)
            rfoo <- as.integer(round(rfoo))
          foo <- fn1(rfoo)          
        }
	foo.length <- length(as.vector(foo))
	if(lexical > 1 && foo.length != lexical) {
	  warning(paste("Function returns a vector of length", foo.length, 
                        "\nbut you specified lexical =", lexical))
	}
        lexical <- foo.length
    }
  if (lexical > 0)
    {
      if(is.null(BFGSfn))
         {
           #All derivative stuff is turned off if we are going to do lexical if BFGSfn is not provided
           BFGS=FALSE
           gradient.check=FALSE
           if(hessian) {
             warning("'hessian' being set to false because of lexical optimization.  See 'BFGSfn' for workaround")
             hessian=FALSE             
           }

           P9 = 0
         } else {
           fn1.bfgs <- function(par, helper = NA) {
             fit <- if(is.null(BFGShelp)) BFGSfn(par, ...) else BFGSfn(par, helper, ...) 
             
             if(is.null(fit))
               fit <- FiniteBadFitValue
             
             if(length(fit)==1)
               if(!is.finite(fit))
                 fit <- FiniteBadFitValue
             
             return(fit)
           }#end of fn1.bfgs

           if(is.null(gr)) {
             gr <- function(par, helper, ...)
               {
                  gr.fn1.bfgs <- function(par, helper, FBFV) {
                    fit <- if(is.null(BFGShelp)) BFGSfn(par, ...) else BFGSfn(par, helper, ...) 
                    
                    if(is.null(fit))
                      fit <- FBFV
                    
                    if(length(fit)==1)
                      if(!is.finite(fit))
                        fit <- FBFV
                    
                    return(fit)
                  }  # end of gr.fn1.bfgs               
                 if(is.na(helper) && !is.null(BFGShelp)) {
                   helper <- do.call(BFGShelp, args = list(initial = par), envir = environment(fn))
                 }
                 genoud.wrapper101.env <- new.env()
                 assign("x", par, env = genoud.wrapper101.env)
                 assign("helper", helper, env = genoud.wrapper101.env)
                 assign("FiniteBadFitValue", FiniteBadFitValue, env = genoud.wrapper101.env)
                 foo <- as.real(attr(numericDeriv(quote(gr.fn1.bfgs(x, helper, FiniteBadFitValue)), theta=c("x"), genoud.wrapper101.env), "gradient"))
                 return(foo)
               } #end of gr
             gr1 <- function(par, helper = NA) gr(par, helper, ...)
	          } # end of if(!is.null(gr))
           gr1func <- gr1
         }# end of else
    }#if lexical > 0
  if (lexical==0)
    lexical <- 1

  #optim st
  if(is.null(BFGSfn))
     {
       genoud.optim.wrapper101 <- function(foo.vals)
         {
           ret <- optim(foo.vals, fn=fn1, gr=gr1, method="BFGS",
                        control=list(fnscale=g.scale));
           return(c(ret$value,ret$par));
         } # end of genoud.optim.wrapper101
     } else {
       genoud.optim.wrapper101 <- function(foo.vals)
         {
	       if(print.level > 2) {
      		 fit <- fn1(foo.vals)
       		 cat("\nPre-BFGS Complete Lexical Fit:\n")
		       print(fit)
         }
         if(is.null(BFGShelp)) {
             ret <- optim(foo.vals, fn=fn1.bfgs, gr=gr1, method="BFGS",
                          control=list(fnscale=g.scale));
         }
         else {
             ret <- optim(foo.vals, fn=fn1.bfgs, gr=gr1, method="BFGS",
                          control=list(fnscale=g.scale),
                          helper = do.call(BFGShelp, args = list(initial = foo.vals), envir = environment(fn)) );
         }

         if(print.level > 2)
             {
               cat("BFGS Fit:",ret$value,"\n")

               fit <- fn1(ret$par)
               cat("Post-BFGS Complete Lexical Fit:\n")
               print(fit)
             }           
           
         foo <- c(ret$value,ret$par)
         return(foo);
      } # end of genoud.optim.wrapper101       
     }

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

      ExpandDots  <- function(...)
        {
          return(match.call())
        }      

      dots  <- ExpandDots(...)
      if ( (length(dots) > 1) & balance==TRUE)
        {
          balance  <- FALSE
          warning("load balancing has been turned off because the function to be optimized requires extra arguments")
        }
      
      if (class(cluster)[1]=="SOCKcluster" | class(cluster)[1]=="PVMcluster" | class(cluster)[1]=="spawnedMPIcluster" | class(cluster)[1]=="MPIcluster") {
        clustertrigger=1
        cl <- cluster
        GENclusterExport(cl, "fn")
        GENclusterExport(cl, "fn1")
      } else {
        clustertrigger=2
        cluster <- as.vector(cluster)
        cat("Initializing Cluster\n")
        cl <- makeSOCKcluster(cluster)
        GENclusterExport(cl, "fn")        
        GENclusterExport(cl, "fn1")
      }

      if (length(cl) < 2 )
        {
          warning("You only have one node. You probably shouldn't be using the cluster option")

          #override snow parRapply to work with only one node
          parRapply.1node  <- function (cl, x, fun, ...)
            {
              splitRows  <- function (x, ncl)
                {
                  lapply(list(1:nrow(x)), function(i) x[i, , drop = F])
                }
              
              docall(c, clusterApply(cl, splitRows(x, length(cl)), apply, 1,
                                     fun, ...))
            } #end of parRapply.1node
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
                cl.out <- clusterApplyLB(cl, cl.in, fn1)
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
                if(length(cl) > 1 )
                  {
                    cl.out = parRapply(cl, in.mat, fn1)
                  } else {
                    cl.out = parRapply.1node(cl, in.mat, fn1)
                  } 

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
                as.integer(UserGradient), as.function(gr1func), as.real(P9mix),
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
      hes <- optim(gout[5:(nvars+4)], fn=if(lexical == 1) fn1 else fn1.bfgs,
                   gr=gr1, method="BFGS", hessian=TRUE, control=list(fnscale=g.scale));
                   
      
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

