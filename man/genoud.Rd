\name{genoud}
\alias{genoud}
\title{
  GENetic Optimization Using Derivatives
}
\description{
  \code{Genoud} is a function that combines evolutionary algorithm
  methods with a derivative-based (quasi-Newton) method to solve
  difficult optimization problems.  \code{Genoud} may also be used for
  optimization problems for which derivatives do not exist.
  \code{Genoud} solves problems that are nonlinear or perhaps even
  discontinuous in the parameters of the function to be optimized.
  When a statistical model's estimating function (for example, a
  log-likelihood) is nonlinear in the model's parameters, the function
  to be optimized will generally not be globally concave and may have
  irregularities such as saddlepoints or discontinuities.
  Optimization methods that rely on derivatives of the objective
  function may be unable to find any optimum at all.  Multiple local
  optima may exist, so that there is no guarantee that a
  derivative-based method will converge to the global optimum. On the
  other hand, algorithms that do not use derivative information (such
  as pure genetic algorithms) are for many problems needlessly poor at
  local hill climbing.  Most statistical problems are regular in a
  neighborhood of the solution.  Therefore, for some portion of the
  search space, derivative information is useful.  \code{Genoud}, via
  the \code{cluster} option, supports the use of multiple computers,
  CPUs or cores to perform parallel computations.
}

\usage{genoud(fn, nvars, max=FALSE, pop.size=1000, max.generations=100, wait.generations=10,
              hard.generation.limit=TRUE, starting.values=NULL, MemoryMatrix=TRUE, 
              Domains=NULL, default.domains=10, solution.tolerance=0.001,
              gr=NULL, boundary.enforcement=0, lexical=FALSE, gradient.check=TRUE, BFGS=TRUE,
              data.type.int=FALSE, hessian=FALSE, unif.seed=812821, int.seed=53058,
              print.level=2, share.type=0, instance.number=0,
              output.path="stdout", output.append=FALSE, project.path=NULL,
              P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
              cluster=FALSE, balance=FALSE, debug=FALSE, ...)
}
\arguments{
  \item{fn}{
    The function to be minimized (or maximized if
    max=\code{TRUE}).  The first argument of the function must be the
    vector of parameters over which minimizing is to
    occur.  The function must return a scalar result.
    
    For example, if we wish to \emph{maximize} the \code{sin()}
    function.  We can simply call genoud by \code{genoud(sin,
      nvars=1,max=TRUE)}.}
  \item{nvars}{The number of parameters to be selected for the function to be minimized (or
    maximized).}
  \item{max}{
    Maximization (\code{TRUE}) or Minimizing (\code{FALSE}).  Determines
    if \code{genoud} minimizes or maximizes the objective function.}
  \item{pop.size}{
    Population Size.  This is the number of individuals \code{genoud} uses to
    solve the optimization problem.  There are several restrictions on
    what the value of this number can be.  No matter what population
    size the user requests, the number is automatically adjusted to make
    certain that the relevant restrictions are satisfied.  These
    restrictions originate
    in what is required by several of the operators.  In particular,
    operators 6 (Simple Crossover) and 8 (Heuristic
    Crossover) require an even number of individuals to work on---i.e., they
    require two parents.  Therefore, the \code{pop.size} variable and the
    operators sets must be such that these three operators have an even
    number of individuals to work with.  If this does not occur, the
    population size is automatically increased until this constraint is
    satisfied.}
  \item{max.generations}{
    Maximum Generations.  This is the maximum number of generations that
    \code{genoud} will run when attempting to optimize a function.  This is a
    \emph{soft} limit.  The maximum generation limit will be binding for
    \code{genoud} only if \code{hard.generation.limit} has
    been set equal to \code{TRUE}.  If it has not been set equal to
    \code{TRUE}, two soft  triggers control when \code{genoud} stops:
    \code{wait.generations} and \code{gradient.check}. \cr
    
    Although the \code{max.generations} variable is not, by default,
    binding, it is nevertheless important because many operators use it
    to adjust
    their behavior.  In essence, many of the operators become less random
    as the generation count gets closer to the \code{max.generations}
    limit.  If
    the limit is hit and \code{genoud} decides to
    continue working, \code{genoud} automatically increases the
    \code{max.generation}
    limit.\cr
    
    Please see \code{MemoryMatrix} for some important interactions
    with memory management.}
  \item{wait.generations}{
    If there is no improvement in the objective function in this number
    of generations, \code{genoud} will think that it has
    found the optimum.  If the
    \code{gradient.check} trigger has been
    turned on, \code{genoud} will only start counting \code{wait.generations}
    if the gradients are within
    \code{solution.tolerance} of zero.  The
    other variables controlling termination are
    \code{max.generations} and \code{hard.generation.limit}.}
  \item{hard.generation.limit}{
    This logical variable determines if the \code{max.generations}
    variable is a binding constraint for \code{genoud}.  If
    \code{hard.generation.limit} is \code{FALSE}, then \code{genoud} may exceed
    the \code{max.generations} count if either the objective function
    has improved within a given number of generations (determined by
    \code{wait.generations}) or if the gradients are not zero
    (determined by \code{gradient.check}). \cr

    Please see \code{MemoryMatrix} for some important interactions
    with memory management.}
  \item{starting.values}{
    This vector contains the starting values which \code{genoud} will use
    at startup.  The \code{starting.values} vector is a way for the user
    to insert \emph{one} individual into the starting population.
    \code{genoud} will randomly create the other individuals.}
  \item{MemoryMatrix}{
    This variable controls if \code{genoud} sets up a memory matrix.  Such a
    matrix ensures that \code{genoud} will request the fitness evaluation of a
    given set of parameters only once. The variable may
    be \code{TRUE} or \code{FALSE}.  If it is \code{FALSE}, \code{genoud} will
    be aggressive in conserving memory.  The most significant negative
    implication of this variable being set to \code{FALSE} is that
    \code{genoud} will no longer maintain a memory matrix of all evaluated
    individuals.  Therefore, \code{genoud} may request evaluations which it has
    already previously requested. \cr

    Note that when \code{nvars} or \code{pop.size} are large, the memory
    matrix consumes a large amount of RAM. \code{Genoud's} memory matrix will
    require somewhat less memory if the user sets
    \code{hard.generation.limit} equal to \code{TRUE}.}
  \item{Domains}{
    This is a \code{nvars} \eqn{\times 2}{*2}
    matrix.  For each variable, in the first column is the lower bound and
    in the second column the upper bound.  None of \code{genoud}'s
    starting population will be
    generated outside of the bounds.  But some of the operators may
    generate children which
    will be outside of the bounds unless the
    \code{boundary.enforcement} flag is
    turned on. \cr
    
    If the user does not provide any values for Domains, \code{genoud} will setup
    default domains using \code{default.domains}. \cr
    
    For linear and nonlinear constraints please see the discussion in
    the \code{Note} section.}
  \item{default.domains}{
    If the user does not want to provide a \code{Domains} matrix,
    domains may nevertheless be set by the user with this easy to use
    scalar option.  \code{Genoud} will create a
    Domains matrix by setting the lower bound for all of the parameters
    equal to -1 \eqn{\times}{*} \code{default.domains} and the upper
    bound equal to \code{default.domains}.}
  \item{solution.tolerance}{
    This is the tolerance level used by \code{genoud}.  Numbers within
    \code{solution.tolerance} are considered to be equal.  This is
    particularly
    important when it comes to evaluating \code{wait.generations} and
    conducting the \code{gradient.check}.}  
  \item{gr}{A function to provide the gradient for the \code{BFGS}
    optimizer.  If it is \code{NULL}, numerical gradients will be used
    instead.}
  \item{boundary.enforcement}{
    This variable determines the degree to which \code{genoud} obeys the
    boundary constraints.  Notwithstanding the value of the variable,
    none of \code{genoud}'s starting population values will be outside
    of the bounds.  But some of the operators may generate children
    which will be outside of the bounds unless the
    \code{boundary.enforcement} flag is turned on. \cr
    
    \code{boundary.enforcement} has three possible values: 0 (anything goes), 1
    (partial), and 2 (no trespassing):
    \itemize{
      \item{\emph{0: Anything Goes}}{
	This option allows any of the operators to
	create out-of-bounds individuals and these individuals will be
	included in the
	population if their fit values are good enough.  The boundaries are only
	important when generating random individuals.}
      \item{\emph{1: partial enforcement}}{
	This allows operators (particularly those operators
	which use the derivative based optimizer, BFGS) to go out-of-bounds
	during the creation of an individual (i.e., out-of-bounds values
	will often be evaluated).  But when the operator has decided
	on an individual, it \emph{must} be in bounds to be acceptable.}
      \item{\emph{2: No Trespassing}}{
	No out-of-bounds evaluations will ever be requested.}}}
  \item{lexical}{This option enables lexical optimization.  This is
    where there are multiple fit criteria and the parameters are chosen
    so as to maximize fitness values in lexical order---i.e., the second fit
    criterion is only considered if the parameters have the same fit for
    the first.  The fit function used with this option should return a
    vector of fits in the order of precedence.  This option can take on
    the values of FALSE, TRUE or an integer equal to the number of fit
    criterion which are returned by \code{fn}.  The \code{par} object
    which is returned by \code{genoud} will include all of the fit
    criterion at the solution.  The \code{\link[Matching]{GenMatch}}
    function makes extensive use of this option.}
  \item{gradient.check}{
    If this variable is \code{TRUE}, \code{genoud} will not start counting
    \code{wait.generations} unless each gradient is
    \code{solution.tolerance} close to zero.  This
    variable has no effect if the \code{max.generations} limit has been
    hit and the \code{hard.generation.limit} option has been set to
    \code{TRUE}.}
  \item{BFGS}{
    This variable denotes whether or not \code{genoud} applies a
    quasi-Newton derivative optimizer (BFGS) to the best individual at
    the end of each generation after the initial one.  Setting BFGS to
    \code{FALSE} does not mean that the BFGS will never be used.  In
    particular, if you want BFGS never to be used, \code{P9} (the
    Local-Minimum Crossover operator) must also be set to zero.}
  \item{data.type.int}{
    This option sets the data type of the parameters of the function to
    be optimized.  If the variable is \code{TRUE},
    \code{genoud} is
    informed that it is dealing with integer values.  \cr
    
    With integer parameters, \code{genoud} never uses derivative
    information.  This implies that the BFGS quasi-Newton optimizer is
    never used---i.e., the \code{BFGS} flag is set to \code{FALSE}.  It
    also implies
    that Operator 9 (Local-Minimum Crossover) is set to zero and that
    gradient checking (as a convergence criterion) is turned off.  No
    matter what other options have been set to,
    \code{data.type.int} takes precedence---i.e., if \code{genoud} is told that
    it is searching over an integer parameter space, gradient
    information is never considered. \cr
    
    There is no option to mix integer and floating point parameters.  If
    one wants
    to mix the two, it is suggested that the user pick integer type and in
    the objective function map a particular integer range into a
    floating point number range.  For example, tell \code{genoud} to search
    from 0 to 100 and divide by 100 to obtain a search grid of 0 to 1.0
    (by .1).}
  \item{hessian}{
    When this flag is set to \code{TRUE}, \code{genoud} will return the
    hessian matrix at the solution as part of its return list.  A user
    can use this matrix to calculate standard errors.}
  \item{unif.seed}{This sets the seed for the floating-point random
    number generator which \code{genoud} uses.  The default value of this seed
    is 81282. \code{genoud} uses its own internal random number generator
    (a Tausworthe-Lewis-Payne) to allow for recursive and parallel calls
    to \code{genoud}.  It does not use the seed set by the
    \code{\link{set.seed}} function.}
  \item{int.seed}{
    This sets the seed for the integer random number generator which
    \code{genoud} uses. The default value of this seed is
    53058. \code{genoud} uses its own internal random number generator
    (a Tausworthe-Lewis-Payne) to allow for recursive and parallel calls
    to \code{genoud}.  It does not use the seed set by the
    \code{\link{set.seed}} function.}
  \item{print.level}{
    This variable controls the level of printing that \code{genoud} does.  There
    are four possible levels: 0 (minimal printing), 1 (normal), 2
    (detailed), and 3 (debug).  If level 2 is selected, \code{genoud} will
    print details about the population at each generation.  The
    \code{print.level} variable also significantly affects how much
    detail is placed in the project file---see \code{project.path}.
    Note that R convention would have
    us at print level 0 (minimal printing).  However, because \code{genoud}
    runs may take a long time, it is important for the user to receive
    feedback.  Hence, print level 2 has been set as the default.}
  \item{share.type}{
    If \code{share.type} is equal to 1, then \code{genoud}, at startup, checks
    to see if there is an existing project file (see
    \code{project.path}).  If such a file exists, it initializes its
    original population using it. \cr
    
    If the project file contains a smaller population than the current
    \code{genoud} run, \code{genoud} will randomly create the necessary individuals.  If
    the project file contains a larger population than the current \code{genoud}
    run, \code{genoud} will kill the necessary individuals using exponential
    selection. \cr
    
    If the number of variables (see \code{nvars})
    reported in the project file is different from the current \code{genoud} run,
    \code{genoud} does not use the project file (regardless of the value of
    \code{share.type}) and \code{genoud} generates the necessary starting
    population at random.}
  \item{instance.number}{
    This number (starting from 0) denotes the number of recursive
    instances of \code{genoud}.  \code{genoud} then sets up its random number
    generators and other such structures so that the multiple instances
    do not interfere with each other.  It is
    up to the user to make certain that the different instances of
    \code{genoud} are not writing to
    the same output file(s): see \code{output.path} and
    \code{project.path}. \cr
    
    For the R version of \code{genoud} this variable is of limited use.
    It is basically there in case a \code{genoud} run is being used to optimize
    the result of another \code{genoud} run (i.e., a recursive
    implementation).}
  \item{output.path}{
    This is the full (relative) path to where \code{genoud}'s output is to go.
    If the value of \code{output.path}= ``stdout'', then \code{genoud}'s output
    will go to standard output in UNIX and to the GUI console in Windows.
    Also see \code{output.append} and \code{project.path}.}
  \item{output.append}{
    If output is being sent to a file (see \code{output.path}), this logical variable
    tells \code{genoud} whether it should append to the file if it already
    exists or if it should overwrite an existing file.}
  \item{project.path}{
    This is the path of the \code{genoud} project file. The project file prints
    one individual per line with the fit value(s) printed first and then
    the parameter values. By default \code{genoud} places its output in a file
    called "genoud.pro" located in the temporary directory provided by
    \code{\link{tempdir}}.  The behavior of the project file depends on
    the \code{print.level} chosen.  If the \code{print.level} variable
    is set to 1, then the project file is rewritten after each
    generation.  Therefore, only the currently fully completed
    generation is included in the file.  If the \code{print.level}
    variable is set to 2, then each new generation is simply appended to
    the project file.  For all other values of \code{print.level}, the
    project file is not created.}
  \item{P1}{
    This is the cloning operator.  \code{genoud} always clones the best
    individual each generation.  But this operator clones others as
    well.  Please see the Operators Section for details about operators
    and how they are weighted.}
  \item{P2}{
    This is the uniform mutation operator.  One parameter of the parent is
    mutated.  Please see the Operators Section for details about operators and how
    they are weighted.}
  \item{P3}{
    This is the boundary mutation operator.  This operator finds a parent
    and mutates one of its parameters towards the boundary.  Please see
    the Operators Section for details about operators and how
    they are weighted.}
  \item{P4}{
    Non-Uniform Mutation.  Please see the Operators Section for details
    about operators and how they are weighted.}
  \item{P5}{
    This is the polytope crossover.  Please see the Operators Section
    for details about operators and how they are weighted.}
  \item{P6}{
    Simple Crossover.  Please see the Operators Section
    for details about operators and how they are weighted.}
  \item{P7}{
    Whole Non-Uniform Mutation.  Please see the Operators Section for
    details about operators and how they are weighted.}  
  \item{P8}{
    Heuristic Crossover.  Please see the Operators Section for details
    about operators and how they are weighted.}
  \item{P9}{
    Local-Minimum Crossover: BFGS.  This is rather CPU intensive,
    and should be generally used less than the other operators.  Please
    see the Operators Section for details about operators and how they
    are weighted.}
  \item{cluster}{This can either be an object of the  'cluster' class
    returned by one of the \code{\link[snow]{makeCluster}} commands in
    the snow package or a vector of machine names so \code{genoud} can
    setup the cluster automatically.  If it is the latter, the vector should
    look like: \cr
    \code{c("localhost","musil","musil","deckard")}.\cr
    This vector would create a cluster with four nodes: one on the
    localhost another on "deckard" and two on the machine named
    "musil".  Two nodes on a given machine make sense if the machine has two or more
    chips/cores.  \code{genoud} will setup a SOCK cluster by a call to
    \code{\link[snow]{makeSOCKcluster}}.  This will require the user
    to type in her password for each node as the cluster is by default
    created via \code{ssh}.  One can add on usernames to the machine
    name if it differs from the current shell: "username@musil".  Other
    cluster types, such as PVM and MPI,
    which do not require passwords can be created by directly calling
    \code{\link[snow]{makeCluster}}, and then passing the returned
    cluster object to \code{genoud}. For an example of how to manually setup up
    a cluster with a direct call to \code{\link[snow]{makeCluster}} see
    \url{http://sekhon.berkeley.edu/rgenoud/R/genoud_cluster_manual.R}.
    For an example of how to get around a firewall by ssh tunneling see:
    \url{http://sekhon.berkeley.edu/rgenoud/R/genoud_cluster_manual_tunnel.R}.}
  \item{balance}{This logical flag controls if load balancing is done
    across the cluster.  Load balancing can result in better cluster
    utilization; however, increased communication can reduce
    performance.  This option is best used if the function being
    optimized takes at least several minutes to calculate or if the
    nodes in the cluster vary significantly in their performance. If
    cluster==FALSE, this option has no effect.}
  \item{debug}{
    This variable turns on some debugging information.  This variable may
    be \code{TRUE} or \code{FALSE}.}
  \item{...}{Further arguments to be passed to \code{fn} and \code{gr}.}
  }
  \value{
  \code{genoud} returns a list with 7 objects.  8 objects are returned if the
  user has requested the hessian to be calculated at the
    solution.  Please see the \code{hessian} option.  The returned
    objects are:

  \item{value}{
    This variable contains the fitness value at the solution. If
    \code{lexical} optimization was requested, it is a vector.}
  \item{par}{
    This vector contains the parameter values found at the solution.}
  \item{gradients}{
    This vector contains the gradients found at the solution.  If no
    gradients were calculated, they are reported to be \code{NA}.}
  \item{generations}{
    This variable contains the number of generations \code{genoud} ran for.}
  \item{peakgeneration}{
    This variable contains the generation number at which \code{genoud} found
    the solution.}
  \item{pop.size}{
    This variable contains the population size that \code{genoud} actually used.
    See \code{pop.size} for why this value may differ from the
    population size the user requested.}
  \item{operators}{
    This vector reports the actual number of operators (of each type)
    \code{genoud} used.  Please see the Operators Section for details.}    
  \item{hessian}{
    If the user has requested the hessian
      matrix to be returned (via the \code{hessian} flag), the hessian
      at the solution will be returned.  The user may use this matrix to calculate standard
      errors.}    
}
\section{Operators}{
  \code{Genoud} has nine operators that it uses.  The integer values which are
  assigned to each of these operators (P1\eqn{\cdots}{...}P9) are
  weights.
  \code{Genoud} calculates the sum of \eqn{s = P1+P2+\cdots+P9}{s =
  P1+P2+...+P9}.  Each operator is
  assigned a weight equal to \eqn{W_{n} = \frac{s}{P_{n}}}{W_n =
  s/(P_n)}. The number of
  times an operator is called usually equals \eqn{c_{n} = W_{n} \times
  pop.size}{c_n = W_n * pop.size}. \cr
  
  Operators 6 (Simple Crossover) and 8 (Heuristic
  Crossover) require an even number of individuals to work on---i.e.,
  they require two parents.  Therefore, the \code{pop.size} variable and
  the operators sets must be such that these three operators have an
  even number of individuals to work with.  If this does not occur,
  \code{genoud} automatically upwardly adjusts the population size to make this
  constraint hold. \cr
  
  Strong uniqueness checks have been built into the operators to help
  ensure that the operators produce offspring different from their
  parents, but this does not always happen. \cr
  
  Note that \code{genoud} always keeps the best individual each generation. \cr
  
  \code{genoud}'s 9 operators are: \cr
  
  \enumerate{
    \item Cloning
    \item Uniform Mutation
    \item Boundary Mutation
    \item Non-Uniform Crossover
    \item Polytope Crossover
    \item Simple Crossover
    \item Whole Non-Uniform Mutation
    \item Heuristic Crossover
    \item Local-Minimum Crossover: BFGS}
  For more information please see Table 1 of the reference article: 
  \url{http://sekhon.berkeley.edu/papers/rgenoudJSS.pdf}.}
\note{
  The most important options affecting performance are those determining
  population size (\code{pop.size}) and the
  number of generations the algorithm runs
  (\code{max.generations}, \code{wait.generations},
  \code{hard.generation.limit} and \code{gradient.check}).  Search
  performance is expected to improve as
  the population size and the number of generations the program runs for
  increase.  These and the other options should be adjusted for the
  problem at hand.  Please pay particular attention to the search
  domains (\code{Domains} and \code{default.domains}). For more information
  please see the reference article. \cr

  Linear and nonlinear constraints among the parameters can be
  introduced by users in their fit function.  For example, if
  the sum of parameters 1 and 2 must be less than 725, the following can
  be placed in the fit function the user is going to have \code{genoud}
  maximize: \code{if ( (parm1 + parm2) >= 725) \{ return(-99999999) \}}.
  In this example, a very bad fit value is returned to \code{genoud} if the
  linear constrain is violated.  \code{genoud} will then attempt to find
  parameter values that satisfy the constraint.}
\references{
  Walter, Mebane R. Jr. and Jasjeet Singh Sekhon.  2007.  ``Genetic
  Optimization Using Derivatives: The rgenoud package for R.''
  \url{http://sekhon.berkeley.edu/papers/rgenoudJSS.pdf} \cr
  
  Sekhon, Jasjeet Singh and Walter R. Mebane, Jr. 1998.  ``Genetic
  Optimization Using Derivatives: Theory and Application to Nonlinear
  Models.'' \emph{Political Analysis}, 7: 187-210.
  \url{http://sekhon.berkeley.edu/genoud/genoud.pdf} \cr

  Mebane, Walter R., Jr. and Jasjeet S. Sekhon. 2004. ``Robust
  Estimation and Outlier Detection for Overdispersed Multinomial Models
  of Count Data.''  \emph{American Journal of Political Science}, 48
  (April): 391-410. \url{http://sekhon.berkeley.edu/multinom.pdf}
}
\author{Walter R. Mebane, Jr., Cornell University,
  \email{wrm1@cornell.edu}, \url{http://macht.arts.cornell.edu/wrm1/} \cr

  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}}
\seealso{
  \code{\link{optim}}.
}
\examples{
#maximize the sin function
 sin1 <- genoud(sin, nvars=1, max=TRUE)

#minimize the sin function
 sin2 <- genoud(sin, nvars=1, max=FALSE)

#maximize a univariate normal mixture which looks like a claw
 claw <- function(xx) {
   x <- xx[1]
   y <- (0.46*(dnorm(x,-1.0,2.0/3.0) + dnorm(x,1.0,2.0/3.0)) +
   (1.0/300.0)*(dnorm(x,-0.5,.01) + dnorm(x,-1.0,.01) + dnorm(x,-1.5,.01)) +
   (7.0/300.0)*(dnorm(x,0.5,.07) + dnorm(x,1.0,.07) + dnorm(x,1.5,.07))) 
   return(y)
 }
 claw1   <- genoud(claw, nvars=1,pop.size=3000,max=TRUE)

\dontrun{
#Plot the previous run
 xx <- seq(-3,3,.05)
 plot(xx,lapply(xx,claw),type="l",xlab="Parameter",ylab="Fit",main="GENOUD: Maximize the Claw Density")
 points(claw1$par,claw1$value,col="red")

# Maximize a bivariate normal mixture which looks like a claw.
 biclaw <- function(xx) {
  mNd2 <- function(x1, x2, mu1, mu2, sigma1, sigma2, rho)
    {
      z1 <- (x1-mu1)/sigma1
      z2 <- (x2-mu2)/sigma2
      w <- (1.0/(2.0*pi*sigma1*sigma2*sqrt(1-rho*rho))) 
      w <- w*exp(-0.5*(z1*z1 - 2*rho*z1*z2 + z2*z2)/(1-rho*rho)) 
      return(w)
    }
  x1 <- xx[1]+1
  x2 <- xx[2]+1
  
  y <- (0.5*mNd2(x1,x2,0.0,0.0,1.0,1.0,0.0) +
	    0.1*(mNd2(x1,x2,-1.0,-1.0,0.1,0.1,0.0) +
		 mNd2(x1,x2,-0.5,-0.5,0.1,0.1,0.0) +
		 mNd2(x1,x2,0.0,0.0,0.1,0.1,0.0) +
		 mNd2(x1,x2,0.5,0.5,0.1,0.1,0.0) +
		 mNd2(x1,x2,1.0,1.0,0.1,0.1,0.0)))

  return(y)
 }
 biclaw1 <- genoud(biclaw, default.domains=20, nvars=2,pop.size=5000,max=TRUE)
}
# For more examples see: http://sekhon.berkeley.edu/rgenoud/R.
}
\keyword{optimize}
\keyword{nonlinear}
%  LocalWords:  genoud GENetic GAs fn nvars max MemoryMatrix BFGS hessian unif
%  LocalWords:  stdout Polytope cr emph eqn GENOUD's peakgeneration dfrac optim
%  LocalWords:  seealso saddlepoints lapply xlab ylab biclaw mNd polytope cdots
%  LocalWords:  frac url jasjeet RGENOUD PVM MPI parm Overdispersed UC dontrun
%  LocalWords:  bivariate CPUs startup Genoud's GenMatch tempdir makeCluster
%  LocalWords:  localhost musil deckard makeSOCKcluster usernames univariate
