\name{genoud}
\alias{genoud}
\title{
  GENetic Optimization Using Derivatives
}
\description{
  GENOUD is a function that combines evolutionary algorithm methods with
  a derivative-based, quasi-Newton method to solve difficult
  unconstrained optimization problems.  GENOUD is made to solve problems
  that are nonlinear or perhaps even discontinuous in the parameters of
  the function to be optimized.  When a statistical model's estimating
  function (for example, a log-likelihood) is nonlinear in the model's
  parameters, the function to be optimized will usually not be globally
  concave and may contain irregularities such as saddlepoints or
  discontinuous jumps.  Optimization methods that rely on derivatives of
  the objective function may be unable to find any optimum at all.
  Multiple local optima may exist, so that there is no guarantee that a
  derivative-based method will converge to the global optimum.  On the
  other hand, algorithms that do not use derivative information (such as
  pure GAs) are for many problems needlessly poor at local hill
  climbing.  Most statistical problems are regular in the neighborhood
  of the solution.  Therefore, for some portion of the search space
  derivative information is useful.
  }

\usage{ genoud
  (fn, nvars, max=FALSE, pop.size=100, max.generations=100, wait.generations=10,
   hard.generation.limit=TRUE, starting.values=NULL, MemoryMatrix=NULL, Debug=FALSE, 
   Domains=NULL, default.domains=10,
   gradient.check=TRUE, boundary.enforcement=2,
   solution.tolerance=0.001, BFGS=TRUE, data.type.int=FALSE, hessian=FALSE,
   unif.seed=812821, int.seed=53058,
   print.level=2, share.type=0, instance.number=0,
   output.path="stdout", output.append=FALSE, project.path="genoud.pro",
   P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0)
}
\arguments{
  \item{fn}{
    The function to be minimized by default (or maximized if max=T), with
    first argument the vector of parameters over which minimizing is to
    take place.  It should return a scalar result.
    
    For example, say we wish to maximize the \code{sin()} function.  We
    can simply call genoud by \code{genoud(sin, nvars=1)}.}
  \item{nvars}{
    This is the number of variables the function to be minimized (or
    maximized) takes.}
  \item{max}{
    Maximization (True) or Minimizing (False).  This variable tells GENOUD
    if it is to minimize or maximize the objective function.}
  \item{pop.size}{
    Population Size.  This is the number of individuals GENOUD uses to
    solve the optimization problem.  There are several restrictions on
    what the value of this number can be.  These constraints originate
    in what is required by several of the operators.  In particular,
    operators 5 (Polytope Crossover), 6 (Multiple Point Simple
    Crossover) and 8 (Heuristic Crossover) require an even number of individuals
    to work on---i.e., they require two parents.  Therefore,
    the \code{pop.size} variable and the operators sets must be such
    that these three operators have
    an even number of individuals to work with.  If this does not occur,
    the population size is increased until this is so.}
  \item{max.generations}{
    Maximum Generations.  This is the maximum number of generations that
    GENOUD will run when attempting to optimize a function.  This is a
    \emph{soft} limit.  The maximum generation limit will be binding for
    GENOUD only if \code{hard.generation.limit} has
    been set equal to TRUE.  If it has not been set equal to TRUE,, two
    soft  triggers control when GENOUD stops:
    \code{wait.generations} and \code{gradient.check}. \cr
    
    Although, the \code{max.generations} variable is not, by default,
    binding, it is nevertheless important because many operators use it
    to adjust
    their behavior.  In essence, many of the operators become less random
    as the generation count gets closer to the max.generations limit.  If
    the limit is hit and GENOUD, because of the soft triggers, decides to
    continue working, GENOUD automatically increases the max.generation
    limit.  GENOUD will print to the screen what it has decided to do. \cr
    
    Please see \code{MemoryMatrix} for some interesting interactions
    with memory management.}
  \item{wait.generations}{
    If there is no improvement in the objective function of
    wait.generations number of generations, GENOUD will think that it has
    found the optimum.  If the
    \code{gradient.check} trigger has been
    turned on, GENOUD will only start counting wait.generations if the
    gradients are zero (within
    \code{solution.tolerance}).  The
    other variable controlling termination is
    \code{max.generations}.}
  \item{hard.generation.limit}{
    By default the \code{max.generations}
    variable is a binding constraint for GENOUD.  But it fails to be
    binding if the hard.generation.limit is set to FALSE. \cr

    Please see \code{MemoryMatrix} for some
    interesting interactions with memory management.}
  \item{starting.values}{
    This is a vector contains the starting values that GENOUD will use
    at startup.  The starting.values vector is a way for the user to
    insert \emph{one} individual into the starting population.
    GENOUD will randomly create the other individuals in the
    population.}
  \item{MemoryMatrix}{
    This variable controls how aggressive GENOUD is in saving memory.  The
    variable my be \eqn{T}{T} or \eqn{F}{F}.  If it is \eqn{F}{F} GENOUD
    will
    be aggressive in
    conserving memory.  The most significant implication of this variable
    being set to \eqn{F}{F} is that GENOUD will no longer maintain a memory
    matrix off all evaluated individuals.  Therefore, GENOUD may request
    evaluations which it has already previously requested.  When
    \code{nvars} is large---i.e., \eqn{> 50}{> 50}, the memory
    matrix consumes a large amount of RAM.  Therefore, it is highly
    recommended that will systems with 512MB of RAM, this variable be set
    to \eqn{F}{F} if the number of variables being optimized is greater
    than 50. \cr
    
    By default GENOUD sets MemoryMatrix equal to \eqn{T}{T} unless the
    number of
    parameters is greater than 20.  In that case, GENOUD sets the default
    value equal to \eqn{F}{F}. \cr
    
    GENOUD's MemoryMatrix will usually require \emph{significantly} less
    memory if the user sets
    \code{hard.generation.limit} equal
    to TRUE.  Doing this is usually a good way of conserving memory while
    still making use of the MemoryMatrix structure.}
  \item{Debug}{
    This variable turns on some debugging information.  This variable my
    be \eqn{T}{T} or \eqn{F}{F}.}
  \item{Domains}{
    This is a \code{nvars} \eqn{\times 2}{*2}
    matrix.  The fist column is the lower bound, and the second column is
    the upper bound.  None of GENOUD's starting population will be outside
    of the bounds.  But some of the operators may generate children which
    will be outside of the bounds unless the
    \code{boundary.enforcement} flag is
    turned on. \cr
    
    If the user does not provide any values for Domains, GENOUD will setup
    default domains using \code{default.domains}.}
  \item{default.domains}{
    If \code{Domains}, GENOUD will create a
    Domains matrix by setting the lower bound for all of the parameters
    equal to -1 \eqn{\times}{*} \code{default.domains} and the upper
    bound equal to \code{default.domains}.}
  \item{gradient.check}{
    If this variable is \eqn{T}{T}, GENOUD will not start counting
    \code{wait.generations} unless the gradients are zero or
    \code{solution.tolerance} close to
    zero.  This variable has not effect if the
    \code{max.generations} limit is hit.}
  \item{boundary.enforcement}{
    This variable determines the degree to which GENOUD obeys the
    \code{boundary constraints}.  Notwithstanding
    the value of the variable, none of GENOUD's starting population will
    be outside of the bounds.  But some of the operators may generate
    children which will be outside of the bounds unless the
    boundary.enforcement flag is turned on. \cr
    
    \code{boundary.enforcement} has three possible values: 0 (anything goes), 1
    (regular), and 2 (no trespassing):
    \itemize{
      \item{\emph{0: Anything Goes}}{
	This option allows any of the operators to
	create out-of-bounds individuals and these will be included in the
	population if their fit is good enough.  No special attention is
	payed to the boundaries after random individual creation (which is
	guaranteed to create in-bounds individuals).}
      \item{\emph{1: Regular}}{
	This allows operators (particularly those operators
	which use the derivative based optimizer, BFGS) to go out-of-bounds
	during the creation of an individual (i.e., out-of-bounds) values
	will often have to be evaluated.  But when the operator has decided
	on an individual, it \emph{must} be in bounds to be accepted. \cr}
      \item{\emph{2: No Trespassing}}{
	No out-of-bounds evaluations will ever be requested.}}}
  \item{solution.tolerance}{
    This is the tolerance level used by GENOUD.  Numbers within
    solution.tolerance are considered to be equal.  This is particularly
    important when it comes to evaluating \code{wait.generations} and
    conducting the \code{gradient.check}.}
  \item{BFGS}{
    This variable denotes whether or not GENOUD applies a quasi-Newton
    derivative optimizer (BFGS) to the best individual at the end of each
    generation after the initial one.  Turning BFGS to false does not mean
    that the BFGS will never be used.  In particular, Operator 9
    (Local-Minimum Crossover) has to be set to zero.}
  \item{data.type.int}{
    This variable informs GENOUD of the data type of the parameters that
    the function to be optimized uses.  If the variable is \eqn{T}{T},
    GENOUD is
    informed that it is dealing with integer values.  \emph{Use of the
      integer data type is supported only as a beta feature.}  Some of the
    included operators will not work well with the integer data type. \cr
    
    When the integer type is used, GENOUD never uses derivative
    information.  This implies that the BFGS quasi-Newton optimizer is
    never used---i.e., the \code{BFGS} flag is false.  It also implies
    that Operator 9 (Local-Minimum Crossover) is set to zero and that
    gradient checking (as a convergence criterion) is turned off.  No
    matter what the values of these other variables is, the value of
    data.type.int takes precedence---i.e., if GENOUD is told that its is
    searching over an integer parameter space all everything that
    requires gradient information is turned off. \cr
    
    There is no option to mix the two types of data.  If one wants to
    mix the two, it is suggested that the user pick integer type and in
    her objective function map a particular integer range into a
    floating point number range.  For example, tell GENOUD to search
    from 0 to 100 and divide by 100 to obtain a search grid of 0 to 1.0
    (by .1).}
  \item{hessian}{
    When this flag is set to TRUE, GENOUD will return the hessian matrix
    at the solution as part of its return list.  A user can use this
    matrix to calculate standard errors.}
  \item{unif.seed}{
    This is the first of the two random number seeds which GENOUD uses.
    The default value of this seed is 81282.}
  \item{int.seed}{
    This is the second of the two random number seeds which GENOUD uses.
    The default value of this seed is 53058.}
  \item{print.level}{
    This variable controls the level of printing that GENOUD does.  There
    are three possible levels: 0 (minimal printing), 1 (normal), and 2
    (detailed).  If level 2 is selected, GENOUD will print details about
    the population at each generation.  The print.level variable also
    significantly affects how much detail is placed in the Project
    File---see \code{project.path}.}
  \item{share.type}{
    If share.type is equal to 1, then GENOUD, at startup, checks to see if
    there is an existing Project File (see \code{project.path}).  If such a file exists,
    it initializes its original population using it. \cr
    
    If the Project File contains a smaller population than the current
    GENOUD run, GENOUD will randomly create the necessary individuals.  If
    the Project File contains a larger population than the current GENOUD
    run, GENOUD will kill the necessary individuals using exponential
    selection. \cr
    
    If the number of variables (see \code{nvars})
    reported in the Project File is different from the current GENOUD run,
    GENOUD does not use the Project File (regardless of the value of
    share.type) and generates the necessary starting population at
    random.}
  \item{instance.number}{
    Number of Concurrent instances of GENOUD.  This number (starting from
    0) denotes the number of parallel running instances of GENOUD.  GENOUD
    then sets up its random number generators and other such structures so
    that the multiple instances do not interfere with each other.  It is
    up to the user to make certain that GENOUD is not always writing to
    the same output file(s): see \code{output.path} and
    \code{project.path}. \cr
    
    For the R version of GENOUD, this variable is of greatly limited use.
    It is basically there in case GENOUD is being used to optimize
    GENOUD.}
  \item{output.path}{
    This is the full (relative) path to where GENOUD's output is to go.
    If the value of output.path= ``stdout'', then GENOUD's output will go
    to standard output.  Also see \code{output.append} and \code{project.path}.}
  \item{output.append}{
    If output is being sent to a file (see \code{output.path}), this logical variable
    tells GENOUD whether it should append to the file if it already
    exists, or if it should overwrite an existing file.}
  \item{project.path}{
    This is the path for the GENOUD project file.  The behavior of the
    Project File depends on the \code{print.level} chosen.  If the
    \code{print.level} variable is set to
    1, then the project file is rewritten after each generation.
    Therefore, only the currently fully completed generation is included
    in the file.  If the \code{print.level} variable is set
    to 2, then each new generation is simply appended to the Project File.
    For all other values of \code{print.level}, the Project File is not created.}
  \item{P1}{
    This is the cloning operator.  GENOUD always clones the best
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
    This is the polytope crossover.  This operator requires two
    parents.  Please see the Operators Section for details about
    operators and how they are weighted.}
  \item{P6}{
    Multiple Point Simple Crossover.  Please see the Operators Section
    for details about operators and how they are weighted.}
  \item{P7}{
    Whole Non-Uniform Mutation.  Please see the Operators Section for
    details about operators and how they are weighted.}  
  \item{P8}{
    Heuristic Crossover.  Please see the Operators Section for details
    about operators and how they are weighted.}
  \item{P9}{
    Local-Minimum Crossover: BFGS.  This is \emph{very} CPU intensive,
    and should be rarely used.  Please see the Operators Section for
    details about operators and how they are weighted.}  
}
\value{
  GENOUD returns a list with 7 objects.  8 objects are returned if the
  user as requested the hessian to be calculated at the
    solution---see the \code{hessian} variable.  They are:

  \item{value}{
    This variable contains the fitness value at the solution.}
  \item{generations}{
    This variables contains the number of generations GENOUD ran for.}
  \item{peakgeneration}{
    This variables contains the generation number at which GENOUD found
    the solution.}
  \item{pop.size}{
    This variable contains the population size that GENOUD actually used.
    See \code{pop.size} for why this value may differ from the
    population size the user requested.}
  \item{par}{
    This vector contains the parameter values found at the solution.}
  \item{gradients}{
    This vector contains the gradients found at the solution.  If no
    gradients were calculated, they are reported to by \eqn{-1.0}{-1.0}.}
  \item{operators}{
    This vector reports the actual number of operators (of each type)
    GENOUD used.  Please see the Operators Section for details.}    
  \item{hessian}{
    If the user has requested the hessian
      matrix to be returned---via the \code{hessian} flag, the hessian
      at the solution will be returned.  The user may use this matrix to calculate standard
      errors.}    
}
\section{Operators}{
  GENOUD has nine operators that it uses.  The integer values which are
  assigned to each of these operators (P1\eqn{\cdots}{...}P9) are
  weights.
  GENOUD calculates the sum of \eqn{s = P1+P2+\cdots+P9}{s =
  P1+P2+...+P9}.  Each operator is
  assigned a weight equal to \eqn{W_{n} = \frac{s}{P_{n}}}{W_n =
  s/(P_n)}. The number of
  times an operator is called usually equals \eqn{c_{n} = W_{n} \times
  pop.size}{c_n = W_n * pop.size}. \cr
  
  Operators 5 (Polytope Crossover), 6 (Multiple Point Simple Crossover)
  and 8 (Heuristic Crossover) require an even number of individuals
  to work on---i.e., they require two parents.  Therefore, the
  \code{pop.size} variable and the operators sets must be such that
  these three operators have an even number of individuals to
  work with.  If this does not occur, GENOUD (upwardly) adjusts the
  population size to make this constraint hold. \cr
  
  Strong uniqueness checks have been built into the operators to help
  ensure that the operators produce offspring different from their
  parents, but this does not always happen! \cr
  
  Note that GENOUD always keeps the best individual each generation. \cr
  
  GENOUD's 9 operators are: \cr
  
  \itemize{
    \item{1}{Cloning}
    \item{2}{Uniform Mutation}
    \item{3}{Boundary Mutation}
    \item{4}{Non-Uniform Crossover}
    \item{5}{Polytope Crossover}
    \item{6}{Multiple Point Simple Crossover}
    \item{7}{Whole Non-Uniform Mutation}
    \item{8}{Heuristic Crossover}
    \item{9}{Local-Minimum Crossover: BFGS}
  }}
\references{
  Sekhon, Jasjeet Singh and Walter R. Mebane, Jr. 1998.  ``Genetic
  Optimization Using Derivatives: Theory and Application to Nonlinear
  Models.'' \emph{Political Analysis} 7: 187-210.
  HTTP://jsekhon.fas.harvard.edu/genoud
}
\seealso{
  \code{\link{optim}}.
}
\examples{
#maximize the sin function
 sin1 <- genoud(sin, nvars=1, max=TRUE);

#minimize the sin function
 sin2 <- genoud(sin, nvars=1, max=FALSE);

#maximimize a univariate normal mixture which looks like a claw
claw <- function(xx) {
  Nd <- function(x, mu, sigma) {
    w <- (1.0/sqrt(2.0*pi*sigma*sigma)) ;
    z <- (x-mu)/sigma;
    w <- w*exp(-0.5*z*z) ;
    as.double(w);
  }
  x <- xx[1];
  y <- (0.46*(Nd(x,-1.0,2.0/3.0) + Nd(x,1.0,2.0/3.0)) +
	   (1.0/300.0)*(Nd(x,-0.5,.01) + Nd(x,-1.0,.01) + Nd(x,-1.5,.01)) +
	   (7.0/300.0)*(Nd(x,0.5,.07) + Nd(x,1.0,.07) + Nd(x,1.5,.07))) ;
  as.double(y);
}
 claw1   <- genoud(claw, nvars=1,P9=100,max=TRUE);

#maximimize a bivariate normal mixture which looks like a claw
biclaw <- function(xx) {
  mNd2 <- function(x1, x2, mu1, mu2, sigma1, sigma2, rho)
    {
      z1 <- (x1-mu1)/sigma1;
      z2 <- (x2-mu2)/sigma2;
      w <- (1.0/(2.0*pi*sigma1*sigma2*sqrt(1-rho*rho))) ;
      w <- w*exp(-0.5*(z1*z1 - 2*rho*z1*z2 + z2*z2)/(1-rho*rho)) ;
      as.double(w);
    }
  x1 <- xx[1]+1;
  x2 <- xx[2]+1;
  
  y <- (0.5*mNd2(x1,x2,0.0,0.0,1.0,1.0,0.0) +
	    0.1*(mNd2(x1,x2,-1.0,-1.0,0.1,0.1,0.0) +
		 mNd2(x1,x2,-0.5,-0.5,0.1,0.1,0.0) +
		 mNd2(x1,x2,0.0,0.0,0.1,0.1,0.0) +
		 mNd2(x1,x2,0.5,0.5,0.1,0.1,0.0) +
		 mNd2(x1,x2,1.0,1.0,0.1,0.1,0.0)));

  as.double(y);
}
 biclaw1 <- genoud(biclaw, nvars=2,P9=100,max=TRUE);
}
\keyword{optimize}
\keyword{nonlinear}
%  LocalWords:  genoud GENetic GAs fn nvars max MemoryMatrix BFGS hessian unif
%  LocalWords:  stdout Polytope cr emph eqn GENOUD's peakgeneration dfrac optim
%  LocalWords:  seealso
