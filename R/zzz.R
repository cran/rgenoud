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
#  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/zzz.R,v 2.15 2005/10/29 06:14:44 jsekhon Exp jsekhon $
#

#.First.lib <- function(lib, pkg) library.dynam("rgenoud", pkg, lib)

# use .onLoad instead of .First.lib for use with NAMESPACE and R(>= 1.7.0)
.onLoad <- function(lib, pkg) {
  library.dynam("rgenoud", pkg, lib)
}#end of .onLoad

.onUnload <- function(libpath) {
   library.dynam.unload("rgenoud", libpath)
}
