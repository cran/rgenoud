#  RGENOUD
#
#  Walter R. Mebane, Jr.                                        
#  Cornell University
#  http://macht.arts.cornell.edu/wrm1
#  wrm1@macht.arts.cornell.edu
#
#  Jasjeet Singh Sekhon 
#  Harvard University
#  http://jsekhon.fas.harvard.edu/
#  jsekhon@fas.harvard.edu
#
#  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/zzz.R,v 1.31 2005/03/01 06:36:36 jsekhon Exp $
#

#.First.lib <- function(lib, pkg) library.dynam("rgenoud", pkg, lib)

# use .onLoad instead of .First.lib for use with NAMESPACE and R(>= 1.7.0)
.onLoad <- function(lib, pkg) {
  library.dynam(pkg, pkg, lib)
}#end of .onLoad

.onUnload <- function(libpath) {
   library.dynam.unload("rgenoud", libpath)
}
