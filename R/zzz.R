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
#  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/zzz.R,v 1.25 2004/03/03 22:56:19 jsekhon Exp $
#

#.First.lib <- function(lib, pkg) library.dynam("rgenoud", pkg, lib)

# use .onLoad instead of .First.lib for use with NAMESPACE and R(>= 1.7.0)
.onLoad <- function(lib, pkg) {
  library.dynam(pkg, pkg, lib)
}#end of .onLoad

.onUnload <- function(libpath) {
   library.dynam.unload("rgenoud", libpath)
}
