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

#.First.lib <- function(lib, pkg) library.dynam("rgenoud", pkg, lib)

# use .onLoad instead of .First.lib for use with NAMESPACE and R(>= 1.7.0)
.onLoad <- function(lib, pkg) {
  library.dynam("rgenoud", pkg, lib)
}#end of .onLoad

.onUnload <- function(libpath) {
   library.dynam.unload("rgenoud", libpath)
}


.onAttach <- function( ... )
{
  rgenoudLib <- dirname(system.file(package = "rgenoud"))
  version <- packageDescription("rgenoud", lib = rgenoudLib)$Version
  BuildDate <- packageDescription("rgenoud", lib = rgenoudLib)$Date
  
  cat(paste("##  rgenoud (Version ", version, ", Build Date: ", BuildDate, ")\n", sep = ""))
  cat("##  See http://sekhon.berkeley.edu/rgenoud for additional documentation.\n")
}


