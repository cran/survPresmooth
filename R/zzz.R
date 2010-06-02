.First.lib <- function(lib, pkg){
	library.dynam("survPresmooth", pkg, lib)
}

.Last.lib <- function (libpath){
  library.dynam.unload("survPresmooth", libpath=libpath) 
}
