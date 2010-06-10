".First.lib" <- function(lib, pkg)
{
    cat("\n")
    cat("-------------------------------------------------------------\n")
    pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION",
                                package="geofd"),
                              fields=c("Title","Version","Date")))
    cat(pkg.info["Title"])
    cat("\n")
    cat("For further information of the library go to http://code.google.com/p/geofd\n")
    cat(paste("geofd version ", pkg.info["Version"],
              " (built on ", pkg.info["Date"], ") is now loaded\n", sep=""))
    cat("-------------------------------------------------------------\n")
    cat("\n")

  return(invisible(0))
}

