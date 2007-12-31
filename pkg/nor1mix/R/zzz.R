## no S4 methodology here; speedup :
.noGenerics <- TRUE

".First.lib" <- function(lib, pkg) {
    if(interactive() || getOption("verbose")) { # not in test scripts
	packageStartupMessage(sprintf("Package %s (%s) loaded.  To cite, see citation(\"%s\")\n",
			pkg, packageDescription(pkg)$Version, pkg))
	warning("'MW.nm2' has been changed (linear transform only) to match\n",
		"the Annals paper. 'MW.nm2.old' is the former version",
		call. = FALSE, immediate. = TRUE)
    }
}
