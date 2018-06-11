## R package `lpridge` -- with lpepa() and lpridge()

- Original S-plus interface and Fortran code by Burkhard Seiffert, Biostat Uni Zurich,
  April 1995.
- The port of the lpepa and lpridge to R has been done by Martin Maechler,
  see `DESCRIPTION`.

- A "link" (to Martin Maechler's files) of the original fortran routines and
  S-plus interface are in the R-forge repository
   https://r-forge.r-project.org/scm/viewvc.php/pkg/lpridge/?root=curves-etc
   (source via  `svn checkout svn://svn.r-forge.r-project.org/svnroot/curves-etc/pkg/lpridge` )

 in sub directory

    Orig/

 but not in the package sources (`lpridge_<ver>.tar.gz`) on CRAN.

## Short History:

- MM did get the original sources in April 1995 as shar files,
  separately for `lpepa` and `lpridge`
- found there were too many communalities to keep these separate.

- one of the relatively first CRAN packages, May 18, 2001:

    lpridge_1.0-0.tar.gz
