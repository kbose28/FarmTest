This is an update fixing C++ problems revealed on Solaris
---------------------------------------------------------

-   sqrt and log functions were previously applied to int. This is
    now corrected.

Test environments
-----------------

-   local OS X install, R 3.3.2
-   Oracle Solaris 10 x86 (using Solaris VMware)
-   Ubuntu 14.04 (on Travis-CI), R 3.3.2
-   win\_builder (devel and release)

R CMD check results
-------------------

There were NO ERRORs, WARNINGs, or NOTEs

Downstream dependencies
-----------------------

There are currently no downstream dependencies for this package
