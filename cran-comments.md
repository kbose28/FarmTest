Test environments
-----------------

-   local OS X install, R 3.3.2
-   Ubuntu 14.04 (on Travis-CI), R 3.3.2
-   win-builder (devel and release)

R CMD check results
-------------------

No NOTES, ERRORs or WARNINGs

Downstream dependencies
-----------------------

I have also run R CMD check on downstream dependencies using
devtools::revdep\_check() followed by
devtools::revdep\_check\_save\_summary(). Results at
<https://github.com/kbose28/farmtest/blob/master/revdep/problems.md>.
There were no problems.
