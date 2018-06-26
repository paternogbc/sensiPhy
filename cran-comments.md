## Resubmission
This is a resubmission. In this version we have:

Fixed issues request by CRAN (see below):

"These fail when checking with env var R_CHECK_LENGTH_1_CONDITION set
to true, which turns using if() with conditions of length greater than
one from a warning into an error."

The package (sensiPhy v0.8.1) was archievied by CRAN because we were not able to fix the above issues in time. The mantainers are really sorry for the delay in fixing these issues but we now believe we have managed to fix all CRAN requests. See details below:

## Test environments
* local Ubuntu 17.10 R 3.4.1
* win-builder (devel and release)
* travis Ubuntu 12.04.5 LTS 

## R CMD check results
There was 3 NOTE and 1 Warning:

## Test results:

### Windows (Status: 1 WARNING, 1 NOTE)
WARNINGS:
1. * checking DESCRIPTION meta-information ... WARNING
Dependence on R version '3.4.1' not with patchlevel 0

_R:_ The following issue was fixed by setting the R Version to 3.4.0 in the description file.

NOTES:
1. New submission
Package was archived on CRAN

_R_: The new submission fixed these issues.

2. Possibly mis-spelled words in DESCRIPTION:
  PCM (13:56)
  clades (14:73)
  intraspecific (16:26)
  topologies (15:43)

_R_: These words are not mis-spelled.

### Debian (Status: 3 NOTEs)

NOTES:
1. * checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Gustavo Paterno <paternogbc@gmail.com>’

New submission

Package was archived on CRAN

Possibly mis-spelled words in DESCRIPTION:
  PCM (13:56)
  clades (14:73)
  intraspecific (16:26)
  topologies (15:43)

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2018-06-01 as check problems were not

_R:_ The new submission fixed these issues.

2. * checking DESCRIPTION meta-information ... NOTE
Dependence on R version ‘3.4.1’ not with patchlevel 0

_R:_ The following issue was fixed by setting the R Version to 3.4.0 in the description file.

3. * checking examples ... [76s/67s] NOTE
Examples with CPU time > 2.5 times elapsed time
               user system elapsed ratio
tree_discrete 2.568    0.7   1.182 2.765

The following example was set to /dontrun{} to avoid long elapsed time (but the code was properly tested in multiple systems).