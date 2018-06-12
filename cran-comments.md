## Resubmission
This is a resubmission. In this version we have:

Fixed issues request by CRAN (see below):

"These fail when checking with env var R_CHECK_LENGTH_1_CONDITION set
to true, which turns using if() with conditions of length greater than
one from a warning into an error."

The package (sensiPhy v0.8.1) was archievied by CRAN because we were not able to fix the above issues in time. The mantainers are really sorry for the delay in fixing these issues but we now believe we have managed to fix all CRAN requests.

## Test environments
* local Ubuntu 17.10 R 3.4.1
* win-builder (devel and release)
* travis Ubuntu 12.04.5 LTS 

## R CMD check results
There were no ERRORs. 

There was 1 NOTE and 1 Warning:

Warning:
CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2018-06-01 as check problems were not
    corrected despite reminders.

__R:__ The mantainers apologize for the delay in fixing the issues raised. We have now fixed all requests and ERRORs reported.

NOTE:
Possibly mis-spelled words in DESCRIPTION:
  PCM (13:56)
  clades (14:73)
  intraspecific (16:26)
  topologies (15:43)

__R:__ These words are not misspelled.