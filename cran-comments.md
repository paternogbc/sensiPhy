## Test environments
* local ubuntu 16.04 install, R 3.2.3
* local Windows XXXXX install, R XXXX
* local OS XXXX install, R XXXXX
* win-builder (devel and release)
* travis Ubuntu 12.04.5 LTS 

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTES:

* New submission  
Possibly mis-spelled words in DESCRIPTION:  
  PGLS (13:20)
  Phylogenetic (15:10)
  clades (14:68)
  intraspecific (16:29)
  phylogenetic (10:59, 11:34)
  
  R: These words are not mis-spelled.
  
* checking DESCRIPTION meta-information ... NOTE
Malformed Description field: should contain one or more complete sentences.

* checking R code for possible problems ... NOTE
  print.data.phy: no visible global function definition for 'str'
  Undefined global functions or variables:
