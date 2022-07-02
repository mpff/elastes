## Resubmission

This is a resubmission. In this version I have:

* Changed T/F to TRUE/FALSE in function's default parameters.
* Documented missing return values in exported functions.
* Removed 'par()' call in 'elastes.Rmd' vignette.


## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Manuel Pfeuffer <mnl.pfeuffer@gmail.com>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    Greven (15:69)
    Steyer (15:51)
    al (16:46)
    cker (15:63, 16:37)
    et (16:42)
    
 -> These words are not misspelled.
  
  Found the following (possibly) invalid URLs:
    URL: https://doi.org/10.1111/biom.13706
      From: README.md
      Status: 503
      Message: Service Unavailable
      
 -> HTTP site https://doi.org/10.1111/biom.13706 works in a browser.
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.1111/biom.13706
      From: DESCRIPTION
      Status: Service Unavailable
      Message: 503
  
 -> HTTP site https://doi.org/10.1111/biom.13706 works in a browser.
    

## Downstream dependencies

There are no downstream dependencies for elastes.
