# Overview
This folder contains the script to reproduce the results of the simulation study, both for Normal and Poisson likelihoods.

To obtain the summaries and figures (that must be stored in a folder called "images") from the output (that must be stored in a folder called "out"), you can run the script _output\_normal.R_ or _output\_poisson.R_ from the command line. This is an example:

`Rscript --vanilla output_normal.R Banded`

where the first and only argument is either Banded or Random
