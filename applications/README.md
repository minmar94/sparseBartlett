# Overview
This folder contains the script to reproduce the results of the application, both for gene expressions data (Normal) and doubs data (Poisson).

To obtain the summaries and figures (that must be stored in a folder called "images") from the output (that must be stored in a folder called "out"), you can run the script _output\_realdata.R_ from the command line. This is an example:

`Rscript --vanilla output_realdata.R Gene 50`

where 
   - the first argument is either Gene or Doubs
   - the second argument -- only needed if the first argument is Gene -- must be either 50 or 100

