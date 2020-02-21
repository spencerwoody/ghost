
# Code for model summary and projections paper

["Model interpretation through lower-dimensional posterior
summarization"][arXiv]  
by Spencer Woody, Carlos M. Carvalho, and Jared
S. Murray

## Contents

- `toy-sigmoid-example.R` contains code for the toy example in Section
  XX which illustrates how our summarization approach can approximate
  the partial effects in nonparametric regression models. 
- `sim-itx-search.R` contains code for the simulation example in
  Section XX which demonstrates our method for discovering the most
  significant interactions present in nonparametric regression models.
  This script creates 20 simulations, while `sim-itx-search-big.R`
  contains 1000 simulations.
- `crime-example.R` contains code for the linear summary of linear
  regression model used the crime example regression analysis in
  Section XX.
- `CA-example-[xxx].R` these scripts replicate the extensive case study
  of Section 4, summarizing a nonparametric regression model for
  housing prices in California at the census tract level. 
    - `CA-example-01dataprep.R` reads in and transforms the data.
    - `CA-example-02GPfit.R` fits a Gaussian process regression model
      to the data. 
    - `CA-example-03GPposterior.R` computes MCMC draws from the
      posterior of the GP regression model and the error variance
      parameter.
    - `CA-example-04globalsummaries.R` computes global summaries
      (linear, additive, additive with bivariate interaction), and the
      interaction search mentioned in Section 4.1. 
    - `CA-example-05localsummaries.R` computes local linear summaries
      for the GP regression model, at various resolutions of locality
      (metropolitan area, county, city, and neighborhood). 
- `R/` contains R functions for the rest of the analysis.
- `src/` contains source C++ functions used in R scripts.
- `data/` contains the California housing price data.
<!-- - `figures` is an empty directory for storing generated graphics. -->
<!-- - `R-output` is an empty directory for storing intermediate R output, -->
<!--   since the California analysis contains multiple steps. -->

Details for these analyses can be found in the original paper. 

## Contact

[Spencer Woody][Personal website]  
Statistics PhD Candidate, UT-Austin  
Email: `spencer.woody@utexas.edu`

[arXiv]: https://arxiv.org/abs/1905.07103 
[Personal website]: https://spencerwoody.github.io/

