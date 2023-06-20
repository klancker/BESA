# BESA
Example code for the BESA bioeconomic stock assessment method for data limited stock assessment. 

Download the following files into the same directory in order to run the BESA example: BESA_Stan_NorShr_Demo.R (main R file), BESAStan.stan (Stan model file) , Data.csv (Exemplary data), CMSYfun.R (sub-code to obtain prior distribution parameters, based on the original description in: Froese, R., Demirel, N., Coro, G., Kleisner, K. M., & Winker, H. (2016). Estimating fisheries reference points from catch and resilience. Fish and Fisheries. http://onlinelibrary.wiley.com/doi/10.1111/faf.12190/full.)

BESA was developed using R version 4.0.4 (2021-02-15).

PLEASE NOTE: Stan is currently having compatibility issues with the current rstan CRAN version and R version. If you have trouble compiling the Stan model with rstan, please use the following (and check whether all dependencies are installed correctly): 

remove.packages(c("StanHeaders", "rstan"))
install.packages("StanHeaders", repos = c(https://mc-stan.org/r-packages/, getOption("repos")))
install.packages("rstan", repos = c(https://mc-stan.org/r-packages/, getOption("repos")))


This content is licensed under the terms of the MIT license. 
