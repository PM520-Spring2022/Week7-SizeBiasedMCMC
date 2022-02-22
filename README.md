# Week7-SizeBiasedMCMC

This is the repo for exercise 2 of week 7. the instructions were as follows:

Code up the Metropolis-Hastings MCMC algorithm for sized-biased sampling from an arbitrary distribution over the range [0,10].

Test your algorithm by using it to construct sized-biased samples from an exp(𝝀) distribution. 

Plot your results and compare them to the curve x𝝀e-𝝀x. 

Recall: Sized-biased samples from a density f(x) have density xf(x), rather than f(x)

Note that your algorithm will work for any density function f.

Pseudocode is found in the file SizeBiasedMCMCPseudocode.R

An exmaple of the completed code, along with the implemnentation of the Gelman-Rubin stationarity diagnostic can be found in the file SizeBiasedMCMCWorking.
