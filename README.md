# DEEPR
Dirichlet-multinomial evolutionary event randomization test

This program performs a randomization test in R to test whether two groups of cophylogenetic reconstructions are the same or not.  It assumes that individual reconstructions within a group are draws/realizations from a parent (prior) Dirichlet distribution.  Technically, the DEEPR test whether the Dirichlet parameters are the same across groups.  The set of possible coevoluionary event types are: cospeciation, duplication, sorting and host switch.  

Data should be in a flat file with a group ID column and a separate column for each coevolutionary event type.  Each row corresponds to an individual cophylogenetic reconstruction and each cell corresponds to the coevolutionary event counts as estimated by some cophylogenetic reconstruction method (e.g. CoRe-PA).

The test proceed by iteratively (and randomly) permuting group labels and calculating the log-likelihood associated with each permutation under a Dirichlet-multinomial model. While the output is a Monte Caralo exact p-value, each time the test is run on the same data the p-value outputted may differ as it depends on the specific permutations.  The R code depends on the dirmult package in R.  Note that there is a DEEPR package in R, but the version here uses parallel computing in order to speed up calculations.
