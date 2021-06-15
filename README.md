# DEEPR
Dirichlet-multinomial evolutionary event randomization test

This program performs a randomization test in R to test whether two groups of cophylogenetic reconstructions are the same or not.  It assumes that individual reconstructions within a group are draws/realizations from a parent (prior) Dirichlet distribution.  Technically, the DEEPR test whether the Dirichlet parameters are the same across groups.  The set of possible coevoluionary event types are: cospeciation, duplication, sorting and host switch.  

Data should be in a flat file with a group ID column and a separate column for each coevolutionary event type.  Each row corresponds to an individual cophylogenetic reconstruction and each cell corresponds to the coevolutionary event counts as estimated by some cophylogenetic reconstruction method (e.g. CoRe-PA).

The test proceed by iteratively (and randomly) permuting group labels and calculating the log-likelihood associated with each permutation under a Dirichlet-multinomial model. While the output is a Monte Caralo exact p-value, each time the test is run on the same data the p-value outputted may differ as it depends on the specific permutations.  The R code depends on the dirmult package in R.  Note that there is a DEEPR package in R, but the version here uses parallel computing in order to speed up calculations.

##################################################################################################
##################################################################################################
#To use parDEEPR on ecto/endo symbiont data in R
#assuming data file is in appropriate folder for 
#uploading.
##################################################################################################

> mydata <- read.csv("EctoEndo_eventcounts.csv", header=T)
> colnames(mydata) <- c("Type",colnames(mydata)[-1])
> head(mydata)
Type            Host       Symbiont eventC eventS eventD eventH                   Source
1 ecto Sphenisciformes   Phthiraptera      9      7      2      3  Banks and Paterson 2004
2 ecto        Primates Sarcoptiformes      4      2      3      1      Bochkov et al. 2011
3 ecto      Chiroptera   Mesostigmata      7      6      1      5  Bruyndonckx et al. 2009
4 ecto   Columbiformes   Phthiraptera      9      7      0      9      Clayton et al. 2003
5 ecto        Rodentia   Phthiraptera      5      4      0      1 Demastes and Hafner 1993
6 ecto        Rodentia   Phthiraptera      9      1      1      6       Hafner et al. 1994

> mydata <- data.frame(mydata)
> grp1 <- mydata[mydata$Type=="ecto", 4:7]
> grp2 <- mydata[mydata$Type=="endo", 4:7]
> test <- parDEEPR(grp1,grp2,number_perms=19)  # just to test, do 19 permutations... takes ~ 6 secs
> test
$p_value
[1] 0.1

$group_1_pi
   eventC     eventS     eventD     eventH 
0.37708553 0.24417008 0.09972235 0.27902204 

$group_2_pi
  eventC    eventS    eventD    eventH 
0.3910896 0.2535684 0.1121499 0.2431922 

$group_1_alphas
[1] 3.3190549 2.1491514 0.8777424 2.4559135

$group_2_alphas
[1] 14.554925  9.436888  4.173810  9.050725

$group_1_theta
[1] 0.1020214

$group_2_theta
[1] 0.02616681

#To run in serial... takes longer

parDEEPR(grp1,grp2,number_perms=19, in.parallel=FALSE)
##################################################################################################
##################################################################################################

