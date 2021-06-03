######################################################
##
##  DEEPR code written in parallel
##  modified from DEEPR package code written by
##  Mark Merilo in 2015.
##
##  written by: Ayesha Ali  (June 12, 2019)
##  updated by: Augustine Wigle
##
##  last update: June 14 2019
########################################################

parDEEPR <- function (group_1, group_2, perm_number) 
{
 ###################################################################
 ## This function takes in coevolutionary counts from 
 ## each group, assuming group counts appear in same
 ## order, and performs a DEEPR test 
 ## (Dirichlet-multinomial evolutionary event profile
 ## randomization test). Null hypothesis is that counts
 ## were generated from multinomial distributions with the 
 ## same parent Dirichlet distribution.
 ##
 ##  group_1 = matrix of counts from reference cophylogenetic group 
 ##  group_2 = matrix of counts from comparison cophylogenetic group
 ##  perm_number = number of permutations to use
 ##
 ##  This function returns: 
 ##  p_value     = p-value of the test
 ##  group*_pi   = expected multinomial proportions for group *
 ##  group*_alpha= estimated parameters of parent Dirichlet 
 ##                  distribution for group *
 ##  group*_theta= dispersion parameter associated with Dirichlet
 ##                  distribution for group *
 ##
 ## DEEPR uses parllel computing to randomly permute labels of the
 ## cophylogenetic groups and calculate the LLRT statistic perm_number
 ## of times. A Monte-Carlo exact p-value is calculated based on the 
 ## proportion of permutation test statistics that are greater than  
 ## or equal to that based on the original observed data/group labels.
 ###################################################################
 
    library(dirmult)    # uses dirmult function in dirmult package
	
    n_group_1  <- nrow(group_1) # Size of group 1
    n_group_2  <- nrow(group_2) # Size of group 2
    group_0    <- rbind(group_1, group_2) # combining groups
    n          <- nrow(group_0) # Total number of cophylogenetic samples
	    
	
	# model observed data and calculate test statistic
    group_1_model <- dirmult(group_1, trace = FALSE) # Estimate dirichlet parameters for group 1
    group_2_model <- dirmult(group_2, trace = FALSE) # Estimate dirichlet parameters for group 2
    obs_teststat  <- group_1_model$loglik + group_2_model$loglik # observed test statistic
	
    do.one.perm <- function(i) { 
	    # do once: Permute group labels and calculate test statistic
	    perm_group_1_rows <- sample(1:n, n_group_1) # Pick out random observations to make up new 
		                                            # permuted group 1, same size as observed group 1
        perm_group_1      <- group_0[ perm_group_1_rows,] # create permuted group 1
        perm_group_2      <- group_0[-perm_group_1_rows,] # create permuted group 2
        perm_group_1_model<- dirmult(perm_group_1, trace = FALSE) # Estimate dirichlet parameters for permuted group 1
        perm_group_2_model<- dirmult(perm_group_2, trace = FALSE) # Estimate dirichlet paramteres for permuted group 2

        # Calculate LLRT statistic for permuted groups
        perm_teststat <- perm_group_1_model$loglik + perm_group_2_model$loglik 
        return(perm_teststat)
    }
    
    #Code to set up parallel computing - requires parallel package in R
    library(parallel)                # load parallel package
    cores <- detectCores()           # detect number of cores
    cl    <- makeCluster(cores)      # create cluster of cores
    clusterEvalQ(cl, library(DEEPR)) # load DEEPR library for all cores
    
	# In parallel, generate perm_number permutations; calculate and return test statistics
    perm_data      <- parLapply(cl, 1:perm_number, do.one.perm) # Compute vector listing LLRT statistics	  
    total_teststat <- c(perm_data, obs_teststat)           # Combine all test statistics
    p_value        <- mean(total_teststat >= obs_teststat) # Calculate p value

    # Cleaning up...
    stopCluster(cl) # Close cluster 

    # return p value and estimates of parameters for groups
    list(p_value = p_value, group_1_pi = group_1_model$pi, group_2_pi = group_2_model$pi, 
        group_1_alphas = group_1_model$gamma, group_2_alphas = group_2_model$gamma, 
        group_1_theta = group_1_model$theta, group_2_theta = group_2_model$theta) 
}
