# Copyright 2012 gcalmettes (gcalmettes@ucla.edu), BSD-style copyright and disclaimer apply
# adapt to multiple dimensional array : Shuping Li
#!/usr/local/bin/python

# needed libraries
import numpy as np
#import random as rd

###################################################
def bootstrap(data, nboot=30, replacement=True, dimension = True):
    """
        Generate n=nboot bootstrap samples of the data with/without replacement,
        and return a 3-d array numpy array of them.

        Input:    data (anything numpy can handle as a numpy array)
        Output:   ndarray numpy array of size (ntime x nboot x nlat x nlon )
    """
    if dimension == True:  # 3D array 
        nt,nlat,nlon = data.shape
    else:
        nt, nlat = data.shape     #2D array

    # Create a ndarray  of bootstrap samples indexes
    if replacement==True: # with replacement (note: 50x-ish faster than without)
        idx = np.random.randint(0,nt,(nt,nboot))
        return data[idx]
        
    elif replacement==False: # without replacement (an event may not occur more than once in a particular sample)
        #print (data[:,3,5])
        data1 = np.array([np.random.permutation(data) for x in np.arange(nboot)])
        data = data1.swapaxes(0,1)
        #print (data[:,1:3,3,5])
        return data


###################################################
def bootci(data, stat=np.mean, nboot=30, replacement=True, alpha=0.05,dimension=True):
    """
        Compute the (1-alpha) confidence interval of a statistic (i.e.: mean, median, etc)
        of the data using bootstrap resampling.
        
        Arguments:
            stat:        statistics we want the confidence interval for (must be a function)
            nboot:       number of bootstrap samples to generate
            replacement: resampling done with (True) or without (False) replacement
            alpha:       level of confidence interval
            method:      type of bootstrap we want to perform
            keepboot:    if True, return the nboot bootstrap statistics from which
                         the confidence intervals are extracted
        
        Methods available:
            - 'pi' = Percentile Interval
            - 'bca' = Bias-Corrected Accelerated Interval (available soon)
    """

    # apply bootstrap to data
    boot = bootstrap(data, nboot=nboot, replacement=replacement, dimension = dimension)

    # calculate the statistics percentile 
    print ('after resampling',boot.shape)
    boot_stat = stat(boot, axis=0)
    print (boot_stat.shape)
    upper_ci = np.percentile (boot_stat, 100-(alpha*100)/2., axis=0)
    lower_ci = np.percentile (boot_stat, (alpha*100)/2., axis=0)
    
    return lower_ci, upper_ci
        

###################################################
def bootci_diff(data1, data2, stat=np.mean, func=np.subtract, nboot=30,
                alpha=0.05, replacement=True):
    """
        Calculate the effect size of a treatment and compute the confidence
        interval of that effect size using bootstrap resampling method.

        Arguments:
            stat:        statistics we want the confidence interval for (must be a function)
            nboot:       number of bootstrap samples to generate
            replacement: resampling done with (True) or without (False) replacement
            alpha:       level of confidence interval
            keepboot:    if True, return the nboot bootstrap statistics from which
                         the confidence intervals are extracted
    """


    # bootstrap samples
    boot_d1 = bootstrap(data1, nboot=nboot, replacement=replacement)
    boot_d2 = bootstrap(data2, nboot=nboot, replacement=replacement)
    
    print ('after resampling',boot_d1.shape, boot_d2.shape)
    
    # stat of bootstrap samples
    boot_d1stat = stat(boot_d1, axis=0)
    boot_d2stat = stat(boot_d2, axis=0)
    
    print ('calculating average', boot_d1stat.shape, boot_d2stat.shape)
  
    
    # ci
    bootdiff = func(boot_d1stat, boot_d2stat)
    #print (bootdiff.shape)
    #sorted_bootdiff = np.sort(bootdiff, axis=0)
    
    upper_ci = np.percentile (bootdiff, 100-(alpha*100)/2., axis=0)
    lower_ci = np.percentile (bootdiff, (alpha*100)/2., axis=0)
    #ci_diff = (sorted_bootdiff[np.round(nboot*alpha/2).astype(int),:,:], 
    #         sorted_bootdiff[np.round(nboot*(1-alpha/2)).astype(int),:,:])
    print ('ci difference')
    print(lower_ci.shape, upper_ci.shape)

   
    return lower_ci, upper_ci


###################################################
def bootpv(data1, data2, stat=np.mean, func=np.subtract, nboot=10000, 
           replacement=True, printout=True, keepboot=False):
    """
        Assuming that there is no significant difference in the means of 
        the two samples, return the one-tailed pvalue to see the probability of
        getting a difference greater than or equal to the observed difference
        in the means by chance alone.
        Permutation (replacement = False) or bootstrapping (default, replacement = True)
        is used to get a distribution to compare to the observed statistic.

        Input:  data1 and data2 (anything numpy can handle as a numpy array)
        Output: one-tailed pvalue

        Arguments:
            stat:        statistics we want the pvalue for (mean, median, etc..)
            func:        the function we want to use to compare data (default algebric diff)
            nboot:       number of bootstrap samples to generate
            replacement: resampling done with (True) or without (False) replacement
            printout:    if True (default), print the result of the analysis
            keepboot:    if True, return the nboot bootstrap statistics from which
                         the confidence intervals are extracted        
    """

    # Ensure that our data are 1D arrays
    data1, data2 = np.ravel(data1), np.ravel(data2)

    # Mesure the difference of the statistic between the two groups
    diff_groups = func(stat(data1), stat(data2))

    # Size of the two samples
    samplesize1, samplesize2 = data1.size, data2.size

    # Put all the values from the two groups in one array
    pool = np.hstack((data1, data2))

    if replacement == False: # resampling without replacement (permutation)
        # Create a matrix of n_samples shuffled indexes of the pool
        pool_idxtable = np.array([np.random.permutation(pool.size) for n in range(nboot)])

        # Create n_samples sets of data from size of the two groups 
        sample1 = pool[pool_idxtable[:, :samplesize1]]
        sample2 = pool[pool_idxtable[:, samplesize1:]]

    elif replacement == True: # bootstrap
        bootpool = bootstrap(pool, nboot=nboot, replacement=True)
        sample1 = bootpool[:, :samplesize1]
        sample2 = bootpool[:, samplesize1:]

    param_sample1 = np.array([stat(x) for x in sample1])
    param_sample2 = np.array([stat(x) for x in sample2])
    
    param_func = func(param_sample1, param_sample2)
    
    # Warning: the criteria set is > or < zero (ok for difference)
    # if the func argument is not np.subtract (algebric difference),
    # this has to be changed (ex: if np.divide is used, criteria will be
    # > or < 1, since a ratio cannot be negative)
    if diff_groups < 0:
        pvalue = np.sum(param_func <= diff_groups)/float(nboot)
    elif diff_groups > 0:
        pvalue = np.sum(param_func >= diff_groups)/float(nboot)
    
    if printout == True:
        sign = "greater than or equal to"
        ndiff = np.sum(param_func >= diff_groups)
        if diff_groups < 0:
            sign = "less than or equal to"
            ndiff = np.sum(param_func <= diff_groups)
        print("Observed difference of two means: {0:.2f}".format(diff_groups))
        print("{0:d} out of {1:d} experiments had a difference of two means {2:s} {3:.2f}"\
            .format(ndiff, nboot, sign, diff_groups))
        print("The chance of getting a difference of two means {0:s} {1:.2f} is {2:.4f}"\
              .format(sign, diff_groups, pvalue))
       
    if keepboot == True:
        return pvalue, param_func
    else:
        return pvalue
