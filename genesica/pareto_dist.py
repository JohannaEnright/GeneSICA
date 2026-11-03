#!/usr/bin/env python3

#from scipy.stats import pareto
from scipy.stats import genpareto, norm, goodness_of_fit
import numpy as np
import warnings

#https://pmc.ncbi.nlm.nih.gov/articles/PMC2687965/


#in the paper, the GPD is fitted to the tail of an F-distribution?
# approximation of the tail is [1-N_exc/N, 1] interval
# test statistic must also be in the tail distribution of the permutation values

### Need to do maximum likelihood to estimate optimal parameters
### Estimate Shape and scale given the exceedences (permutation values - threshold)
### only applies for permutation values that are beyond the threshold
### alternatives to maximum likelihood include 'method of moments' and 'probability weighted moments'
### Method must be: 'MLE' of 'MM'
### x values are esceedences: perm_values-t
def optimize_parameters(y_list, t, method="MLE"):

    c, loc, scale = genpareto.fit(y_list, floc=t, method=method) #c is the shape parameter

    #print(f"Optimized params: c={round(c, 4)}, loc={round(loc, 4)}, scale={round(scale, 4)}")
    return (c, loc, scale)  #shape, shift and scale



## determine exceedence threshold
## order permutation values from smallest to largest...
#pluck the smallest one that is just above the exceedence threshold and the largest one that is just below the threshodef exceedence_threshold(z_list, y_list):
def exceedence_threshold_calc(y_list_rev, Nexc):
    
    #y_list is sorted from high to low
    y_Nexc = y_list_rev[Nexc-1] # adjust for 0 indexing
    y_Nexc_1 = y_list_rev[Nexc]

    t = (y_Nexc + y_Nexc_1)/2
    #print(f"t={round(t, 4)}")
    return t


##k is the shape param, a is the scale param, z are the exceedences
## cumulative distribution function
#def gpd_cdf(z, k, a):
#
#    if k != 0:
#        y = 1 - (1-k*z/a)**(1/k)
#    else:
#        y = 1-np.exp(-1*z/a)
#    return y
#
#def gpd_pdf(z, k, a):
#
#    if k != 0:
#       y = a**(-1)*(1-k*z/a)**(1/(k-1))
#    else:
#       y = a**(-1)*np.exp(-1*z/a)
#    return y


## what is the threshold t to determine the exceedences?
## x is the permutation values, t is the threshold
## return the list of exceedences
### paper states using N_exc = 250 (where N_exc = len(z_list))
#def compute_exceedences(y_list, t):
#
#    z_list = []
#    for y in y_list:
#        if y > t:  #passing a list of ys that should exceed t anyways
#            z = y-t
#            z_list.append(z)
#        else:
#            raise TypeError(f"Permutation value y: {y}, does not exceed threshold t: {t}\nBut it should.\n")
#    return z_list

# N_exc of permutation statistics that exceed threshold
# N is number of permutation statistics in total
#x_0 is the original statistic
# t is the exceedence threshold (based on where the extreme tail end begins I beleive?)
#def tail_approximation(N_exc, N, x_0, t):
#
#    p_gpd = (1 - F(compute_exceedences(x_0, t))*N_exc/N)


# Employ goodness of fit tests to determine if extreme tail permutation values follow a GPD.
# start at N_exc = 250, iteratively decreased by 10, until a GPD goodness of fit is reached (p > 0.05)
def exceedence_goodness_of_fit(y_list, c, t, scale, alpha=0.05):

    ## insert goodness of fit function here, to compute p value
    #print(f"Performing Goodness of Fit Test: c={c}, scale={scale}")
    result = goodness_of_fit(
        dist=genpareto,
        data=y_list,
        known_params={
            "loc":t
            },
        fit_params={
            "c":c,
            "scale":scale
            },
        #guessed_params=None,
        statistic='ad',        #{“ad”, “ks”, “cvm”, “filliben”} #ad='anderson-darling'
        n_mc_samples=1000,   #10 000 monte carlo samples drawn for null dist
        random_state=10
        )
    
    if result.pvalue <= alpha:
        reject_null = True
    else:
        reject_null = False
    return reject_null

# input is all permutation values
# determine the permutation values that will be used to compute the exceedences
def determine_Y_exc(y_list, init_chop=500, interval=25, secondary_int=10, low_bound=50):

    #print("Determining Y values for Exceedences")
    y_list = np.sort(np.array(y_list), descending=True)  #sort in descending order (needed to t calc)

    #transform permutations and test  statistic to the power of 3 (heavier tails (c<0) converge faster with fewer perms
    y_list_new = np.sort(np.array(y_list), descending=True) #make a copy
    reject_null = True
    n_iter = 0

    #1. Determine Exceedences
    while reject_null:
        L = len(y_list_new)
        if n_iter == 0:
            if init_chop < low_bound:
                break
            elif L > init_chop:
                y_list_new = y_list_new[:init_chop]  #chose the most extreme 250 values to start
        else:
            if L <= low_bound:
                break
            if L > 100:
                new_L = max(low_bound, L - interval)
            elif L <= 100:
                new_L = max(low_bound, L - secondary_int) #decrease slower once you get below 100
            y_list_new = y_list_new[:new_L]  #decrease by interval until it follows GPD distribution

        #2. Calculate exceedence threshold: t and exceedence values
        t = exceedence_threshold_calc(y_list, len(y_list_new)) #
        #z_list = y_list_new-t
        #3. Optimize parameters 
        c, loc, scale = optimize_parameters(y_list_new, t)
        #4. Calculate Goodness of Fit
        reject_null = exceedence_goodness_of_fit(y_list_new, c, t, scale)

        n_iter += 1

    ## If it doesn't follow a gpd.
    if reject_null:
        return False  #False as in it does not follow at gpd

    ## If it does follow a gpd
    fit_info = {
            "t": t,
            "c": c,
            "scale": scale,
            "y_exc": y_list_new
            }
    print(f"Final GPD params: t (loc)={round(t,4)}, shape (c)={round(c,4)}, scale={round(scale,4)}, Nexc={len(y_list_new)}")
    return fit_info



def check_z_range(c, scale, y_list, t):

    if c <= 0:
        #print(f"c (shape param) is less than 0: {c}, this should provide more stable low p-val estimates")
        in_range=True
    else:
        upper_bound = scale/c
        n_beyond_bound = sum((y_list-t) > upper_bound)
        if n_beyond_bound > 0:
            warnings.warn(f"Less reliable low p-val estimates. c (shape param) is greater than 0: {c}.\n{n_beyond_bound} exceedences exceed the upper limit: {upper_bound}.")
            in_range=False
        else:
            #print(f"c (shape param) is greater than 0: {c}, however all exceedence values are within the upper limit.")
            in_range=True

    return in_range
            


## Let M be the number of permutation values exceeding the test statistic
## Helps determine if the test statistic is in the tail end of the distribution or not
def calc_M(test_stat, perm_values):

    perm_values = np.array(perm_values)
    M = np.sum(perm_values > test_stat)
    return M

## p-value is the CDF evaluated at the value of the test statistic
def approximate_p_val(M, ts_exc=False, c=False, t=False, scale=False, N_perm=False, N_exc=False):
    
    ## If more than 10 exceed the test statistic, it is likely not in the tail of the distribution?
    ## Also if M < 10 but it does not follow a GPD
    if (M >= 10) or (not c):
        #print(f"Falling back to permutation approximation of pval.  M={M}")
        y = (M + 1)/(N_perm + 1)  #add 1 pseudocount to avoid p=0
        #se = np.sqrt(p_prime * (1 - p_prime) / N_perm)
        #P'ecdf approximation (can follow normal approximation according to the central limit theorem)
        # mu = p'ecdf  
        # sigma = p'ecdf*(1-p'ecdf)/N
        #conf = norm.pdf(ts_exc, mu=M, sigma=M*(1-M)/N_perm)
    else:
        #x, c (shape parameter), loc=0 (loc shifts the distribution), scale=1 
        #y_alt1 = (N_exc/N_perm)*(1-genpareto.cdf(ts_exc, c=c, scale=scale, loc=t))  #multiply by N_exc/N_perm
        y = (N_exc/N_perm) * genpareto.sf(ts_exc, c=c, scale=scale, loc=t)  #sf = 1-cdf (but more stable)
    return y


def main():

    ## Read in intermediate File
    with open("pareto_sample_1000.txt", "r") as f:
        samp = f.readlines()
        samp = np.array([float(i.strip()) for i in samp]) #permuted values (numpy array)
    # test statistic
    test_stat = 30.0
    transformation_exp = 0.01  #can adjust trasformation (>0) to get a negative k value

    #t = exceedence_threshold_calc(samp, len(y_exc_list))
    #z_list = compute_exceedences(y_exc_list, t)
    M = calc_M(test_stat, samp)
    #print("M: ", M)
    #c, loc, scale = optimize_parameters(z_list)
    if M < 10:
        fit_info = determine_Y_exc(samp**transformation_exp)  #list of y values needed to calculate exceedences

        test_stat = test_stat**transformation_exp  #can use any strictly increasing function to get a heavier tail (can use 1 if no transformation is desired
        c = fit_info["c"], 
        t = fit_info["t"], 
        scale = fit_info["scale"]
        in_range = check_z_range(c, scale, samp, t) #if in_range == False, could maybe adjust the transformation_exp
    
        pval = approximate_p_val(
               M, 
               test_stat, 
               c, 
               t, 
              scale, 
              N_perm=len(samp), 
              N_exc=len(fit_info["y_exc"])
              )
        #print("pval:  ", pval)
    else:
        pass
        #print("M is greater than 10")

#if __name__ == "__main__":
#    main()

    

#Your fitted params were c≈0.49, loc≈6.07, scale≈3.23. Since c > 0, the GPD has finite support: t ≤ x ≤ t + scale/c ≈ 6.07 + 6.59 ≈ 12.66. Your test statistic is 30.0 — far beyond that upper bound. Once the code is correct, sf(30, ...) will likely evaluate to exactly 0, giving you a p-value of zero. This isn't a bug in your code — it's the exact instability the paper describes in §3.1.1/§3.2.1 for light-tailed distributions (k>0) under large extrapolation. If that happens, the paper's suggested fix is an order-preserving power transform (e.g., cube both the permutation values and the test statistic) to push the shape parameter negative before refitting.#
