using KissMCMC
using Random
using Distributions
using GLMakie

# Define the forward model
function glacier_mb(P, T, m)
    return P - T * m
end

# Define prior
function prior(theta)
    P, T = theta
    - (P - 1000)^2 / (2*100^2) - (T - 2)^2 / (2*0.5^2)
end

# Define likelihood
function likelihood(theta, data)
    P, T = theta
    mb_d, mb_sigma = data
    mb = glacier_mb(P, T, 200)
    return - (mb - mb_d)^2 / (2*mb_sigma^2)    
end

# define posterior
function posterior(theta, data)
    prior(theta) + likelihood(theta, data)
end

posterior([879, 1.5], [400, 56])
