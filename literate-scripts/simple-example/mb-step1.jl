#src # This is needed to make this run as normal Julia file:
using Markdown #src

using KissMCMC, Random, Distributions, CairoMakie

# Define the forward model
"""
    glacier_mb(P, T, m)

Calculate the mass balance of a glacier.

Arguments
- `P::Float64`: Precipitation amount.
- `T::Float64`: Temperature.
- `m::Float64`: Melt factor.

Returns
- `Float64`: The mass balance of the glacier calculated as `P - T * m`.
"""
function glacier_mb(P, T, m)
    #hint ...
    return P - T * m #sol
end

# Define prior
"""
   prior(theta)

   Log-prior.  Assumes a normal distribution of P and T around their
   means (1000 and 2) with standard deviations of 100 and 0.5, respectively.
"""
function prior(theta)
    #hint ...
    P, T = theta #sol
    return - (P - 1000)^2 / (2*100^2) - (T - 2)^2 / (2*0.5^2) #sol
end

# Define likelihood
"""
    likelihood(theta::Tuple, data::Tuple) -> Float64

Calculate the log-likelihood of the observed data given the parameters. The
function computes the likelihood based on the difference between the modeled
glacier mass balance (`mb`) and the observed mass balance (`mb_d`), normalized
by the variance (`mb_sigma^2`).

Arguments
- `theta::Tuple`: A tuple containing the parameters (P, T).
- `data::Tuple`: A tuple containing the observed data ().

Returns
- `Float64`: The likelihood value.
"""
function likelihood(theta, data)
    #hint ...
    P, T = theta #sol
    mb_d, mb_sigma = data #sol
    mb = glacier_mb(P, T, 200) #sol
    return - (mb - mb_d)^2 / (2*mb_sigma^2) #sol
end

# define posterior
"""
    posterior(theta, data)

Calculate the posterior probability of the parameter `theta` given the `data`.

Arguments
- `theta`: The parameter for which the posterior probability is calculated.
- `data`: The observed data used to update the prior probability.

Returns
- The posterior probability of `theta` given the `data`.
"""
function posterior(theta, data)
    #hint ...
    prior(theta) + likelihood(theta, data) #sol
end

posterior([879, 1.5], [400, 56])
