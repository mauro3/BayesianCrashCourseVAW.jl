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

#############################

# Metropolis-Hastings sampler
function metropolis_hastings(target_pdf,   # Target PDF (unnormalized)
                             proposal,   # Proposal distribution (e.g., Normal step)
                             initial_state, # Starting point
                             n_samples   # Number of samples to generate
                            )
    current_state = initial_state
    samples = Vector{Vector{Float64}}(undef, n_samples)
    for i in 1:n_samples
        # Propose a new state
        proposed_state = proposal(current_state)
        
        # Compute acceptance probability
        acceptance_ratio = target_pdf(proposed_state) - target_pdf(current_state)
        
        # Decide whether to accept or reject
        u = log(rand())
        if u < acceptance_ratio
            current_state = proposed_state
        end
        
        # Store the current state
        samples[i] = current_state
    end
    
    return samples
end

# Helper function: Metropolis-Hastings with normal proposal
function mh_with_normal_proposal(target_pdf, initial_state, n_samples, proposal_sd)
    proposal = x -> x .+ rand(Normal(0, proposal_sd), length(x))
    return metropolis_hastings(target_pdf, proposal, initial_state, n_samples)
end

## Use MH to sample from posterior
# --> adjust the example_2d from mh-example.jl file:


