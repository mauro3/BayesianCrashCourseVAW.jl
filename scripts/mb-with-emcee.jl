
using KissMCMC
using Random
using Distributions
using GLMakie

# Define the forward model
function glacier_mb(P, T, m)
    return P - T * m
end

# Define the log-posterior for the model
function bayesian_glacier_model_logpdf(θ)
    # θ[1]: Precipitation (P), θ[2]: Temperature (T)
    P, T = θ[1], θ[2]
    
    # Fixed "data" and model parameters
    mb_measured = 400.0           # Observed mass balance
    mb_uncertainty = 100.0        # Observation uncertainty
    m = 200.0                     # Fixed melt factor

    # Priors (log form of Normal distributions)
    logprior_P = - (P-1000)^2/(2 * 100^2)
    logprior_T = - (T-2)^2/(2 * 0.5^2)

    # Likelihood
    mb = glacier_mb(P, T, m)
    loglikelihood = - (mb - mb_measured)^2/(2 * mb_uncertainty^2)

    # Return log-posterior (log-prior + log-likelihood)
    return logprior_P + logprior_T + loglikelihood
end

# Create initial walkers
function make_theta0s(theta0, ballradius, nwalkers)
    n_dim = length(theta0)
    return [theta0 .+ ballradius .* randn(n_dim) for _ in 1:nwalkers]
end

# Define starting conditions
n_dim = 2                           # Number of parameters (P and T)
n_walkers = 10                      # Number of walkers
initial_theta = [1000.0, 2.0]       # Initial guess: [P, T]
theta0s = make_theta0s(initial_theta, 1.0, n_walkers)  # Generate initial walker positions

# Run the affine invariant sampler
samples, accept_ratio, logdensities = emcee(
    bayesian_glacier_model_logpdf,
    theta0s;
    niter = 20_000,   # Total number of iterations
    nburnin = 5_000,  # Burn-in period
    a_scale = 2.0,    # Stretch move scale
    use_progress_meter = true
)

# Extract the samples and combine all walker chains into one
flattened_samples = squash_walkers(samples, accept_ratio, logdensities)[1]  # Combine walker samples into a single chain

# Extract the posterior samples for P and T
P_samples = [s[1] for s in flattened_samples]
T_samples = [s[2] for s in flattened_samples]

# Generate prior samples for comparison
P_samples_prior = rand(Normal(1000, 100), 10_000)
T_samples_prior = rand(Normal(2, 0.5), 10_000)

# Plot the results: Prior vs Posterior distributions
fig = Figure()

# Plot Precipitation (P)
ax1 = Axis(fig[1, 1], title = "Posterior Distribution: Precipitation (P)", xlabel = "P (mm)", ylabel = "Density")
hist!(ax1, P_samples_prior; bins = 50, normalization = :pdf, color = :green, strokewidth = 0.0, label = "Prior")
hist!(ax1, P_samples; bins = 50, normalization = :pdf, color = :blue, strokewidth = 0.0, label = "Posterior")
axislegend(ax1)

# Plot Temperature (T)
ax2 = Axis(fig[1, 2], title = "Posterior Distribution: Temperature (T)", xlabel = "T (°C)", ylabel = "Density")
hist!(ax2, T_samples_prior; bins = 50, normalization = :pdf, color = :green, strokewidth = 0.0, label = "Prior")
hist!(ax2, T_samples; bins = 50, normalization = :pdf, color = :blue, strokewidth = 0.0, label = "Posterior")
axislegend(ax2)

display(fig)
