using Turing, GLMakie # using a MCMC package, and plotting

# Define the forward model
glacier_mb(P, T, m) = P - T*m

@model function bayesian_glacier_model(mb_measured, mb_uncertainty)
    # Priors
    P ~ Normal(1000, 100)  # Precipitation
    T ~ Normal(2, 0.5)     # Temperature
    # Fixed melt factor
    m = 200.0
    
    # Likelihood
    mb = glacier_mb(P, T, m)
    mb_measured ~ Normal(mb, mb_uncertainty)
end

# Define model and draw samples from the posterior
model = bayesian_glacier_model(400.0, 100.0)
chain = sample(model, NUTS(0.65), 1000, progress=false);
summarystats(chain)

# Plot prior vs posterior distributions
# Extract samples for precipitation (P) and temperature (T)
P_samples, T_samples = chain[:P][:], chain[:T][:] # Samples of P and T
# make samples from the prior to compare
P_samples_prior,T_samples_prior = rand(Normal(1000, 100), 10000), rand(Normal(2, 0.5), 10000)

# Plot precipitation posterior
fig = Figure(); ax1 = Axis(fig[1, 1], title="Posterior Distribution: Precipitation (P)", xlabel="P (mm)", ylabel="Density")
hist!(ax1, P_samples_prior, bins=50, normalization=:pdf, color=:green, strokewidth=0.0)
hist!(ax1, P_samples, bins=50, normalization=:pdf, color=:blue, strokewidth=0.0)
# Plot temperature posterior
ax2 = Axis(fig[1, 2], title="Posterior Distribution: Temperature (T)", xlabel="T (°C)", ylabel="Density")
hist!(ax2, T_samples_prior, bins=50, normalization=:pdf, color=:green, strokewidth=0.0)
hist!(ax2, T_samples, bins=50, normalization=:pdf, color=:blue, strokewidth=0.0)
fig

#################

using Colors
include("forward-propagation.jl")
mb_samples = glacier_mb.(P_samples, T_samples, m)[:]
fig = Figure(); ax = Axis(fig[1, 1], title="Mass Balance Posterior Distribution",
          xlabel="Mass Balance (mm)", ylabel="Density")
hist!(ax, mb_samples_MC, bins=50, normalization=:pdf, color=:green, strokewidth=0.0)
hist!(ax, mb_samples, bins=50, normalization=:pdf, color=:blue, strokewidth=0.0); fig

################

# Define labels and pairwise combinations
variables = ["P (Precipitation)", "T (Temperature)", "mb (Mass Balance)"]
samples = (P_samples, T_samples, mb_samples)

# Compute correlation matrix
correlations = [cor(samples[i], samples[j]) for i in 1:3, j in 1:3]  # Pairwise correlations

# Create a grid figure for a 3x3 matrix
fig = Figure()

# Loop over each plot position in the 3x3 grid
for i in 1:3
    for j in 1:3
        local ax = fig[i, j] = Axis(fig[i, j], 
                        xlabel=i == 3 ? variables[j] : "",  # Add xlabel to the bottom row
                        ylabel=j == 1 ? variables[i] : "",  # Add ylabel to the left column
                        xticklabelsize=8,
                        yticklabelsize=8)

        if i == j
            # Diagonal: Marginal distribution (histogram)
            hist!(ax, samples[i], bins=50, normalization=:pdf, color=:blue, strokewidth=0.0)
        elseif i > j
            # Lower triangle: Pairwise joint distribution (scatter plot)
            scatter!(ax, samples[j], samples[i], color=:blue, markersize=2)
        else
            # Upper triangle: Correlation coefficient
            corr_val = round(correlations[i, j], digits=2)
            text!(ax, 0.5, 0.5, text=string(corr_val), fontsize=20, color=:black)

           # Add colored background
           if corr_val > 0  # Positive correlation: Green
              poly!(ax, [0, 1, 1, 0], [0, 0, 1, 1], color=RGBA(0.8, 1.0, 0.8, abs(corr_val)))
           else  # Negative correlation: Red
              poly!(ax, [0, 1, 1, 0], [0, 0, 1, 1], color=RGBA(1.0, 0.8, 0.8, abs(corr_val)))
           end

            hidedecorations!(ax)  # Hide ticks and labels in the upper triangle
        end
    end
end

fig
