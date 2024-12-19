using Distributions, Random, Statistics, GLMakie # for stats and plotting

# Define the forward model
glacier_mb(P, T, m) = P - T*m

# Define input distributions
P_dist = Normal(1000, 100)  # Snow precipitation distribution (mean=1000mm, std=100mm)
T_dist = Normal(2, 0.5)     # Annual temperature distribution (mean=2°C, std=0.5°C)
# Fixed melt factor
m = 200.0

# Monte Carlo sampling
num_samples = 10000
P_samples = rand(P_dist, num_samples)  # Sample precipitation
T_samples = rand(T_dist, num_samples)  # Sample temperature

# Compute mb for each sample
mb_samples_MC = glacier_mb.(P_samples, T_samples, m)

# Analyze results
mean_mb, std_mb  = mean(mb_samples_MC), std(mb_samples_MC)
# Print results
println("\nMean annual mass balance: $(mean_mb) mm")
println("Standard deviation of MB: $(std_mb) mm")

# Plot results
fig = Figure()
ax = Axis(fig[1, 1], title="Monte Carlo Mass Balance Distribution", xlabel="Mass balance (mm)", ylabel="Frequency")
hist!(ax, mb_samples_MC, bins=50, color=:blue, strokewidth=0.0,)
fig
