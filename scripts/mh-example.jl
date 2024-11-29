using Random
using Distributions
using GLMakie

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
        acceptance_ratio = target_pdf(proposed_state) / target_pdf(current_state)
        alpha = min(1.0, acceptance_ratio)
        
        # Decide whether to accept or reject
        u = rand()
        if u < alpha
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

# Example 1: 1D Target PDF (Gaussian Mixture Model)
function example_1d()
    # Define the target PDF (a mixture of Gaussians, unnormalized)
    target_pdf(x) = exp(-0.5 * (x[1] - 2)^2) + exp(-0.5 * (x[1] + 2)^2)
    
    initial_state = [0.0]  # Starting point
    n_samples = 10_000     # Number of samples
    proposal_sd = 1.0      # Proposal standard deviation

    # Perform the sampling
    samples = mh_with_normal_proposal(target_pdf, initial_state, n_samples, proposal_sd)
    
    # Extract 1D samples for plotting
    samples_1d = [s[1] for s in samples]

    # Normalize the PDF for display
    xs = -5.0:0.1:5.0
    pdf_vals = [target_pdf([x]) for x in xs]
    pdf_vals /= sum(pdf_vals) * (xs[2] - xs[1]) # Normalize against the histogram

    # Plot results
    fig = Figure()
    ax = Axis(fig[1, 1], title = "1D Metropolis-Hastings", xlabel = "x", ylabel = "Density")
    hist!(ax, samples_1d; bins = 50, normalization=:pdf, color = :blue, label = "Samples")
    lines!(ax, xs, pdf_vals, color = :orange, linewidth = 2, label = "Target PDF")
    axislegend(ax)
    return fig
end

# Example 2: 2D Target PDF (Correlated Gaussian)
function example_2d()
    # Define the target PDF (unnormalized correlated Gaussian)
    target_pdf(x) = exp(-0.5 * (x[1]^2 + x[2]^2 + 2.0 * x[1] * x[2]))
    
    initial_state = [0.0, 0.0]  # Starting point
    n_samples = 10_000          # Number of samples
    proposal_sd = 0.5           # Proposal standard deviation

    # Perform the sampling
    samples = mh_with_normal_proposal(target_pdf, initial_state, n_samples, proposal_sd)

    # Extract 2D samples for plotting
    samples_x = [s[1] for s in samples]
    samples_y = [s[2] for s in samples]

    # Plot the 2D results
    fig = Figure()
    ax = Axis(fig[1, 1], title = "2D Metropolis-Hastings", xlabel = "x1", ylabel = "x2")
    scatter!(ax, samples_x, samples_y, markersize = 2, color = :blue, label = "Samples", alpha = 0.3)
    return fig
end

# Run the examples (uncomment to execute one at a time)
fig_1d = example_1d() # 1D example
display(fig_1d)

fig_2d = example_2d() # 2D example
display(fig_2d)

