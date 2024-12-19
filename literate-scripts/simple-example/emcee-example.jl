using Random
using LinearAlgebra
using GLMakie

# Affine Invariant Ensemble Sampler
function affine_invariant_sampler(
    target_pdf::Function,   # Target PDF (unnormalized)
    walkers::Matrix{Float64}, # Initial positions of walkers, size (n_dim, n_walkers)
    n_steps::Int,           # Total number of steps to run
    a::Float64 = 2.0        # Stretch move scaling parameter (usually 2.0)
)
    n_dim, n_walkers = size(walkers)
    samples = Vector{Matrix{Float64}}(undef, n_steps) # Store all walker positions at each step
    
    for step in 1:n_steps
        for i in 1:n_walkers
            # Choose a random companion walker
            companion_idx = rand(1:n_walkers)
            while companion_idx == i
                companion_idx = rand(1:n_walkers)
            end
            companion = walkers[:, companion_idx]

            # Propose a new position using the stretch move
            z = ((a - 1.0) * rand() + 1.0)^2 / a
            proposal = companion + z * (walkers[:, i] - companion)

            # Acceptance probability
            acceptance_ratio = (z^(n_dim - 1)) *
                               (target_pdf(proposal) / target_pdf(walkers[:, i]))
            α = min(1.0, acceptance_ratio)

            # Accept or reject the proposal
            if rand() < α
                walkers[:, i] = proposal
            end
        end

        # Store current positions of all walkers
        samples[step] = copy(walkers)
    end
    
    return samples
end

# Helper function to flatten walker samples into a list
function flatten_samples(samples)
    n_dim, n_walkers = size(samples[1])
    flattened = Vector{Vector{Float64}}()
    for step in samples
        for walker in axes(step, 2)
            push!(flattened, step[:, walker])
        end
    end
    return flattened
end

function example_1d_affine()
    # 1D target PDF: Mixture of two Gaussians
    target_pdf(x) = exp(-0.5 * (x[1] - 2)^2) + exp(-0.5 * (x[1] + 2)^2)
    
    # Initialize walkers
    n_walkers = 10
    walkers = [randn() for _ in 1:n_walkers] |> x -> hcat(x...)
    walkers = reshape(walkers, :, n_walkers)  # Reshape into (n_dim, n_walkers)

    # Perform sampling
    n_steps = 500
    samples = affine_invariant_sampler(target_pdf, walkers, n_steps)
    
    # Extract and plot samples
    final_samples = flatten_samples(samples)
    flattened_1d = [s[1] for s in final_samples]

    # Plot
    fig = Figure()
    ax = Axis(fig[1, 1], title = "1D Affine Invariant Sampler", xlabel = "x", ylabel = "Density")
    hist!(ax, flattened_1d; bins = 50, normalization=:pdf, color = :blue, label = "Samples")
    xs = -5.0:0.1:5.0
    pdf_vals = [target_pdf([x]) for x in xs]
    pdf_vals /= sum(pdf_vals) * (xs[2] - xs[1]) # Normalize to match histogram
    lines!(ax, xs, pdf_vals, color = :orange, linewidth = 2, label = "Target PDF")
    axislegend(ax)
    return fig
end

fig_1d_affine = example_1d_affine()
display(fig_1d_affine)
function example_2d_affine()
    # Define 2D target PDF: Correlated Gaussian
    target_pdf(x) = exp(-0.5 * (x[1]^2 + x[2]^2 + 2.0 * x[1] * x[2]))
    
    # Initialize walkers (matrix where each column is a walker)
    n_dim = 2
    n_walkers = 10
    walkers = randn(n_dim, n_walkers)  # Shape (2, 10)

    # Perform sampling using the Affine Invariant Sampler
    n_steps = 500
    samples = affine_invariant_sampler(target_pdf, walkers, n_steps)
    
    # Extract samples and flatten
    final_samples = flatten_samples(samples)
    samples_x = [s[1] for s in final_samples]
    samples_y = [s[2] for s in final_samples]

    # Scatter plot results
    fig = Figure()
    ax = Axis(fig[1, 1], title = "2D Affine Invariant Sampler", xlabel = "x1", ylabel = "x2")
    scatter!(ax, samples_x, samples_y, markersize = 2, color = :blue, alpha = 0.5, label = "Samples")
    axislegend(ax)
    return fig
end

fig_2d_affine = example_2d_affine()
