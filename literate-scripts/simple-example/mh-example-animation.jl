
using Random
using Distributions
using GLMakie

# Metropolis-Hastings sampler with animation support
function metropolis_hastings_animate(
    target_pdf::Function,  # Target PDF (unnormalized)
    proposal::Function,    # Proposal distribution
    initial_state::Vector, # Starting point
    n_steps::Int,          # Number of steps in the walk
    update_callback::Function # Callback to update the visualization
)
    current_state = initial_state
    samples = Vector{Vector{Float64}}(undef, n_steps)
    
    for i in 1:n_steps
        # Generate new proposal
        proposed_state = proposal(current_state)

        # Compute acceptance probability
        acceptance_ratio = target_pdf(proposed_state) / target_pdf(current_state)
        alpha = min(1.0, acceptance_ratio)
        
        # Decide whether to accept the proposal
        u = rand()
        if u < alpha
            current_state = proposed_state
        end

        # Store current state
        samples[i] = current_state
        
        # Call the update function to visualize the current state
        update_callback(current_state, i, samples[1:i])
    end
    
    return samples
end

# Define 2D example with animation
function example_2d_animated()
    # Define the target unnormalized PDF (2D correlated Gaussian)
    target_pdf(x) = exp(-0.5 * (x[1]^2 + x[2]^2 + 2.0 * x[1] * x[2]))
    
    # Setup parameters
    initial_state = [0.0, 0.0]  # Starting point
    total_steps = 500           # Total number of steps to animate
    proposal_sd = 0.5           # Proposal standard deviation
    proposal = x -> x .+ rand(Normal(0, proposal_sd), length(x))

    # Prepare the figure and axis
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], title = "2D Metropolis-Hastings Random Walk", xlabel = "x1", ylabel = "x2")
    # scatter!(ax, [], [], markersize = 4, color = :blue, label = "Samples", alpha = 0.5) # Initialize empty scatter
    current_point = scatter!(ax, [0.0], [0.0], color = :red, markersize = 6, label = "Current State") # Current step
    axislegend(ax)

    # Helper to update the visualization
    function update_callback(current_state, step, samples)
        x_samples = [s[1] for s in samples]
        y_samples = [s[2] for s in samples]

        # Update the scatter plot for all samples
        scatter!(ax, x_samples, y_samples, markersize = 4, color = :blue, alpha = 0.5; overdraw = true)
        
        # Update current point dynamically
        current_point[1] = (current_state[1], current_state[2]) # Update the red dot position
        sleep(0.01) # Slow down to make the process viewable
    end

    # Run the Metropolis-Hastings sampling with animation support
    metropolis_hastings_animate(target_pdf, proposal, initial_state, total_steps, update_callback)

    # Display the final figure
    return fig
end
# Run the animated 2D walk
fig = example_2d_animated()
display(fig)

