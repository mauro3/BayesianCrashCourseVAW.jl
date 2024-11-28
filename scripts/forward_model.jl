# TODO
# 1. Implement likelyhood
# 2. Filter the data to remove borders of the glacier as per Huss et al 2010

using DelimitedFiles, GLMakie

include("utils.jl")

glacier_dh(hr, a, b, c, g) = (hr + a)^g + b * (hr + a) + c

# mm / a
glacier_mb(P, T, m) = P - max(T, 0) * m

glacier_T(z, zref, Tref, dT_dz) = Tref + (z - zref) * dT_dz
glacier_P(z, zref, Pref, dP_dz) = Pref + (z - zref) * dP_dz

function likelihood()
    # TODO
end

function forward_model(H0, B, A, Pref, Tref, zref, dP_dz, dT_dz)
    # known physics
    rhow = 1000.0 # water density (kg m⁻³)
    rhoi = 920.0  # ice density (kg m⁻³)

    # fixed surface mass balance parameters
    m = 1000.0 # melt rate factor (mm a⁻¹ C⁻¹)

    # fixed dh-model parameters for large valley glaciers
    a = -0.02
    b = 0.12
    c = 0.0
    g = 6

    # numerics
    nt = 8 # number of time steps
    dt = 1.0 # time step (a)

    # time loop
    H = copy(H0)
    time = 0
    for it in 1:nt
        time = it*dt
        # precipitation rate
        P = glacier_P.(H, zref, Pref, dP_dz) # mm / a

        # temperature
        T = glacier_T.(H, zref, Tref, dT_dz) # C

        icemask = H.>B
        @show sum(icemask)
        sum(icemask)==0 && break
        # annual mass balance
        MB = rhow .* A .* icemask .* glacier_mb.(P, T, m) .* 1e-3 # kg / a
        # integrate over glacier and over time step
        Ba = sum(MB) * dt
        @show it, Ba

        hmin, hmax = extrema(H[icemask])
        hr = (hmax .- H) ./ (hmax - hmin)
        dh = glacier_dh.(hr, a, b, c, g)

        fs = Ba / (rhoi * sum(dh .* A .* icemask)) # m

        dhh = H .- max.(H .+ fs .* dh .* icemask, B)
        a, b = sum(fs .* dh .* A .* icemask) * rhoi, sum(dhh.*A.*rhoi) 
        @show (Ba-a)/Ba, (Ba-b)/Ba

        H .= max.(H .+ fs .* dh .* icemask, B)
    end

    return H, time
end

preprocess_flowline("scripts/flowline_rhone_2007.txt", "scripts/rhone_2007.txt")

length, surface, bed, area = read_flowline("scripts/rhone_2007.txt")

Pref = 1000.0 # mm / a
Tref = 0.0    # C

zref  = 3100.0 # ELA
dP_dz = 0.0
dT_dz = -0.65/100 # C/m

a = -0.02
b = 0.12
c = 0.0
g = 6

icemask = H.>B
hmin, hmax = extrema(surface[icemask])
hr = (hmax .- surface) ./ (hmax - hmin)
dh = glacier_dh.(hr, a, b, c, g)

x = LinRange(0, 1, 100)
y = glacier_dh.(x, a, b, c, g)

        P = glacier_P.(surface, zref, Pref, dP_dz) # mm / a

        # temperature
        T = glacier_T.(surface, zref, Tref, dT_dz) # C
glacier_mb.(P, T, m)

@time surface1, tt = forward_model(surface, bed, area, Pref, Tref, zref, dP_dz, dT_dz)

fig = Figure(size=(600, 300))
axs = (Axis(fig[1, 1]; yreversed=true),
       Axis(fig[1, 2]))

plt = (lines!(axs[1], x, y; color=:red, linewidth=1),
       scatter!(axs[1], hr, dh; markersize=5),
       lines!(axs[2], length, bed),
       lines!(axs[2], length, surface),
       lines!(axs[2], length, surface1, linewidth=3))

display(fig)

