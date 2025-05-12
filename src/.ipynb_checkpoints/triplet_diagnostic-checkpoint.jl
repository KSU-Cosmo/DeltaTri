using HDF5

function triplet_diagnostic(x::Vector{Float32}, y::Vector{Float32}, z::Vector{Float32},
                                 bin_edges::Vector{Float64}, n_angle_bins::Int;
                                 periodic::Bool=false, boxsize::Float64=0.0)
    N = length(x)
    Nr = length(bin_edges) - 1
    Na = n_angle_bins

    if periodic && boxsize <= 0.0
        error("If periodic=true, boxsize must be > 0")
    end

    # Precompute bin pairs
    histograms = Dict{Tuple{Int, Int}, Matrix{Int}}()
    for r1 in 1:Nr
        for r2 in r1:Nr
            histograms[(r1, r2)] = zeros(Int, Na, Na)
        end
    end

    # Loop over each object
    for i in 1:N
        counts = [zeros(Int, Na) for _ in 1:Nr]
        xi, yi, zi = x[i], y[i], z[i]

        for j in 1:N
            if i == j
                continue
            end
            dx = x[j] - xi
            dy = y[j] - yi
            dz = z[j] - zi

            if periodic
                dx = mod(Float64(dx) + boxsize/2, boxsize) - boxsize/2
                dy = mod(Float64(dy) + boxsize/2, boxsize) - boxsize/2
                dz = mod(Float64(dz) + boxsize/2, boxsize) - boxsize/2
            end

            r = sqrt(dx^2 + dy^2 + dz^2)
            if r >= bin_edges[end] - 1e-6
                continue
            end

            r_bin = Int(searchsortedfirst(bin_edges, r))
            if r_bin == 1 || r_bin > length(bin_edges)
                continue  # out of bounds
            end
            r_bin -= 1

            mu = abs(dz / r)  # |cos(theta)|
            a_bin = clamp(Int(floor(mu * Na)) + 1, 1, Na)

            counts[r_bin][a_bin] += 1
        end

        # Cross products of angular bins for each radial pair
        for r1 in 1:Nr
            for r2 in r1:Nr
                c1 = counts[r1]
                c2 = counts[r2]

                for a1 in 1:Na
                    for a2 in 1:Na
                        histograms[(r1, r2)][a1, a2] += c1[a1] * c2[a2]
                    end
                end
            end
        end
    end

    return histograms
end

function save_histograms_to_hdf5(filename::String, histograms::Dict{Tuple{Int, Int}, Matrix{Int}})
    h5open(filename, "w") do file
        for ((r1, r2), mat) in histograms
            dset_name = "r$(r1)_r$(r2)"
            write(file, dset_name, mat)
        end
    end
end

function load_histograms_from_hdf5(filename::String)
    histograms = Dict{Tuple{Int, Int}, Matrix{Int}}()
    h5open(filename, "r") do file
        for name in keys(file)
            if startswith(name, "r")
                parts = split(name, "_")
                r1 = parse(Int, parts[1][2:end])
                r2 = parse(Int, parts[2][2:end])
                histograms[(r1, r2)] = read(file[name])
            end
        end
    end
    return histograms
end
