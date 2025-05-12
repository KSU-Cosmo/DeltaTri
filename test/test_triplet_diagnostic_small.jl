include("../src/triplet_diagnostic.jl")
using Test

# Small test dataset
x = Float32[-999.0, -998.0, -999.0, 999.5]
y = Float32[-999.0, -999.0, -998.0, 999.5]
z = Float32[-999.0, -999.0, -999.0, 999.5]

bin_edges = [0.0, 1.1, 2.0]
n_angle_bins = 2
boxsize = 1.0

expected_nonper = Dict{Tuple{Int,Int}, Matrix{Int}}()
expected_nonper[(1,1)] = [0 2; 0 0]
expected_nonper[(1,2)] = [2 2; 0 2]
expected_nonper[(2,2)] = [0 2; 0 0]

expected_per = Dict{Tuple{Int,Int}, Matrix{Int}}()
expected_per[(1,1)] = [0 2; 0 0]
expected_per[(1,2)] = [2 2; 0 2]
expected_per[(2,2)] = [0 2; 0 0]


@testset "triplet_diagnostic full bin test" begin
    hist_nonper = triplet_diagnostic(x, y, z, bin_edges, n_angle_bins, periodic=false)
    hist_per    = triplet_diagnostic(x, y, z, bin_edges, n_angle_bins, periodic=true, boxsize=boxsize)

    for bin_key in keys(expected_nonper)
        @test hist_nonper[bin_key] == expected_nonper[bin_key]
    end
    for bin_key in keys(expected_per)
        @test hist_per[bin_key] == expected_per[bin_key]
    end
end
