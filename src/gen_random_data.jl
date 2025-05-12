using Random, HDF5

Random.seed!(1234)

N = 100_000
range_min = -1000.0
range_max = 1000.0

x = rand(Float32, N) .* (range_max - range_min) .+ range_min
y = rand(Float32, N) .* (range_max - range_min) .+ range_min
z = rand(Float32, N) .* (range_max - range_min) .+ range_min

# Save to the data/ folder
output_file = "data/random_points.h5"
h5write(output_file, "x", x)
h5write(output_file, "y", y)
h5write(output_file, "z", z)