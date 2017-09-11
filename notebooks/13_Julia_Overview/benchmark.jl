# number of runs for each test
runs = 10

# Julia benchmark
include("benchmark_julia.jl")

# R benchmark
using RCall

@rput runs
R"""
source('R-benchmark-25.R')
"""
@rget times_R

# Markdown table

using DataFrames

results = DataFrame(
    Test = ["Matrix creation, trans., deform. (2500 x 2500)",
    "Power of matrix (2500 x 2500, `A.^1000`)",
    "Quick sort (\$n = 7 \\times 10^6\$)",
    "Cross product (2800 x 2800, \$A^TA\$)",
    "LS solution (\$n = p = 2000\$)",
    "FFT (\$n = 2,400,000\$)",
    "Eigen-values (\$600 \\times 600\$)",
    "Determinant (\$2500 \\times 2500\$)",
    "Cholesky (\$3000 \\times 3000\$)",
    "Matrix inverse (\$1600 \\times 1600\$)",
    "Fibonacci (vector calculation)",
    "Hilbert (matrix calculation)",
    "GCD (recursion)",
    "Toeplitz matrix (loops)",
    "Escoufiers (mixed)"],
    R_v3_4_1 = vec(times_R),
    Julia_v0_6_0 = vec(times_julia),
    Speedup = vec(times_R) ./ vec(times_julia)
    )

println("| Test | R 3.4.1 | Julia 0.6.0 | Speedup |  ")
println("|:-------- |:-------:|:-------:|:-------:|  ")
for r in 1:size(results, 1)
    @printf("| %s | %4.2f | %4.2f | %4.2f |  \n", results[r, :Test], results[r, :R_v3_4_1], results[r, :Julia_v0_6_0], results[r, :Speedup])
end
