# Julia port of the R Benchmark 2.5
# http://r.research.att.com/benchmarks/R-benchmark-25.R
versioninfo()

using BenchmarkTools

srand(123) # seed

# Number of times the tests are executed
#runs = 3
BenchmarkTools.DEFAULT_PARAMETERS.samples = runs

times_julia = zeros(5, 3)

println("\n\n   Julia Benchmark")
println("   ===============")
println("Number of times each test is run__________________________: $runs")
println("\n")

#************************ I. Matrix calculation ********************#

println("   I. Matrix calculation")
println("   ---------------------")

# (1)
function test01()
    a = randn(2500, 2500) / 10
    b = a.'
    b = reshape(b, 1250, 5000)
    a = b.'
end

t = @benchmark test01()
times_julia[1, 1] = mean(t).time / 1e9
println("Creation, transp., deformation of a 2500x2500 matrix (sec): $(times_julia[1, 1])")

# (2)
t = @benchmark a.^1000 setup=(a = abs.(randn(2500, 2500) / 2))
times_julia[2, 1] = mean(t).time / 1e9
println("2500x2500 normal distributed random matrix ^1000____ (sec): $(times_julia[2, 1])")

# (3)
t = @benchmark sort(a, alg = QuickSort) setup=(a = randn(7000000))
times_julia[3, 1] = mean(t).time / 1e9
println("Sorting of 7,000,000 random values__________________ (sec): $(times_julia[3, 1])")

# (4)
a = randn(2800, 2800)
t = @benchmark a' * a setup=(a = randn(2800, 2800))
times_julia[4, 1] = mean(t).time / 1e9
println("2800x2800 cross-product matrix (b = a' * a)_________ (sec): $(times_julia[4, 1])")

# (5)
t = @benchmark (a' * a) \ (a' * b) setup=(a = randn(2000, 2000); b = collect(1.0:2000.0))
times_julia[5, 1] = mean(t).time / 1e9
println("Linear regr. over a 2000x2000 matrix (c = a \\ b')___ (sec): $(times_julia[5, 1])")


#************************ II. Matrix functions ********************#

println("   II. Matrix functions");
println("   ---------------------");

# (1)
t = @benchmark fft(a) setup = (a = randn(2400000))
times_julia[1, 2] = mean(t).time / 1e9
println("FFT over 2,400,000 random values____________________ (sec): $(times_julia[2, 1])")

# (2)
t = @benchmark eigvals(a) setup = (a = randn(600, 600))
times_julia[2, 2] = mean(t).time / 1e9
println("Eigenvalues of a 600x600 random matrix______________ (sec): $(times_julia[2, 2])")

# (3)
t = @benchmark det(a) setup = (a = randn(2500, 2500))
times_julia[3, 2] = mean(t).time / 1e9
println("Determinant of a 2500x2500 random matrix____________ (sec): $(times_julia[3, 2])")

# (4)
t = @benchmark chol(a) setup = (tmp = randn(3000, 3000); a = tmp' * tmp)
times_julia[4, 2] = mean(t).time / 1e9
println("Cholesky decomposition of a 3000x3000 matrix________ (sec): $(times_julia[4, 2])")

# (5)
t = @benchmark inv(a) setup = (a = randn(1600, 1600))
times_julia[5, 2] = mean(t).time / 1e9
println("Inverse of a 1600x1600 random matrix________________ (sec): $(times_julia[5, 2])")


#************************ III. Programmation ********************#

println("   III. Programmation");
println("   ---------------------");

# (1)
function fibfun(a)
  const ϕ = 1.6180339887498949
  return @. (ϕ ^ a - (-ϕ) ^ (-a)) / sqrt(5)
end
t = @benchmark fibfun(a) setup = (a = floor.(rand(3500000) * 1000))
times_julia[1, 3] = mean(t).time / 1e9
println("3,500,000 Fibonacci numbers calculation (vector calc)(sec): $(times_julia[1, 3])")

# (2)
function hilbertfun(n)
    [1 / (i + j - 1) for i in 1:n, j in 1:n]
end
t = @benchmark hilbertfun(3000)
times_julia[2, 3] = mean(t).time / 1e9
println("Creation of a 3000x3000 Hilbert matrix (matrix calc) (sec): $(times_julia[2, 3])")

# (3)
function gcd2(x, y)
  if sum(y .> 1.0e-4) == 0
    return x
  else
    for i in eachindex(y)
      y[i] == 0 && (y[i] = x[i])
    end
    gcd2(y, x .% y)
  end
end
t = @benchmark gcd2(a, b) setup = (a = ceil.(rand(400000) * 1000); b = ceil.(rand(400000) * 1000))
times_julia[3, 3] = mean(t).time / 1e9
println("Grand common divisors of 400,000 pairs (recursion)__ (sec): $(times_julia[3, 3])")

# (4)
function toeplitzfun(n)
  b = zeros(n, n)
  for j = 1:500
    for k = 1:500
      b[k, j] = abs(j - k) + 1
    end
  end
  b
end
t = @benchmark toeplitzfun(500)
times_julia[4, 3] = mean(t).time / 1e9
println("Creation of a 500x500 Toeplitz matrix (loops)_______ (sec): $(times_julia[4, 3])")

# (5)
function escoufier(x)
  p = size(x, 1)
  vt = collect(1:p)
  vr = Int[]
  RV = collect(1.0:p)
  vrt = 0
  RxytRxy = zeros(p, p)
  RyyRyy = zeros(p, p)
  for j = 1:p
    Rvmax = 0
    for k = 1:(p-j+1)
      x2 = [x x[:, vr] x[:, vt[k]]]
      R = cor(x2)
      @views Ryy = R[1:p, 1:p]
      @views Rxx = R[(p+1):(p+j), (p+1):(p+j)]
      @views Rxy = R[(p+1):(p+j), 1:p]
      At_mul_B!(RxytRxy, Rxy, Rxy)
      A_mul_B!(RyyRyy, Ryy, Ryy)
      rvt = trace(RxytRxy) / sqrt(trace(RyyRyy) * trace(Rxx * Rxx))
      if rvt > Rvmax
        Rvmax = rvt
        vrt = vt[k]
      end
    end
    append!(vr, [vrt])
    RV[j] = Rvmax
    vt = vt[vt .≠ vr[j]]
  end
end
t = @benchmark escoufier(x) setup = (x = abs.(randn(45, 45)))
times_julia[5, 3] = mean(t).time / 1e9
println("Escoufier's method on a 45x45 matrix (mixed)________ (sec): $(times_julia[5, 3])")
