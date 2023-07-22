# Optimization and parallelization
# !!!
using BenchmarkTools

@benchmark include("sod_shock.jl")

#= !!!
Devs say that using global causes poor performance, and instead to use const
!global benchmark 
BenchmarkTools.Trial: 
  memory estimate:  5.29 GiB
  allocs estimate:  51848546
  --------------
  minimum time:     6.283 s (4.64% GC)
  median time:      6.283 s (4.64% GC)
  mean time:        6.283 s (4.64% GC)
  maximum time:     6.283 s (4.64% GC)
  --------------
  samples:          1
  evals/sample:     1

!const benchmark
BenchmarkTools.Trial: 
  memory estimate:  5.28 GiB
  allocs estimate:  51671465
  --------------
  minimum time:     4.704 s (5.62% GC)
  median time:      4.729 s (5.56% GC)
  mean time:        4.729 s (5.56% GC)
  maximum time:     4.754 s (5.50% GC)
  --------------
  samples:          2
  evals/sample:     1

Sure looks like a performance improvement to me!
=#

#= !!!
Explicit type declarations is also supposed to improve performance
Previous benchmark values are with implicit declarations
=#

#= !!!
Using @. instead of .* .+ ./ etc is supposed to be significantly faster

BenchmarkTools.Trial: 
  memory estimate:  5.15 GiB
  allocs estimate:  52712081
  --------------
  minimum time:     4.476 s (7.09% GC)
  median time:      4.588 s (6.65% GC)
  mean time:        4.588 s (6.65% GC)
  maximum time:     4.699 s (6.23% GC)
  --------------
  samples:          2
  evals/sample:     1

Not seeing much of a difference here
With functional tests I'm not really seeing 
    any performance difference between the two methods, 
    but the @. macro is much easier to write
=#