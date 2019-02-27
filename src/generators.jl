module BilevelGenerators

import Random
using Distributions: Uniform
import ..BilevelLP

export UniformBilevelGenerator

"""
Size and bound specifications for generating random `BilevelLP` instances
Each bound is stored in a tuple `(lower, upper)`
"""
struct UniformBilevelGenerator
    ml::Int
    mu::Int
    nl::Int
    nu::Int
    cx_bounds::Tuple{Float64,Float64}
    cy_bounds::Tuple{Float64,Float64}
    G_bounds::Tuple{Float64,Float64}
    H_bounds::Tuple{Float64,Float64}
    q_bounds::Tuple{Float64,Float64}
    d_bounds::Tuple{Float64,Float64}
    A_bounds::Tuple{Float64,Float64}
    B_bounds::Tuple{Float64,Float64}
    b_bounds::Tuple{Float64,Float64}
end

"""
Default constructor only specifies sizes for each element of the problem and
the same lower and upper bounds for each of them, passed as keyword arguments `lb` and `ub`
"""
function UniformBilevelGenerator(ml,mu,nl,nu; lb = 0., ub = 1.)
    UniformBilevelGenerator(ml,mu,nl,nu, (lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub))
end

"""
Generate a new random instance from a `UniformBilevelGenerator`
"""
function Base.rand(rng::Random.AbstractRNG, ubg::UniformBilevelGenerator)
    BilevelLP(
        rand(rng, Uniform(ubg.cx_bounds...), ubg.nu),
        rand(rng, Uniform(ubg.cy_bounds...), ubg.nl),
        rand(rng, Uniform(ubg.G_bounds...), ubg.mu, ubg.nu),
        rand(rng, Uniform(ubg.H_bounds...), ubg.mu, ubg.nl),
        rand(rng, Uniform(ubg.q_bounds...), ubg.mu),
        rand(rng, Uniform(ubg.d_bounds...), ubg.nl),
        rand(rng, Uniform(ubg.A_bounds...), ubg.ml, ubg.nu),
        rand(rng, Uniform(ubg.B_bounds...), ubg.ml, ubg.nl),
        rand(rng, Uniform(ubg.b_bounds...), ubg.ml)
    )
end

Base.rand(ubg::UniformBilevelGenerator) = Base.rand(Random.GLOBAL_RNG, ubg)

end
