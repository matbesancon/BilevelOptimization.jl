module BilevelGenerators

import Random
using Distributions: Uniform
import ..BilevelLP

export UniformBilevelGenerator

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

function UniformBilevelGenerator(ml,mu,nl,nu; lb = 0., ub = 1.)
    UniformBilevelGenerator(ml,mu,nl,nu, (lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub),(lb,ub))
end

function Base.rand(rng::Random.AbstractRNG, ubg::UniformBilevelGenerator)
    BilevelLP(
        rand(Uniform(ubg.cx_bounds...), ubg.nu),
        rand(Uniform(ubg.cy_bounds...), ubg.nl),
        rand(Uniform(ubg.G_bounds...), ubg.mu, ubg.nu),
        rand(Uniform(ubg.H_bounds...), ubg.mu, ubg.nl),
        rand(Uniform(ubg.q_bounds...), ubg.mu),
        rand(Uniform(ubg.d_bounds...), ubg.nl),
        rand(Uniform(ubg.A_bounds...), ubg.ml, ubg.nu),
        rand(Uniform(ubg.B_bounds...), ubg.ml, ubg.nl),
        rand(Uniform(ubg.b_bounds...), ubg.ml)
    )
end

Base.rand(ubg::UniformBilevelGenerator) = Base.rand(Random.GLOBAL_RNG, ubg)

end
