# BilevelOptimization

[![Build Status](https://travis-ci.org/matbesancon/BilevelOptimization.jl.svg?branch=master)](https://travis-ci.org/matbesancon/BilevelOptimization.jl)
[![codecov.io](http://codecov.io/github/matbesancon/BilevelOptimization.jl/coverage.svg?branch=master)](http://codecov.io/github/matbesancon/BilevelOptimization.jl?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/matbesancon/BilevelOptimization.jl/badge.svg?branch=master)](https://coveralls.io/github/matbesancon/BilevelOptimization.jl?branch=master)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3249187.svg)](https://doi.org/10.5281/zenodo.3249187)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01278/status.svg)](https://doi.org/10.21105/joss.01278)

This package is a Julia toolbox based on JuMP.jl for solving
[bilevel optimization problems](https://en.wikipedia.org/wiki/Bilevel_optimization).
These are encountered in various applications, including power grids, security games,
market equilibria or chemical reaction optimization.

## Generic bilevel linear problems (BLP)

The generic problem can be written:
```julia
min_{x} F(x,y)
such that
  G * x + H * y ⩽ q
  y ∈ arg min_y {d^T y + x^T F * y
                 such that
                   A * x + B * y ⩽ b
                   y ⩾ 0
  }
  x_j integer ∀ j ∈ Jx
```

`x` represents the upper-level decision variable and `y` the lower-level one.
`y` is thus the solution to a parametric optimization sub-problem, depending
on the value of `x`.
The required data describing this problem are
the feasibility domains of the upper and lower level and the coefficients
of the objective functions. All these are regrouped within the `BilevelLP`
type of this package.

The formulation is made as general as possible
for the problem to remain approachable with plain Mixed-Integer Solvers
(CBC, GLPK, SCIP, Gurobi, CPLEX). For a simple linear-linear problem,
the user can set `Jx = ∅` and `F` as a zero matrix of appropriate dimension.
The problem can be made as complex as wanted **at the upper level**,
as long as JuMP and the solver used support the constraints and objective.

### Usage:

The main function is `build_bilevel_lp`, which will build the JuMP model
or modify it for the bilevel problem.

The signature:  
`build_blp_model(bp::BilevelOptimization.BilevelLP, solver; comp_method)`
builds the model from scratch. It will return a `Tuple`: `(m, x, y, λ, s)` with:
- `m` the JuMP model
- `x` the upper-level variable vector
- `y` the lower-level variable vector
- `λ` the dual of the lower-level problem
- `s` the lower-level slack variable vector

The function can also be called with a model already built:  
`build_blp_model(m::JuMP.Model, bp::BilevelOptimization.BilevelLP, x, y; comp_method)`
In which case it will add the lower-level optimality constraints, it returns the same tuple.

If the user is not willing to describe the whole problem using a `BilevelLP`,
the following signature can be used:
`build_blp_model(m::JuMP.Model, B::M, d, s; comp_method)`

With `B` the lower-level constraint matrix, d the lower-level objective and `s` the lower-level
slack variable. Only the KKT conditions are added to the model in that case (not the lower-level
feasibility constraints).

## Installation

The package can be installed using Julia `Pkg` tool:
```julia
julia> ]
(v1.0) pkg> add BilevelOptimization
```

You will also need an optimization solver up and running with [JuMP](https://github.com/juliaopt/JuMP.jl).

## Testing

Tests can be performed using `Pkg`:
```julia
julia> ]
(v1.0) pkg> test BilevelOptimization
```

## API documentation

From the Julia REPL, type `?` to show the help prompt, then type the
identifier you want the documentation for.

```julia
julia> import BilevelOptimization

help?> BilevelOptimization.BilevelLP
  A bilevel linear optimization problem of the form:

  min cx^T * x + cy^T * y
  s.t. G x + H y <= q
       x_j ∈ [xl_j,xu_j]
       x_j ∈ ℤ ∀ j ∈ Jx
       y ∈ arg min {
          d^T * y + x^T * F * y
          s.t. A x + B y <= b
               y_j ∈ [yl_j,yu_j]
          }

  Note that integer variables are allowed at the upper level. 
```

## Resolution method

The "hard" part of the reduction of a bilevel problem is the set of
complementarity constraints of the form:

```julia
λ ⋅ (b - Ax - By) = 0
```

These constraints cannot be handled directly, different methods have been
developed in the literature and implemented in this package.
The standard way is to give a different algorithm in `build_blp_model`:

```julia
build_blp_model(args..., comp_method::ComplementarityMethod = my_method)
```


### Special ordered sets 1

Special-ordered Sets of type 1 or [SOS1](https://en.wikipedia.org/wiki/Special_ordered_set)
are used by default for complementarity constraints.
The option to pass is:

```julia
build_blp_model(args..., comp_method = SOS1Complementarity())
```

### Dual and primal bounds

The most common technique for these constraints is the linearization of the
constraint with a formulation developed in Fortuny-Amat and McCarl, 1981,
using so-called big-M constraints.

```julia
build_blp_model(args..., comp_method = BoundComplementarity(MD, MP))
```

`MD`, `MP` are primal and dual bounds, both can be either a scalar
for one bound per variable type or an abstract vector for one bound per
variable.

### Custom method

Users can create a custom method for the complementarity constraint,
by creating a type `T <: ComplementarityMethod` (sub-typing is optional but
helps for clarity). They also need to implement a method:

```julia
add_complementarity_constraint(m, cm::T, s, λ)
```

Within which the complementarity constraints are added to the JuMP model `m`.

## The toll-setting problem

As a special application of the above model, the module `BilevelFlowProblems`
offers the following problem:
* The upper-level, acting as a leader of the [Stackelberg game](https://en.wikipedia.org/wiki/Stackelberg_competition), chooses taxes to set on some arcs of a directed graph.
* The lower-level, acting as the follower, makes a
[minimum-cost flow](https://en.wikipedia.org/wiki/Minimum-cost_flow_problem) with
a given minimum amount from the source to the sink.
* Each arc has an invariant base cost and a tax level decided upon by the leader.

This has been investigated in the literature as the "toll-setting problem".
The required data include:
- the initial cost of each arc for all `i,j`
- which edges can be taxed by the leader for all `i,j`
- the tax options (at which level can each edge be taxed) for all `i,j,k`
- flow capacities of each edge
- the minimum flow the follower has to pass from source to sink

### Example:

```julia
init_cost = [
	0. 1. 1. 4.
	0. 0. 0. 1.
	0. 0. 0. 1.
	0. 0. 0. 0.
]
taxable_edges = [
	false true true false
	false false false true
	false false false true
	false false false false
]
tax_options = zeros(4,4,5)
for i in 1:4, j in 1:4
	if taxable_edges[i,j]
    	tax_options[i,j,:] .= (0.0, 0.5, 1.0, 1.5, 2.0)
    end
end

capacities = [
	0. 3. 2. 3.
    0. 0. 0. 2.
    0. 0. 0. 1.
    0. 0. 0. 0.
]

minflow = 3.

BilevelFlowProblem(init_cost,taxable_edges,capacities,tax_options, minflow)

(m, r, y, f, λ) = build_blp_model(bfp, CbcSolver())
st = JuMP.solve(m)

# st === :Optimal
# getobjectivevalue(m) ≈ 6.
#     for j in 1:size(r)[2]
#        for i in 1:size(r)[1]
#            @test getvalue(r[i,j]) ≈ sum(getvalue(y[i,j,:]).*bfp.tax_options[i,j,:]) * getvalue(f[i,j])
#        end
#    end
```

## Questions, issues, contributions

Problems with the package and its usage can be explained through Github issues,
ideally with a minimal working example showing the problem.
Pull requests (PR) are welcome.

Please read detailed information in **CONTRIBUTING.md**.

## Related packages

* [Complementarity.jl](https://github.com/chkwon/Complementarity.jl) solving a generic class
including bilevel problems using non-linear techniques
* [MibS](https://github.com/coin-or/MibS) for problems where the lower-level also includes integer variables. KKT conditions can therefore not be used and other branching and cutting plane techniques are leveraged.
* [YALMIP](https://yalmip.github.io/tutorial/bilevelprogramming/) includes a bilevel solver and offers roughly the same features (and a bit more) as BilevelOptimization.jl

## Citing

See *CITATION.bib*, prefer citing the paper published in the Journal of Open-Source Software.
