# BilevelOptimization

[![Build Status](https://travis-ci.org/matbesancon/BilevelOptimization.jl.svg?branch=master)](https://travis-ci.org/matbesancon/BilevelOptimization.jl)
[![codecov.io](http://codecov.io/github/matbesancon/BilevelOptimization.jl/coverage.svg?branch=master)](http://codecov.io/github/matbesancon/BilevelOptimization.jl?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/matbesancon/BilevelOptimization.jl/badge.svg?branch=master)](https://coveralls.io/github/matbesancon/BilevelOptimization.jl?branch=master)

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

## Installation

The package can be installed using Julia `Pkg` tool:
```julia
julia> ]
(v1.0) pkg> add https://github.com/matbesancon/BilevelOptimization.jl
```

You will also need an optimization solver up and running with [JuMP](https://github.com/juliaopt/JuMP.jl).

## Testing
The tests require the [cbc solver](https://projects.coin-or.org/Cbc) installed on the system, along with its corresponding `julia` package.

```bash
# Grab the tests (after the cbc solver is installed)
mkdir bilvlopt_test && cd bilvlopt_test
wget https://github.com/matbesancon/BilevelOptimization.jl/raw/master/test/runtests.jl
julia
```
```julia
julia> import Pkg; Pkg.add("Cbc")
julia> exit()
```
```bash
# In the same folder (bilvlopt_test)
julia runtests.jl
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

## Questions, issues, contributions

Problems with the package and its usage can be explained through Github issues,
ideally with a minimal working example showing the problem.
Pull requests (PR) are welcome.

Please read detailed information in **CONTRIBUTING.md**.

## Related packages

* [Complementarity.jl](https://github.com/chkwon/Complementarity.jl solving a generic class
including bilevel problems using non-linear techniques
* [MibS](https://github.com/coin-or/MibS) for problems where the lower-level also includes integer variables. KKT conditions can therefore not be used and other branching and cutting plane techniques are leveraged.
* [YALMIP](https://yalmip.github.io/tutorial/bilevelprogramming/) includes a bilevel solver and offers roughly the same features (and a bit more) as BilevelOptimization.jl

## Citing

```
@misc{besancon19,
    author = {{Mathieu Besançon}},
    title  = "BilevelOptimization.jl, a JuMP-based toolbox for bilevel optimization",
    url = {https://github.com/matbesancon/BilevelOptimization.jl},
    version = {0.1},
    year = {2019}
}
```

A software paper may be written.
