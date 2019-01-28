# BilevelOptimization

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

## Resolution method

`BilevelOptimization.jl` uses Special-ordered Sets of type 1 or [SOS1](https://en.wikipedia.org/wiki/Special_ordered_set) for complementarity constraints.
This avoids the bound estimation phase which is often tricky for dual variables.
This avoids solving sub-problems to estimate primal and dual bound and
still allows the solver to branch on either `(λ,s)` efficiently.

## The toll-setting problem

As a special application of the above model, the module `BilevelFlowProblems`
offers the following problem:
* The upper-level, acting as a leader of the [Stackelberg game](https://en.wikipedia.org/wiki/Stackelberg_competition), chooses taxes to set on some arcs of a directed graph.
* The lower-level, acting as the follower, makes a
[minimum-cost flow](https://en.wikipedia.org/wiki/Minimum-cost_flow_problem) with
a given minimum amount from the source to the sink.
* Each arc has an invariant base cost and a tax level decided upon by the leader.

This has been investigated in the literature as the "toll-setting problem".

## Related packages

* [Complementarity.jl](https://github.com/chkwon/Complementarity.jl solving a generic class
including bilevel problems using non-linear techniques
* [MibS](https://github.com/coin-or/MibS) for problems where the lower-level also includes integer variables. KKT conditions can therefore not be used and other branching and cutting plane techniques are leveraged.
* [YALMIP](https://yalmip.github.io/tutorial/bilevelprogramming/) includes a bilevel solver and offers roughly the same features as BilevelOptimization.jl

## Citing

```
@misc{bilevel19,
    author = {{Mathieu Besançon}},
    title  = "BilevelOptimization.jl, a JuMP-based toolbox for bilevel optimization",
    url = {https://github.com/matbesancon/BilevelOptimization.jl},
    version = {0.1},
    year = {2019}
}
```

A software paper may be written.
