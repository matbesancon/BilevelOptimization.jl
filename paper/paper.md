---
title: 'A Julia package for bilevel optimization problems'
tags:
  - Julia
  - JuMP
  - optimization
  - game theory
  - bilevel optimization
authors:
  - name: Mathieu Besançon
    orcid: 0000-0002-6284-3033
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
affiliations:
 - name: École Polytechnique de Montréal, QC, Canada
   index: 1
 - name: GERAD, QC, Canada
   index: 2
 - name: INOCS, INRIA Lille Nord-Europe, France
   index: 3
date: 29 January 2019
bibliography: paper.bib
---

# Summary

Mathematical optimization is the discipline dealing with
the determination of the best
(or almost best) decision with respect to a specific cost function and to
a set of constraints on the decision.
Bilevel optimization is a class of mathematical optimization problems
with the optimality conditions of a lower-level problem embedded in the
constraints. BilevelOptimization.jl is a toolbox built on top of the JuMP.jl
modeling package [@dunning2017jump].
Bilevel optimization is used to tackle various problems in areas such as
power systems, security applications, network design or market equilibria.
See [@dempe2018bilevel] for an overview of applications and recent
formulations and theoretical progress.  

The computation of an optimal solution to a bilevel problem is in general hard.
Even with all the constraints and the objectives at the two level being linear,
the resulting problem is non-convex and NP-hard, with possibly a disjoint
feasible set. Optimization practitioners often rely on problem-specific
properties and modeling techniques or heuristics, the goal of this package
is to offer a both flexible model of a general class of bilevel problems
and a solving method which is compatible with the JuMP workflow.    

# Bilevel optimization

A generic formulation for a bilevel problem is:

$$\min_{x} F(x,y)
\text{s.t.}
G_i(x,y) \leq 0 \forall i \in \{1..m_u\}
y \in arg \min_y \{ f(x,y),
                   g_i(x,y) \leq 0 \forall i \in \{1..m_l\}
                 \}.$$

If the lower-level problem is convex, i.e. if the functions $f(x,y)$ and
$g_i(x,\cdot)$ are convex and if Slater's qualification constraints hold,
the Karush-Kuhn-Tucker conditions can be used to characterize the optimality
at the lower-level.  

For most lower-level problems, there are several optimal solutions
(different solutions yielding the same optimal value of the objective).
Several methodologies have been developed for such case, the two principal
being the optimistic and pessimistic bilevel formulations, turning the
set-valued problem into a regular one. The approach used for now in
BilevelOptimization.jl is the optimistic one, allowing for more
easily reformulated problems.  

This package is initially designed for a restricted form:
$$\min_{x} c_x^T x + c_y^T y
\text{s.t.}
G x + H y \leq q
x \geq 0
x_j \in \mathcal{Z}_+ \forall j \in Jx
y \in arg \min_y \{ d^T y + x^T F y,
                   A x + B y \leq b
                   y \geq 0
                 \}.$$

The single-level reduction of the optimistic version of this problem is:
$$\min_{x} c_x^T x + c_y^T y
\text{s.t.}
G x + H y \leq q
A x + B y \leq b
x \geq 0
y \geq 0
\lambda \geq 0
x_j \in \mathcal{Z}_+ \forall j \in J_x
d + F^T x + B^T \lambda = 0
\lambda_i \cdot (b_i - Ax_i - By_i) = 0 \forall i \in \{1..m_l\}.$$

The last equation is a complementarity constraint, corresponding
to the fact that at least one of $(\lambda_i, s_i)$ has to be equal
to zero. This non-convex, non-linear constraint cannot be tackled
efficiently by common optimization solvers and needs to be linearized.
The two common approaches are linearization using a binary variable and
"big-M" primal and dual upper bounds [@pineda19] and Special Ordered Sets
of type 1 (SOS1). Due to their greater flexibility, the latter option is used
in BilevelOptimization.jl, thus not requiring users to provide additional
parameters. A special ordered set 1 is a group of two or more variables,
of which at most one can be non-zero. Mixed-Integer Linear solvers use this
information for branching directly on the two variables.
JuMP supports the modeling of Special Ordered Sets 1 with the following syntax:

```julia
for i in 1:ml
    JuMP.addSOS1(m, [λ[i], s[i]])
end
```

and passes on the information to the solver. In the case of the bilevel
problem presented above, the sets contain the slack variable and dual variable
associated with each lower-level constraint, forcing at least one of them to 0.  

# Application to toll-setting problems

The toll-setting problem is a class of bilevel optimization where the two
levels of decision are taken on a graph [@brotcorne2001bilevel].
It belongs to the more general framework of network pricing problems with
applications in road management [@harks2018toll] or telecommunication
network reliability [@Hayrapetyan2007].  

In this problem, the upper level decides on a toll to apply on some arcs
of a directed graph. Each arc has an initial cost, the lower-level then
finds the minimum-cost flow from a source to a sink with a minimum circulating
flow. This problem can entirely be modeled using the framework
presented above, a composite type holding all required data is defined
in the package, allowing users to bypass the re-formulation of the model
from its algebraic JuMP form to a standard form.

# More general problem formulations

Even though BilevelOptimization.jl is designed initially for linear-linear
bilevel problems, the design allows users to bypass the upper-level problem
specification by providing a JuMP `Model` with the upper-level objective
and constraints already set, for instance for quadratic or conic upper level
formulations. The only requirement is that the solver must support
the type of constraints specified and special ordered sets 1.
This flexibility allows users to leverage some recent advances on
mixed-integer convex optimization and solvers tackling these problems
[@LubinYamangilBentVielma2016]. As of the current state of BilevelOptimization,
the only restricted part of the model is the linear-quadratic lower-level.

# References
