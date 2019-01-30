# Contributing to BilevelOptimization.jl

Contributions to BilevelOptimization.jl are welcome both for direct
improvement and issue reports.

## Issues / bug reports

Issues should be raised on GitHub at:
https://github.com/matbesancon/BilevelOptimization.jl/issues  

If possible, an issue should come with a Minimum Working Example (MWE),
the expected behavior and observed error.
Since the package depends on optimization solvers, the example should
preferably be given either:
- in a solver-independent example
- with an open-source solver (CBC, GLPK) and corresponding JuMP wrapper.
So that any contributor without access to commercial software can
reproduce the issue.

## Questions

If you have any questions about BilevelOptimization.jl,
feel encouraged to ping **@matbesancon** on the
[JuliaLang Slack](https://slackinvite.julialang.org/) in the
`#maths-optimization` channel, or on the
[JuliaLang Discourse forum](https://discourse.julialang.org/c/domain/opt)
again in the optimization channel.  

It is also permissible to ask questions by opening issues on this GitHub
repository. Ideally, such questions be phrased as requests for documentation.
If your question is not answered by the docs and README, then it is a sign the
documentation could be improved.  

For help on the JuMP ecosystem itself, the Discourse and Slack channels are
also good starting points. The [JuMP.jl](https://github.com/juliaopt/jump.jl)
also contains tutorial notebooks.

## Contributing Code

Code contributions are welcome and encouraged, including for example
adding support for other classes of bilevel problems.  

Code contributions should be made in the normal way of making a pull request.
In general they should try to match the style of the code already present.
This follows the
[JuliaLang/julia repo conventions](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md#code-formatting-guidelines).
