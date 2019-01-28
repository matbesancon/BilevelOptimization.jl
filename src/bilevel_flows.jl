"""
Data structure and resolution of special types of bilevel flow problems
"""
module BilevelFlowProblems

using SparseArrays: spzeros
using JuMP: @variable, @constraint, @objective
import JuMP

const MT = AbstractMatrix{<:Real}
const VT = AbstractVector{<:Real}

export BilevelFlowProblem

import ..build_blp_model

"""
Holding the information of a bilevel flow problem:
* Upper-level picks arcs to tax and a discrete tax level
* Lower-level solves a mincost flow problem
"""
struct BilevelFlowProblem{M1<:MT,M2<:AbstractMatrix{Bool},M3<:MT,M4<:AbstractArray{<:Real,3}}
    init_cost::M1
    taxable_edges::M2
    capacities::M3
    tax_options::M4
    nv::Int
    ne::Int
    minflow::Float64
    function BilevelFlowProblem(init_cost::M1,taxable_edges::M2,capacities::M3,tax_options::M4, minflow) where {M1<:MT,M2<:AbstractMatrix{Bool},M3<:MT,M4<:AbstractArray{<:Real,3}}
        size(init_cost) == size(taxable_edges) == size(capacities) == size(tax_options)[1:end-1] || DimensionMismatch("Matrix dimensions")
        nv = size(init_cost)[1]
        nv == size(init_cost)[2] || DimensionMismatch("Cost matrix should be square")
        ne = sum(1 for i in Base.OneTo(nv) for j in Base.OneTo(nv) if capacities[i,j] > 0. && i!=j)
        new{M1,M2,M3,M4}(init_cost, taxable_edges, capacities, tax_options, nv, ne, minflow)
    end
end

function build_blp_model(bfp::BilevelFlowProblem, solver)
    m = JuMP.Model(solver = solver)
    (nv, _, nopt) = size(bfp.tax_options)
    @variable(m, y[i=1:nv,j=1:nv,k=1:nopt], Bin)
    @variable(m, x[i=1:nv,j=1:nv] >= 0.)
    @variable(m, f[i=1:nv,j=1:nv] >= 0.)
    @constraint(m, limit2[i=1:nv,j=1:nv,k=1:nopt], x[i,j] <= bfp.tax_options[i,j,k] * f[i,j] + maximum(bfp.tax_options[i,j,:]) * bfp.capacities[i,j]*(1-y[i,j,k]))
    @constraint(m, unique_opt[i=1:nv,j=1:nv], sum(y[i,j,k] for k in 1:nopt) == 1.)
    @objective(m, Max,
        sum(x[i,j] for i in 1:nv for j in 1:nv if bfp.taxable_edges[i,j])
    )
    flat_flow = [f[i] for i in eachindex(f)]
    lin_cost  = [bfp.init_cost[i,j] + sum(y[i,j,k] * bfp.tax_options[i,j,k] for k in Base.OneTo(nopt)) for j in Base.OneTo(bfp.nv) for i in Base.OneTo(bfp.nv)]
    (B, b) = flow_constraint_standard(bfp)
    # sc = maximum(bfp.capacities) # upper bound on capacity slack
    # # sbound = map(eachindex(b)) do i
    # #     if i == 1
    # #         1. # is never greater than 0. anyway
    # #     elseif i < bfp.nv
    # #         0. # slack at 0 for flow conservation constraints
    # #     else
    # #         sc
    # #     end
    # # end
    @variable(m,
        s[i=1:length(b)] >= 0. # slack at 0 if flow constraints at nodes
    )
    @constraint(m, B*flat_flow .+ s .== b)
    (_, λ) = build_blp_model(m, B, lin_cost, s)
    return (m, x, y, f, λ)
end

function flow_constraint_standard(bfp::BilevelFlowProblem)
    B = spzeros( # constraint matrix
        1 +        # min flow from source
        bfp.nv-2 + # flow conservation at all nodes but source & sink
        bfp.nv*bfp.nv,    # capacity at all edges
        bfp.nv*bfp.nv
    )
    b = zeros(size(B)[1])
    b[1] = -bfp.minflow
    if bfp.capacities[1,bfp.nv] > 0.
        B[1,end-bfp.nv+1] = -1
    end
    for k in 2:bfp.nv-1
        if bfp.capacities[1,k] > 0.
            B[1,bfp.nv*(k-1)+1] = -1. # arc s -> k
        end
        for i in Base.OneTo(bfp.nv)
            if i != k
                if bfp.capacities[i,k] > 0.
                    B[k,bfp.nv*(k-1)+i] = 1.
                end
                if bfp.capacities[k,i] > 0.
                    B[k,bfp.nv*(i-1)+k] = -1.
                end
            end
        end
    end
    for i in Base.OneTo(bfp.nv), j in Base.OneTo(bfp.nv)
        B[bfp.nv-1+i+bfp.nv*(j-1),i+bfp.nv*(j-1)] = 1.
        b[bfp.nv-1+i+bfp.nv*(j-1)] = bfp.capacities[i,j]
    end
    sc = sum(bfp.capacities) # upper bound on slack variable

    return (B, b)
end

end
