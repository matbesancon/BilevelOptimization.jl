"""
Data structure and resolution of special types of bilevel flow problems
"""
module BilevelFlowProblems

using SparseArrays: spzeros
using JuMP: @variable, @constraint, @objective
import JuMP

const MT = AbstractMatrix{<:Real}
const VT = AbstractVector{<:Real}

"""
Holding the information of a bilevel flow problem:
* Upper-level picks arcs to tax and a discrete tax level
* Lower-level solves a mincost flow problem
"""
struct BilevelFlowProblem{M1<:MT,M2<:AbstractMatrix{Bool},M3<:MT,M4<:AbstractArray{3,<:Real}}
    init_cost::M1
    taxable_edges::M2
    capacities::M3
    tax_options::M4
    nv::Int
    na::Int
    minflow::Float64
    function BilevelFlowProblem(init_cost::M1,taxable_edges::M2,capacities::M3,tax_options::M4) where {M1<:MT,M2<:AbstractMatrix{Bool},M3<:MT,M4<:MT}
        size(init_cost) == size(taxable_edges) == size(capacities) == size(tax_options)[1:end-1] || DimensionMismatch("Matrices dimensions")
        nv = size(init_cost)[1]
        nv == size(init_cost)[2] || DimensionMismatch("Matrices should be squared")
        na = sum(1 for i in Base.OneTo(nv) for j in Base.OneTo(nv) if capacities[i,j] > 0. && i!=j)
        new{M1,M2,M3,M4}(init_cost, taxable_edges, capacities, tax_options, nv, na)
    end
end

function bilevel_flow(bfp::BilevelFlowProblem, solver)
    m = JuMP.Model(solver = solver)
    (nv, _, nopt) = size(bfp.tax_options)
    @variable(m, y[i=1:nv,j=1:nv,k=1:nopt], Bin)
    @variable(m, x[i=1:nv,j=1:nv] >= 0.)
    @variable(m, f[i=1:nv,j=1:nv] >= 0.)
    @constraint(m, limit1[i=1:nv,j=1:nv], x[i,j] <= bfp.capacities[i,j] * sum(y[i,j,k] * bfp.tax_options[i,j,k] for k in 1:nopt))
    @constraint(m, limit2[i=1:nv,j=1:nv,k=1:nopt], x[i,j] <= f[i,j] * bfp.tax_options[i,j,k])
    @constraint(m, unique_opt[i=1:nv,j=1:nv], sum(y[i,j,k] for k in 1:nopt) <= 1.)
    @objective(m, Max,
        sum(x[i,j] for i in 1:nv, j in nv if bfp.taxable_edges[i,j])
    )
    sc = sum(bfp.capacities) # upper bound on slack variable
    @variable(m,
        0. <= s[i=1:length(b)] <= ifelse(i>1 & i < bfp.nv, 0., sc) # slack at 0 if flow constraints at nodes
    )
    flat_flow = Iterators.flatten(f)
    flat_bin  = Iterators.flatten(y)
    base_cost = Iterators.flatten(bfp.init_cost)
    (B, b, F) = flow_constraint_standard(bfp)
    @constraint(m, B*f + s .== b)
    build_blp_model(m, B, base_cost, y, f, s, F)
end

function flow_constraint_standard(bfp::BilevelFlowProblem)
    B = spzeros( # constraint matrix
        1 +        # min flow from source
        bfp.nv-2 + # flow conservation at all nodes but source & sink 
        bfp.ne,    # capacity at all edges
        bfp.nv*bfp.nv
    )
    b = zeros(size(B)[1])
    b[1] = -bfp.fmin
    if bfp.capacities[1,nv] > 0.
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
    edge_idx = bfp.nv
    b[edge_idx:end] .= bfp.capacities
    for i in Base.OneTo(bfp.nv), j in Base.OneTo(bfp.nv)
        if i != j & bfp.capacities[i,j] > 0.
            B[edge_idx,bfp.nv*(j-1)+i] = 1.
            edge_idx += 1
        end
    end
    sc = sum(bfp.capacities) # upper bound on slack variable

    # Quadratic cost matrix F <=> second level objective c^T f + y^T F f
    F = spzeros(bfp.nv^2*nopt,bfp.nv^2)
    for i in Base.OneTo(bfp.nv)
        for j in Base.OneTo(bfp.nv)
            if bfp.taxable_edges[i,j]
                col = i + bfp.nv*(j-1)
                for k in Base.OneTo(nopt)
                    F[ncol + (k-1)*bfp.nv*bfp.nv,ncol] = bfp.tax_options[i,j,k]
                end
            end
        end
    end
    return (B, b, F)
end
end
