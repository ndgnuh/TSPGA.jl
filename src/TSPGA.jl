module TSPGA

using Random
using Statistics
using Evolutionary
using LightGraphs
using GraphPlot
include("Example.jl")
using .Examples
include("plot.jl")
include("other.jl")


function p2edge(p)
    #= p′ = [p[2:end]; first(p)] =#
    n = length(p)
    p′ = Iterators.map(i -> i % n + 1, p)
    zip(p, p′)
end


function firstequal(x, X)
    for i in eachindex(X)
        if X[i] == x
            return i
        end
    end
end

function TSPObject(W::AbstractMatrix{T}) where {T}
    function (x)
        (W[last(x), first(x)] + sum(W[i, j] for (i, j) in zip(x, Iterators.drop(x, 1))))
    end
end
TSPObject(G::AbstractGraph) = TSPObject(adjacency_matrix(G))

function first_few(k, W, p, n)
    sel = first(sortperm([W[i, j] for (i, j) in p2edge(p)]), k)
    append!(sel, sel .% n .+ 1)
    unique!(sel)
    p[sel]
end

function full_crossover(W::AbstractMatrix{S}, T = Int) where {S}
    n::T = size(W, 1)
    m::T = max(n ÷ 4, one(T))
    function (p1, p2)
        k1, k2 = rand(1:m, 2)
        c1 = first_few(k1, W, p1, n)
        c2 = first_few(k2, W, p2, n)
        append!(c1, setdiff(p2, c1))
        append!(c2, setdiff(p1, c2))
        [c1, c2]
    end
end

function TSPCrossover(W::AbstractMatrix{S}) where {S}
    n::Int = size(W, 1)
    cycleind(i) = i % n + 1
    #= m = max(n÷2, 1) =#
    m = n - 1
    function breed(p1::AbstractVector{T}, p2::AbstractVector{T}) where {T<:Integer}
        w = [W[i, j] for (i, j) in p2edge(p1)]
        c = sort(p1, by=i -> w[i])[1:rand(1:m)]
        append!(c, setdiff(p2, c))
        [c, p1[p2]]
    end
end
TSPCrossover(G::AbstractGraph) = TSPCrossover(adjacency_matrix(G))

function mutate(p)
    i, j = rand(eachindex(p), 2)
    p[i], p[j] = p[j], p[i]
    return p
end

function initialize_objective(f, x::AbstractArray{T}, S = typeof(f(x))) where {T}
    NonDifferentiable{S,typeof(x)}(f, zero(S), zeros(T, length(x)), [0])
end

function solve_tsp(
    W::AbstractMatrix{T};
    sortresult = false,
    ga_kwargs = Dict(),
    options_kwargs = Dict(),
) where {T}
    n = size(W, 1)

    options = Evolutionary.Options(;
        successive_f_tol = 30,
        store_trace = true,
        options_kwargs...,
    )

    popsize = get(ga_kwargs, :populationSize, min(n, 100))
    individual = [randperm(n) for _ = 1:popsize]

    obj = initialize_objective(TSPObject(W), first(individual))

    alg = GA(;
        crossover = TSPCrossover(W),
        selection = rouletteinv,
        mutation = mutate,
        mutationRate = 0.5,
        crossoverRate = 0.5,
        epsilon = 0.2,
        populationSize = popsize,
        ga_kwargs...,
    )

    result =
        Evolutionary.optimize(obj, Evolutionary.NoConstraints(), individual, alg, options)
    if result.converged && sortresult
        let r = result.minimizer
            _, i = findmin(r)
            if i == 1
                r
            elseif i == lastindex(r)
                reverse(r)
            end
            result.minimizer = [r[i:end]; r[1:(i - 1)]]
        end
    end
    return result
end
solve_tsp(G::AbstractGraph; kwargs...) = solve_tsp(adjacency_matrix(G); kwargs...)

function bsolve_tsp(G, n; kwargs...)
    results = map(1:n) do
        solve_Tsp(G; kwargs)
    end
    findmin(x -> x.minimum, results)[2]
end

export Evolutionary, TSPCrossover, TSPObject, solve_tsp, mutate, wgplot

end # module
