module TSPGA

using Random
using Evolutionary
using LightGraphs
using GraphPlot
include("Example.jl"); using .Examples
include("plot.jl")

struct TSPObject{T} <: Function
	W::Matrix{T}
	TSPObject(W::AbstractMatrix{T}) where T = new{T}(W)
	TSPObject(G::AbstractGraph) = TSPObject(adjacency_matrix(G))
end
struct TSPCrossover{T} <: Function
	W::Matrix{T}
	TSPCrossover(W::AbstractMatrix{T}) where T = new{T}(W)
	TSPCrossover(G::AbstractGraph) = TSPCrossover(adjacency_matrix(G))
end

function p2edge(p)
	c = Iterators.take(Iterators.cycle(p), length(p) + 1)
	zip(c, Iterators.drop(c, 1))
end

function (f::TSPObject)(p)
	sum(f.W[I[1], I[2]] for I in p2edge(p))
end

function (f::TSPCrossover)(p1, p2)
	n = size(f.W, 1)
	c1 = [[i, j] for (i, j) in p2edge(p1)]
	c2 = [[i, j] for (i, j) in p2edge(p2)]
	w1 = [f.W[i, j] for (i, j) in c1]
	w2 = [f.W[i, j] for (i, j) in c2]
	perm1 = sortperm(w1)[1:rand(1:n÷4)]
	perm2 = sortperm(w2)[1:rand(1:n÷4)]
	cities1 = unique(reduce(vcat, c1[perm1]))
	cities2 = unique(reduce(vcat, c2[perm2]))
	child1 = [cities1; setdiff(p2, cities1)]
	child2 = [cities2; setdiff(p1, cities2)]
	[child1, child2]
end

function mutate(p)
	p = copy(p)
	i, j = rand(p, 2)
	p[i], p[j] = p[j], p[i]
	return p
end

function initialize_objective(f, x::AbstractArray{T}, S = typeof(f(x))) where T
	NonDifferentiable{S,typeof(x)}(f, zero(S), zeros(T, length(x)), [0,])
end

function solve_tsp(W::AbstractMatrix{T}; sortresult = false, ga_kwargs = Dict(), options_kwargs = Dict()) where T
	n = size(W ,1)

	options = Evolutionary.Options(;successive_f_tol = 30,
								   iterations = 2000,
								    store_trace = true,
								    options_kwargs...
								   )

	alg = GA(;crossover = TSPCrossover(W),
	   selection = rouletteinv,
	   mutation = mutate,
	   mutationRate = 0.5,
	   crossoverRate = 0.5,
	   epsilon = 0.2,
	   populationSize = n,
	   ga_kwargs...)

	individual = [randperm(n) for _ in 1:alg.populationSize]

	obj = initialize_objective(TSPObject(W), first(individual))
	result = Evolutionary.optimize(obj, Evolutionary.NoConstraints(), individual, alg, options)
	if result.converged && sortresult
		let r = result.minimizer
			_, i = findmin(r)
			if i == 1
				r
			elseif i == lastindex(r)
				reverse(r)
			end
			result.minimizer = [r[i:end]; r[1:i-1]]
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
