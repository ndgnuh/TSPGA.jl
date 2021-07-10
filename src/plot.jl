using GraphPlot
using LightGraphs
using SparseArrays

function wgplot(g)
	g = copy(g)
	T = eltype(g.weights)
	g.weights = map(g.weights) do w
		if w == typemax(T)
			zero(T)
		else
			w
		end
	end
	dropzeros!(g.weights)
	n = nv(g)
	edgelabels = [g.weights[i, j] for (i, j) in Iterators.product(1:n, 1:n) if i < j]
	gplot(g; nodelabel=1:n, edgelabel=edgelabels)
end
