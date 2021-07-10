module Examples
using LightGraphs
using SimpleWeightedGraphs

const all_examples = SimpleWeightedGraph[]

function wgraph(E, n)
	eweights = [e[3] for e in E]
	T = eltype(eweights)
	g = SimpleWeightedGraph{T, Int}(n)
	for e in E
		add_edge!(g, e...)
	end
	for i in 1:nv(g)
		add_edge!(g, i, i, typemax(T))
	end
	push!(all_examples, g)
	return g
end

const g1 = let n = 5
	E = [(1, 2, 12)
		 (1, 3, 10) 
		 (1, 4, 19)
		 (1, 5, 8)
		 (2, 3, 3)
		 (2, 4, 7) 
		 (2, 5, 2)
		 (3, 4, 6)
		 (3, 5, 20)
		 (4, 5, 4)
		 ]
	wgraph(E, n)
end

const g2 = let n = 5
	E = [(1, 2, 14)
		 (1, 3, 15) 
		 (1, 4, 4)
		 (1, 5, 9)
		 (2, 3, 18)
		 (2, 4, 5) 
		 (2, 5, 13)
		 (3, 4, 19)
		 (3, 5, 10)
		 (4, 5, 12)
		 ]
	wgraph(E, n)
end

const g3 = let n = 5
	E = [(1, 2, 5)
		 (1, 3, 7)
		 (1, 4, 7)
		 (1, 5, 10)
		 (2, 3, 2)
		 (2, 4, 3)
		 (2, 5, 3)
		 (3, 4, 5)
		 (3, 5, 10)
		 (4, 5, 4)
		 ]
	wgraph(E, n)
end

const g4 = let n = 6
	E = [(1, 2, 12)
		 (1, 3, 29)
		 (1, 4, 22)
		 (1, 5, 13)
		 (1, 6, 24)
		 (2, 3, 19)
		 (2, 4, 3)
		 (2, 5, 25)
		 (2, 6, 6)
		 (3, 4, 21)
		 (3, 5, 23)
		 (3, 6, 28)
		 (4, 5, 4)
		 (4, 6, 5)
		 (5, 6, 16)
		 ]
	wgraph(E, n)
end

const g5 = let n = 4
	E = [(1, 2, 12)
		 (1, 3, 14)
		 (1, 4, 17)
		 (2, 3, 15)
		 (2, 4, 18)
		 (3, 4, 29)
		 ]
	wgraph(E, n)
end

end
