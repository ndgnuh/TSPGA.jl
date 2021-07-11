struct Permutations{T}
	a::T
	t::Int
end

Base.eltype(::Type{Permutations{T}}) where {T} = Vector{eltype(T)}

Base.length(p::Permutations) = (0 <= p.t <= length(p.a)) ? factorial(length(p.a), length(p.a)-p.t) : 0

"""
permutations(a)

Generate all permutations of an indexable object `a` in lexicographic order. Because the number of permutations
can be very large, this function returns an iterator object.
Use `collect(permutations(a))` to get an array of all permutations.
"""
permutations(a) = Permutations(a, length(a))

"""
permutations(a, t)

Generate all size `t` permutations of an indexable object `a`.
"""
function permutations(a, t::Integer)
	if t < 0
		t = length(a) + 1
	end
	Permutations(a, t)
end

function Base.iterate(p::Permutations, s = collect(1:length(p.a)))
	(!isempty(s) && max(s[1], p.t) > length(p.a) || (isempty(s) && p.t > 0)) && return
	nextpermutation(p.a, p.t ,s)
end

function nextpermutation(m, t, state)
	perm = [m[state[i]] for i in 1:t]
	n = length(state)
	if t <= 0
		return(perm, [n+1])
	end
	s = copy(state)
	if t < n
		j = t + 1
		while j <= n &&  s[t] >= s[j]; j+=1; end
	end
	if t < n && j <= n
		s[t], s[j] = s[j], s[t]
	else
		if t < n
			reverse!(s, t+1)
		end
		i = t - 1
		while i>=1 && s[i] >= s[i+1]; i -= 1; end
		if i > 0
			j = n
			while j>i && s[i] >= s[j]; j -= 1; end
			s[i], s[j] = s[j], s[i]
			reverse!(s, i+1)
		else
			s[1] = n+1
		end
	end
	return (perm, s)
end

function bruteforce_tsp(W::AbstractMatrix{T}, s = 1) where T
	n = size(W, 1)
	non_source = filter(!isequal(s), 1:n)
	result = collect(1:n)
	minpath = typemax(T)
	for p in permutations(non_source)
		path = [s; p]
		obj = sum(W[i, j] for (i, j) in p2edge(path))
		if obj < minpath
			minpath = obj
			result .= path
		end
	end
	(result, minpath)
end
bruteforce_tsp(g::AbstractGraph, s = 1) = bruteforce_tsp(adjacency_matrix(g), s)

function nna_tsp(g::AbstractMatrix{T}, u = rand(1:size(g, 1))) where T
	n = size(g, 1)
	unvisited = ones(Bool, n)
	unvisited[u] = false
	path = [u]
	sizehint!(path, n)
	weight = zero(eltype(g))
	while any(unvisited)
		min = typemax(T)
		v = 0
		for i in 1:n
			if unvisited[i] && g[u, i] < min
				min = g[u, i]
				v = i
			end
		end
		unvisited[v] = false
		push!(path, v)
		u = v
		weight = weight + min
	end
	weight = weight + g[last(path), first(path)]
	(weight, path)
end
nna_tsp(g::AbstractGraph) = nna_tsp(adjacency_matrix(g))
