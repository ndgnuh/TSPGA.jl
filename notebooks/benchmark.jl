### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 2966421a-fd83-4aef-8197-1771eb01d8cc
begin
	using CairoMakie
	using BenchmarkTools
	using Pkg
	Pkg.activate("..")
	using TSPGA
	using SimpleWeightedGraphs
	using LightGraphs
end

# ╔═╡ c3539b19-e69d-46c6-8144-7aeceda26139
using Distributions

# ╔═╡ c31600f1-436f-4f9d-9f10-c26ff45f77ef
using LinearAlgebra

# ╔═╡ 5f531589-280c-480e-8ebe-df5290328d2a


# ╔═╡ 5d7329f2-5600-4fc5-83f5-7e03300027fb
const d = Normal(Float16(0), Float16(1))

# ╔═╡ 5671a4a2-e1bb-4122-8bd1-4a4dc71d9fd1
const N = range(25, step=25, stop=1000)

# ╔═╡ 9355ce0b-c94f-4500-9576-841da3732d44
function run_bench_mark(N)
	benches = map(reverse(N)) do n
		@info "Start $n"
		println()
		g = TSPGA.Examples.generate_random_graph(n, Float16(0.1), Float16(5), Float16(0.01))
		GC.gc()
		b = @benchmark $solve_tsp($g) samples=2 evals=2 seconds=Inf gctrial=true gcsample=true
		@info "Done $n"
		println()
		display(b)
		b
	end
	reverse(benches)
end

# ╔═╡ 08edbc9b-ba81-4f67-984f-3cb98f91f0b5
function crf(W::AbstractMatrix{S}) where S
	function cr1(p1::AbstractVector{T}, p2::AbstractVector{T}) where T
    	n = size(W, 1)
		c1 = p2edge(p1)
		c2 = p2edge(p2)
		w1 = [@inbounds W[i, j] for (i, j) in c1]
		w2 = [@inbounds W[i, j] for (i, j) in c2]
		k1 = rand(1:(n÷4))
		k2 = rand(1:(n÷4))
		perm1 = sortperm(w1)[1:k1]
		perm2 = sortperm(w2)[1:k2]
		cities1 = T[]
		cities2 = T[]
		sizehint!(cities1, n)
		sizehint!(cities2, n)
		for i in perm1
			p = p1[i]
			push!(cities1, p)
			push!(cities1, p == n ? 1 : (p + 1))
		end
		for i in perm2
			p = p2[i]
			push!(cities2, p)
			push!(cities2, p == n ? 1 : (p + 1))
		end
		for i in p1
			if !(i in cities2)
				push!(cities2, i)
			end
		end
		for i in p2
			if !(i in cities1)
				push!(cities1, i)
			end
		end
		[cities1, cities2]
	end
end

# ╔═╡ 2eae3f53-d92e-42b4-8db9-e645a416bc8a
benches = run_bench_mark(N)

# ╔═╡ 2549c63a-bae1-49ee-be54-3acd2f132fbf
let g = TSPGA.Examples.generate_random_graph(601, 1, 20)
	GC.gc()
	@benchmark solve_tsp($g) evals=5 samples=5 seconds=Inf gctrial=true gcsample=true
end

# ╔═╡ 1800f400-a2c2-4856-9dc2-fe0f61ab59e6
N[12]

# ╔═╡ 98b2702b-bfff-400e-a463-0287006da4d2
memory(benches[end]) / 1024 / 1024

# ╔═╡ f81055d9-1689-4fe6-a8ce-5c2e50d736d1
scatterlines(sizeof(Int) * N.^2 / 1024^2)

# ╔═╡ 52b58522-02ee-4821-84e9-a4a974929949
times = time.(benches) / 10^9

# ╔═╡ c5afe161-4c71-4363-8edb-f88a6d5a7b83
memo = (memory.(benches) .+ sizeof(Int) * N.^2) ./ 1024^2

# ╔═╡ 5934ac7d-de44-44d2-8223-51e490b86d96
scatterlines(N, times)

# ╔═╡ 5c2e5b65-dd93-4fb4-8287-b224557337bc
scatterlines(N, memo)

# ╔═╡ 4a9b2028-8fb0-43fd-bf2e-908761ab77da
begin
	function quick_poly_ls(N, X, n)
		A = reduce(hcat, [reshape(N.^i, :, 1) for i in 0:n])
		pinv(A) * X
	end
	function quick_poly_ls(N, X)
		results = map(1:4) do n
			b = quick_poly_ls(N, X, n)
			p(x) = sum(x^(i-1) * b[i] for i in eachindex(b))
			err = sqrt(sum((p(i) - X[i])^2 for i in eachindex(N)))
			(b, err)
		end
		findmin(x -> x[2], results)[2]
	end
end

# ╔═╡ d7025a21-bccf-430a-8ba5-97829dccd6c8
quick_poly_ls(N, memo)

# ╔═╡ d215faad-e1a3-4b52-be5e-6512107c8a14
quick_poly_ls(N, times)

# ╔═╡ 5963c89c-585f-4d9a-b2a8-608032f927a3
let n = length(memo)
	A = [ones(n) N N.^2 N.^3]
	pinv(A) * memo
end

# ╔═╡ Cell order:
# ╠═2966421a-fd83-4aef-8197-1771eb01d8cc
# ╠═c3539b19-e69d-46c6-8144-7aeceda26139
# ╠═5f531589-280c-480e-8ebe-df5290328d2a
# ╠═5d7329f2-5600-4fc5-83f5-7e03300027fb
# ╠═5671a4a2-e1bb-4122-8bd1-4a4dc71d9fd1
# ╠═9355ce0b-c94f-4500-9576-841da3732d44
# ╠═08edbc9b-ba81-4f67-984f-3cb98f91f0b5
# ╠═2eae3f53-d92e-42b4-8db9-e645a416bc8a
# ╠═2549c63a-bae1-49ee-be54-3acd2f132fbf
# ╠═1800f400-a2c2-4856-9dc2-fe0f61ab59e6
# ╠═98b2702b-bfff-400e-a463-0287006da4d2
# ╠═f81055d9-1689-4fe6-a8ce-5c2e50d736d1
# ╠═52b58522-02ee-4821-84e9-a4a974929949
# ╠═c5afe161-4c71-4363-8edb-f88a6d5a7b83
# ╠═5934ac7d-de44-44d2-8223-51e490b86d96
# ╠═5c2e5b65-dd93-4fb4-8287-b224557337bc
# ╠═c31600f1-436f-4f9d-9f10-c26ff45f77ef
# ╠═4a9b2028-8fb0-43fd-bf2e-908761ab77da
# ╠═d7025a21-bccf-430a-8ba5-97829dccd6c8
# ╠═d215faad-e1a3-4b52-be5e-6512107c8a14
# ╠═5963c89c-585f-4d9a-b2a8-608032f927a3
