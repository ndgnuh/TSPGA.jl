using TSPGA
using TSPGA.Examples
using BenchmarkTools
using Test
using Random

@testset "Solve" begin 
	for G in Examples.examples
		solve_tsp(G)
		@test true
	end
	for tsp in Examples.tsp100_examples
		solve_tsp(tsp.W)
		@test true
	end
end

@testset "Plot" begin
	for G in Examples.examples
		wgplot(G)
		@test true
	end
end

@testset "Other solver" begin
	r1 = TSPGA.nna_tsp(Examples.g5)
	@test true
	display(r1)
	r2 = TSPGA.bruteforce_tsp(Examples.g4)
	@test true
	display(r2)
end

GC.gc()
Random.seed!(0)
let n = 1000,
	g = repeat(Examples.tsp100_examples[1].W, 10, 10)
	#= g = Examples.generate_random_graph(n, 0.1, 3.0, 0.1) =#
	@info size(g)
	bench = (@benchmark solve_tsp($g) samples=2 evals=2 seconds=Inf gctrial=true)
	display(bench)
end

@info "Solve sample"
let r = solve_tsp(Examples.g4)
	display(r)
	display(r.trace)
	display(r.minimizer)
end
let rs = map(Examples.tsp100_examples) do tsp
		r = solve_tsp(tsp.W, ga_kwargs = (populationSize = 100,))
		r.minimum, tsp.F
	end
	display(rs)
end


