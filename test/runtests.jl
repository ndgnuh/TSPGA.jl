using TSPGA
using TSPGA.Examples
using Test

@testset "Solve" begin 
	for G in Examples.all_examples
		solve_tsp(G)
		@test true
	end
end

@testset "Plot" begin
	for G in Examples.all_examples
		wgplot(G)
		@test true
	end
end
