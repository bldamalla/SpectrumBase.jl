using SpectrumBase
using Test

@testset "SpectrumBase.jl" begin
    # Write your tests here.
    include("types.jl")
    include("views.jl")
    include("maths.jl")
end
