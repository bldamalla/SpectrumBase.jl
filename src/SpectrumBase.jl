module SpectrumBase

# Write your package code here.
include("types.jl")
include("lineshapes.jl")
include("views.jl")
include("manip.jl")
include("maths.jl")

# include the extension
include(joinpath(@__DIR__, "../ext/SpectrumBaseOptimExt/SpectrumBaseOptimExt.jl"))

end
