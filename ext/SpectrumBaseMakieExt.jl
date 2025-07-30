# ext/SpectrumBaseMakieExt.jl --- extension for plotting recipes

module SpectrumBaseMakieExt

using Makie
using SpectrumBase

# one dimensional views etc
function Makie.convert_arguments(P::PointBased, spec::AbstractSpectrum{1})
    rg = only(range(spec))
    ints = intensities(spec)
    convert_arguments(P, rg, ints)
end

function Makie.convert_arguments(P::PointBased, 
                spec::SpectrumView{<:AbstractSpectrum{1}})
    rg = only(range(spec))
    ints = intensities(spec)
    convert_arguments(P, rg, ints)
end

# two dimensional views etc
function Makie.convert_arguments(P::GridBased, spec::AbstractSpectrum{2})
    rA, rB = range(spec)
    convert_arguments(P, rA, rB, intensities(spec))
end

function Makie.convert_arguments(P::GridBased,
                spec::SpectrumView{<:AbstractSpectrum{2}})
    rA, rB = range(spec)
    convert_arguments(P, rA, rB, intensities(spec))
end

end # module

