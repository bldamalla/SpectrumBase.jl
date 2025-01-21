# src/ext/SpectrumBaseIntegralsExt.jl === extension for numerical integration
# the original implementations for integrations will still remain

module SpectrumBaseIntegralsExt

using Integrals
using SpectrumBase

function Integrals.SampledIntegralProblem(spec)
    ndims(spec) == 1 || throw(ArgumentError("currently only one-dimensional arrays are supported"))
    xs = range(spec)
    ys = intensities(spec)
    return SampledIntegralProblem(ys, xs)
end

end # module

