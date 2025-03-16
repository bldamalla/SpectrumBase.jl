# src/ext/SpectrumBaseIntegralsExt.jl === extension for numerical integration
# the original implementations for integrations will still remain

module SpectrumBaseIntegralsExt

using Integrals
using SpectrumBase

# sampled integrals
function Integrals.SampledIntegralProblem(spec, f=identity)
    ndims(spec) == 1 || throw(ArgumentError("currently only one-dimensional arrays are supported"))
    xs = range(spec)
    ys = intensities(spec) .|> f
    return SampledIntegralProblem(ys, xs)
end

# lineshape composite shape integrals
function Integrals.IntegralProblem(ls::Union{LineShape,CompositeShape}, domain, f=identity; kwargs...)
    evalfunc(x, _) = evaluate(ls, x) |> f
    return IntegralProblem(evalfunc, domain; kwargs...)
end

end # module

