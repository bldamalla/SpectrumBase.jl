# maths.jl --- mathematical operations on spectra

export LeftRiemann, RightRiemann, Midpoint
export nominalmin, nominalmax
export integrate, moment
export firstderivative

################
#
#   INTEGRATION
#
################
abstract type IntegrationScheme end

struct LeftRiemann <: IntegrationScheme end
struct RightRiemann <: IntegrationScheme end
struct Midpoint <: IntegrationScheme end

_integspeccheck(spec) = isevenspaced(spec) ? true : throw(ArgumentError("only even spaced spectra are supported for integration."))

@generated function integrate(scheme, spec::AbstractSpecOrView)
    idcs = if scheme === LeftRiemann
        :(CartesianIndices(sz_targ))
    elseif scheme === RightRiemann
        :(CartesianIndices(sz_targ) .+ CartesianIndex(1,1))
    end

    return quote
        _integspeccheck(spec)
        yvals = intensities(spec)
        sz_targ = size(yvals) .- 1
        return @views sum(yvals[$idcs]) * prod(step(spec))
    end
end
function integrate(::Midpoint, spec::AbstractSpecOrView)
    _integspeccheck(spec)
    yvals = intensities(spec)
    sz = size(spec); N = ndims(spec)
    return prod(step(spec)) * sum(CartesianIndices(yvals)) do ix
        nflags = count(zip(sz, Tuple(ix))) do (mx, i)
            1 < i < mx
        end
        scalefactor = 2^(N-nflags)
        @inbounds yvals[ix] / scalefactor
    end
end

################
#
#   SPECTRAL MOMENTS / NOMINAL EXTREMA
#
################
"""
    nominalmax(spectrum::Union{SpectrumView,AbstractSpectrum})

Return the nominal (global) maximum intensity of the spectrum and its corresponding coordinate.
"""
function nominalmax(spectrum::AbstractSpecOrView)
    xvals = range(spectrum)
    yvals = intensities(spectrum)

    maxint, idx = findmax(yvals); tidx = Tuple(idx)
    rvalues = ntuple(length(tidx)) do i
        @inbounds xvals[i][tidx[i]]
    end
    return rvalues, maxint
end

"""
    nominalmin(spectrum::Union{SpectrumView,AbstractSpectrum})

Return the nominal (global) minimum intensity of the spectrum and its corresponding coordinate.
"""
function nominalmin(spectrum::AbstractSpecOrView)
    xvals = range(spectrum)
    yvals = intensities(spectrum)

    minint, idx = findmin(yvals); tidx = Tuple(idx)
    rvalues = ntuple(length(tidx)) do i
        @inbounds xvals[i][tidx[i]]
    end
    return rvalues, minint
end

"""
    moment(scheme, spectrum, degree=1)

Calculate the `degree`th spectral moment of `spectrum`, where integrations are performed using
the integration scheme `scheme`. See `integrate` for a list of valid schemes.
"""
moment(scheme, spectrum, degree=1, center=zero(eltype(range(spectrum)))) = _unscaledmoment(scheme, spectrum, degree, center) / integrate(scheme, spectrum)

@generated function _unscaledmoment(scheme, spectrum, degree, center)
    sm_ex = _sumexpression(scheme)          # calculate sum expression according to scheme

    return quote
        _integspeccheck(spectrum)           # check if spectrum is even spaced
        pre_yvals = intensities(spectrum)   # "pre" because we still need to multiply xvals^degree
        rg = range(spectrum)
        yvals = @. pre_yvals * (rg - center) ^ degree
        Δ = step(spectrum)
        # len = length(spectrum) - 1  # one point is lost in all above schemes

        sm = $(sm_ex)
        return sm * Δ
    end
end

################
#
#   DERIVATIVES
#
################

"""
    firstderivative(spectrum)

Calculate the first derivative of the `spectrum` using the central difference method.
"""
function firstderivative(spectrum)
    _integspeccheck(spectrum)   # check if the spectrum is even spaced
    yvals = intensities(spectrum)
    Δ = step(spectrum)
    len = length(spectrum) - 2  # length of vector containing actual derivatives

    # initialize the output vector; first, last elements =0 to indicate loss of info
    ysimilar = similar(yvals)
    ysimilar[begin] = ysimilar[end] = zero(eltype(ysimilar))

    # get the derivatives y'[j] = (y[j+1] - y[j-1])/2Δ
    for j in Iterators.take(Iterators.drop(eachindex(ysimilar), 1), len)
        @inbounds ysimilar[j] = (yvals[j+1] - yvals[j-1]) / (2*Δ)
    end

    return ysimilar
end

