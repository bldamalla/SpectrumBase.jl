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
# util function for determining sum expression of intensities inside integrals
function _sumexpression(scheme)
    if scheme === LeftRiemann
        return :(sum(yvals) - last(yvals))
    elseif scheme === RightRiemann
        return :(sum(yvals) - first(yvals))
    elseif scheme === Midpoint
        return :(sum(yvals) - (first(yvals) + last(yvals)) / 2)
    else
        throw(ArgumentError("Integration scheme not internally supported."))
    end
end

"""
    integrate(R::IntegrationScheme, spectrum)

Numerically integrate `spectrum` using the integration scheme `R`.

# Examples

    # integration using LeftRiemann scheme
    integrate(LeftRiemann(), spectrum)

    # integration using RightRiemann scheme
    integrate(RightRiemann(), spectrum)

    # integration using Midpoint scheme (midpoint rule)
    integrate(Midpoint(), spectrum)
"""
@generated function integrate(scheme, spectrum)
    sm_ex = _sumexpression(scheme)      # calculate sum expression according to scheme

    return quote
        _integspeccheck(spectrum)
        yvals = intensities(spectrum)
        Δ = step(spectrum)
        len = length(spectrum) - 1  # one point is lost in all above schemes

        sm = $(sm_ex)
        return sm * Δ * len
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
function nominalmax(spectrum)
    xvals = range(spectrum)
    yvals = intensities(spectrum)

    maxint, idx = findmax(yvals)
    zippedidx = idx - firstindex(yvals) + firstindex(xvals)
    return xvals[zippedidx], maxint
end

"""
    nominalmin(spectrum::Union{SpectrumView,AbstractSpectrum})

Return the nominal (global) minimum intensity of the spectrum and its corresponding coordinate.
"""
function nominalmin(spectrum)
    xvals = range(spectrum)
    yvals = intensities(spectrum)

    minint, idx = findmin(yvals)
    zippedidx = idx - firstindex(yvals) + firstindex(xvals)
    return xvals[zippedidx], minint
end

"""
    moment(scheme, spectrum, degree=1)

Calculate the `degree`th spectral moment of `spectrum`, where integrations are performed using
the integration scheme `scheme`. See `integrate` for a list of valid schemes.
"""
moment(scheme, spectrum, degree=1) = _unscaledmoment(scheme, spectrum, degree) / integrate(scheme, spectrum)

@generated function _unscaledmoment(scheme, spectrum, degree)
    sm_ex = _sumexpression(scheme)          # calculate sum expression according to scheme

    return quote
        pre_yvals = intensities(spectrum)   # "pre" because we still need to multiply xvals^degree
        yvals = @. pre_yvals * range(spectrum) ^ degree
        Δ = step(spectrum)
        len = length(spectrum) - 1  # one point is lost in all above schemes

        sm = $(sm_ex)
        return sm * Δ * len
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
    len = length(spectrum) - 1  # length of vector containing actual derivatives

    # initialize the output vector; first, last elements =0 to indicate loss of info
    ysimilar = similar(yvals)
    ysimilar[firstindex(ysimilar)] = ysimilar[lastindex(ysimilar)] = zero(eltype(ysimilar))

    # get the derivatives y'[j] = (y[j+1] - y[j-1])/2Δ
    for j in Iterators.take(Iterators.drop(eachindex(ysimilar), 1), len)
        @inbounds ysimilar[j] = (yvals[j+1] - yvals[j-1]) / (2*Δ)
    end

    return ysimilar
end

