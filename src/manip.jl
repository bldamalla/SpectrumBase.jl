# manip.jl --- basic spectrum manipulations

abstract type AbstractFilter end

export SGFilter
export vshift, vscale, hshift

"""
    vshift(spectrum, y)

Vertically shift the intensities of the `spectrum` by `y`. Returns only shifted intensities.
Does not construct spectra.
"""
function vshift(spec::AbstractSpectrum, y)
    yvals = intensities(spec)
    yT = eltype(yvals)

    y_ = convert(yT, y)
    return yvals .+ y_
end

"""
    vscale(spectrum, α)

Vertically scale the intensities of the `spectrum` by `α`. Returns only scaled intensities.
Does not construct spectra.
"""
function vscale(spec::AbstractSpectrum, α)
    yvals = intensities(spec)
    return yvals .* α
end

"""
    hshift(spectrum, x)

Horizontally shift the spectrum coordinates by `x`. Returns only shifted ranges.
Does not construct spectra.
"""
function hshift(spec::AbstractSpectrum, x)
    rng = range(spec)
    xT = eltype(rng)

    x_ = convert(xT, x)
    return rng .+ x_
end

using Memoization

"""
    SGFilter(windowsize, fitdegree, derivativedegree)
    (ft::SGFilter)(::AbstractVector)
    (ft::SGFilter)(xdata, ydata) --> unimplemented
    (ft::)(::Union{<:SpectrumView,<:EvenSpacedSpectrum}) --> unimplemented

Savitzky--Golay filter for noise reduction and derivative estimation.
Convolution coefficients are calculated using https://pubs.acs.org/doi/epdf/10.1021/ac00205a007
"""
struct SGFilter{wT<:Integer,fT} <: AbstractFilter
    windowsize::wT
    fitdegree::wT
    derv::wT
    coeffs::Vector{fT}
    function SGFilter{wT,fT}(w::wT, f::wT, d::wT, coeffs::Vector{fT}) where {wT<:Integer,fT}
        # set some constraints on the window size, derv degree, etc
        isodd(w) & (w > 1) || throw(ArgumentError("window size should be odd and greater than 1, got $w"))
        w - f > 1 || throw(ArgumentError("difference between window size and polynomial order 
            should be at least 2, got $(w-f)"))
        length(coeffs) == w || throw(ArgumentError("number of coefficients should match 
            window size, got $(length(coeffs)) and $w"))
            d ≥ 0 || throw(ArgumentError("derivative order should be nonnegative, got $d"))
        return new{wT,fT}(w, f, d, coeffs)
    end
end
function SGFilter(w, degree, derv)
    w_, degree_, d_ = promote(w, degree, derv); T = typeof(w_)
    return SGFilter{T,Rational{Int}}(w_, degree_, d_, _filtercoeffs(w_, degree_, d_))
end
SGFilter(w, d) = SGFilter(w, d, 0)

function (ft::SGFilter)(ydata::AbstractVector)
    wts = ft.coeffs; w = ft.windowsize; hw = _halfwindow(w)
    ylen = length(ydata); rT = promote_type(eltype(ydata), eltype(wts))
    rvec = zeros(rT, ylen)
    @inbounds for i in 1+hw:ylen-hw, j in -hw:hw
        rvec[i] += wts[j+hw+1]*ydata[i+j]
    end
    _fixendpoints!(rvec, ft, ydata)
    return rvec
end

function _fixendpoints!(rvec, ft::SGFilter, ydata::AbstractVector)
    w = ft.windowsize; hw = _halfwindow(w)
    deriv = ft.derv; order = ft.fitdegree
    @views @inbounds for i in 1:hw, j in -hw:hw
        coeff  = _precoeff( j,  i, order, hw, deriv)
        coeff2 = _precoeff(-j, -i, order, hw, deriv)
        rvec[i] += coeff*ydata[i]
        rvec[end-i+1] += coeff2*ydata[end-i+1]
    end
    return rvec
end

_halfwindow(w::Int) = div(w-1, 2)

# what follows are all integer operations, so it might be better to preserve
# the coefficients to full precision using Rational

function _filtercoeffs(w, order, deriv)
    # do something
    hwidth = _halfwindow(w)
    return map(-hwidth:hwidth) do i
        _precoeff(i, 0, order, hwidth, deriv)
    end
end

function _grampoly(i, k, m, s)
    if k > 0
        prefac1 = (4k-2) // (k*(2m-k+1))
        prefac2 = (1-k)*(2m+k) // (k*(2m-k+1))
        return prefac1 * (i*_grampoly(i,k-1,m,s) + s*_grampoly(i,k-1,m,s-1)) +
            prefac2 * _grampoly(i,k-2,m,s)
    elseif (k == 0 && s == 0) return 1
    else return 0
    end
end

function _precoeff(i, t, o, m, s)
    return sum(0:o) do k
        prefac = (2k+1) * _genfac(2m, k) // _genfac(2m+k+1,k+1)
        return prefac * _grampoly(i,k,m,0) * _grampoly(t,k,m,s)
    end
end

_genfac(a, b) = ifelse(b>0, prod(a-b+1:a),  1)

