# manip.jl --- basic spectrum manipulations

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

