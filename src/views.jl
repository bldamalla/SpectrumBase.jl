# views.jl --- handling views of spectra

export getview, section

"""
    getview(spec, inner, method=:binsearch)

Return the smallest `SpectrumView` of a spectrum such that the range encompasses the interval indicated
by the `inner` argument. Note that the types of the interval endpoints should be promotable to
`eltype(range(spec))`.

Implemented methods for searching smallest intervals are:
+ `:binsearch` (default) - naive binary search
+ `:linear` - naive linear search
"""
function getview(spec::AbstractSpectrum{xT}, inner, method=:binsearch) where xT
    rng = range(spec)
    la, ra = first(rng), last(rng)
    lb, rb = convert.(xT, inner)

    @assert all(inner) do ix
        (la <= ix <= ra) | (la >= ix >= ra)
    end "Endpoints in second argument should be within the `range` of the first argument."

    @assert (la - ra) * (lb - rb) > zero(xT) "Ordering of the endpoints should be maintained."

    # get the sections
    targrange = envelope(rng, inner, method)
    return SpectrumView(spec, targrange)
end

"""
    section(specview::SpectrumView)

Return the `range` and `intensities` of the `SpectrumView` from `getview(spec)`. This may be useful when
generating separate copies of portions of the spectra for further manipulation.

See also: [`getview`](@ref)
"""
function section(specview::SpectrumView)
    range(specview), copy(intensities(specview))
end

function envelope(parent, inner, method)
    if method === :linear
        _envelope_linear(parent, inner)
    elseif method === :binsearch
        _envelope_binsearch(parent, inner)
    else
        throw(ArgumentError("See `?getview` for implemented methods."))
    end
end

function _envelope_linear(parent, inner)
    cpstart, cpstop = last(parent)-first(parent) > 0 ? (<, >) : (>, <)

    istart, istop = inner
    T = eltype(eachindex(parent))
    envstart = zero(T); envstop = zero(T)

    for (ix, p) in pairs(parent)
        if cpstart(p, istart)
            envstart = ix
            continue
        end
        if cpstop(p, istop)
            envstop = ix
            break
        end
    end

    return (envstart, envstop)
end

function _envelope_binsearch(parent, inner)
    istart, istop = inner
    isrev = istart > istop

    # get left end
    envstart = searchsortedlast(parent, istart, rev=isrev)
    # get right end
    envstop = searchsortedfirst(parent, istop, rev=isrev)

    return (envstart, envstop)
end

