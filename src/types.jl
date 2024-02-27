# types.jl --- types for spectra keeping and manip

export AbstractSpectrum, EvenSpacedSpectrum, SpectrumView
export endpoints, intensities

"""
    abstract type AbstractSpectrum

Umbrella type that encompasses spectrum types. The package focuses on one-dimensional spectra, _i.e._
those that can be represented in two dimensions. Certain examples of primary applications are absorption/fluorescence
spectra and one-dimensional nuclear magnetic resonance (NMR) spectra.

# Implementation

Derived types should implement the following methods:
+ `Base.range` - values for which the spectrum is obtained
+ `intensities` - values for each ``x`` value
"""
abstract type AbstractSpectrum{N,xT<:Number,yT<:Number} end
Base.ndims(::Type{<:AbstractSpectrum{N}}) where N = N
Base.ndims(obj::AbstractSpectrum) = ndims(typeof(obj))

"""
    abstract type EvenSpacedSpectrum <: AbstractSpectrum

Spectra that are recorded at evenly spaced intervals, _e.g._ chemical shifts or wavelengths. Have definite
`step`s between points along the spectrum.

# Implementation

Derived types should, in addition to those for `SpectrumBase.intensities`, implement the following methods:
+ `endpoints` (instead of `Base.range`) - endpoints for ``x`` values used to define the spectrum

*Note*: `Base.range` is calculated using `endpoints` for this spectrum type.
"""
abstract type EvenSpacedSpectrum{xT,yT} <: AbstractSpectrum{1,xT,yT} end
const ESS = EvenSpacedSpectrum

abstract type GridSpectrum{xT,yT} <: AbstractSpectrum{2,xT,yT} end
const GS = GridSpectrum

# from here below, definitions are mostly based on properties of evenly spaced spectra

function endpoints end
function intensities end

"""
    range(ev::EvenSpacedSpectrum)

Return the range for which the spectrum is obtained.
"""
Base.range(ev::ESS) = (range(endpoints(ev)...; length=length(ev)),)
Base.size(ev::AbstractSpectrum) = size(intensities(ev))

function Base.range(gr::GS)
    eA, eB = endpoints(gr)
    sA, sB = size(gr)
    return (range(eA...; length=sA), range(eB...; length=sB))
end
Base.step(gr::Union{ESS,GS}) = map(step, range(gr))

#################
#
#   ViewFrame
#
#################

"""
    ViewFrame{iT}

A two-field container for the start and stop indices for a one-dimensional view.
This represents an iterator for a unit range (of `Integer`s). Identical to a 
`UnitRange`. Only for "type hygiene" within the project.
"""
struct ViewFrame{iT<:Integer}
    start::iT
    stop::iT
    function ViewFrame(sta::iT, sto::iT) where iT
        @assert sta < sto
        return new{iT}(sta, sto)
    end
end
Base.first(vf::ViewFrame) = vf.start
Base.last(vf::ViewFrame) = vf.stop
Base.length(vf::ViewFrame) = vf.stop - vf.start + 1
asrange(vf::ViewFrame) = range(first(vf), last(vf))
function Base.iterate(vf::ViewFrame, state=first(vf))
    state > last(vf) && return nothing
    return (state, state+1)
end
function viewshift(vf::ViewFrame{iT}, i::iT) where iT
    sta, sto = first(vf), last(vf)
    return ViewFrame(sta+i, sto+i)
end
function ViewFrame(tpl::Tuple{iT,jT}) where {iT,jT}
    cT = promote_type(iT,jT)
    a, b = convert.(cT, tpl)
    return ViewFrame(a, b)
end
function ViewFrame(len::T) where T
    return ViewFrame(one(T), len)
end

#################
#
#   SpectrumView
#
#################

"""
    SpectrumView

A lazy view of an `AbstractSpectrum`. Stores the `parent` integer values indicating
indices where the slices will be obtained.
"""
struct SpectrumView{sT<:AbstractSpectrum,iT}
    parent::sT
    viewframe::NTuple{N,ViewFrame{iT}} where N
    function SpectrumView(p::AbstractSpectrum{N}, vrange::NTuple{N}) where N
        st = typeof(p)
        iT = eltype(vrange[begin])
        return new{st,iT}(p, vrange)
    end
end
Base.ndims(::SpectrumView{T}) where T = ndims(T)

function Base.range(esv::SpectrumView)
    return ntuple(ndims(esv)) do i
        vrange = @inbounds esv.viewframe[i] |> asrange
        return @inbounds rgs[i][vrange]
    end
end
Base.range(esv::SpectrumView, dim::Integer) = @inbounds range(esv)[dim]
function Base.step(esv::SpectrumView)
    @assert isevenspaced(esv) "parent spectrum should be even spaced"
    return step.(range(esv))
end
# intensities(esv::SpectrumView) = intensities(esv.parent)[vrange(esv)]
intensities(esv::SpectrumView) = view(intensities(esv.parent), range(esv)...)

isevenspaced(::Type{<:EvenSpacedSpectrum}) = true
isevenspaced(::Type{<:GridSpectrum}) = true
isevenspaced(::Type{<:AbstractSpectrum}) = false
isevenspaced(::Type{SpectrumView{sT}}) where sT = isevenspaced(sT)
isevenspaced(spec) = isevenspaced(typeof(spec))

################
#
#   Plotting recipe
#
###############
import MakieCore

function MakieCore.convert_arguments(P::MakieCore.PointBased, spec::AbstractSpectrum{1})
    MakieCore.convert_arguments(P, range(spec), intensities(spec))
end

function MakieCore.convert_arguments(P::MakieCore.GridBased, spec::GridSpectrum)
    rA, rB = range(spec)
    MakieCore.convert_arguments(P, rA, rB, intensities(spec))
end

