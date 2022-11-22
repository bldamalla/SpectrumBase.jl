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
abstract type AbstractSpectrum{xT<:Number,yT<:Number} end

"""
    abstract type EvenSpacedSpectrum <: AbstractSpectrum

Spectra that are recorded at evenly spaced intervals, _e.g._ chemical shifts or wavelengths. Have definite
`step`s between points along the spectrum.

# Implementation

Derived types should, in addition to those for `SpectrumBase.intensities`, implement the following methods:
+ `endpoints` (instead of `Base.range`) - endpoints for ``x`` values used to define the spectrum

*Note*: `Base.range` is calculated using `endpoints` for this spectrum type.
"""
abstract type EvenSpacedSpectrum{xT,yT} <: AbstractSpectrum{xT,yT} end
const ESS = EvenSpacedSpectrum

# from here below, definitions are mostly based on properties of evenly spaced spectra

function endpoints end
function intensities end

"""
    range(ev::EvenSpacedSpectrum)

Return the range for which the spectrum is obtained.
"""
Base.range(ev::ESS) = range(endpoints(ev)...; length=length(ev))
Base.step(ev::ESS) = step(range(ev))
Base.length(ev::AbstractSpectrum) = length(intensities(ev))

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
struct SpectrumView{sT<:AbstractSpectrum}
    parent::sT
    viewrange::Tuple{Int,Int}
    function SpectrumView(parent, viewrange)
        st, sto = viewrange
        @assert st ∈ eachindex(intensities(parent)) && sto ∈ eachindex(intensities(parent))
        return new{typeof(parent)}(parent, viewrange)
    end
end

vrange(esv) = esv.viewrange[1]:esv.viewrange[2]
Base.range(esv::SpectrumView) = @inbounds range(esv.parent)[vrange(esv)]
Base.step(esv::SpectrumView{T}) where {T<:ESS} = step(esv.parent)
Base.length(esv::SpectrumView) = esv.viewrange[2] - esv.viewrange[1] + 1
intensities(esv::SpectrumView) = intensities(esv.parent)[vrange(esv)]

isevenspaced(::Type{<:EvenSpacedSpectrum}) = true
isevenspaced(::Type{<:AbstractSpectrum}) = false
isevenspaced(::Type{SpectrumView{sT}}) where sT = isevenspaced(sT)
isevenspaced(spec) = isevenspaced(typeof(spec))

################
#
#   Plotting recipe
#
###############
import MakieCore

function MakieCore.convert_arguments(P::MakieCore.PointBased, spec)
    MakieCore.convert_arguments(P, range(spec), intensities(spec))
end

