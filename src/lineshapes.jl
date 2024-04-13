# lineshapes.jl --- basic lineshapes for testing functionality and fitting

export LineShape, Gaussian, Cauchy, RaisedCosine
export evaluate, lscenter, scale

abstract type LineShape{xT,yT} end

"""
    evaluate(ls<:LineShape, x)

Evaluate the line shape `ls` at the point `x`. Supported line shapes as of writing are:
+ `Gaussian`
+ `Lorentzian`
+ `RaisedCosine`


"""
function evaluate end

function lscenter end
function scale end

_nparams(::T) where T<:LineShape = _nparams(T)

struct Gaussian{xT,yT} <: LineShape{xT,yT}
    μ::xT
    σ::xT
    h::yT
    function Gaussian(μ::xT, σ::xT, h::yT) where {xT,yT}
        σ > 0 || throw(ArgumentError("scale ``σ`` should be positive"))
        return new{xT,yT}(μ,σ,h)
    end
end
lscenter(gs::Gaussian) = gs.μ
scale(gs::Gaussian) = gs.σ
_nparams(::Type{<:Gaussian}) = 3
function evaluate(gs::Gaussian, x)
    scaled = (x - gs.μ) / gs.σ
    expd = exp(scaled^2 * -1//2)
    return expd * gs.h / gs.σ
end

struct Cauchy{xT,yT} <: LineShape{xT,yT}
    x₀::xT
    γ::xT
    h::yT
    function Cauchy(x₀::xT, γ::xT, h::yT) where {xT,yT}
        γ > 0 || throw(ArgumentError("scale ``γ`` should be positive"))
        return new{xT,yT}(x₀,γ,h)
    end
end
lscenter(lt::Cauchy) = lt.x₀
scale(lt::Cauchy) = lt.γ
_nparams(::Type{<:Cauchy}) = 3
function evaluate(lt::Cauchy, x)
    scaled = (x - lt.x₀) / lt.γ
    dnm = lt.γ * (1 + scaled^2)
    return inv(dnm) * lt.h
end

struct RaisedCosine{xT,yT} <: LineShape{xT,yT}
    μ::xT
    s::xT
    h::yT
    function RaisedCosine(μ::xT, s::xT, h::yT) where {xT,yT}
        s > 0 || throw(ArgumentError("scale ``s`` should be positive"))
        return new{xT,yT}(μ,s,h)
    end
end
lscenter(rc::RaisedCosine) = rc.μ
scale(rc::RaisedCosine) = rc.s
_nparams(::Type{<:RaisedCosine}) = 3
function evaluate(rc::RaisedCosine, x)
    scaled = (x - rc.μ) / rc.s
    hvc = 1 + cospi(scaled)
    return hvc * rc.h
end
