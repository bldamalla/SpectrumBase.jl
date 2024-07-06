# lineshapes.jl --- basic lineshapes for testing functionality and fitting

export LineShape, Gaussian, Cauchy, RaisedCosine
export CompositeShape
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
_nparams(T::Type{<:LineShape}) = fieldcount(T)
_outtype(::T) where T<:LineShape = _outtype(T)
function _lsparenttype(::LineShape) end

struct Gaussian{xT,yT} <: LineShape{xT,yT}
    μ::xT
    σ::xT
    h::yT
    function Gaussian{xT,yT}(μ::xT, σ::xT, h::yT) where {xT,yT}
        σ > 0 || throw(ArgumentError("scale ``σ`` should be positive"))
        return new{xT,yT}(μ,σ,h)
    end
end
function Gaussian(μ::Number, σ::Number, h::hT) where hT <: Number
    μ_, σ_ = promote(μ, σ); cT = typeof(μ_)
    return Gaussian{cT,hT}(μ_, σ_, h)
end
lscenter(gs::Gaussian) = gs.μ
scale(gs::Gaussian) = gs.σ
function evaluate(gs::Gaussian, x)
    scaled = (x - gs.μ) / gs.σ
    expd = exp(scaled^2 * -0.5)
    return expd * gs.h / gs.σ
end
_outtype(::Type{Gaussian{xT,yT}}) where {xT,yT} = yT

struct Cauchy{xT,yT} <: LineShape{xT,yT}
    x₀::xT
    γ::xT
    h::yT
    function Cauchy{xT,yT}(x₀::xT, γ::xT, h::yT) where {xT,yT}
        γ > 0 || throw(ArgumentError("scale ``γ`` should be positive"))
        return new{xT,yT}(x₀,γ,h)
    end
end
function Cauchy(x₀::Number, γ::Number, h::hT) where hT<:Number
    x₀_, γ_ = promote(x₀, γ); cT = typeof(γ_)
    return Cauchy{cT,hT}(x₀_, γ_, h)
end
lscenter(lt::Cauchy) = lt.x₀
scale(lt::Cauchy) = lt.γ
function evaluate(lt::Cauchy, x)
    scaled = (x - lt.x₀) / lt.γ
    dnm = lt.γ * (1 + scaled^2)
    return inv(dnm) * lt.h
end
_outtype(::Type{Cauchy{xT,yT}}) where {xT,yT} = yT

struct RaisedCosine{xT,yT} <: LineShape{xT,yT}
    μ::xT
    s::xT
    h::yT
    function RaisedCosine{xT,yT}(μ::xT, s::xT, h::yT) where {xT,yT}
        s > 0 || throw(ArgumentError("scale ``s`` should be positive"))
        return new{xT,yT}(μ,s,h)
    end
end
function RaisedCosine(μ::Number, s::Number, h::hT) where hT<:Number
    μ_, s_ = promote(μ, s); cT = typeof(s_)
    return RaisedCosine{cT,hT}(μ_, s_, h)
end
lscenter(rc::RaisedCosine) = rc.μ
scale(rc::RaisedCosine) = rc.s
function evaluate(rc::RaisedCosine, x)
    scaled = (x - rc.μ) / rc.s
    hvc = 1 + cospi(scaled)
    return hvc * rc.h
end
_outtype(::Type{RaisedCosine{xT,yT}}) where {xT,yT} = yT

using StructArrays
"""
    CompositeShape{base<:LineShape}
"""
struct CompositeShape{T<:LineShape,A}
    tparams::A
end
function CompositeShape(T::Type{<:LineShape}, N::Int)
    m = _nparams(T); filler = zeros(m, N); A = typeof(filler)
    return CompositeShape{T,A}(filler)
end
function CompositeShape(T::Type{<:LineShape}, tparams::AbstractVecOrMat)
    return CompositeShape{T,typeof(tparams)}(tparams)
end
function evaluate(cs::CompositeShape{T}, x::Number) where T
    # create a structvec
    strvec = StructVector{T}(cs.tparams, dims=1)
    xeval(str) = evaluate(str, x)
    return sum(xeval.(strvec))
end

############
##
##  FITTING PROBLEMS
##
############

export FittingProblem, FittingSolution
export build_objective

# define the fitting problem
# make a solution struct for keeping the solution and
# define functions for showing the fitting statistics of the solution
# important ones are residues, sum square error, degrees of freedom, mse (rmsd)

"""
    FittingProblem

Defines a problem for subsequent fitting using a solve function.
"""
struct FittingProblem{LS<:LineShape,iP,xT,yT}
    init_params::iP
    xdata::xT
    ydata::yT
    function FittingProblem(ls::Type{<:LineShape}, ipars, xdata, ydata)
        isbitstype(ls) || error("lineshape is not bits type")
        rem(length(ipars), _nparams(ls)) == 0 || error("improper params vector")
        length(xdata) == length(ydata) || DimensionMismatch("data don't have the same dimensions: got $(length(xdata)) and $(length(ydata))")
        return new{ls,typeof(ipars),
            typeof(xdata),typeof(ydata)}(ipars, xdata, ydata)
    end
end

## TODO: define the gradient for least squares objective

"""
    FittingSolution

Represents a solution to a `FittingProblem`. Points to the original 
`FittingProblem` solved along with the solution object from `Optim`.
"""
struct FittingSolution{fP<:FittingProblem,sT}
    problem::fP
    optimized::sT
end

"""
    solve(::FittingProblem)

Solve a `FittingProblem` and return a `FittingSolution`
"""
function solve end

"""
    build_objective(::FittingProblem)

Return the least squares objective function to be minimized for the
given `FittingProblem`.
"""
function build_objective(prob::FittingProblem{LS}) where LS
    (; xdata, ydata) = prob
    return function _obj(params)
        # get y estimates from xdata
        # get the mean of squares of residuals
        cs = CompositeShape(LS, reshape(params, _nparams(LS), :))
        l = length(ydata)
        cseval(x) = evaluate(cs, x)
        sds = @. abs2(ydata - cseval(xdata))
        return sum(sds) / l
    end
end
