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
function _lsparenttype(::LineShape) end

struct Gaussian{xT,yT} <: LineShape{xT,yT}
    μ::xT
    σ::xT
    h::yT
    function Gaussian(μ::xT, σ::xT, h::yT) where {xT,yT}
        #σ > 0 || throw(ArgumentError("scale ``σ`` should be positive"))
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
_lsparenttype(::Type{<:Gaussian}) = Gaussian

struct Cauchy{xT,yT} <: LineShape{xT,yT}
    x₀::xT
    γ::xT
    h::yT
    function Cauchy(x₀::xT, γ::xT, h::yT) where {xT,yT}
        #γ > 0 || throw(ArgumentError("scale ``γ`` should be positive"))
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
_lsparenttype(::Type{<:Cauchy}) = Cauchy

struct RaisedCosine{xT,yT} <: LineShape{xT,yT}
    μ::xT
    s::xT
    h::yT
    function RaisedCosine(μ::xT, s::xT, h::yT) where {xT,yT}
        #s > 0 || throw(ArgumentError("scale ``s`` should be positive"))
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
_lsparenttype(::Type{<:RaisedCosine}) = RaisedCosine

############
##
##  FITTING PROBLEMS
##
############

export FittingProblem, FittingSolution
export build_objective, solve

# define the fitting problem
# make a solution struct for keeping the solution and
# define functions for showing the fitting statistics of the solution
# important ones are residues, sum square error, degrees of freedom, mse (rmsd)

"""
    FittingProblem

Defines a problem for subsequent fitting using a solve function.
"""
struct FittingProblem{LS<:LineShape,iP,xT,yT}
    targ_ls::Type{LS}
    init_params::iP
    xdata::xT
    ydata::yT
    function FittingProblem(ls::Type{<:LineShape}, ipars, xdata, ydata)
        isbitstype(ls) || error("lineshape is not bits type")
        rem(length(ipars), _nparams(ls)) == 0 || error("improper params vector")
        length(xdata) == length(ydata) || DimensionMismatch("data don't have the same dimensions: got $(length(xdata)) and $(length(ydata))")
        return new{ls,typeof(ipars),
            typeof(xdata),typeof(ydata)}(ls, ipars, xdata, ydata)
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

function _evaluator(LS, params)
    # throw an error if the length of params in the compound is not divisible by the
    # number of params in a single lineshape
    M = _nparams(LS); LS_ = _lsparenttype(LS)
    return function _tp(x)
        partitioned = Iterators.partition(params, M)
        return sum(partitioned) do part
            evaluate(LS_(part...), x)
        end
    end
end

"""
    build_objective(::FittingProblem)

Return the least squares objective function to be minimized for the
given `FittingProblem`.
"""
function build_objective(prob::FittingProblem{LS}) where LS
    (; xdata, ydata) = prob
    return function _obj(params)
        # get y estimates from xdata
        # get the sum of squares of residuals
        evtor_ = _evaluator(LS, params)
        return sum(zip(xdata, ydata)) do (x, y)
            Δ = y - evtor_(x)
            abs2(Δ)
        end
    end
end
