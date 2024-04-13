# Optimization extension for SpectrumBase

module SpectrumBaseOptimExt

using SpectrumBase
using StructArrays, Optim

export FittingProblem
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
    init_params::iP
    xdata::xT
    ydata::yT
end
function _evaluator(lineshape, params)
    npars = SpectrumBase._nparams(lineshape)
    pars_reshaped = reshape(params, npars, :)
    structvec = StructArray{lineshape}(pars_reshaped; dims=1)
    return function _tp(x)
        return sum(structvec) do shape
            evaluate(x, shape)
        end
    end
end

function build_objective(prob::FittingProblem{LS}) where LS
    (; xdata, ydata) = prob
    return function _obj(params)
        # get y estimates from xdata
        # get the sum of squares of residuals
        evtor_ = _evaluator(LS, params)
        return sum(zip(xdata, ydata)) do (x, y)
            Δ = y - evtor_(x, params)
            abs2(Δ)
        end
    end
end

"""
    FittingSolution

Represents a solution to a `FittingProblem`. Points to the original 
`FittingProblem` solved along with the solution object from `Optim`.
"""
struct FittingSolution{fP<:FittingProblem,sT}
    problem::fP
    optimized::sT
end

# define the solve function
"""
    solve(prob::FittingProblem)

Solve a `FittingProblem` prob and return its optimized solution.
"""
function solve(prob::FittingProblem)
    initparams = prob.init_params
    objfunc = build_objective(prob)
    soln = optimize(objfunc, initparams, LBFGS(); autodiff=:forward)
    return FittingSolution(prob, soln)
end

end
