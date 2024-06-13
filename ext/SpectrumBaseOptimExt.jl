# Optimization extension for SpectrumBase

module SpectrumBaseOptimExt

using SpectrumBase, Optim

# define the solve function
"""
    solve(prob::FittingProblem)

Solve a `FittingProblem` prob and return its optimized solution.
"""
function SpectrumBase.solve(prob::FittingProblem; autodiff=:forward)
    initparams = prob.init_params
    objfunc = build_objective(prob)
    soln = optimize(objfunc, initparams, LBFGS(), autodiff=autodiff)
    return FittingSolution(prob, soln)
end
Optim.optimize(prob::FittingProblem) = SpectrumBase.solve(prob)

# TODO: extend optim functions to work on types defined here
# especially fitting solution
# define for common Optim methods
# summary, minimizer, minimum, iterations, iteration_limit_reached,
# trace, x_trace, f_trace, f_calls, converged
Optim.summary(sol::FittingSolution) = Optim.summary(sol.optimized)
Optim.minimizer(sol::FittingSolution) = Optim.minimizer(sol.optimized)
Optim.minimum(sol::FittingSolution) = Optim.minimum(sol.optimized)
Optim.iterations(sol::FittingSolution) = Optim.iterations(sol.optimized)
Optim.iteration_limit_reached(sol::FittingSolution) = Optim.iteration_limit_reached(sol.optimized)
Optim.trace(sol::FittingSolution) = Optim.trace(sol.optimized)
Optim.x_trace(sol::FittingSolution) = Optim.x_trace(sol.optimized)
Optim.f_trace(sol::FittingSolution) = Optim.f_trace(sol.optimized)
Optim.f_calls(sol::FittingSolution) = Optim.f_calls(sol.optimized)
Optim.converged(sol::FittingSolution) = Optim.converged(sol.optimized)

end
