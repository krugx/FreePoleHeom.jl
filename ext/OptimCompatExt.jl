module OptimCompatExt

using Optim
import GRAPE: run_optimizer, step_width, search_direction, _apply_convergence_check!
using GRAPE: GrapeWrk, update_result!

# Compatibility override for Optim.jl result fields (ls_success removed in newer Optim)
function run_optimizer(
    optimizer::Optim.LBFGS,
    wrk,
    fg!,
    callback,
    check_convergence
)

    tol_options = Optim.Options(
        x_tol = get(wrk.kwargs, :x_tol, 0.0),
        f_tol = get(wrk.kwargs, :f_tol, 0.0),
        g_tol = get(wrk.kwargs, :g_tol, 1e-8),
    )
    if any(wrk.lower_bounds .> -Inf) || any(wrk.upper_bounds .< Inf)
        error("bounds are not implemented for Optim.jl optimization")
    end
    initial_x = wrk.pulsevals
    method = optimizer
    objective = Optim.promote_objtype(method, initial_x, :finite, true, Optim.only_fg!(fg!))
    wrk.optimizer_state = Optim.initial_state(method, tol_options, objective, initial_x)

    function optim_callback(optimization_state::Optim.OptimizationState)
        iter = wrk.result.iter + 1
        update_result!(wrk, iter)
        info_tuple = callback(wrk, wrk.result.iter)
        if hasproperty(objective, :DF)
            wrk.gradient .= objective.DF
        elseif (optimization_state.iteration == 1)
            @error "Cannot determine guess gradient"
        end
        copyto!(wrk.pulsevals_guess, wrk.pulsevals)
        wrk.fg_count .= 0
        if !(isnothing(info_tuple) || isempty(info_tuple))
            push!(wrk.result.records, info_tuple)
        end
        _apply_convergence_check!(wrk.result, check_convergence)
        return wrk.result.converged
    end

    options = Optim.Options(
        callback = optim_callback,
        iterations = (wrk.result.iter_stop - wrk.result.iter_start),
        x_tol = get(wrk.kwargs, :x_tol, 0.0),
        f_tol = get(wrk.kwargs, :f_tol, 0.0),
        g_tol = get(wrk.kwargs, :g_tol, 1e-8),
        show_trace = get(wrk.kwargs, :show_trace, false),
        extended_trace = get(wrk.kwargs, :extended_trace, false),
        store_trace = get(wrk.kwargs, :store_trace, false),
        show_every = get(wrk.kwargs, :show_every, 1),
        allow_f_increases = get(wrk.kwargs, :allow_f_increases, false),
    )

    res = Optim.optimize(objective, initial_x, method, options, wrk.optimizer_state)

    if hasproperty(res, :ls_success)
        if !res.ls_success
            @error "optimization failed (linesearch)"
            wrk.result.message = "Failed linesearch"
        end
    end
    if hasproperty(res, :stopped_by)
        if res.stopped_by.f_increased
            @error "loss of monotonic convergence (try allow_f_increases=true)"
            wrk.result.message = "Loss of monotonic convergence"
        end
    end
    if !wrk.result.converged
        @warn "Optimization failed to converge"
    end

    return nothing

end

function run_optimizer(
    optimizer::Optim.GradientDescent,
    wrk,
    fg!,
    callback,
    check_convergence
)

    tol_options = Optim.Options(
        x_tol = get(wrk.kwargs, :x_tol, 0.0),
        f_tol = get(wrk.kwargs, :f_tol, 0.0),
        g_tol = get(wrk.kwargs, :g_tol, 1e-8),
    )
    if any(wrk.lower_bounds .> -Inf) || any(wrk.upper_bounds .< Inf)
        error("bounds are not implemented for Optim.jl optimization")
    end
    initial_x = wrk.pulsevals
    method = optimizer
    objective = Optim.promote_objtype(method, initial_x, :finite, true, Optim.only_fg!(fg!))
    wrk.optimizer_state = Optim.initial_state(method, tol_options, objective, initial_x)

    function optim_callback(optimization_state::Optim.OptimizationState)
        iter = wrk.result.iter + 1
        update_result!(wrk, iter)
        info_tuple = callback(wrk, wrk.result.iter)
        if hasproperty(objective, :DF)
            wrk.gradient .= objective.DF
        elseif (optimization_state.iteration == 1)
            @error "Cannot determine guess gradient"
        end
        copyto!(wrk.pulsevals_guess, wrk.pulsevals)
        wrk.fg_count .= 0
        if !(isnothing(info_tuple) || isempty(info_tuple))
            push!(wrk.result.records, info_tuple)
        end
        _apply_convergence_check!(wrk.result, check_convergence)
        return wrk.result.converged
    end

    options = Optim.Options(
        callback = optim_callback,
        iterations = (wrk.result.iter_stop - wrk.result.iter_start),
        x_tol = get(wrk.kwargs, :x_tol, 0.0),
        f_tol = get(wrk.kwargs, :f_tol, 0.0),
        g_tol = get(wrk.kwargs, :g_tol, 1e-8),
        show_trace = get(wrk.kwargs, :show_trace, false),
        extended_trace = get(wrk.kwargs, :extended_trace, false),
        store_trace = get(wrk.kwargs, :store_trace, false),
        show_every = get(wrk.kwargs, :show_every, 1),
        allow_f_increases = get(wrk.kwargs, :allow_f_increases, false),
    )

    res = Optim.optimize(objective, initial_x, method, options, wrk.optimizer_state)

    if hasproperty(res, :ls_success)
        if !res.ls_success
            @error "optimization failed (linesearch)"
            wrk.result.message = "Failed linesearch"
        end
    end
    if hasproperty(res, :stopped_by)
        if res.stopped_by.f_increased
            @error "loss of monotonic convergence (try allow_f_increases=true)"
            wrk.result.message = "Loss of monotonic convergence"
        end
    end
    if !wrk.result.converged
        @warn "Optimization failed to converge"
    end

    return nothing

end

function step_width(wrk::GrapeWrk{Optim.LBFGS})
    return wrk.optimizer_state.alpha
end

function search_direction(wrk::GrapeWrk{Optim.LBFGS})
    return wrk.optimizer_state.s
end

function step_width(wrk::GrapeWrk{Optim.GradientDescent})
    return wrk.optimizer_state.alpha
end

function search_direction(wrk::GrapeWrk{Optim.GradientDescent})
    return wrk.optimizer_state.s
end

end # module
