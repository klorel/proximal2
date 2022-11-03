function update_center(ctr_proximal, nVariables, center)
    for i in 1:nVariables
        JuMP.setRHS(ctr_proximal[i], center[i])
    end
end
function launch_proximal(X, xMin, xMax)
    master = Model(solver=CplexSolver(CPX_PARAM_THREADS=1, CPX_PARAM_QPMETHOD=2, CPXPARAM_ScreenOutput=0));
    nVariables, nCassures = size(X)
    x, α = build_master_cutting_plane(X, master, xMin, xMax)
    var_proximal = @variable(master, var_proximal[1:nVariables])
    ctr_proximal = @constraint(master, ctr_proximal[i in 1:nVariables], x[i]==var_proximal[i])

    stop = false
    lb = -1e20
    ub = +1e20
    n_ite=0

    best_ub = ub
    prediction = 0
    nb_ss=0
    nb_ns=0
    nb_update = 3
    # nb_update=3
    step="NONE"
    weight=1e-0
    tol=1e-1
    @objective(master, :Min, α +sum(weight*var_proximal[i]^2 for i in 1:nVariables))

    center=Base.zeros(nVariables)
    update_center(ctr_proximal, nVariables, center)
    while ! stop
        n_ite+=1
        solve(master)
        lb = getvalue(α)

        x_value, rhs, sub_gradient = build_cut(X, x)
        ub = rhs

        prediction = lb-best_ub
        if ub - best_ub <= tol*prediction || n_ite==1
            step="SS"
            center = x_value
            best_ub = ub
            update_center(ctr_proximal, nVariables, center)
            nb_ss+=1
            # print(master)
        else
            step="NS"
            nb_ns+=1

        end

        if nb_ss>nb_update
            weight*=0.5
            nb_ss=0
            @objective(master, :Min, α +sum(weight*var_proximal[i]^2 for i in 1:nVariables))
        end
        if nb_ns>nb_update
            weight*=2
            nb_ns=0
            @objective(master, :Min, α +sum(weight*var_proximal[i]^2 for i in 1:nVariables))
        end

        if n_ite>1 && -prediction<=1e-6
            stop = true
        else
            @constraint(master, α>=rhs+sum( sub_gradient[i] * (x[i]-x_value[i]) for i in 1:nVariables))
        end
        @printf("%10d%6s%20.10E%20.10E%20.10E%20.10E\n", n_ite, step, weight, lb, ub, prediction)
    end
end
