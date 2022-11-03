function update_center(ctr_proximal, nVariables, center)
    for i in 1:nVariables
        JuMP.set_normalized_rhs(ctr_proximal[i], center[i])
    end
end
function launch_proximal(X, xMin, xMax)
    master = get_my_model();
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
    weight=0
    tol=1e-1
    @objective(master, Min, α +sum(weight*var_proximal[i]^2 for i in 1:nVariables))

    center=Base.zeros(nVariables)
    update_center(ctr_proximal, nVariables, center)
    while ! stop
        n_ite+=1
        optimize!(master)
        lb = value(α)

        x_value, rhs, sub_gradient = build_cut(X, x)
        ub = rhs

        prediction = lb-best_ub

        best_ub = min(ub, best_ub)
        if n_ite>1 && lb >= best_ub-EPS
            stop = true
        else
            @constraint(master, α>=rhs+sum( sub_gradient[i] * (x[i]-x_value[i]) for i in 1:nVariables))
        end
        @printf("%10d%6s%20.10E%20.10E%20.10E%20.10E\n", n_ite, step, weight, lb, ub, prediction)
    end
end
