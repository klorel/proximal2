include("common.jl")
include("cutting_plane.jl")
include("lp.jl")
include("proximal.jl")

###
# general setup
###
nVariables = 50;
nCassures = 50;
xMin, xMax = build_bounds(nVariables)
X = build_X(nVariables, nCassures, xMin, xMax)

# println("-----------------------------------")
# println("-------- CUTTING PLANE ------------")
# println("-----------------------------------")
# launch_cutting_plane(X, xMin, xMax)

println("-----------------------------------")
println("-------- PROXIMAL -----------------")
println("-----------------------------------")

launch_proximal(X, xMin, xMax)

println("-------- CHECK         ------------")
problem = Model(solver=CplexSolver(CPX_PARAM_THREADS=1, CPX_PARAM_LPMETHOD=4));
x_full, x_left, x_right = full_jump_model(X, problem, xMin, xMax)

status = solve(problem)
#
# x_value, rhs, sub_gradient = build_cut(X, x_full)
# # print(problem);
# println("obj : ", getvalue(getobjective(problem)))
# println("rhs : ", rhs)
# println("X = ", X);
# println("x_left = ", getvalue(x_left));
# println("x_right = ", getvalue(x_right));
# println("x_full = ", getvalue(x_full));
