using OrdinaryDiffEq
using Trixi

equations = CompressibleEulerEquations2D(1.4) # set gas gamma = 1.4

# initial condition
initial_condition = initial_condition_convergence_test
source_term       = source_terms_convergence_test

# boundary condition types
boundary_condition_dirichlet= BoundaryConditionDirichlet(initial_condition)

# boundary condition dictionary
boundary_conditions = Dict( :OuterCircle => boundary_condition_dirichlet,
                            :LeftSlant => boundary_condition_dirichlet,
                            :RightSlant => boundary_condition_dirichlet,
                            :IceCream => boundary_condition_dirichlet)

# DGSEM solver.
#    1) polydeg must be >= the polynomial order set in the HOHQMesh control file to guarantee
#       freestream preservation. As a extra task try setting polydeg=3
#    2) VolumeIntegralFluxDifferencing with central volume flux is activated
#       for dealiasing
volume_flux = flux_ranocha
solver = DGSEM(polydeg=6, surface_flux=flux_hll,
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

# create the unstructured from your mesh file
mesh_file = joinpath(@__DIR__, "out/mesh_for_live_tutorial.mesh")
mesh = UnstructuredMesh2D(mesh_file)

# Create semidiscretization with all spatial discretization-related components
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions=boundary_conditions,
                                    source_terms=source_term)

# Create ODE problem from semidiscretization with time span from 0.0 to 2.0
tspan = (0.0, 5.0)
ode = semidiscretize(semi, tspan)


# Create the callbacks to output timing information, solution files, and adapt the time step
summary_callback = SummaryCallback()
save_solution = SaveSolutionCallback(interval=10,
                                     save_initial_solution=true,
                                     save_final_solution=true)
stepsize_callback = StepsizeCallback(cfl=1.0)

callbacks = CallbackSet(summary_callback, save_solution, stepsize_callback)

# Evolve ODE problem in time using `solve` from OrdinaryDiffEq
sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
# print the timer summary
summary_callback()