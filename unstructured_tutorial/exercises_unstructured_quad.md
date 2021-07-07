# Tutorial on Trixi.jl at ICOSAHOM 2021

## Exercises: Unstructured curvilinear quadrilateral solver

### Andrew Winters

## Authors and license

This material is distributed by Michael Schlottke-Lakemper, Hendrik Ranocha and Andrew Winters under the MIT license.

# Before you begin these exercises

Follow the instructions to install Julia, Trixi, and Trixi2Vtk provided in the Jupyter notebook `introduction_to_trixi.ipynb`.

Copy the files `box_with_object.control` and `tutorial_unstructured_exercise_3.jl` into your Trixi.jl directory.

The commands and problem setups in this exercise set were set up and tested with Julia v1.6.1 but may also work with other (newer) versions.

# Exercise 1: Run and visualize your first unstructured simulation

Trixi supports solving hyperbolic problems on several mesh types. Unstructured curvilinear quadrilateral meshes is
one such option. There is a default example for this mesh type that can be executed by
```julia
julia> trixi_include(default_example_unstructured())
```
This will compute a smooth, manufactured solution test case for the 2D compressible Euler equations
on the curved quadrilateral mesh described in the
[Trixi documentation](https://trixi-framework.github.io/Trixi.jl/stable/meshes/unstructured_quad_mesh/).

Apart from the usual error and timing output provided by the Trixi run, it is useful to visualize and inspect
the solution. Currently, for solutions on unstructured quadrilateral meshes, this requires post-processing the
Trixi output file(s) using the `Trixi2Vtk` tool and plotting them with [Paraview](https://www.paraview.org/download/).

To convert the `.h5` file(s) into VTK format execute the following
```julia
julia> using Trixi2Vtk

julia> trixi2vtk("out/solution_000180.h5", output_directory="out")
```
Note this step takes 15-30 seconds as the package `Trixi2Vtk` must be precompiled and executed for the first time
in your REPL session. The above `trixi2vtk` command will convert the solution file at the final time into a `.vtu`
which can be readin and visualize with Paraview. A required argument for `trixi2vtk` is to point to the `output_directory`
where the new files will be saved. An optional argument that can be set with `trixi2vtk` is to specify the number of
visualization nodes. For instance, if we want to use 12 uniformly spaced nodes for visualization we can execute
```julia
julia> trixi2vtk("out/solution_000180.h5", output_directory="out", nvisnodes=12)
```
By default `trixi2vtk` will use the same number of nodes as specified in the `elixir` file used to
run the simulation.

Finally, if you want to convert all the solution files to VTK execute
```julia
julia> trixi2vtk("out/solution_000*", output_directory="out", nvisnodes=12)
```
then it is possible to open the `.pvd` file with Paraview and create a video of the simulation.

# Exercise 2: Generate an unstructured quadrilateral mesh for use in Trixi

Where did the mesh used in the previous exercise come from? For the elixir files present in `examples/unstructured_2d_dgsem/` certain
mesh files are automatically downloaded and used for the Trixi execution. For instance, after executing *Exercise 1*
you will now see the file `mesh_trixi_unstructured_mesh_docs.mesh` in the `examples/unstructured_2d_dgsem/` folder.

This exercise, broken into three components, provides some background detail on the mesh generator
of the Trixi toolchain and how it can be called directly in Julia to create mesh files locally for a Trixi
simulation.

## Exercise 2a: Obtain the mesh generator

To obtain unstructured curvilinear quadrilateral meshes in the format required by Trixi we use the
[*High-Order Hex-Quad (HOHQ) Mesh*](https://github.com/trixi-framework/HOHQMesh) generator created and developed by David Kopriva.
HOHQMesh is a mesh generator specifically designed for spectral element methods where elements can be larger (due to the high accuracy
of the spatial approximation) and provides high-order boundary curve information (needed to accurate set boundary conditions).
For more information about the design and features of HOQHMesh you can refer to its
[official documentation](https://trixi-framework.github.io/HOHQMesh/).

HOHQMesh is incorporated in the Trixi framework via the registered package [HOHQMesh.jl](https://github.com/trixi-framework/HOHQMesh.jl).
This package provides a Julia wrapper for the HOHQMesh generator that allows users to easily create mesh files without the need to build
HOHQMesh from source. To install the HOHQMesh package execute
```julia
julia> import Pkg; Pkg.add("HOHQMesh")
```
Now we are ready to generate an unstructured quadrilateral mesh that can be used by Trixi.

## Exercise 2b: Explanation of a HOHQMesh control file

The creation of a mesh using the HOHQMesh generator is driven by a **control file**. Is this file the user dictates
the domain to be meshed, prescribes any desired boundary curvature, the polynomial order of said boundaries, etc.
In this tutorial we only cover several basic features of the possible control inputs. For a complete discussion
on the control file see the [HOHQMesh documentation](https://trixi-framework.github.io/HOHQMesh/).

Open the file `box_with_object.control` provided in this tutorial. To begin we note that blank space or anything after a `%` is ignored
by HOHQMesh at readin. The first three block of information are wrapped within a `CONTROL_INPUT` environment block as they define the
core components of the quadrilateral mesh that will be generated.

The first block of information in `RUN_PARAMETERS` is
```
\begin{RUN_PARAMETERS}
   mesh file name   = box_with_object.mesh
   plot file name   = box_with_object.tec
   stats file name  = none
   mesh file format = ISM-v2
   polynomial order = 4
   plot file format = skeleton
\end{RUN_PARAMETERS}
```
The mesh and plot file names will be the files created by HOHQMesh once successfully executed. The stats file name is
available if you wish to also save a collection of mesh statistics. For this example it is deactivated.
These file names given within `RUN_PARAMETERS` **should match** that of the control file, although this is not required by
HOHQMesh it is a useful style convention.
The mesh file format `ISM-v2` is the format currently required by Trixi. The `polynomial order` prescribes the order
of an interpolant constructed on the Chebyshev-Gauss-Lobatto nodes that is used to represent any curved boundaries on a particular element.
The plot file format of `skeleton` means that visualizing the plot file will only draw the element boundaries (and no internal nodes).
Alternatively, the format can be set to `sem` to visualize the interior nodes of the approximation as well.

The second block of information in `BACKGOUND_GRID`is
```
\begin{BACKGROUND_GRID}
   x0 = [-3.0, -3.0, 0.0]
   dx = [1.0, 1.0, 0.0]
   N  = [6,6,1]
\end{BACKGROUND_GRID}
```
This lays a grid of Cartesian elements for the domain beginning at the point `x0`. The value of `dx`, which could differ in each direction
if desired, controls the step size taken in each Cartesian direction. The values in `N` set how many Cartesian box elements
are set in each coordinate direction. The above parameters define a $`6\times 6`$ element square mesh on $`[-3,3]^2`$. Further, this sets up four outer boundaries of the domain that are given the default names: `Top, Left, Bottom, Right`.

The third block of information in `SPRING_SMOOTHER` is
```
\begin{SPRING_SMOOTHER}
   smoothing            = ON
   smoothing type       = LinearAndCrossBarSpring
   spring constant      = 1.0
   mass                 = 1.0
   rest length          = 0.0
   damping coefficient  = 5.0
   number of iterations = 25
   time step            = 0.1
\end{SPRING_SMOOTHER}
```
Once HOHQMesh generates the mesh a spring-mass-dashpot model is created to smooth the mesh and create "nicer" quadrilateral elements.
The parameters for this spring-mass-dashpot model have been selected after a fair amount of experimentation across many meshes
and typically will not need to be altered. However, if you ever wish to deactivate this feature you can set `smoothing = OFF`
(or remove this block from the control file).

After the `CONTROL_INPUT` environment block comes the `MODEL` environment block. It is here where the user can prescribe curved boundary information either an `OUTER_BOUNDARY` (not covered in this tutorial) or `INNER_BOUNDARIES`. There are several options to
describe the boundary curve data to HOHQMesh like splines or parametric curves.

For the example `box_with_object.control` we define a single internal boundary using a parametric equation for a circle
of radius $`r`$ centered at the point $`(x_c, y_c)`$, i.e.,
```math
x(t) = x_c + r\cos(2\pi t), y(t) = y_c + r\sin(2\pi t)
```
where we select the radius to be $`r=0.3`$ and the center to be the origin.
Within the HOHQMesh control input each curve must be assigned to a `CHAIN` as shown below in the complete
`INNER_BOUNDARIES` block
```
\begin{INNER_BOUNDARIES}

   \begin{CHAIN}
       name = InnerCircle1
       \begin{PARAMETRIC_EQUATION_CURVE}
          name = Circle
          xEqn = f(t) = 0.0 + 0.3*cos(2*pi*t)
          yEqn = f(t) = 0.0 + 0.3*sin(2*pi*t)
          zEqn = z(t) = 0.0
      \end{PARAMETRIC_EQUATION_CURVE}
   \end{CHAIN}

\end{INNER_BOUNDARIES}
```
It is important to note there are two `name` quantities one for the `CHAIN` and one for the `PARAMETRIC_EQUATION_CURVE`.
The name for the `CHAIN` is used internally by HOHQMesh, so if you have multiple `CHAIN`s they **must be given a unique name**.
The name for the `PARAMETRIC_EQUATION_CURVE` with be printed to the appropriate boundaries within the `.mesh` file produced by
HOHQMesh. Trixi uses this boundary name to assign boundary conditions in an elixir file as done next in *Exercise 3*.

## Exercise 2c: Generate an unstructured quadrilateral mesh

To generate the mesh with a box around a circular object execute
```julia
julia> control_file = joinpath(@__DIR__, "box_with_object.control");
julia> output = generate_mesh(control_file);
julia> println(output)
 2D Mesh Statistics:
    Total time         =    1.6957000000000000E-002
    Number of nodes    =           79
    Number of Edges    =          144
    Number of Elements =           65

 Mesh Quality:
         Measure         Minimum         Maximum         Average  Acceptable Low Acceptable High       Reference
     Signed Area      0.06278469      1.00000000      0.38069910      0.00000000    999.99900000      1.00000000
    Aspect Ratio      1.00000000      1.80473785      1.42299887      1.00000000    999.99900000      1.00000000
       Condition      1.00000000      1.50000000      1.24448991      1.00000000      4.00000000      1.00000000
      Edge Ratio      1.00000000      2.17880680      1.59557964      1.00000000      4.00000000      1.00000000
        Jacobian      0.03911184      1.00000000      0.32191112      0.00000000    999.99900000      1.00000000
   Minimum Angle     45.00000000     90.00000000     67.00242604     40.00000000     90.00000000     90.00000000
   Maximum Angle     90.00000000    135.00000000    111.67181399     90.00000000    135.00000000     90.00000000
       Area Sign      1.00000000      1.00000000      1.00000000      1.00000000      1.00000000      1.00000000
```
The third command that prints the mesh statistics to the screen is optional. The `box_with_object.mesh` and `box_with_object.tec` files
are placed into the `out/` by default. You can visualize the mesh that was just generated also using Paraview simply
select "Tecplot Reader" when prompted after opening the `box_with_object.tec` file.
From such a visualization it appears that the mesh does not have a curved interior boundary, but this is an artifact of plotting software
combined with using `plot file format = skeleton` in the `RUN_PARAMETERS`.

Regenerate the mesh but change the control file to use `plot file format = sem`, execute
```julia
julia> output = generate_mesh(control_file);
```
and re-load the `box_with_object.tec` file in Paraview to also visualize that internal mesh points/boundary curvature.

# Exercise 3: Run and visualize on `box_with_object.mesh`

With the new mesh generated from *Exercise 3* we are ready to run another Trixi simulation on an unstructured quadrilateral mesh.
For this we must create a new elixir file. As in *Exercise 1* we will solve the 2D compressible Euler equations.

This elixir file already creates a new initial condition for a uniform background flow state with a free stream Mach number of 0.3.
An exercise dedicated to modifying the initial conditions is provided in `exercises_linear_advection.ipynb` for
the linear advection equations.

The focus of this exercise is to specify the boundary conditions and the construct the new mesh from the
file that was generated in the previous exercise. It is straightforward to set the different boundary
condition types in an elixir by assigning a particular function to a boundary name inside a
Julia dictionary, `Dict`, variable. Observe that the names of these boundaries match those that were either
default, e.g. `Bottom`, or user assigned, e.g. `Circle`, within HOHQMesh. For this problem setup use
* Freestream boundary conditions on the four box edges
* Free slip wall boundary condition on the interior circular boundary

Construct the boundary condition and correctly load your mesh file by completing the code provided in
`tutorial_unstructured_exercise_3.jl` (reproduced below).
```julia
using OrdinaryDiffEq
using Trixi

equations = CompressibleEulerEquations2D(1.4) # set gas gamma = 1.4

# freestream flow state with Ma_inf = 0.3
@inline function uniform_flow_state(x, t, equations::CompressibleEulerEquations2D)

  # set the freestream flow parameters
  rho_freestream = 1.0
  u_freestream = 0.3
  p_freestream = inv(equations.gamma)

  theta = 0.0 # zero angle of attack
  si, co = sincos(theta)
  v1 = u_freestream * co
  v2 = u_freestream * si

  prim = SVector(rho_freestream, v1, v2, p_freestream)
  return prim2cons(prim, equations)
end

# initial condition
initial_condition = uniform_flow_state

# boundary condition types
boundary_condition_uniform_flow = BoundaryConditionDirichlet(uniform_flow_state)
boundary_condition_slip_wall = BoundaryConditionWall(boundary_state_slip_wall)

# boundary condition dictionary
boundary_conditions = Dict( :Bottom => ,# Your code can be written here
                            :Top    => ,# Your code can be written here
                            :Right  => ,# Your code can be written here
                            :Left   => ,# Your code can be written here
                            :Circle =>  # Your code can be written here
                            )

# DGSEM solver.
#    1) polydeg must be >= the polynomial order set in the HOHQMesh control file to guarantee
#       freestream preservation. As a extra task try setting poyldeg=3
#    2) VolumeIntegralFluxDifferencing with central volume flux is activated
#       for dealiasing
volume_flux = flux_ranocha
solver = DGSEM(polydeg=4, surface_flux=flux_hll,
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

# create the unstructured from your mesh file
mesh_file = # Your code can be written here
mesh = UnstructuredMesh2D(mesh_file)

# Create semidiscretization with all spatial discretization-related components
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions=boundary_conditions)

# Create ODE problem from semidiscretization with time span from 0.0 to 2.0
tspan = (0.0, 2.0)
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
```

Before we run the next simulation it is good to "clean" the output directory of old plotting files with
```
rm out/*.h5 out/*.vtu
```

Once you have modified the elixir file `tutorial_unstructured_exercise_3.jl` appropriately you are ready to execute
your next unstructured Trixi simulation with
```julia
julia> trixi_include("tutorial_unstructured_exercise_3.jl")
```
The simulation should use 131 time steps. You can convert the solution files to VTK on ten visualization nodes with
```julia
trixi2vtk("out/solution_000*", output_directory="out", nvisnodes=10)
```
Then open the `.pvd` file in Paraview and watch the solution video.

# Exercise 4: Generate a mesh with two objects and run again

The final exercise is to demonstrate the ease with which the user can update a new mesh and
modify an existing elixir file to run a new simulation.

To begin "clean" the output directory of old plotting files with
```
rm out/*.h5 out/*.vtu out/*.pvd
```
Next, create a new HOHQMesh control file
```
cp box_with_object.control box_with_two_objects.control
```
Update the `mesh file name` and `plot file name` appropriately in your new control file.
Next, modify the new `box_with_two_objects.control` file to include an additional `CHAIN` inside
the `INNER_BOUNDARIES` control block to include a second inner boundary that is an ellipse with the
parametric equation
```math
x(t) = 0.6\cos(2\pi t), y(t) = 1.0 + 0.3\sin(2\pi t)
```
using the skeleton
```
   \begin{CHAIN}
       name =
       \begin{PARAMETRIC_EQUATION_CURVE}
          name =
          xEqn = f(t) =
          yEqn = f(t) =
          zEqn = z(t) = 0.0
      \end{PARAMETRIC_EQUATION_CURVE}
   \end{CHAIN}
```
Remember to give this second `CHAIN` and `PARAMETRIC_EQUATION_CURVE` unique names!
Generate this new mesh with a box around a two objects by executing the following
```julia
julia> control_file = joinpath(@__DIR__, "box_with_two_objects.control");
julia> output = generate_mesh(control_file);
julia> println(output)
 2D Mesh Statistics:
    Total time         =    5.1221999999999990E-002
    Number of nodes    =          293
    Number of Edges    =          552
    Number of Elements =          258

 Mesh Quality:
         Measure         Minimum         Maximum         Average  Acceptable Low Acceptable High       Reference
     Signed Area      0.00455048      1.08825939      0.13629291      0.00000000    999.99900000      1.00000000
    Aspect Ratio      1.03122978      2.53677069      1.32813030      1.00000000    999.99900000      1.00000000
       Condition      1.00153121      4.54755364      1.22372928      1.00000000      4.00000000      1.00000000
      Edge Ratio      1.05206342      8.96626016      1.67425779      1.00000000      4.00000000      1.00000000
        Jacobian      0.00084206      1.00000000      0.09816746      0.00000000    999.99900000      1.00000000
   Minimum Angle     35.35965372     88.73787332     68.90605888     40.00000000     90.00000000     90.00000000
   Maximum Angle     91.44219947    135.33499769    112.90953573     90.00000000    135.00000000     90.00000000
       Area Sign      1.00000000      1.00000000      1.00000000      1.00000000      1.00000000      1.00000000
```
We provide the mesh statistics output as a reference. There should now appear new
`box_with_two_objects.mesh` and `box_with_two_objects.tec` files in your `out/` folder.

Next, modify your elixir file from *Exercise 3* to set a free slip wall boundary condition
on the new elliptical object and load your new mesh file into the `UnstructuredMesh2D`.
Execute this modified elixir file
```julia
julia> trixi_include("tutorial_unstructured_exercise_3.jl")
```
which will run for 1362 time steps. Similar as before, you can convert the solution files to VTK on
ten visualization nodes with
```julia
trixi2vtk("out/solution_00*", output_directory="out", nvisnodes=10)
```
Then open the `.pvd` file in Paraview and watch the solution video on your two object simulation.
