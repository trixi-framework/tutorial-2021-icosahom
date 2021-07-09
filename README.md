# Trixi.jl: High-Order Numerical Simulations of Hyperbolic PDEs in Julia

This is the companion repository for the tutorial on
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl) at
[ICOSAHOM 2021](https://www.icosahom2020.org):

**Tutorial 2: Trixi.jl**<br />
*Venue:* gather.town (link to be added later)<br />
*Date & time:* Wednesday, 14th July 2021, 6:15pm - 8:15pm (CEST)<br />
*Chaired by:* Hendrik Ranocha, Michael Schlottke-Lakemper<br />
*Link:* [conference agenda](https://www.conftool.com/icosahom2020/index.php?page=browseSessions&form_session=71&presentations=show)

**Note: This repository is still work-in-progress and will receive further
updates until the day before the tutorial session.**

In case of questions before the beginning of the tutorial, please get in touch with
[Hendrik](https://ranocha.de) or
[Michael](https://www.mi.uni-koeln.de/NumSim/schlottke-lakemper),
[create an issue](https://github.com/trixi-framework/tutorial-2021-icosahom/issues/new),
or
[join the Trixi.jl Slack workspace](https://join.slack.com/t/trixi-framework/shared_invite/zt-sgkc6ppw-6OXJqZAD5SPjBYqLd8MU~g).

**Table of contents**
1. [Schedule](#schedule)
2. [Floor plan](#floor-plan)
3. [Tutorial files](#tutorial-files)
4. [Abstract](#abstract)
5. [Getting started](#getting-started)
   1. [Using mybinder.org](#using-mybinderorg)
   2. [Setting up a local Julia/Jupyter installation](#setting-up-a-local-juliajupyter-installation)
6. [Authors](#authors)
7. [License](#license)


## Schedule

The [Auditorium](#auditorium) will be used to give a few guided sessions on how
to use Trixi.jl and Julia. The [Computer Lab](#computer-lab) has work spaces to
try out Trixi.jl on your own and with the help of our instructors. Please check
out the [floor plan](#floor-plan) to help with finding out where everything is
located.

*Note:* All times are in Central European Summer Time (CEST).

**Note: The schedule is still work-in-progress and will receive further updates.**

### Auditorium
1. **6:15pm: What is Trixi.jl and how can I use it for my own projects?**<br />
   *Speaker:* Hendrik Ranocha<br />
   *Duration:* approx. 20 minutes + Q&A<br />

   * Overview of Trixi's capabilities
   * Learn how to set up and run simulations
   * Extend Trixi for your own research
   * Q&A

   *Hands-on:* the Jupyter notebook `introduction_to_trixi.ipynb` (see
   [below](#tutorial-files)) can be used to simultaneously try out the
   examples in the presentation

2. **6:45pm: Getting started with Julia**<br />
   *Speaker:* Hendrik Ranocha<br />
   *Duration:* approx. 10 minutes + Q&A

   * Brief introduction to Julia for newcomers
   * Learn about the syntax, the type model, and multiple dispatch
   * Q&A
 <br />
   *Hands-on:* the Jupyter notebook `introduction_to_julia.ipynb` (see
   [below](#tutorial-files)) can be used to simultaneously try out the
   examples in the presentation

3. **7:00pm: From grid generation to visualization: Working with unstructured curved meshes in Trixi.jl**<br />
   *Speaker:* Andrew R. Winters<br />
   *Duration:* approx. 15 minutes + Q&A

   * Introduction to using HOHQMesh.jl to generate high-order meshes
   * Understand how to set up an unstructured simulation in Trixi
   * Learn about the visualization pipeline with Trixi2Vtk.jl and ParaView
   * Q&A

### Computer Lab (hands-on sessions)
The Computer Lab can be used any time to try out Julia and Trixi.jl either in
small groups or alone. There are 18 lab tables available (see [floor plan](#floor-plan)),
each being a private space where you can talk to all people at the same table.
Throughout the room and/or at the organizer table near the south wall, you will
find a number of instructors who you can ask for assistance with using Trixi.jl
or Julia (just walk up to us and ask away!):

* Christof Czernik
* Erik Faulhaber
* Hendrik Ranocha
* Michael Schlottke-Lakemper
* Andrew R. Winters


## Floor plan
**tbd.**


## Tutorial files
| Item | [nbviewer](https://nbviewer.jupyter.org/) | [mybinder](https://mybinder.org/) |
|:-|:-:|:-:|
| [`introduction_to_trixi.ipynb`](introduction_to_trixi.ipynb) | [![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/trixi-framework/tutorial-2021-icosahom/blob/main/introduction_to_trixi.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/trixi-framework/tutorial-2021-icosahom/HEAD?filepath=introduction_to_trixi.ipynb) |
| [`introduction_to_julia.ipynb`](introduction_to_julia.ipynb) | [![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/trixi-framework/tutorial-2021-icosahom/blob/main/introduction_to_julia.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/trixi-framework/tutorial-2021-icosahom/HEAD?filepath=introduction_to_julia.ipynb) |
| [`exercises_linear_advection.ipynb`](exercises_linear_advection.ipynb) | [![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/trixi-framework/tutorial-2021-icosahom/blob/main/exercises_linear_advection.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/trixi-framework/tutorial-2021-icosahom/HEAD?filepath=exercises_linear_advection.ipynb) |
| [`exercises_euler.ipynb`](exercises_euler.ipynb) | [![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/trixi-framework/tutorial-2021-icosahom/blob/main/exercises_euler.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/trixi-framework/tutorial-2021-icosahom/HEAD?filepath=exercises_euler.ipynb) |
| [`unstructured_tutorial/exercises_unstructured_quad.md`](unstructured_tutorial/exercises_unstructured_quad.md) | - | - |

Additional tutorials are available in the
[documentation of Trixi.jl](https://trixi-framework.github.io/Trixi.jl/stable/).


## Abstract

Trixi.jl is a numerical simulation framework for adaptive, high-order
discretizations of conservation laws. It has a modular architecture that
allows users to easily extend its functionality and was designed to be
useful to experienced researchers and new users alike.
In this tutorial we will demonstrate what you can do with Trixi.jl and
how you can use it (and extend it) for your own research. Furthermore,
there will be a hands-on session to try out Trixi for yourself, and time
to talk to the Trixi developers. We will also include a short
introduction to Julia for those who have not used the language before.


## Getting started

You can view a static version of the Jupyter notebooks `*.ipynb`

- directly on GitHub (select the notebook; this may fail sometimes)
- or on [nbviewer.jupyter.org](https://nbviewer.jupyter.org/)
  (select the "render" badges in the table of contents above)

These static versions do not contain output of the code cells.

### Using mybinder.org
The easiest way to get started is to click on the *Launch Binder* badges
in the table of contents above.
This launches the notebook for interactive use in your browser without the need
to download or install anything locally.

In this case, you can skip the rest of this *Getting started* section. A
Jupyter instance will be started automagically in the cloud via
[mybinder.org](https://mybinder.org), and the notebook will loaded directly from
this repository.

*Note:*  Depending on current usage and available resources, it typically takes
a few minutes to launch a notebook with [mybinder.org](https://mybinder.org)
(sometimes a little longer), so try to remain patient. Similarly, the first two
cells of the notebook take much longer to execute than usual (around 1.5 minutes
for the first Trixi simulation and about 1 minute for the first plot), since
Julia compiles all methods "just-ahead-of-time" at first use. Subsequent runs
will be much faster.

### Setting up a local Julia/Jupyter installation
Alternatively, you can also clone this repository and open the notebook on your
local machine. This is recommended if you already have a Julia + Jupyter setup
or if you plan to try out Julia anyways.

#### Installing Julia and IJulia
To obtain Julia, go to https://julialang.org/downloads/ and download the latest
stable release (v1.6.1 as of 2021-07-07; neither use the LTS release nor
Julia Pro). Then, follow the
[platform-specific instructions](https://julialang.org/downloads/platform/)
to install Julia on your machine. Note that there is no need to compile anything
if you are using Linux, MacOS, or Windows.

After the installation, open a terminal and start the Julia *REPL*
(i.e., the interactive prompt) with
```shell
julia
```
To use the notebook, you also need to get the
[IJulia](https://github.com/JuliaLang/IJulia.jl) package, which provides a Julia
backend for Jupyter. In the REPL, execute
```julia
using Pkg
Pkg.add("IJulia")
```
to install IJulia. For more details, especially on how to use an existing Jupyter
installation, please refer to the
[IJulia documentation](https://julialang.github.io/IJulia.jl/stable/).
From here on, we assume that you have a working installation of Julia, Jupyter,
and the Julia kernel for Jupyter.

#### Installing the required Julia packages
To make the notebook fully reproducible, we have used Julia's package manager
to pin all packages to a fixed release. This ensures that you always have a
Julia environment in which all examples in this notebook work. Later you can
always install the latest versions of Trixi and its dependencies by following
the instructions in the Trixi
[documentation](https://trixi-framework.github.io/Trixi.jl/stable/).

If you have not done it yet, clone the repository where this notebook is stored:
```shell
git clone https://github.com/trixi-framework/tutorial-2021-icosahom.git
```
Then, navigate to your repository folder and install the required packages:
```shell
cd tutorial-2021-icosahom
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
This will download and build all required packages, including the ODE package
[OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl), the visualization
package [Plots](https://github.com/JuliaPlots/Plots.jl), and of course
[Trixi](https://github.com/trixi-framework/Trixi.jl).
The `--project=.` argument tells Julia to use the `Project.toml`
and `Manifest.toml` files from this repository to figure out which packages to install.

As an alternative to running the examples in the notebook directly, you may
also just view the notebook *statically* by opening it within
[Jupyter NBViewer](https://nbviewer.jupyter.org/github/trixi-framework/talk-2021-Introduction_to_Julia_and_Trixi/blob/main/Talk.ipynb?flush_cache=true).

*General note:* Make sure that you execute the examples (either in the notebook
or in the REPL) *in order*, at least for the first time. Both the notebook and
the Julia REPL maintain an internal state and and some snippets depend on
earlier statements having been executed.

#### Displaying the presentation

To display the presentation as in the talk (skipping some cells/slides that
provide further information), you need the
[Jupyter extension RISE](https://rise.readthedocs.io/en/stable),
that you can install via
```shell
pip3 install --user RISE
```
After opening the Jupyter notebook, you can enter the RISE presentation mode
with `Alt + R`.


## Authors
This repository was initiated by
[Hendrik Ranocha](https://ranocha.de),
[Michael Schlottke-Lakemper](https://www.mi.uni-koeln.de/NumSim/schlottke-lakemper)
and [Andrew R. Winters](https://liu.se/en/employee/andwi94).


## License
The contents of this repository are licensed under the MIT license
(see [LICENSE.md](LICENSE.md)).
