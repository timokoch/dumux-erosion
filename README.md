Erosion model
===================

Similar to Derr et al. (2020) Phys. Rev. Lett., [10.1103/PhysRevLett.125.158002](https://doi.org/10.1103/PhysRevLett.125.158002).

Governing equations are

```math
\begin{align}
\kappa \frac{\partial p}{\partial t} - \nabla \cdot \left( K(\phi) \nabla p \right) &= 0, \\
\frac{\partial \phi}{\partial t} - \nabla \cdot \left( D \nabla \phi \right) &= -\phi \mathrm{max}\left\lbrace 0, \nabla p : \nabla p - \sigma^2(g) \right\rbrace, \\
R \frac{\partial g}{\partial t} - \nabla \cdot \left( B \nabla g \right) &= \phi - g,
\end{align}
```

where $\phi$ is the solid volume fraction, $K(\phi) = \frac{(1-\phi)^3}{\phi^2}$,
$\sigma^2(g) = \frac{H(g)-H(0)}{H(1)-H(0)}$ with $H(g) = \frac{1}{2}\left[1 -  \mathrm{tanh}(\omega (g^* - g))\right]$.
We erode when $\nabla p : \nabla p > \sigma^2(g)$, where $\sigma^2(g)$ is a yield stress that depends on the solidity $g$.


Installation
=========================

Requirements
--------------

* cmake >= 3.13
* C++17 compliant compiler (e.g g++ >= 8 or clang >= 6)
* pkg-config

Setup
--------------

Make a new folder (e.g. `dumux`) that will contain all modules.
Inside of this folder, clone this repo.
Folder structure would look like this

```
dumux
└───dumux-granular
```

Then run from the top folder (`dumux`):

```
./dumux-granular/setup.sh
```

to download Dune/DuMux dependencies and configure and build the project.
The script takes care of this but if you are manually cloning the dependencies
make sure to use that branch.

After that folder structure will look like this:

```
dumux
├───dumux-granular
│    ├───build-cmake/app
│    ├───app
│    ├───CMakeLists.txt
│    ....
├───dune-common
├───dune-geometry
├───dune-grid
├───dune-istl
├───dune-localfunctions
├───dune-alugrid
└───dumux
```

Applications
--------------

Programs are build and run inside the build folder `build-cmake/app`. The build folder
mimicks the structure of the source folder `app` where you find the header files (source files).

2D example
------------------------------------

Go to `dumux-granular/build-cmake/app`
and run

* `make test_erosion_2d`
* `DUMUX_NUM_THREADS=1 mpirun -np 8 ./test_erosion_2d params.input`

Runtime parameters can be set in `params.input`.
View result files with ParaView.
