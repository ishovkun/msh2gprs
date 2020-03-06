# Easy Reservoir Simulation Preprocessor (ERSP)
ERSP is a tool to easily create models for reservoir simulation.
It currently only supports output to AD-GRPS but can easily be extended
to other file formats.
It currently supports mesh files created with GMsh tool.
The properties can be assigned YAML human-readable formats.

## Features
- Supports arbitrary keywords
- Properties can be assigned as expressions as opposed to simply numbers
- Wells
- Supports discrete and embedded fractures for flow and mechanics
- Currently only TPFA is available

## Build

msh2gprs requires a C++-17-compatible compiler and minimum CMake 3.7
(build was tested on GCC 8.2 and clang 7.0).

### External optional libraries
- Optional boost improved the performance
- The DFEM feature requires linking with Gmsh and Eigen.

#### Build Gmsh for DFEM
The DFEM feature requires Gmsh.
A simple tarball download does not cut it since it sometimes causes conflicts
with the system metis (used for multiscale partitioning).
Threfore, the correct way is to build Gmsh from source.
``` cmake
git clone https://gitlab.onelab.info/gmsh/gmsh ./gmsh_git
cd gmsh_git;
mkdir build; cd build;
cmake -DCMAKE_INSTALL_PREFIX=$HOME/build/gmsh-git-install \
    -DENABLE_BUILD_DYNAMIC=ON \
    -DENABLE_BUILD_SHARED=ON \
    -DENABLE_BLAS_LAPACK=ON \
    -DENABLE_PLUGINS=ON \
    -DENABLE_QUADTRI=ON \
    -DENABLE_MSH=ON ..
```
If you're getting errors while linking the main code with gmsh due to metis
issues, try to disable METIS at this step.

#### Eigen
I was writing the code using a system-wide Eigen 3.3.7 installation.
Currently CMake does not support using Eigen from a custom location.

### Boost (optional)
CMake will automatically detect whether boost is available and use it.
If Boost is not available, CMake will stick with using a custom library
for 256-bit integers (used for hashing by angem library).

### msh2gprs
To build mshgprs use the following commands.
``` cmake
git clone --recursive https://github.com/ishovkun/msh2gprs
cd msh2gprs
mkdir build; cd build;
cmake ..
make
```

If you'd like to use a custom Gmsh install, during the cmake run issue

``` cmake
cmake -DGMSH_INSTALL_PATH=$HOME/build/gmsh-install-git ..
```

Metis and boost libraries are picked up automatically.

## Examples
The example models are located in examples directory.
Check out the Wiki of the project to get a handle on the usage.
