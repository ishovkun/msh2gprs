# MSH2GRPS
msh2grps is a tool to easily create models for AD-GRPS (preprocessor).
It currently supports mesh files created with GMsh tool.
The properties are be assigned with a choice of JSON or YAML human-readable formats.

## Features
- Supports arbitrary keywords
- Properties are be evaluated as expressions as opposed to simply numbers
- Supports discrete and embedded fractures
- Currently only TPFA is available but MPFA is planned out for the future releases
- Wells are not yet supported but will be soon (only need to finish output)

## Build

### Gmsh (optional)
The DFEM feature requires Gmsh.
A simple tarball download does not cut it since it sometimes causes conflicts
with the system metis (used for multiscale partitioning).
Threfore, the correct way is to build Gmsh from source.
``` cmake
git clone https://gitlab.onelab.info/gmsh/gmsh ./gmsh_git
cd gmsh_git;
mkdir build; cd build;
cmake -DCMAKE_INSTALL_PREFIX=$HOME/build/gmsh-git-install ..
```
Essentially, you can tick off all flags except ENABLE_BUILD_DYNAMIC .
If that does not work, try ticking on ENABLE_ONELAB.
msh2gprs requires a C++-17-compatible compiler and minimum CMake 3.7
(build was tested on GCC 8.2 and clang 7.0).
There is also a Boost optional dependecy (using boost improves the performance
by a lot).

To build mshgprs use the following commands.
```
git clone --recursive https://github.com/ishovkun/msh2gprs
cd msh2gprs
mkdir build; cd build
cmake ..
```

### Boost (optional)
CMake will automatically detect whether boost is available and use it.
If Boost is not available, CMake will stick with using a custom library
for 256-bit integers (used for hashing by angem library).

### msh2gprs

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

## Examples
The example models are located in examples directory.
Check out the Wiki of the project to get a handle on the usage.
