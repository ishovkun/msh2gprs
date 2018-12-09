# MSH2GRPS
msh2grps is a tool to easily create models for AD-GRPS (preprocessor).
It currently supports mesh files created with GMsh tool.
The properties are be assigned with a choice of JSON or YAML formats human-readable.

## Features
- Supports arbitrary keywords
- Properties are be evaluated as expressions as opposed to simply numbers
- Supports discrete and embedded fractures
- Currently only TPFA is available but MPFA is planned out for the future releases

## Build
msh2gprs requires a C++-17-compatible compiler and minimum CMake 3.7.
There is also a Boost optional dependecy (using boost improves the performance
a lot).

To build mshgprs use the following commands.
```
git clone https://github.com/ishovkun/msh2gprs
cd msh2gprs
mkdir build; cd build
cmake ..
```
CMake will automatically detect whether boost is available and use it.
If Boost is not available, CMake will stick with using a custom library
for 256-bit integers (used for hashing).

## Examples
The example models are located in examples/ directory.
