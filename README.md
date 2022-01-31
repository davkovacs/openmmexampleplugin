OpenMM evaluator for Atomic Cluster Expansion models
====================================================

This project implements an evaluator of [ACE1](https://acesuit.github.io/ACE1docs.jl/dev/) models in [OpenMM](https://openmm.org). It includes a CPU implementation, serialization support and a Python API. 

Building The Plugin
===================

0. Make sure you have ACE1.jl installed and running by following the instructions [HERE](https://acesuit.github.io/ACE1docs.jl/dev/gettingstarted/installation/)

If you don't have OpenMM:

1. Create a new conda environment `conda create -n openmm python=3.9` , `conda activate openmm`

2. Install OpenMM: `conda install -c conda-forge openmm`

From here the same if you have OpenMM already. 

3. This project uses [CMake](http://www.cmake.org) for its build system. If you don't have CMake installed add it by `conda install -c anaconda cmake`.

4. This project uses swig to generate the Python interface. If you don't have swig installed, and want to create the python interface run `conda install -c anaconda swig`

5. Add the julia folder to the LD_LIBRARY_PATH: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/julia-1.7.1/lib/` 

6. Than clone this repo, and create a dictionary in which to build the plugin (inside this repo):

```
mkdir build
cd build
```

7. Run CMake by `ccmake ..`

8. Press c ("Configure").

9. Set OPENMM_DIR to point to the directory where OpenMM is installed.  This is needed to locate
the OpenMM header files and libraries. For example the directory of your conde environment eg. `~/miniconda3/envs/openmm`

10. Set CMAKE_INSTALL_PREFIX to the directory where the plugin should be installed.  Usually,
this will be the same as OPENMM_DIR, so the plugin will be added to your OpenMM installation.

11. Set JULIA_DIR to the directory of your julia installation eg. `~/Applications/julia-1.7.1`

12. Press c ("Configure") again if necessary, then press g ("Generate").

13. Use the build system you selected to build and install the plugin.  For example, if you
selected Unix Makefiles, type `make install`.

Python API
==========

In the `build` directory run `make PythonInstall`

Details:

OpenMM uses [SWIG](http://www.swig.org) to generate its Python API.  SWIG takes an "interface
file", which is essentially a C++ header file with some extra annotations added, as its input.
It then generates a Python extension module exposing the C++ API in Python.

When building OpenMM's Python API, the interface file is generated automatically from the C++
API.  That guarantees the C++ and Python APIs are always synchronized with each other and avoids
the potential bugs that would come from having duplicate definitions. 

To build and install the Python API, build the "PythonInstall" target, for example by typing
"make PythonInstall".  (If you are installing into the system Python, you may need to use sudo.)
This runs SWIG to generate the C++ and Python files for the extension module
(ExamplePluginWrapper.cpp and exampleplugin.py), then runs a setup.py script to build and
install the module. 


License
=======

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2014-2021 Stanford University and the Authors.

Authors: David P. Kovacs, Peter Eastman

Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.

