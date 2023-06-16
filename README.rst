
Unnamed Project
===============

Building
--------

Use CMake to create the build instructions.
To build the applications you will need CMake >= 3.18, and a compiler toolchain.
To build the Python bindings you will need in addition to a CMake and a compiler toolchain SWIG and Python with development headers.
The following instructions are for building from the command line.
Using a CMake GUI application will be different but the process is essentially the same.

Setup
+++++

Get the code, should be a clone of a repository somewhere but maybe an unzip command??

::

  git clone <some repository somewhere>

or::

  unzip unnamed-project.zip

Create a directory for building the project::

  mkdir build-unnamed-project-release

Configure with CMake
++++++++++++++++++++

::

  cd build-unnamed-project-release
  cmake -DCMAKE_BUILD_TYPE=Release ../unnamed-project

Build
+++++

Execute the compiler toolchains build command for the CMake build generator the project was configured with.
For the Makefile generator on Linux or macOS this would be::

  make

Usage
+++++

The applications are meant to be run in-situ.
Many of the applicaitons require additional inputs from the source directory to be run from the build directory.
Copy of any files that application requires, alternatively run the applications from the source directory.
The Pennati_02 model requires no input files and can be run as is.

To use the Python bindings make the bindings build directory available to your Python virtual environment.
The easiest (but not only, or necessarily the best) way is to set the environment variable PYTHONPATH::

  export PYTHONPATH=<absolute-path-build-directory>/src/bindings/python

Then, start a Python interactive session::

  python

The pennati_02 application can then be run with the following python commands::

  >>> from placenta import pennati_02
  >>> pennati_02.run()

Notes
-----

C, V, Rvec -- variables from the Matlab workspace placenta.mat saved in txt format.
templinks & tempnode_elts -- the variables links and node_elts from the Matlab workspace placenta.mat saved in txt format, *flattened* column wise before exporting!!

The above 5 files with '2' at the end: same, but from placenta2.mat
The above 5 files with '3' at the end: same, but from placenta3.mat
The above 5 files with '3B' at the end: same, but from placenta3B.mat
The above 5 files with '4' at the end: same, but from placenta4.mat

plac_flow.c: equivalent of the Matlab code of the same name. The C version proved to be considerably faster.

pennati_02.c: Pennati model *without* the O2 layer, implemented in C. Same as the Python and Matlab codes of the same name, except doesn't bother with the O2 concentrations, only hemodynamics.

combie.c: Pennati combined with the detailed placenta model. Hemodynamics only, no O2. Equivalent to Matlab code of the same name, but faster. Starts to evolve from t=0.

continue.c: continues the evolution started by combie.c from the last point saved.

plot_combie_fetal/placenta.c: codes that process the data saved by combie & continue and prepare vectors for plotting. These are written out and can be plotted with the help of a code in teh Matlab folder. For placenta variables, this is much faster than processing in Matlab.

Only the pennati_02 has been ported to be available from Python.
To create the other applications replicate what has beenn done for that module.
