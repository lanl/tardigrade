# tardigrade

C20048 Tardigrade

© 2021. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.

An implementation of micromorphic continuum mechanics that uses functions and 
values defined in the tardigrade-micromorphic-element repository available using SSH 
from
> `git@github.com:lanl/tardigrade-micromorphic-element.git`

All required repositories should be located in a common directory which must 
be identified in the Makefile as the variable `MICROMORPHIC_DIR`. All the 
default repostitory names should be used except for Voro++ which should be
located in the directory `voro++` in `MICROMORPHIC_DIR`.

The user also must define the location of the C++ library for the compiler
used to build libmicromat and libmicrobalance in Tardigrade's makefile. This
is done by specifying the path in `MICROMORPHIC_COMPILER_PATH`.

In order to generate the material models which are used by MicromorphicMaterial 
one should run

> `make -f makefile_libmicromat`

> `make -f makefile_libmicrobalance`

 in `micromorphic_element/src/cpp`
which will generate the shared libraries `libmicromat.so.1` and `libmicrobalance.so.1` 
which will need to be linked to by setting the environment variable

> `LD_LIBRARY_PATH=/path/to/libs:$LD_LIBRARY_PATH`

The user is expected to have built the xdmf library in the root directory of tardigrade
( i.e. from the location of this Readme the xmdf directory is ../xmdf ) using the typical
nomenclature of, "build," as the build directory name ( ../xmdf/build ). xmdf requires
the boost header files ( which should be referenced by the environment variable BOOST_ROOT )
and libxml2 which is included in the moose conda environment. Xmdf should be built against
those libraries. The location of libxml2 needs to be updated in the Makefile along with 
the root micromorphic directory.

The documentation for this implementation detailing the theory is available 
in the micromorphic_element repository as the balance equations and other 
materials are located there.

Questions should be addressed to Nathan Miller who can be reached at 
nathanm@lanl.gov.

=====

"Fork tardigrade" to create a new MOOSE-based application.

For more information see: [http://mooseframework.org/create-an-app/](http://mooseframework.org/create-an-app/)
