# tardigrade

An implementation of micromorphic continuum mechanics that uses functions and 
values defined in the micromorphic_element repository available using SSH 
from
> `git@bitbucket.org:NateAM2/micromorphic_element.git`

TODO: UPDATE GIT LOCATIONS

All required repositories should be located in a common directory which must 
be identified in the Makefile as the variable `MICROMORPHIC_DIR`. All the 
default repostitory names should be used except for Voro++ which should be
located in the directory `voro++` in `MICROMORPHIC_DIR`.

In order to generate the material models which are used by MicromorphicMaterial 
one should run

> `make -f makefile_libmicromat`

> `make -f makefile_libmicrobalance`

 in `micromorphic_element/src/cpp`
which will generate the shared libraries `libmicromat.so.1` and `libmicrobalance.so.1` 
which will need to be linked to by setting the environment variable

> `LD_LIBRARY_PATH=/path/to/libs:$LD_LIBRARY_PATH`

The documentation for this implementation detailing the theory is available 
in the micromorphic_element repository as the balance equations and other 
materials are located there.

Questions should be addressed to Nathan Miller who can be reached at 
nathanm@lanl.gov.

=====

"Fork tardigrade" to create a new MOOSE-based application.

For more information see: [http://mooseframework.org/create-an-app/](http://mooseframework.org/create-an-app/)
