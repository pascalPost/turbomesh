#!/bin/bash

gfortran -std=legacy -shared -O3 -fPIC -fdefault-real-8 -o libdierckx.so *.f

# TODO move this into the build system
# TODO remove the *.so from the repo!!
