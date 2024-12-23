#!/bin/bash

gfortran -std=legacy -shared -O3 -fPIC -fdefault-real-8 -o libdierckx.so *.f
