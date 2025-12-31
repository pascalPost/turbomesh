# Turbomesh

An open-source mesh generator for the CFD simulation of turbomachinery applications.

Reimplementation of [turbomesh](https://github.com/pascalPost/turbomesh).


## Development Roadmap

- [x] 2D linear tfi
- [x] elliptic block smoothing
- [x] elliptic mesh smoothing (multiple blocks including inter-block boundaries)
- [x] boundary layer meshes
- [ ] automated blocking for turbomachinery (ongoing)
- [ ] gui (ongoing)
- [ ] WASM to run in browser
- [ ] extension to 3D based on meshed 2d cuts
- [ ] possible extension to radial configurations

Planned features:
- [ ] axial turbines
  - [ ] w/o tip clearance
  - [ ] w/ tip clearance
- [ ] axial compressors
- [ ] cascade configurations
- [ ] rotational configurations
- [ ] high order meshes
  - [ ] equidistant
  - [ ] other distributions
- output formats
  - vtk
    - [ ] lagacy vtk
    - [ ] modern vtk formats
  - cgns
    - [x] structured
    - [ ] unstructured


## Developers

### Dependencies

- Zig (0.15.2)
- PETSc (and all its dependencies, e.g. mpi), UMFPACK, MUMPS
- CGNS

### PETSc

For an updated version run:
```
zig translate-c -lc -I/opt/petsc/linux-c-opt/include /opt/petsc/linux-c-opt/include/petsc.h
```

TODO: put this into the build (e.g. as an option).

```
pub extern fn get_petsc_comm_self() callconv(.c) MPI_Comm;
```
