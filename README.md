# turbomesh

An automated mesher for turbomachinery applications (written in rust).

Development roadmap:
- [x] 2D linear tfi
- [ ] elliptic block smoothing
- [ ] elliptic mesh smoothing (multiple blocks including inter-block boundaries)
- [ ] automated blocking for turbomachinery
- [ ] extension to 3D based on meshed 2d cuts
- [ ] possible extension to radial configurations
- [ ] gui


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
    - [x] lagacy vtk
    - [ ] modern vtk formats
  - cgns
    - [ ] structured
    - [ ] unstructured