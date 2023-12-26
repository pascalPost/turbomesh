# turbomesh

An automated mesher for turbomachinery applications (written in rust).

Development roadmap:

- [x] 2D linear tfi
- [x] elliptic block smoothing
- [x] elliptic mesh smoothing (multiple blocks including inter-block boundaries)
- [ ] boundary layer meshes
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
        - [ ] lagacy vtk
        - [ ] modern vtk formats
    - cgns
        - [x] structured
        - [ ] unstructured

## Debugging w/ lldb

To enhance the debugging experience, the `lldb/formatter.py` script can be used
to provide a better view of the `ndarray` objects.

- [ ] remove the unnecessary parts of `lldb/formatter.py`

### Debugging with VSCode

1. Install
   the [CodeLLDB](https://marketplace.visualstudio.com/items?itemName=vadimcn.vscode-lldb)
   extension
2. Add the following to your `settings.json`:

```json
  "lldb.launch.initCommands": [
"command script import ${workspaceFolder}/lldb/formatter.py"
],
```

### Debugging with RustRover (bundled renderes)

When using the bundled renderers, the lldb pretty printers are located
in `rustrover/plugins/intellij-rust/prettyPrinters` and can be modified,
see https://youtrack.jetbrains.com/issue/RUST-10219/LLDB-debugger-add-pretty-printer-for-NaiveDate.

You may cop all (except the `__lldb_init_module` function)
of `lldb/formatter.py` into `lldb_formatters/__init__.py` and
add the following at the end of the `__lldb_init_module` function:

```python
initialize_category(debugger)
```

- [ ] enhance the formatter for RustRover