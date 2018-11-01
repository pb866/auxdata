auxdata
=======

> Notice: 
> `auxdata` is now deprecated and will slowly be replaced by individual repositories for each Julia package.

Auxiliary data and scripts used mainly for Julia scripting. Some auxiliary python
scripts to derive parameterisations for MCM photolysis processes are provided,
but no longer maintained.

Content
-------

- Folder `py.fcn`: Auxiliary python scripts for the derivation of MCM photolysis processes
  (no longer maintained)
- `gnu.head` in folder `py.dat`: Obsolete header for a gnuplot script (as include file)
- `jl.mod`: Auxiliary Julia modules
  - `fhandle`: Functions for file handling, I/O reading/writing
  - `rddata`: Functions to read _j_ values from TUV output files
  - `fitfcn`: Functions to derive parameterisations for MCM photolysis processes from TUV output data
  - `pltfcn`: Functions to plot MCM photolysis parameterisations derived from TUV output
  - `ChemPhysConst`: Important chemical and physical constants in natural sciences
  - `numerics`: Functions dealing with manipulating numbers
  - `VectorMaths`: Functions to derive the mean and sum of a set of vectors
  - `pyp`: Functions to load and plot data using `PyPlot`


Version history
===============

Version 0.8
-----------
- Removed module `make_plots` for plotting of DSMACCanalysis output (see https://github.com/pb866/DSMACCanalysis.git)
- New module `ChemPhysConst` with important chemical and physical constants in natural sciences

Version 0.7 – 0.7.2
-------------------
- Updated Julia plotting scripts
- Mean/sum of multiple vectors in `VectorMaths`; extrapolation of mean to minimum/maximum of overall data by scaling or offset
- Bug fixes

Version 0.6
-----------
- Initial auxiliary files with obsolete python scripts and Julia scipts for the derivation of MCM photolysis parameterisations and plotting
