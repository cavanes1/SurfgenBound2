# SurfgenBound 2
This is an updated version of SurfgenBound. This program generates a least-squares polynomial fit to ab initio points via diabatization by ansatz.
The generated surface can be analyzed using [these tools](https://github.com/cavanes1/PES-analysis).
The initial commit version was provided by Yafu Guan.

The provided makefile runs succesfully on ARCH.

## Notes
* Energies of higher-energy states not used in the fit can be included in energy.all without issue
* Does not work for only 1 electronic state (dat.x cannot read Hd.CheckPoint)

## Python tools

### surfextractor.pl
Extracts the necessary data for Surfgen

### points.py
Identifies quasidegeneracies from energy.all to generate points.in

### autosort.py
Reorders the geometries used by Surfgen according to an input list of names.
After running, one must run points.py.

## Input and output

### Required input files
* basis.in
* intcfl
* Hd.CheckPoint.old (if using existing surface as starting point)
* fit.in
* points.in
* energy.all
* geom.all
* cartgrd.drt1.state$.all
* cartgrd.nad.drt1.state$.drt1.state$.all

### Generated output files
* Hd.CheckPoint
* surfgen.log

### Errors
Examine output surface quickly using `grep Epoch surfgen.log`


| Symbol  | Variable | Meaning |
| ------------- | ------------- | ------------- |
| d[E] | rmsee | RMS absolute energy error |
| <d[E]> | mue | average absolute value of absolute energy error |
| d[g] | rmseg | RMS or relative errors of 2-norms of absolute vector difference |
| <d[g]> | mueg | Average of 2-norms of absolute vector difference errors |
| d[hij] | rmsec | RMS of 2-norms of relative vector difference errors |
