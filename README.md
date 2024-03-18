# SurfgenBound 2
Property of the Yarkony Group at Johns Hopkins University.

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

### fit.in parameters

| Keyword  | Meaning |
| ------------- | ------------- |
| npoints | Number of points to fit |
| enfDiab | Point number of reference point |
| epmax | Maximum number of epochs |
| w_energy | Energy weights |
| w_grad | Gradient weights |
| w_fij | Derivative coupling weights |
| gradcutoff | when <=0, value of gradcutoff is used instead |
| cpcutoff |  |
| deggrdbinding |  |
| deg_cap |  |
| lambda | Has to do with regularization |
| eshift | Uniform shift on ab initio energies |
| energyT | Threshold energy above which gradients/energies will be scaled down |
| highEScale |  The scaling factor for these high energy data |
| nrmediff |  |
| ediffcutoff |  |
| fixref | Fix the coeffs of zeroth and first order if .true. |
| natoms | Number of atoms |
| nstates | Number of electronic states |

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
