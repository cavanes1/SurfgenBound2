# SurfgenBound 2
This program was used for and is cited by [this research paper](https://doi.org/10.1063/5.0214635) published in the Journal of Chemical Physics.
This program generates a least-squares polynomial fit to _ab initio_ energies, energy gradients, and derivative couplings _via_ diabatization by _ansatz_ for molecules with a bounded potential energy surface.
The generated matrix surface can be analyzed using [these tools](https://github.com/cavanes1/PES-analysis).
This is an improved version of [SurfgenBound](https://github.com/YifanShenSZ/SurfGenBound).
The provided makefile is designed for [ARCH](https://www.arch.jhu.edu/).

### A couple of important notes
* Energies of higher-energy states not used in the fit can be included in energy.all without issue.
* This program does not work for only 1 electronic state (dat.x cannot read Hd.CheckPoint). A quick workaround is to fit a two-by-two matrix where the (equal) diagonal components are the electronic energy and the off-diagonals are zero.

## Python tools
| File  | Description |
| ------------- | ------------- |
| surfextractor.pl | Extracts the necessary data for Surfgen |
| points.py | Identifies quasidegeneracies from energy.all to generate points.in |
| autosort.py | Reorders the geometries used by Surfgen according to an input list of names. After running, one must run points.py. |

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

### Errors
Tip: You can examine the output surface quickly with `grep Epoch surfgen.log`


| Symbol  | Variable | Meaning |
| ------------- | ------------- | ------------- |
| d[E] | rmsee | RMS (root mean square) absolute energy error |
| <d[E]> | mue | Average absolute value of absolute energy error |
| d[g] | rmseg | RMS or relative errors of 2-norms of absolute vector difference |
| <d[g]> | mueg | Average of 2-norms of absolute vector difference errors |
| d[hij] | rmsec | RMS of 2-norms of relative vector difference errors |

## Acknowledgments
The initial commit version was provided by Yafu Guan.
