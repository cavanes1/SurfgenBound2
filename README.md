# SurfgenBound 2
This is an updated version of SurfgenBound. This program generates a polynomial fit to ab initio points via diabatization by ansatz. The initial version was provided by Yafu Guan.

The provided makefile runs succesfully on ARCH.

## Notes
* Energies of higher-energy states not used in the fit can be included in energy.all without issue
* Does not work for only 1 electronic state (dat.x cannot read Hd.CheckPoint)

## Errors

| Symbol  | Variable | Meaning |
| ------------- | ------------- | ------------- |
| d[E] | rmsee | RMS absolute energy error |
| <d[E]> | mue | average absolute value of absolute energy error |
| d[g] | rmseg | RMS or relative errors of 2-norms of absolute vector difference |
| <d[g]> | mueg | Average of 2-norms of absolute vector difference errors |
| d[hij] | rmsec | RMS of 2-norms of relative vector difference errors |
