Julia code to produce and compress a pz-like wannier function for Graphene. The basis functions for the compression
are linear combinations of gaussian-polynomials, whose centers are located somewhere between the center of the wannier
function and the middle of a pi-bond. The compression is adapted, specifically for graphene wannier functions with D3 symmetry,
from the general method of [WannierCompression][^1].


# Requirements: 
Julia 1.8. and above.

# Installing all dependancies
Open a Julia shell with `julia --project` in your local copy of this repository and call
```
using Pkg; Pkg.instantiate(".")
``` 
to install all the needed dependancies.

# Usage

All that concerns the compression routine has not been designed to be launched by the user.
It is not optimized and has only been run once to produce the compressed wannier data.
For the sake of completeness, the file `tests/launch_compression.jl` launches the
compression routine from a precomputed wannier function stored in `wannier_pz.json`.

Only the part of the code that is related to the construction of tight-binding quantities from
the compressed wannier is to be used. # This part is still to be coded..

# Contact
This is research code, not user-friendly, actively maintened, extremely robust or optimized.
It is desinged to be called only once to produce a set of coefficient. 
If you have questions contact me at: Laurent(dot)vidal(at)enpc(dot)fr

[^1]: Athmane, B., Eric, C., Paul, C., Shiang, F., & Efthimios, K. (2017). Compression of Wannier functions into Gaussian-type orbitals. (arXiv:1712.02996v1)
