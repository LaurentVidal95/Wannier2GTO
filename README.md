Here are some conventions of the code:

 - The fundamental building blocs of the code are gaussian-polynomial functions.
   Their polynomial part is described by its decomposition in the basis of 3 dimensional monoms, hence as a list of powers `(n_x, n_y, n_z)` and corresponding  coefficients. Their are denoted `Χ` throughout the code. A family of `Χ` is indexed by greek letters.
 - A linear combination of Gaussian Polynomials is called an AO and denoted Φ. A family of AOs is indexed with latin letters.
 - Wannier functions are to be projected on a basis of AOs, hence is ultimately described as a list of coefficient and corresponding Φs.
 - All scalar products are by default associated to the L2 norm. However any H^s norm can be chosen instead.
 - All functions (GTOs, AOs and Wannier) are stored by their Fourier coefficients in the plane-wave basis associated to scf computations that generated the Bloch waves. All operations utlimately acts on these arrays. If needed (e.g. for ploting) switching between real and Fourier space is done via FFT routines already implemented in DFTK (`G_to_r` and `r_to_G`)
