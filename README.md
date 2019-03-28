### Redshift Reduction

This is a suite of codes written in IDL and Python that have been used to determine redshifts and optical classifications for ~3000 spectra taking with the VIMOS instrument on the VLT.

The bulk of the IDL code that handles the redshift determination comes from [DEIMOS](https://github.com/MCTwo/DEIMOS) and has been optimized for VIMOS spectra.

[`idl_wrapper.pro`](idl_wrapper.pro) can be modified to read in a set of spectra, fit it against models of your choosing (currently PCA reconstructions using SDSS 2012 eigenspectra of galaxies, AGN, stars, and cataclysmic variables), and output redshifts, fit statistics, and best fitting models.
