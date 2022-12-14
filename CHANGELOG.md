# [FortEPiaNO](https://bitbucket.org/ahep_cosmo/fortepiano_public/)
FORTran-Evolved PrImordiAl Neutrino Oscillations  
by S. Gariazzo (gariazzo@ific.uv.es), P.F. de Salas (pablo.fernandez@fysik.su.se) and S. Pastor (pastor@ific.uv.es)

## v1.0.0 (19/03/2021)
* first public version

## v1.1.0 (19/05/2021)
* neutrino-neutrino collision terms take into account the presence of sterile neutrinos
* muons are now disabled by default, enable their contributions compiling with DO_MUONS=1
* first (very incomplete) python wrapper for the fortran code
* possibility to store intermediate quantities for later checks (the overall normalization, Y, Ydot, Heff, the commutator, collision terms as functions of x)
* Makefile now defaults to gfortran
* several other internal improvements and restructuring

## v1.2.0 (14/10/2021)
* use FULL_F_AB=1 when compiling in order to use non-standard interactions between neutrinos and electrons (nsi_GL_ij, nsi_GR_ij input parameters)
* updated references in README
* fixes to previous problems
