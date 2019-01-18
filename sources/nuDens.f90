program nuDens
	use precision
	use ndConfig
	use ndErrors
	use ndEquations
!	use omp_lib
	implicit none

	call initConfig

!	call init_interp_FD		!interpolation slower than function!
!	call init_interp_d123	!interpolation slower than function!
	call init_interp_jkyg12
	call init_interp_ElDensity
	call init_interp_dme2_e
	call allocate_interpNuDens

	call solver

	print *,"Neff=",Neff_from_rho_z(nuDensVec(ntot))

	call addToLog("Finished. closing log file.")
	call closeLogFile
end program nuDens
