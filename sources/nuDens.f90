program nuDens
	use precision
	use ndConfig
	use ndErrors
	use ndEquations
!	use omp_lib
	implicit none

	call initConfig

	call init_interp_jkyg12
	call init_interp_d123
	call loadElDensity
	call dme2_e_load

	call solver
	
	call addToLog("Finished. closing log file.")
	call closeLogFile
end program nuDens
