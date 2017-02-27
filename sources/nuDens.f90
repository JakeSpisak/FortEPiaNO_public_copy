program nuDens
	use precision
	use ndConfig
	use ndErrors
	use ndEquations
!	use omp_lib
	implicit none

	call initConfig
	
	J_func   => J_funcFull 
	K_func   => K_funcFull 
	Jprime   => JprimeFull 
	Kprime   => KprimeFull 
	Y_func   => Y_funcFull 
	G12_func => G12_funcFull
	electronDensity => electronDensityFull
	dme2_electron => dme2_electronFull
	
	call init_interp_jkyg12
	call loadElDensity
	call dme2_e_load

	call solver
	
	call addToLog("Finished. closing log file.")
	call closeLogFile
end program nuDens
