program nuDens
	use precision
	use ndConfig
	use ndErrors
	use ndEquations
!	use icsFlux
!	use icsMC
!	use icsMinim
!	use icsLLH
	
	implicit none

	call openLogFile
	call initConfig
	call init_interp_jkyg12
	
	electronDensity => electronDensityFull
	call loadElDensity
	dme2_electron => dme2_electronFull
	call dme2_e_load

!	call initData
!	call initmc
!	call initminim(.true.)
	
	call solver
	
	call addToLog("Finished. closing log file.")
	call closeLogFile
end program nuDens
