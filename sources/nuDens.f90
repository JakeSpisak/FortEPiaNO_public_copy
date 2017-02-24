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
!	call initData
!	call initmc
!	call initminim(.true.)
	
	call solver
!	call printVec(nuDensVec,flavorNumber**2)
	
	call addToLog("Finished. closing log file.")
	call closeLogFile
end program nuDens
