program icsteriles
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
	
!	if (doLLHEn)      call minimizer(0,1) !Energy bins only
!	if (doLLHCosTh)   call minimizer(1,1) !CosTh bins only
!	if (doLLHEnCosTH) call minimizer(0,2) !Energy-cosTh bins
	
	call addToLog("finished")
	call closeLogFile
end program icsteriles
