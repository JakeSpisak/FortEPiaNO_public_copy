program icsteriles
	use precision
	use ndConfig
!	use icsData
!	use icsFlux
!	use icsMC
!	use icsMinim
!	use icsLLH
	
	implicit none

	call initConfig
!	call initData
!	call initmc
!	call initminim(.true.)
	
!	if (doLLHEn)      call minimizer(0,1) !Energy bins only
!	if (doLLHCosTh)   call minimizer(1,1) !CosTh bins only
!	if (doLLHEnCosTH) call minimizer(0,2) !Energy-cosTh bins
	print *,"finished"
end program icsteriles
