module ndCosmology
	use precision
	use constants
	use ndErrors
	use ndInteractions
	implicit none
	
	contains
	
	function radDensity(x,y,z)
		real(dl) :: radDensity,x,y,z
		
		radDensity = photonDensity(z) + electronDensity(x,z) 
	end function
	
	function photonDensity(z)
		real(dl) :: photonDensity, z
		
		photonDensity = PISQD15 * z**4
	end function photonDensity
	
	function integr_rho_e(vec,y)
		real(dl) :: integr_rho_e, y,a
		real(dl), dimension(2) :: vec
		a = E_k_m(y,vec(1))
		integr_rho_e = y*y*a*fermiDirac_massless(a,vec(2))
	end function integr_rho_e
	
	function electronDensity(x,z)
		real(dl) :: electronDensity, x,z, rombint_obj
		real(dl), dimension(2) :: vec
		external rombint_obj
		
		vec(1)=x
		vec(2)=z
		electronDensity = rombint_obj(vec, integr_rho_e, 0., 60., 1d-3)
		electronDensity = electronDensity / PISQD2
	end function electronDensity

end module
