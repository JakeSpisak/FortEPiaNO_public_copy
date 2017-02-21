module ndCosmology
	use precision
	use constants
	use ndErrors
	use ndInteractions
	implicit none
	
	type nuDensArgs
		real(dl) :: x,z
		integer iFl
	end type nuDensArgs
	
	contains
	
	function radDensity(x,y,z)
		real(dl) :: radDensity,x,y,z
		
		radDensity = photonDensity(z) + &
			electronDensity(x,z) + &
			allNuDensity(x,z)
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

	function integr_rho_nu(vec,y)
		real(dl) :: integr_rho_nu, y,a
		type(nuDensArgs) :: vec
		type(cmplxMatNN) :: mat
		
		mat = interp_nuDensIJ(y, vec%iFl, vec%iFl)
		
		integr_rho_nu = y*y*y*mat%re(vec%iFl, vec%iFl)
	end function integr_rho_nu
	
	function nuDensity(x,z, iFl)
		real(dl) :: nuDensity, x,z, rombint_obj
		type(nuDensArgs) :: vec
		integer :: iFl
		external rombint_obj
		
		vec%x=x
		vec%z=z
		vec%iFl=iFl
		nuDensity = rombint_obj(vec, integr_rho_nu, y_min, y_max, 1d-3) / PISQD2
	end function nuDensity
	
	function allNuDensity(x,z)
		real(dl) :: allNuDensity, x,z
		integer :: ix
		
		allNuDensity = 0.d0
		do ix=1, flavorNumber
			allNuDensity = allNuDensity + nuDensity(x,z,ix)
		end do
		allNuDensity = allNuDensity / PISQD2
	end function allNuDensity

end module
