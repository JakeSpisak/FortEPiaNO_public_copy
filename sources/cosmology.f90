module ndCosmology
	use precision
	use constants
	use utilities
	use ndErrors
	use ndInteractions
	use bspline_module
	implicit none
	
	logical :: loadedEDens = .false.
	type(bspline_2d) :: elDens
	
	type nuDensArgs
		real(dl) :: x,z
		integer iFl
	end type nuDensArgs

	procedure (funcXY), pointer :: electronDensity => null ()
	
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
	
	function electronDensityFull(x,z)
		real(dl) :: electronDensityFull, x,z, rombint_obj
		real(dl), dimension(2) :: vec
		external rombint_obj
		
!		if (loadedEDens) then
!			electronDensityFull = elDens%Value(log10(x),z)
!		else
			vec(1)=x
			vec(2)=z
		
			electronDensityFull = rombint_obj(vec, integr_rho_e, 0.d0, 60.d0, 1d-3, maxiter)
			electronDensityFull = electronDensityFull / PISQD2
!		end if
	end function electronDensityFull
	
	function electronDensityInterp(x,z)
		real(dl) :: electronDensityInterp, x,z
		integer :: iflag
		
		call elDens%evaluate(x,z,0,0,electronDensityInterp,iflag)
	end function electronDensityInterp

	subroutine loadElDensity
		real(dl), dimension(:,:), allocatable :: ed_vec
		integer :: ix, iz, iflag
		real(dl) :: x,z
		
		call addToLog("[cosmo] Initializing interpolation for electron density...")
		allocate(ed_vec(interp_nx,interp_nz))
		do ix=1, interp_nx
			do iz=1, interp_nz
				ed_vec(ix,iz) = electronDensityFull(interp_xvec(ix),interp_zvec(iz))
			end do
		end do
		call elDens%initialize(interp_xvec,interp_zvec,ed_vec,4,4,iflag)
		electronDensity => electronDensityInterp
		
		call random_seed()
		call random_number(x)!=0.13151
		call random_number(z)!=0.13151
		x=(x_fin-x_in)*x + x_in
		z=0.4d0*z + z_in
		write(*,"(' [cosmo] test electronDensityInterp in ',*(E12.5))") x,z
		write(*,"(' [cosmo] comparison (true vs interp): ',*(E17.10))") electronDensityFull(x,z), electronDensity(x,z)
		
		loadedEDens = .true.
		call addToLog("[cosmo] ...done!")
	end subroutine loadElDensity
	
	function integr_rho_nu(vec,y)
		real(dl) :: integr_rho_nu, y,a
		type(nuDensArgs) :: vec
		type(cmplxMatNN) :: mat
		
		call allocateCmplxMat(mat)
		
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
		nuDensity = rombint_obj(vec, integr_rho_nu, y_min, y_max, 1d-3, maxiter) / PISQD2
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
