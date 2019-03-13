module ndCosmology
	use precision
	use constants
	use utilities
	use ndErrors
	use ndInteractions
	use linear_interpolation_module
	implicit none

	type(linear_interp_2d) :: elDens, muDens

	contains

	function radDensity(x,z)
		real(dl) :: radDensity
		real(dl), intent(in) :: x,z

		radDensity = photonDensity(z) + &
			electronDensity(x,z) + &
			muonDensity(x,z) + &
			allNuDensity(z)
	end function

	elemental function photonDensity(z)
		real(dl) :: photonDensity
		real(dl), intent(in) :: z

		photonDensity = PISQD15 * z**4
	end function photonDensity

	function integr_rho_e(vec,y)
		real(dl) :: integr_rho_e, Emk
		real(dl), intent(in) :: y
		real(dl), dimension(3), intent(in) :: vec
		Emk = Ebare_i_dme(y, vec(1), vec(3))
		integr_rho_e = y*y*Emk * fermiDirac(Emk/vec(2))
	end function integr_rho_e

	function electronDensityFull(x,z)!electron + positron!
		real(dl) :: electronDensityFull
		real(dl), intent(in) :: x,z
		real(dl), dimension(3) :: vec

		vec(1)=x
		vec(2)=z
		vec(3)=dme2_electronFull(x, 0.d0, z)

		electronDensityFull = rombint_vec(vec, integr_rho_e, fe_l, fe_u, toler_ed, maxiter)
		electronDensityFull = electronDensityFull / PISQD2 !the factor is given by g = 2(elicity) * 2(e+e-)
	end function electronDensityFull

	function electronDensity(x,z)
		real(dl) :: electronDensity
		real(dl), intent(in) :: x,z

		call elDens%evaluate(x,z,electronDensity)
	end function electronDensity

	function integr_rho_mu(vec,y)
		real(dl) :: integr_rho_mu, Emk
		real(dl), intent(in) :: y
		real(dl), dimension(3), intent(in) :: vec
		Emk = E_k_m(y, vec(1)*m_mu_o_m_e)
		integr_rho_mu = y*y*Emk * fermiDirac(Emk/vec(2))
	end function integr_rho_mu

	function muonDensityFull(x,z)!electron + positron!
		real(dl) :: muonDensityFull
		real(dl), intent(in) :: x,z
		real(dl), dimension(3) :: vec

		vec(1)=x
		vec(2)=z
		vec(3)=0.

		if (x.gt.0.5d0) then
			muonDensityFull = 0.d0
		else
			muonDensityFull = rombint_vec(vec, integr_rho_mu, fe_l, fe_u, toler_ed, maxiter)
			muonDensityFull = muonDensityFull / PISQD2 !the factor is given by g = 2(elicity) * 2(e+e-)
		end if
	end function muonDensityFull

	function muonDensity(x,z)
		real(dl) :: muonDensity
		real(dl), intent(in) :: x,z

		call muDens%evaluate(x,z,muonDensity)
	end function muonDensity

	subroutine init_interp_ElDensity
		real(dl), dimension(:,:), allocatable :: ed_vec, md_vec
		integer :: ix, iz, iflag
		real(dl) :: x,z, t1,t2
		real(8) :: timer1

		call addToLog("[cosmo] Initializing interpolation for electron and muon density...")
		allocate(ed_vec(interp_nx, interp_nz), md_vec(interp_nx, interp_nz))
		!$omp parallel do default(shared) private(ix, iz) schedule(dynamic)
		do ix=1, interp_nx
			do iz=1, interp_nz
				ed_vec(ix,iz) = electronDensityFull(interp_xvec(ix), interp_zvec(iz))
				md_vec(ix,iz) = muonDensityFull(interp_xvec(ix), interp_zvec(iz))
			end do
		end do
		!$omp end parallel do
		call elDens%initialize(interp_xvec, interp_zvec, ed_vec, iflag)!linear
		call muDens%initialize(interp_xvec, interp_zvec, md_vec, iflag)!linear

		call random_seed()
		call random_number(x)
		call random_number(z)
		x=(x_fin-x_in)*x + x_in
		z=0.4d0*z + z_in
		write(*,"(' [cosmo] test electronDensityInterp in ',*(E12.5))") x,z
		t1 = electronDensityFull(x,z)
		t2 = electronDensity(x,z)
		write(*,"(' [cosmo] comparison electron density (true vs interp): ',*(E17.10))") t1,t2
		t1 = muonDensityFull(x,z)
		t2 = muonDensity(x,z)
		write(*,"(' [cosmo] comparison muon density (true vs interp): ',*(E17.10))") t1,t2

		if (timing_tests) then
			call tic(timer1)
			write (*,*) "[interactions] now doing some timing..."
			call tic(timer1)
			do ix=1, 1000000
				call random_number(x)
				call random_number(z)
				x=(x_fin-x_in)*x + x_in
				z=0.4d0*z + z_in
				t1 = electronDensity(x,z)
			end do
			call toc(timer1, "<interpolated>")

			call tic(timer1)
			do ix=1, 1000000
				call random_number(x)
				call random_number(z)
				x=(x_fin-x_in)*x + x_in
				z=0.4d0*z + z_in
				t1 = electronDensityFull(x,z)
			end do
			call toc(timer1, "<full>")
		end if

		call addToLog("[cosmo] ...done!")
	end subroutine init_interp_ElDensity

	function nuDensityLin(z, iFl)
		real(dl) :: nuDensityLin, y
		real(dl), intent(in) :: z
		integer, intent(in) :: iFl
		integer :: ix

		do ix=1, Ny
			y = y_arr(ix)
			fy_arr(ix) = y*y*y * nuDensMatVecFD(ix)%re(iFl, iFl)
		end do
		nuDensityLin = integral_linearized_1d(Ny, dy_arr, fy_arr) / PISQ
	end function nuDensityLin

	function allNuDensity(z)
		real(dl) :: allNuDensity
		real(dl), intent(in) :: z
		integer :: ix

		allNuDensity = 0.d0
		do ix=1, flavorNumber
			allNuDensity = allNuDensity + nuDensityLin(z,ix)*nuFactor(ix)
		end do
		allNuDensity = allNuDensity
	end function allNuDensity

	function nuDensityLinEq(z)
		real(dl) :: nuDensityLinEq, y
		real(dl), intent(in) :: z
		integer :: ix

		do ix=1, Ny
			y = y_arr(ix)
			fy_arr(ix) = y*y*y * fermiDirac(y)
		end do
		nuDensityLinEq = integral_linearized_1d(Ny, dy_arr, fy_arr) / PISQ
	end function nuDensityLinEq
end module
