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

	function totalRadiationDensity(x,z)
		real(dl) :: totalRadiationDensity
		real(dl), intent(in) :: x,z

		totalRadiationDensity = photonDensity(z) + &
			electronDensity(x,z) + &
			muonDensity(x,z) + &
			allNuDensity()
	end function

	elemental function photonDensity(z)
		real(dl) :: photonDensity
		real(dl), intent(in) :: z

		photonDensity = PISQD15 * z**4
	end function photonDensity

	pure function integr_rho_e(x, z, dme2, y)
		real(dl) :: integr_rho_e, Emk
		real(dl), intent(in) :: x, z, dme2, y
		Emk = Ebare_i_dme(y, x, dme2)
		integr_rho_e = Emk * fermiDirac(Emk/z)
	end function integr_rho_e

	pure function electronDensityFull(x, z)!electron + positron!
		real(dl) :: electronDensityFull, dme2
		real(dl), intent(in) :: x,z
		integer :: i

		dme2 = dme2_electronFull(x, 0.d0, z)

		electronDensityFull = 0.d0
		do i=1, N_opt_y
			electronDensityFull = electronDensityFull &
				+ opt_y_w(i)*integr_rho_e(x, z, dme2, opt_y(i))
		end do
		electronDensityFull = electronDensityFull / PISQD2 !the factor is given by g = 2(elicity) * 2(e+e-)
	end function electronDensityFull

	function electronDensity(x,z)
		real(dl) :: electronDensity
		real(dl), intent(in) :: x,z

		call elDens%evaluate(x,z,electronDensity)
	end function electronDensity

	pure function integr_rho_mu(x, z, y)
		real(dl) :: integr_rho_mu, Emk
		real(dl), intent(in) :: x, z, y
		Emk = E_k_m(y, x*m_mu_o_m_e)
		integr_rho_mu = Emk * fermiDirac(Emk/z)
	end function integr_rho_mu

	pure function muonDensityFull(x,z)!electron + positron!
		real(dl) :: muonDensityFull
		real(dl), intent(in) :: x,z
		integer :: i

		muonDensityFull = 0.d0
		if (x .lt. x_muon_cut) then
			do i=1, N_opt_y
				muonDensityFull = muonDensityFull &
					+ opt_y_w(i)*integr_rho_mu(x, z, opt_y(i))
			end do
			muonDensityFull = muonDensityFull / PISQD2
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

	function nuDensityNC(i1, i2, reim)
		real(dl) :: nuDensityNC, y
		integer, intent(in) :: i1, i2
		logical, intent(in), optional :: reim
		integer :: ix

		if (present(reim) .and. .not.reim) then
			do ix=1, Ny
				y = y_arr(ix)
				fy_arr(ix) = y*y*y * nuDensMatVecFD(ix)%im(i1, i2)
			end do
		else
			do ix=1, Ny
				y = y_arr(ix)
				fy_arr(ix) = y*y*y * nuDensMatVecFD(ix)%re(i1, i2)
			end do
		end if
		nuDensityNC = integral_NC_1d(Ny, dy_arr, fy_arr) / PISQ
	end function nuDensityNC

	function nuDensityGL(i1, i2, reim)
		real(dl) :: nuDensityGL
		integer, intent(in) :: i1, i2
		logical, intent(in), optional :: reim
		integer :: ix

		if (present(reim) .and. .not.reim) then
			do ix=1, Ny
				fy_arr(ix) = nuDensMatVecFD(ix)%im(i1, i2)
			end do
		else
			do ix=1, Ny
				fy_arr(ix) = nuDensMatVecFD(ix)%re(i1, i2)
			end do
		end if
		nuDensityGL = integral_GL_1d(w_gl_arr, fy_arr) / PISQ
	end function nuDensityGL

	function allNuDensity()
		use ndInterfaces1
		real(dl) :: allNuDensity
		integer :: ix
		procedure (nuDensity_integrator), pointer :: nuDensityInt

		if (use_gauss_laguerre) then
			nuDensityInt => nuDensityGL
		else
			nuDensityInt => nuDensityNC
		end if

		allNuDensity = 0.d0
		do ix=1, flavorNumber
			allNuDensity = allNuDensity + nuDensityInt(ix, ix)*nuFactor(ix)
		end do
		allNuDensity = allNuDensity
	end function allNuDensity

	function nuDensityEq(w)
		real(dl) :: nuDensityEq, y
		real(dl), intent(in) :: w
		integer :: ix

		if (use_gauss_laguerre) then
			do ix=1, Ny
				fy_arr(ix) = fermiDirac(y_arr(ix)/w)
			end do
			nuDensityEq = integral_GL_1d(w_gl_arr, fy_arr) / PISQ
		else
			do ix=1, Ny
				y = y_arr(ix)
				fy_arr(ix) = y*y*y * fermiDirac(y/w)
			end do
			nuDensityEq = integral_NC_1d(Ny, dy_arr, fy_arr) / PISQ
		end if
	end function nuDensityEq

	pure function integrand_K_uX(u, w, n)
		real(dl) :: integrand_K_uX
		real(dl), intent(in) :: u, w
		integer, intent(in) :: n
		real(dl) :: squ2w2, expsqu2w2
		squ2w2 = sqrt(u**2+w**2)
		expsqu2w2 = exp(squ2w2)
		integrand_K_uX = u**n/(squ2w2*(1.d0+expsqu2w2))
	end function integrand_K_uX

	!check GL method
	function int_K_uX(w, n)
		!computes the integral of u**n/(sqrt(u**2+w**2)(1+exp(sqrt(u**2+w**2))))
		!useful to get the pressure for non-relativistic species
		real(dl) :: int_K_uX
		real(dl), intent(in) :: w
		integer, intent(in) :: n
		integer k, Ninte
		real(dl) :: u_ini,u_fin,du
		Ninte = 10000
		u_ini = 0.d0
		u_fin =500.d0 ! 50.d0 is OK except for muons
		du = u_fin/Ninte
		int_K_uX = 0.d0
		do k = 1,Ninte
			int_K_uX = int_K_uX + integrand_K_uX(u_ini+du*(k-1),w,n)+integrand_K_uX(u_ini+du*k,w,n)
		enddo
		int_K_uX = int_K_uX*du*0.5d0/PISQ
	end function int_K_uX

	function photonEntropy(z) ! comoving
		real(dl) :: photonEntropy
		real(dl), intent(in) :: z
		photonEntropy = 4.d0/45.d0*PISQ*z**3
	end function photonEntropy

	function nonRelativistic_Entropy(x, z, m)
		real(dl) :: nonRelativistic_Entropy
		real(dl), intent(in) :: x, z, m
		nonRelativistic_Entropy = 2.d0*z**3*(4.d0/3.d0*int_K_uX(m*x/z, 4) + (x/z)**2*int_K_uX(m*x/z, 2))
	end function nonRelativistic_Entropy

	function electronEntropy(x,z)
		real(dl) :: electronEntropy
		real(dl), intent(in) :: x,z
		electronEntropy = nonRelativistic_Entropy(x, z, 1.d0)
	end function electronEntropy

	function muonEntropy(x,z)
		real(dl) :: muonEntropy
		real(dl), intent(in) :: x,z
		if (x .lt. x_muon_cut) then
			muonEntropy = nonRelativistic_Entropy(x, z, m_mu_o_m_e)
		else
			muonEntropy = 0.d0
		end if
	end function muonEntropy

end module
