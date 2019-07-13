module ndCosmology
	use precision
	use constants
	use utilities
	use ndErrors
	use ndInteractions
	use linear_interpolation_module
	implicit none

	real(dl), parameter :: upper = 1.d2

#ifndef NOINTERPOLATION
	type(linear_interp_2d) :: elDens, muDens
#endif

	type nonRelativistic_fermion
		character(len=20) :: fermionName
		logical :: isElectron !if true, compute the electromagnetic corrections to the mass
		real(dl) :: mass_factor !mass of the particle divided by electron mass
		real(dl) :: x_enDens_cut !for x larger than this, do not consider this particle (fix energy density and other things to 0)
		type(linear_interp_2d) :: enDensInterp !information for interpolating the energy density
		contains
		procedure :: initialize => nonRelativistic_initialize !save properties and prepare energy density interpolation
		procedure :: dzodx_terms => nonRelativistic_dzodx_terms !contributions to the dz/dx numerator and denominator
		procedure :: energyDensityFull => nonRelativistic_energyDensity_full !full integral of the energy density
		procedure :: energyDensity => nonRelativistic_energyDensity !interpolated energy density
		procedure :: entropy => nonRelativistic_entropy !entropy
		procedure :: pressure => nonRelativistic_pressure !pressure
	end type nonRelativistic_fermion
	
	integer, parameter :: fermions_number = 2
	type(nonRelativistic_fermion), dimension(fermions_number), target :: fermions
	!define these only for easier reference in updateMatterDensities, output and tests:
	type(nonRelativistic_fermion), pointer :: electrons, muons

	contains

	!total energy density
	function totalRadiationDensity(x,z) result(rho_r)
		real(dl) :: rho_r
		real(dl), intent(in) :: x,z
		integer :: ix

		rho_r = photonDensity(z) &
			+ allNuDensity()
		do ix=1, fermions_number
			rho_r = rho_r + fermions(ix)%energyDensity(x,z)
		end do
		rho_r = rho_r + deltaRhoTot_em(x, z)
	end function

	!photons
	elemental function photonDensity(z)
		real(dl) :: photonDensity
		real(dl), intent(in) :: z

		photonDensity = PISQD15 * z**4
	end function photonDensity

	elemental function photonEntropy(z)
		real(dl) :: photonEntropy
		real(dl), intent(in) :: z
		photonEntropy = four_thirds*PISQD15*z**3
	end function photonEntropy

	!functions for non relativistic particles
	pure function integrand_rho_nonRel(x, z, y)
		real(dl) :: integrand_rho_nonRel, Emk
		real(dl), intent(in) :: x, z, y
		Emk = E_k_m(y, x)
		integrand_rho_nonRel = Emk * fermiDirac(Emk/z)
	end function integrand_rho_nonRel

	pure function integrand_uX_Ek_nonRel(u, w, n)
		real(dl) :: integrand_uX_Ek_nonRel
		real(dl), intent(in) :: u, w
		integer, intent(in) :: n
		real(dl) :: squ2w2, expsqu2w2
		squ2w2 = E_k_m(u, w)
		expsqu2w2 = exp(squ2w2)
		integrand_uX_Ek_nonRel = u**n/(squ2w2*(1.d0+expsqu2w2))
	end function integrand_uX_Ek_nonRel

	pure function nonRelativistic_energyDensity_full(cls, x, z) result (nredf)!fermion + antifermion
		class(nonRelativistic_fermion), intent(in) :: cls
		real(dl) :: nredf
		real(dl), intent(in) :: x,z
		integer :: i

		nredf = 0.d0
		if (x .lt. cls%x_enDens_cut) then
			do i=1, N_opt_y
				nredf = nredf + opt_y_w(i)* integrand_rho_nonRel(x*cls%mass_factor, z, opt_y(i))
			end do
			nredf = nredf / PISQD2
			!the factor is given by g = 2(elicity) * 2(f+\bar f)
		end if
	end function nonRelativistic_energyDensity_full

	function nonRelativistic_energyDensity(cls, x, z) result (ed)
		class(nonRelativistic_fermion) :: cls
		real(dl) :: ed
		real(dl), intent(in) :: x, z

#ifdef NOINTERPOLATION
		ed = nonRelativistic_energyDensity_full(cls, x, z)
#else
		call cls%enDensInterp%evaluate(x, z, ed)
#endif
	end function nonRelativistic_energyDensity

	pure function integrate_uX_Ek_nonRel(w, n)
		!computes the integral of u**n/(sqrt(u**2+w**2)(1+exp(sqrt(u**2+w**2))))
		!useful to get the pressure for non-relativistic species
		real(dl) :: integrate_uX_Ek_nonRel
		real(dl), intent(in) :: w
		integer, intent(in) :: n
		integer :: i
		integrate_uX_Ek_nonRel = 0.d0
		do i = 1, N_opt_xoz
			integrate_uX_Ek_nonRel = integrate_uX_Ek_nonRel &
				+ opt_y_w(i)*integrand_uX_Ek_nonRel(opt_y(i), w, n-2)!n-2 because the weights already take into account y^2
		end do
		integrate_uX_Ek_nonRel = integrate_uX_Ek_nonRel/PISQx2
	end function integrate_uX_Ek_nonRel

	pure function nonRelativistic_pressure(cls, x, z) result(press)
		real(dl) :: press
		class(nonRelativistic_fermion), intent(in) :: cls
		real(dl), intent(in) :: x, z
		real(dl) :: t

		if (x .lt. cls%x_enDens_cut) then
			t = cls%mass_factor*x/z
			press = 4.d0*z**4/3.d0*(integrate_uX_Ek_nonRel(t, 4))
		else
			press = 0.d0
		end if
	end function nonRelativistic_pressure

	function nonRelativistic_entropy(cls, x, z) result(entropy)
		real(dl) :: entropy
		class(nonRelativistic_fermion), intent(in) :: cls
		real(dl), intent(in) :: x, z

		if (x .lt. cls%x_enDens_cut) then
			entropy = (cls%energyDensity(x, z) + cls%pressure(x, z))/z
		else
			entropy = 0.d0
		end if
	end function nonRelativistic_entropy

	subroutine nonRelativistic_initialize(cls, fermionName, isElectron, mass_factor, xcut)
		class(nonRelativistic_fermion) :: cls
		character(len=*), intent(in) :: fermionName
		logical, intent(in) :: isElectron
		real(dl), intent(in) :: mass_factor, xcut
		real(dl), dimension(:,:), allocatable :: ed_vec
		integer :: ix, iz, iflag
		real(dl) :: x,z, t1,t2
		real(8) :: timer1

		call addToLog("[cosmo] Initializing "//fermionName//"...")
		cls%fermionName = fermionName
		cls%isElectron = isElectron
		cls%mass_factor = mass_factor
		cls%x_enDens_cut = xcut

#ifndef NOINTERPOLATION
		allocate(ed_vec(interp_nx, interp_nz))
		!$omp parallel do default(shared) private(ix, iz) schedule(dynamic)
		do ix=1, interp_nx
			do iz=1, interp_nz
				ed_vec(ix,iz) = cls%energyDensityFull(interp_xvec(ix), interp_zvec(iz))
			end do
		end do
		!$omp end parallel do
		call cls%enDensInterp%initialize(interp_xvec, interp_zvec, ed_vec, iflag)!linear

		call random_seed()
		call random_number(x)
		call random_number(z)
		x=(x_fin-x_in)*x + x_in
		z=0.4d0*z + z_in
		write(*,"(' [cosmo] test energyDensity interpolation in ',*(E12.5))") x, z
		t1 = cls%energyDensityFull(x, z)
		t2 = cls%energyDensity(x, z)
		write(*,"(' [cosmo] comparison electron density (true vs interp): ',*(E17.10))") t1, t2
#endif
		call addToLog("[cosmo] ...done!")
	end subroutine nonRelativistic_initialize

	!functions derived from rho and drho/dx, for photon temperature evolution
	pure function J_int(o, u, a)
		real(dl) :: J_int, esuo
		real(dl), intent(in) :: o, u
		integer, intent(in) :: a

		esuo=exp(sqrt(u*u+o*o))
		J_int = u**a * esuo / ((1.d0+esuo)**2)
	end function J_int

	pure function J_funcFull(o, a)
		real(dl) :: J_funcFull
		real(dl), intent(in) :: o
		integer, intent(in) :: a
		integer :: i

		if (a.lt.2) then
			J_funcFull = rombint_ri(o, a, J_int, zero, upper, toler_jkyg, maxiter) / PISQ
		else
			J_funcFull = 0.d0
			do i=1, N_opt_xoz
				J_funcFull = J_funcFull + opt_xoz_w(i)*J_int(o, opt_xoz(i), a-2)
			end do
			J_funcFull = J_funcFull / PISQ
		end if
	end function J_funcFull

	pure function Jprime_int(o, u, a)
		real(dl) :: Jprime_int
		real(dl), intent(in) :: o, u
		integer, intent(in) :: a
		real(dl) :: uuoo, sqrtuuoo, expsqrtuuoo

		uuoo=u*u+o*o
		sqrtuuoo = sqrt(uuoo)
		expsqrtuuoo = exp(sqrtuuoo)

		Jprime_int = u**a * expsqrtuuoo * (1.d0-expsqrtuuoo) /(sqrtuuoo*(expsqrtuuoo+1)**3)
	end function Jprime_int

	pure function JprimeFull(o, a)
		real(dl) :: JprimeFull
		real(dl), intent(in) :: o
		integer, intent(in) :: a
		integer :: i

		if (a.lt.2) then
			JprimeFull = rombint_ri(o, a, Jprime_int, zero, upper, toler_jkyg, maxiter) * o / PISQ
		else
			JprimeFull = 0.d0
			do i=1, N_opt_xoz
				JprimeFull = JprimeFull + opt_xoz_w(i)*Jprime_int(o, opt_xoz(i), a-2)
			end do
			JprimeFull = JprimeFull * o / PISQ
		end if
	end function JprimeFull

	pure function K_int(o, u, a)
		real(dl) :: K_int,suo
		real(dl), intent(in) :: o, u
		integer, intent(in) :: a

		suo=sqrt(u*u+o*o)
		K_int = u**a / (suo * (1.d0+exp(suo)))
	end function K_int

	pure function K_funcFull(o, a)
		real(dl) :: K_funcFull
		real(dl), intent(in) :: o
		integer, intent(in) :: a
		integer :: i

		if (a.lt.2 .or. o.lt.5d-2) then
			K_funcFull = rombint_ri(o, a, K_int, zero, upper, toler_jkyg, maxiter) / PISQ
		else
			K_funcFull = 0.d0
			do i=1, N_opt_xoz
				K_funcFull = K_funcFull + opt_xoz_w(i)*K_int(o, opt_xoz(i), a-2)
			end do
			K_funcFull = K_funcFull / PISQ
		end if
	end function K_funcFull

	pure function Kprime_int(o, u, a)
		real(dl) :: Kprime_int
		real(dl), intent(in) :: o, u
		integer, intent(in) :: a
		real(dl) :: uuoo, sqrtuuoo, expsqrtuuoo

		uuoo=u*u+o*o
		sqrtuuoo = sqrt(uuoo)
		expsqrtuuoo = exp(sqrtuuoo)

		Kprime_int = u**a / (uuoo*sqrtuuoo*(expsqrtuuoo+1)**2) * &
			(1.d0 + expsqrtuuoo*(sqrtuuoo+1.d0))
	end function Kprime_int

	pure function KprimeFull(o, a)
		real(dl) :: KprimeFull
		real(dl), intent(in) :: o
		integer, intent(in) :: a

		KprimeFull = - rombint_ri(o, a, Kprime_int, zero, upper, toler_jkyg, maxiter) * o / PISQ
	end function KprimeFull

	pure function G12_funcFull(o)
		real(dl), dimension(2) :: G12_funcFull
		real(dl), intent(in) :: o
		real(dl) :: k2, j2, k2p, j2p, ga
		real(dl) :: k0, j0, k0p, j0p, j4p, gb, gc, osq

		if (ftqed_temperature_corr) then
			k2=k_funcFull(o, 2)
			j2=j_funcFull(o, 2)
			k2p=KprimeFull(o, 2)
			j2p=JprimeFull(o, 2)
			ga = (k2p/6.d0 - k2*k2p + j2p/6.d0 + j2p*k2 + j2*k2p)
			G12_funcFull(1) = PIx2*alpha_fine *(&
				(k2/3.d0 + 2.d0*k2*k2 - j2/6.d0 - k2*j2)/o + &
				ga )
			G12_funcFull(2) = PIx2*alpha_fine*( &
				o * ga &
				- 4.d0*((k2+j2)/6.d0 + k2*j2 - k2*k2/2.d0))

			if (ftqed_log_term) then
				!to edit
				G12_funcFull(1) = G12_funcFull(1) + &
					0.d0
				G12_funcFull(2) = G12_funcFull(2) + &
					0.d0
			end if

			if (ftqed_ord3) then
				osq = o*o
				k0=k_funcFull(o, 0)
				j0=j_funcFull(o, 0)
				k0p=KprimeFull(o, 0)
				j0p=JprimeFull(o, 0)
				j4p=JprimeFull(o, 4)
				gb = sqrt(k2 + osq*k0/2.d0)
				gc = 0.5*(2*j2 + osq*j0)/(2*k2 + osq*k0)
				G12_funcFull(1) = G12_funcFull(1) + &
					alpha_fine*electron_charge * gb * ( &
						(2.d0*j2-4.d0*k2)/o &
						- 2*j2p &
						- osq* j0p &
						- o*(2*k0+j0) &
						- gc * (o*(k0-j0) + k2p) &
					)
				G12_funcFull(2) = G12_funcFull(2) + &
					alpha_fine*electron_charge * gb * ( &
						gc*(2*j2 + osq*j0) &
						- 2*j4p/o &
						- o*(3*j2p+osq*j0p) &
					)
			end if
		else
			G12_funcFull = 0.d0
		end if
	end function G12_funcFull

	pure function nonRelativistic_dzodx_terms(cls, xoz) result(oj)
		real(dl), dimension(2) :: oj
		class(nonRelativistic_fermion), intent(in) :: cls
		real(dl), intent(in) :: xoz
		real(dl) :: j, y, o
		o = xoz * cls%mass_factor
		if (xoz .lt. cls%x_enDens_cut) then
			j = J_funcFull(o, 2)
			y = J_funcFull(o, 4)
		else
			j = 0.d0
			y = 0.d0
		end if
		oj(1) = cls%mass_factor * o * j
		oj(2) = o**2 * j + y
	end function

	!corrections to total energy and pressure
	pure function deltaRhoTot_em(x, z)
		real(dl) :: deltaRhoTot_em, z4, o, k2, j2, k0, j0, osqd
		real(dl), intent(in) :: x, z
		deltaRhoTot_em = 0.d0
		
		if (ftqed_temperature_corr) then
			z4 = z**4
			o = x/z
			k2 = k_funcFull(o, 2)
			j2 = j_funcFull(o, 2)
			deltaRhoTot_em = electron_charge_sq * z4 * ( &
				k2**2/2.d0 - (k2+j2)/6.d0 - k2*j2
			)

			if (ftqed_log_term) then
				deltaRhoTot_em = deltaRhoTot_em + &
					0.d0 !to edit
			end if

			if (ftqed_ord3) then
				osq = o*o/2.d0
				k0 = k_funcFull(o, 0)
				j0 = j_funcFull(o, 0)
				deltaRhoTot_em = deltaRhoTot_em &
					+ electron_charge_cub * z4 * &
					* sqrt(k2 + osqd*k0) &
					* (j2 + osqd*j0)
			end if
		end if
	end function deltaRhoTot_em

	pure function deltaPTot_em(x, z)
		real(dl) :: deltaPTot_em
		real(dl), intent(in) :: x, z
		deltaPTot_em = 0.d0
		
		if (ftqed_temperature_corr) then
			deltaPTot_em = 0.d0!to edit
			if (ftqed_log_term) then
				deltaPTot_em = 0.d0!to edit
			end if
		end if
	end function deltaPTot_em

	pure function deltaEntropyTot_em(x, z)
		real(dl) :: deltaEntropyTot_em
		real(dl), intent(in) :: x, z
		deltaEntropyTot_em = (deltaRhoTot_em(x,z) + deltaPTot_em(x,z)) / z
	end function deltaEntropyTot_em

	!functions for neutrino energy density
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
end module
