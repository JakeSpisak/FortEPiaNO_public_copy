module fpCosmology
	use precision
	use constants
	use utilities
	use fpErrors
	use ftqed
	use linear_interpolation_module
	use sgTestUtils
	implicit none

	type nonRelativistic_fermion
		character(len=20) :: fermionName
		logical :: isElectron !if true, compute the electromagnetic corrections to the mass
		real(dl) :: mass_factor !mass of the particle divided by electron mass
		real(dl) :: x_enDens_cut !for x larger than this, do not consider this particle (fix energy density and other things to 0)
		type(linear_interp_2d) :: enDensNoThMassInterp !information for interpolating the energy density without thermal mass corrections (for total energy density)
		type(linear_interp_2d) :: enDensWThMassInterp !information for interpolating the energy density with thermal mass corrections (for lepton matter potentials)
		type(linear_interp_2d) :: pressInterp !information for interpolating the energy density
		contains
		procedure :: initialize => nonRelativistic_initialize !save properties and prepare energy density interpolation
		procedure :: dzodx_terms => nonRelativistic_dzodx_terms !contributions to the dz/dx numerator and denominator
		procedure :: energyDensityFull => nonRelativistic_energyDensity_full !full integral of the energy density
		procedure :: energyDensity => nonRelativistic_energyDensity !interpolated energy density
		procedure :: entropy => nonRelativistic_entropy !entropy
		procedure :: numberDensity => nonRelativistic_numberDensity_full !number density
		procedure :: pressureFull => nonRelativistic_pressure_full !pressure
		procedure :: pressure => nonRelativistic_pressure !interpolated pressure
	end type nonRelativistic_fermion

#ifdef DO_MUONS
	integer, parameter :: fermions_number = 2
#else
	integer, parameter :: fermions_number = 1
#endif
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
			rho_r = rho_r + fermions(ix)%energyDensity(x, z, .false.)
		end do
		rho_r = rho_r + deltaRhoTot_em(x, z)
	end function

	!photons
	elemental function photonNumberDensity(z)
		real(dl) :: photonNumberDensity
		real(dl), intent(in) :: z

		photonNumberDensity = zetaOfThreex2DPISQ * z**3
	end function photonNumberDensity

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
	pure function integrand_rho_nonRel(x, z, y, elTherMass)
		real(dl) :: integrand_rho_nonRel, Emk
		real(dl), intent(in) :: x, z, y
		logical, intent(in) :: elTherMass
		if (elTherMass) then
			Emk = Ebare_i_dme(x, y, dme2_electronFull(x, y, z))
		else
			Emk = E_k_m(y, x)
		end if
		integrand_rho_nonRel = Emk * fermiDirac(Emk/z)
	end function integrand_rho_nonRel

	pure function integrand_uX_Ek_nonRel(u, x, z, mf, elTherMass, n)
		real(dl) :: integrand_uX_Ek_nonRel
		real(dl), intent(in) :: u, x, z, mf
		logical, intent(in) :: elTherMass
		integer, intent(in) :: n
		real(dl) :: w, squ2w2, expsqu2w2
		w = x*mf/z
		if (elTherMass) then
			squ2w2 = Ebare_i_dme(u, w, dme2_electronFull(x, u, z))
		else
			squ2w2 = E_k_m(u, w)
		end if
		expsqu2w2 = exp(squ2w2)
		integrand_uX_Ek_nonRel = u**n/(squ2w2*(1.d0+expsqu2w2))
	end function integrand_uX_Ek_nonRel

	pure function nonRelativistic_numberDensity_full(cls, x, z, elTherMass) result (nredf)!fermion + antifermion
		class(nonRelativistic_fermion), intent(in) :: cls
		real(dl) :: nredf
		real(dl), intent(in) :: x,z
		logical, intent(in) :: elTherMass
		real(dl) :: nx, ny, Emk
		integer :: i

		nredf = 0.d0
		if (x .lt. cls%x_enDens_cut) then
			nx = x * cls%mass_factor
			do i=1, N_opt_y
			ny = opt_y(i)
			if (elTherMass) then
				Emk = Ebare_i_dme(nx, ny, dme2_electronFull(nx, ny, z))
			else
				Emk = E_k_m(ny, nx)
			end if
				nredf = nredf + opt_y_w(i) * fermiDirac(Emk/z)
			end do
			nredf = nredf / PISQD2
			!the factor is given by g = 2(elicity) * 2(f+\bar f)
		end if
		if (abs(nredf) .lt. 1d-99) &
			nredf = 0.0d0
	end function nonRelativistic_numberDensity_full

	pure function nonRelativistic_energyDensity_full(cls, x, z, elTherMass) result (nredf)!fermion + antifermion
		class(nonRelativistic_fermion), intent(in) :: cls
		real(dl) :: nredf
		real(dl), intent(in) :: x,z
		logical, intent(in) :: elTherMass
		integer :: i

		nredf = 0.d0
		if (x .lt. cls%x_enDens_cut) then
			do i=1, N_opt_y
				nredf = nredf + opt_y_w(i)* integrand_rho_nonRel(x*cls%mass_factor, z, opt_y(i), elTherMass)
			end do
			nredf = nredf / PISQD2
			!the factor is given by g = 2(elicity) * 2(f+\bar f)
		end if
	end function nonRelativistic_energyDensity_full

	function nonRelativistic_energyDensity(cls, x, z, elTherMass) result (ed)
		class(nonRelativistic_fermion) :: cls
		real(dl) :: ed
		real(dl), intent(in) :: x, z
		logical, intent(in) :: elTherMass

#ifdef NO_INTERPOLATION
		ed = nonRelativistic_energyDensity_full(cls, x, z, elTherMass)
#else
		if (cls%isElectron .and. elTherMass) then
			call cls%enDensWThMassInterp%evaluate(x, z, ed)
		else
			call cls%enDensNoThMassInterp%evaluate(x, z, ed)
		end if
#endif
		if (abs(ed) .lt. 1d-99) &
			ed = 0.0d0
	end function nonRelativistic_energyDensity

	pure function integrate_uX_Ek_nonRel(x, z, mf, elTherMass, n)
		!computes the integral of u**n/(sqrt(u**2+w**2)(1+exp(sqrt(u**2+w**2))))
		!useful to get the pressure for non-relativistic species
		real(dl) :: integrate_uX_Ek_nonRel
		real(dl), intent(in) :: x, z, mf
		integer, intent(in) :: n
		logical, intent(in) :: elTherMass
		integer :: i
		integrate_uX_Ek_nonRel = 0.d0
		do i = 1, N_opt_xoz
			integrate_uX_Ek_nonRel = integrate_uX_Ek_nonRel &
				+ opt_y_w(i)*integrand_uX_Ek_nonRel(opt_y(i), x, z, mf, elTherMass, n-2)!n-2 because the weights already take into account y^2
		end do
		integrate_uX_Ek_nonRel = integrate_uX_Ek_nonRel/PISQx2
	end function integrate_uX_Ek_nonRel

	pure function nonRelativistic_pressure_full(cls, x, z, elTherMass) result(press)
		real(dl) :: press
		class(nonRelativistic_fermion), intent(in) :: cls
		real(dl), intent(in) :: x, z
		logical, intent(in) :: elTherMass

		if (x .lt. cls%x_enDens_cut) then
			press = 4.d0*z**4/3.d0*(integrate_uX_Ek_nonRel(x, z, cls%mass_factor, elTherMass, 4))
		else
			press = 0.d0
		end if
	end function nonRelativistic_pressure_full

	function nonRelativistic_pressure(cls, x, z, elTherMass) result (ed)
		class(nonRelativistic_fermion) :: cls
		real(dl) :: ed
		real(dl), intent(in) :: x, z
		logical, intent(in) :: elTherMass

#ifdef NO_INTERPOLATION
		ed = nonRelativistic_pressure_full(cls, x, z, elTherMass)
#else
		call cls%pressInterp%evaluate(x, z, ed)
#endif
	end function nonRelativistic_pressure

	function nonRelativistic_entropy(cls, x, z) result(entropy)
		real(dl) :: entropy
		class(nonRelativistic_fermion), intent(in) :: cls
		real(dl), intent(in) :: x, z

		if (x .lt. cls%x_enDens_cut) then
			entropy = (cls%energyDensity(x, z, .false.) + cls%pressure(x, z, .false.))/z
		else
			entropy = 0.d0
		end if
		if (abs(entropy) .lt. 1d-99) &
			entropy = 0.0d0
	end function nonRelativistic_entropy

	subroutine nonRelativistic_initialize(cls, fermionName, isElectron, mass_factor, xcut)
		class(nonRelativistic_fermion) :: cls
		character(len=*), intent(in) :: fermionName
		logical, intent(in) :: isElectron
		real(dl), intent(in) :: mass_factor, xcut
		real(dl), dimension(:,:), allocatable :: ednt_vec, edwt_vec, p_vec
		logical :: thmass
		integer :: ix, iz, iflag
		real(dl) :: x,z, t1,t2
		real(8) :: timer1
		character(len=300) :: tmpstr
		integer, parameter :: uid = 8324
		logical :: exists

		call addToLog("[cosmo] Initializing "//fermionName//"...")
		cls%fermionName = fermionName
		cls%isElectron = isElectron
		cls%mass_factor = mass_factor
		cls%x_enDens_cut = xcut

#ifndef NO_INTERPOLATION
		allocate(ednt_vec(interp_nx, interp_nz))
		allocate(edwt_vec(interp_nx, interp_nz))
		allocate(p_vec(interp_nx, interp_nz))
		thmass = isElectron .and. ftqed_e_mth_leptondens
		write(tmpstr, "(A,'cosmo_',A,'_',L,'_',L,'_',L,'.dat')") trim(get_interpolation_folder()), fermionName, thmass, ftqed_temperature_corr, ftqed_log_term
		inquire(file=trim(tmpstr), exist=exists)
		if (exists) then
			call addToLog("[cosmo] read values from file: "//trim(tmpstr))
			open(file=trim(tmpstr), unit=uid, form="unformatted")
			do ix=1, interp_nx
				do iz=1, interp_nz
					if (cls%isElectron) then
						read(uid) edwt_vec(ix,iz), ednt_vec(ix,iz), p_vec(ix,iz)
					else
						read(uid) ednt_vec(ix,iz), p_vec(ix,iz)
					end if
				end do
			end do
			close(uid)
			call addToLog("[cosmo] check if few saved values are correct: ")
			ix=123
			iz=89
			if (cls%isElectron) &
				call assert_double_rel_safe( &
					"check saved "//fermionName//" quantities edwt A", &
					edwt_vec(ix,iz), &
					cls%energyDensityFull(interp_xvec(ix), interp_zvec(iz), .true.), &
					1d-7, 1d-6 &
				)
			call assert_double_rel_safe( &
				"check saved "//fermionName//" quantities ednt A", &
				ednt_vec(ix,iz), &
				cls%energyDensityFull(interp_xvec(ix), interp_zvec(iz), .false.), &
				1d-7, 1d-6 &
			)
			call assert_double_rel_safe( &
				"check saved "//fermionName//" quantities p A", &
				p_vec(ix,iz), &
				cls%pressureFull(interp_xvec(ix), interp_zvec(iz), thmass), &
				1d-7, 1d-6 &
			)
			ix=611
			iz=189
			if (cls%isElectron) &
				call assert_double_rel_safe( &
					"check saved "//fermionName//" quantities edwt A", &
					edwt_vec(ix,iz), &
					cls%energyDensityFull(interp_xvec(ix), interp_zvec(iz), .true.), &
					1d-7, 1d-6 &
				)
			call assert_double_rel_safe( &
				"check saved "//fermionName//" quantities ednt A", &
				ednt_vec(ix,iz), &
				cls%energyDensityFull(interp_xvec(ix), interp_zvec(iz), .false.), &
				1d-7, 1d-6 &
			)
			call assert_double_rel_safe( &
				"check saved "//fermionName//" quantities p A", &
				p_vec(ix,iz), &
				cls%pressureFull(interp_xvec(ix), interp_zvec(iz), thmass), &
				1d-7, 1d-6 &
			)
			call addToLog("everything works!")
		else
			!$omp parallel do default(shared) private(ix, iz) schedule(dynamic)
			do ix=1, interp_nx
				do iz=1, interp_nz
					if (cls%isElectron) &
						edwt_vec(ix,iz) = cls%energyDensityFull(interp_xvec(ix), interp_zvec(iz), .true.)
					ednt_vec(ix,iz) = cls%energyDensityFull(interp_xvec(ix), interp_zvec(iz), .false.)
					p_vec(ix,iz) = cls%pressureFull(interp_xvec(ix), interp_zvec(iz), thmass)
				end do
			end do
			!$omp end parallel do
			open(file=trim(tmpstr), unit=uid, status="unknown", form="unformatted")
			do ix=1, interp_nx
				do iz=1, interp_nz
					if (cls%isElectron) then
						write(uid) edwt_vec(ix,iz), ednt_vec(ix,iz), p_vec(ix,iz)
					else
						write(uid) ednt_vec(ix,iz), p_vec(ix,iz)
					end if
				end do
			end do
			close(uid)
			call addToLog("[cosmo] values saved to file: "//trim(tmpstr))
		end if
		if (cls%isElectron) &
			call cls%enDensWThMassInterp%initialize(interp_xvec, interp_zvec, edwt_vec, iflag)!linear
		call cls%enDensNoThMassInterp%initialize(interp_xvec, interp_zvec, ednt_vec, iflag)!linear
		call cls%pressInterp%initialize(interp_xvec, interp_zvec, p_vec, iflag)!linear

		if (tests_interpolations) then
			call random_seed()
			call random_number(x)
			call random_number(z)
			x=(x_fin-x_in)*x + x_in
			z=0.4d0*z + z_in
			write(*,"(' [cosmo] test energyDensity/pressure interpolation in x,z=',*(E12.5))") x, z
			t1 = cls%energyDensityFull(x, z, thmass)
			t2 = cls%energyDensity(x, z, thmass)
			write(*,"(' [cosmo] comparing energy density (true vs interp): ',*(E17.10))") t1, t2
			t1 = cls%pressureFull(x, z, thmass)
			t2 = cls%pressure(x, z, thmass)
			write(*,"(' [cosmo] comparing pressure (true vs interp): ',*(E17.10))") t1, t2
		endif
#endif
		call addToLog("[cosmo] ...done!")
	end subroutine nonRelativistic_initialize

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

	!functions for neutrino number and energy density
	function nuNumberDensityNC(i1, i2, reim)
		real(dl) :: nuNumberDensityNC, y
		integer, intent(in) :: i1, i2
		logical, intent(in), optional :: reim
		integer :: ix
		logical :: useim

		useim=.false.
		if (present(reim)) then
			if (.not.reim) &
				useim=.true.
		end if
		if (useim) then
			do ix=1, Ny
				y = y_arr(ix)
				fy_arr(ix) = y*y * nuDensMatVecFD(ix)%im(i1, i2)
			end do
		else
			do ix=1, Ny
				y = y_arr(ix)
				fy_arr(ix) = y*y * nuDensMatVecFD(ix)%re(i1, i2)
			end do
		end if
		nuNumberDensityNC = integral_NC_1d(Ny, dy_arr, fy_arr) / PISQ
	end function nuNumberDensityNC

	function nuNumberDensityGL(i1, i2, reim)
		real(dl) :: nuNumberDensityGL
		integer, intent(in) :: i1, i2
		logical, intent(in), optional :: reim
		integer :: ix
		logical :: useim

		useim=.false.
		if (present(reim)) then
			if (.not.reim) &
				useim=.true.
		end if
		if (useim) then
			do ix=1, Ny
				fy_arr(ix) = nuDensMatVecFD(ix)%im(i1, i2) / y_arr(ix)
			end do
		else
			do ix=1, Ny
				fy_arr(ix) = nuDensMatVecFD(ix)%re(i1, i2) / y_arr(ix)
			end do
		end if
		nuNumberDensityGL = integral_GL_1d(w_gl_arr, fy_arr) / PISQ
	end function nuNumberDensityGL

	function allNuNumberDensity()
		use fpInterfaces1
		real(dl) :: allNuNumberDensity
		integer :: ix
		procedure (nuDensity_integrator), pointer :: nuNumberDensityInt

		if (use_gauss_laguerre) then
			nuNumberDensityInt => nuNumberDensityGL
		else
			nuNumberDensityInt => nuNumberDensityNC
		end if

		allNuNumberDensity = 0.d0
		do ix=1, flavorNumber
			allNuNumberDensity = allNuNumberDensity + nuNumberDensityInt(ix, ix)*nuFactor(ix)
		end do
		allNuNumberDensity = allNuNumberDensity
	end function allNuNumberDensity

	function nuDensityNC(i1, i2, reim)
		real(dl) :: nuDensityNC, y
		integer, intent(in) :: i1, i2
		logical, intent(in), optional :: reim
		integer :: ix
		logical :: useim

		useim=.false.
		if (present(reim)) then
			if (.not.reim) &
				useim=.true.
		end if
		if (useim) then
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
		logical :: useim

		useim=.false.
		if (present(reim)) then
			if (.not.reim) &
				useim=.true.
		end if
		if (useim) then
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
		use fpInterfaces1
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

	function Neff_from_rho_z(z)
		real(dl) :: Neff_from_rho_z
		real(dl), intent(in) :: z

		Neff_from_rho_z = (zid)**4 * allNuDensity()/photonDensity(z) / 0.875d0
	end function Neff_from_rho_z
end module
