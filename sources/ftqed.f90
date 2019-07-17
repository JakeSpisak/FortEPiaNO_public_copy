module ftqed
	use precision
	use constants
	use utilities
	use fpErrors
	use linear_interpolation_module
	implicit none

	real(dl), parameter :: upper = 1.d2

#ifndef NOINTERPOLATION
	type(linear_interp_2d) :: dmeNLCorr
	type(linear_interp_3d) :: dmeCorr
#endif

	contains

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

	pure function G1_ln_integrand(x, z, y, k)
		real(dl) :: G1_ln_integrand
		real(dl), intent(in) :: x, y, z, k
		real(dl) :: Eksq, Eysq, Ek, Ey, Ekoz, Eyoz, Eekoz, Eeyoz, Nfk, Nfy
		real(dl) :: pxNk, pzNk, pxpzNk, pxNy

		G1_ln_integrand = 0.d0
		if (abs(y-k) .gt. 1d-7) then
			Eksq = x*x + k*k
			Ek = sqrt(Eksq)
			Ekoz = Ek/z
			Eekoz = exp(Ekoz)
			Nfk = 2.*fermiDirac(Ekoz)

			Eysq = x*x + y*y
			Ey = sqrt(Eysq)
			Eyoz = Ey/z
			Eeyoz = exp(Eyoz)
			Nfy = 2.*fermiDirac(Eyoz)

			pxNk = - x*Eekoz * Nfk**2 / (Ek*z)
			pzNk = Eekoz * Nfk**2 * Ek / (z**2)
			pxpzNk = x * Eekoz * Nfk**2 / (z**3) * (1 - 2*Eekoz*Nfk + z/Ek)
			pxNy = - x*Eeyoz * Nfy**2 / (Ey*z)

			G1_ln_integrand = &
				y * k / (Ey * Ek) * log( abs((y+k) / (y-k)) ) * &
				( &
					! -x(z( ... ) - ... )
					-x * ( &
						z * ( &
							pxNy*pzNk + Nfy*pxpzNk &
						) &
						- Nfy*pxNk &
					) &
					!N_y N_k
					- Nfy*Nfk &
					!z N_y \partial_z N_k
					- z * Nfy * pzNk &
					! - x^2 ...
					- x**2 * (Eksq+Eysq)/(2*Eksq*Eysq) * ( &
						2*z*Nfy*pzNk - Nfy*Nfk &
					) &
				)
		end if
	end function G1_ln_integrand

	pure function G2_ln_integrand(x, z, y, k)
		real(dl) :: G2_ln_integrand
		real(dl), intent(in) :: x, y, z, k
		real(dl) :: Eksq, Eysq, Ek, Ey, Ekoz, Eyoz, Nfk, Nfy, NfEek

		G2_ln_integrand = 0.d0
		if (abs(y-k) .gt. 1d-7) then
			Eksq = x*x + k*k
			Ek = sqrt(Eksq)
			Ekoz = Ek/z
			Nfk = 2.*fermiDirac(Ekoz)
			NfEek = Nfk * exp(Ekoz)

			Eysq = x*x + y*y
			Ey = sqrt(Eysq)
			Eyoz = Ey/z
			Nfy = 2.*fermiDirac(Eyoz)

			G2_ln_integrand = &
				y * k / Ey * log( abs((y+k) / (y-k)) ) * &
				Nfk * Nfy * NfEek / z**4 * &
				(2.d0*Ek*NfEek - Ek + Nfy*exp(Eyoz)*Ey - 2.d0*z)
		end if
	end function G2_ln_integrand

	pure function G12_funcFull(x, z)
		real(dl), dimension(2) :: G12_funcFull
		real(dl), intent(in) :: x, z
		real(dl) :: o
		real(dl) :: k2, j2, k2p, j2p, ga
		real(dl) :: k0, j0, k0p, j0p, j4p, gb, gc, osq
		real(dl) :: coeff, integ
		integer :: ix, iy
		real(dl), dimension(:,:), allocatable :: fval

		if (ftqed_temperature_corr) then
			o = x/z
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
				coeff = electron_charge**2 * x / (PIx2**4*z**3)
				allocate(fval(ln_2dint_Npts, ln_2dint_Npts))
				!G1
				do ix=1,ln_2dint_Npts
					do iy=1,ln_2dint_Npts
						fval(ix, iy) = &
							G1_ln_integrand(x, z, ln_2dint_y(ix), ln_2dint_y(iy))
					end do
				end do
				integ = integral_NC_2d(&
						ln_2dint_Npts, &
						ln_2dint_Npts, &
						ln_2dint_dy, &
						ln_2dint_dy, &
						fval &
						)
				G12_funcFull(1) = G12_funcFull(1) + coeff * integ
				!G2
				do ix=1,ln_2dint_Npts
					do iy=1,ln_2dint_Npts
						fval(ix, iy) = &
							G2_ln_integrand(x, z, ln_2dint_y(ix), ln_2dint_y(iy))
					end do
				end do
				integ = integral_NC_2d(&
						ln_2dint_Npts, &
						ln_2dint_Npts, &
						ln_2dint_dy, &
						ln_2dint_dy, &
						fval &
						)
				G12_funcFull(2) = G12_funcFull(2) + x*z*coeff * integ
				deallocate(fval)
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

	!corrections to total energy and pressure
	pure function deltaRho_ln_integrand(x, z, y, k)
		real(dl) :: deltaRho_ln_integrand
		real(dl), intent(in) :: x, y, z, k
		real(dl) :: Eksq, Eysq, Ek, Ey, Ekoz, Eyoz, Nfk, Nfy

		deltaRho_ln_integrand = 0.d0
		if (abs(y-k) .gt. 1d-7) then
			Eksq = x*x + k*k
			Ek = sqrt(Eksq)
			Ekoz = Ek/z
			Nfk = 2.*fermiDirac(Ekoz)

			Eysq = x*x + y*y
			Ey = sqrt(Eysq)
			Eyoz = Ey/z
			Nfy = 2.*fermiDirac(Eyoz)

			deltaRho_ln_integrand = &
				y * k / (Ek * Ey) * log( abs((y+k) / (y-k)) ) * &
				Nfy * Nfk * ( &
					2. * Nfk * exp(Ekoz) * Ek / z - 1.d0 &
				)
		end if
	end function deltaRho_ln_integrand

	pure function deltaP_ln_integrand(x, z, y, k)
		real(dl) :: deltaP_ln_integrand
		real(dl), intent(in) :: x, y, z, k
		real(dl) :: Eksq, Eysq, Ek, Ey, Ekoz, Eyoz, Nfk, Nfy

		deltaP_ln_integrand = 0.d0
		if (abs(y-k) .gt. 1d-7) then
			Eksq = x*x + k*k
			Ek = sqrt(Eksq)
			Ekoz = Ek/z
			Nfk = 2.*fermiDirac(Ekoz)

			Eysq = x*x + y*y
			Ey = sqrt(Eysq)
			Eyoz = Ey/z
			Nfy = 2.*fermiDirac(Eyoz)

			deltaP_ln_integrand = &
				y * k / (Ek * Ey) * log( abs((y+k) / (y-k)) ) * &
				Nfk * Nfy
		end if
	end function deltaP_ln_integrand

	pure function deltaRhoTot_em(x, z)
		real(dl) :: deltaRhoTot_em, z4, o, k2, j2, k0, j0, osqd
		real(dl), intent(in) :: x, z
		integer :: ix, iy
		real(dl), dimension(:,:), allocatable :: fval

		deltaRhoTot_em = 0.d0
		
		if (ftqed_temperature_corr) then
			z4 = z**4
			o = x/z
			k2 = k_funcFull(o, 2)
			j2 = j_funcFull(o, 2)
			deltaRhoTot_em = electron_charge_sq * z4 * ( &
				k2**2/2.d0 - (k2+j2)/6.d0 - k2*j2 &
			)

			if (ftqed_log_term) then
				allocate(fval(ln_2dint_Npts, ln_2dint_Npts))
				do ix=1,ln_2dint_Npts
					do iy=1,ln_2dint_Npts
						fval(ix, iy) = &
							deltaRho_ln_integrand(x, z, ln_2dint_y(ix), ln_2dint_y(iy))
					end do
				end do
				deltaRhoTot_em = deltaRhoTot_em + &
					electron_charge_sq * x**2 /(16.d0*PISQ*PISQ) &
					* integral_NC_2d(&
						ln_2dint_Npts, &
						ln_2dint_Npts, &
						ln_2dint_dy, &
						ln_2dint_dy, &
						fval &
						)
				deallocate(fval)
			end if

			if (ftqed_ord3) then
				osqd = o*o/2.d0
				k0 = k_funcFull(o, 0)
				j0 = j_funcFull(o, 0)
				deltaRhoTot_em = deltaRhoTot_em &
					+ electron_charge_cub * z4 / PI &
					* sqrt(k2 + osqd*k0) &
					* (j2 + osqd*j0)
			end if
		end if
	end function deltaRhoTot_em

	pure function deltaPTot_em(x, z)
		real(dl) :: deltaPTot_em
		real(dl) :: z4, o, k2, k2rk0
		real(dl), intent(in) :: x, z
		integer :: ix, iy
		real(dl), dimension(:,:), allocatable :: fval

		deltaPTot_em = 0.d0
		
		if (ftqed_temperature_corr) then
			z4 = z**4
			o = x/z
			k2 = k_funcFull(o, 2)
			deltaPTot_em = - electron_charge_sq * z4 * &
				k2 * (one_sixth + 0.5d0 * k2)

			if (ftqed_log_term) then
				allocate(fval(ln_2dint_Npts, ln_2dint_Npts))
				do ix=1,ln_2dint_Npts
					do iy=1,ln_2dint_Npts
						fval(ix, iy) = &
							deltaP_ln_integrand(x, z, ln_2dint_y(ix), ln_2dint_y(iy))
					end do
				end do
				deltaPTot_em = deltaPTot_em + &
					electron_charge_sq * x**2 /(16.d0*PISQ*PISQ) &
					* integral_NC_2d(&
						ln_2dint_Npts, &
						ln_2dint_Npts, &
						ln_2dint_dy, &
						ln_2dint_dy, &
						fval &
						)
				deallocate(fval)
			end if

			if (ftqed_ord3) then
				k2rk0 = k2 + o*o/2.d0 * k_funcFull(o, 0)
				deltaPTot_em = deltaPTot_em + &
					2.d0 * electron_charge_cub * z4 / (3.d0*PI) * &
					sqrt(k2rk0) * k2rk0
			end if
		end if
	end function deltaPTot_em

	pure function deltaEntropyTot_em(x, z)
		real(dl) :: deltaEntropyTot_em
		real(dl), intent(in) :: x, z
		deltaEntropyTot_em = (deltaRhoTot_em(x,z) + deltaPTot_em(x,z)) / z
	end function deltaEntropyTot_em

	!electron mass corrections for collision integrals
#ifndef NOINTERPOLATION
	subroutine init_interp_dme2_e
		real(dl), dimension(:,:), allocatable :: dmg_vec
		real(dl), dimension(:,:,:), allocatable :: dme_vec
		integer :: ix, iy, iz, iflag
		real(dl) :: x, z, t1, t2
		real(8) :: timer1
		logical :: initial

		call addToLog("[interactions] Initializing interpolation for electron mass corrections...")
		allocate(dme_vec(interp_nx, interp_ny, interp_nz))
		!$omp parallel do default(shared) private(ix, iy, iz) schedule(dynamic)
		do ix=1, interp_nx
			do iy=1, interp_ny
				do iz=1, interp_nz
					dme_vec(ix,iy,iz) = dme2_electronFull(interp_xvec(ix),interp_yvec(iy),interp_zvec(iz))
				end do
			end do
		end do
		!$omp end parallel do
		call dmeCorr%initialize(interp_xvec, interp_yvec, interp_zvec, dme_vec, iflag)!linear

		!store dme2 without log term for use in collision terms
		initial = ftqed_log_term
		ftqed_log_term = .false.
		allocate(dmg_vec(interp_nx, interp_nz))
		!$omp parallel do default(shared) private(ix, iz) schedule(dynamic)
		do ix=1, interp_nx
			do iz=1, interp_nz
				dmg_vec(ix,iz) = dme2_electronFull(interp_xvec(ix), 0.d0, interp_zvec(iz))
			end do
		end do
		!$omp end parallel do
		call dmeNLCorr%initialize(interp_xvec, interp_zvec, dmg_vec, iflag)!linear

		call random_seed()
		if (timing_tests) then
			call tic(timer1)
			write (*,*) "[interactions] now doing some timing..."
			call tic(timer1)
			do ix=1, 1000000
				call random_number(x)
				call random_number(z)
				x=(x_fin-x_in)*x + x_in
				z=0.4d0*z + z_in
				t1 = dme2_nolog(x,z)
			end do
			call toc(timer1, "<interpolated>")

			call tic(timer1)
			do ix=1, 1000000
				call random_number(x)
				call random_number(z)
				x=(x_fin-x_in)*x + x_in
				z=0.4d0*z + z_in
				t1 = dme2_electronFull(x,0.d0,z)
			end do
			call toc(timer1, "<full>")
		end if
		ftqed_log_term = initial
		call random_number(x)
		call random_number(z)
		x=(x_fin-x_in)*x + x_in
		z=0.4d0*z + z_in
		write(*,"(' [interactions] test dme2_electronInterp in ',*(E12.5))") x, 0.01d0, z
		t1 = dme2_electronFull(x, 0.01d0, z)
		t2 = dme2_electron(x, 0.01d0, z)
		write(*,"(' [interactions] comparison (true vs interp): ',*(E17.10))") t1,t2

		deallocate(dmg_vec)
		deallocate(dme_vec)
		call addToLog("[interactions] ...done!")
	end subroutine init_interp_dme2_e
#endif

	pure function dme2_e_i1(x, z, y)
	!doi:10.1016/S0370-2693(02)01622-2 eq.12 first integral
	!equal to integral in eq.13
		real(dl) :: dme2_e_i1
		real(dl), intent(in) :: x, z, y
		real(dl) :: Ekm
		Ekm = E_k_m(y, x)
		dme2_e_i1 = 1.d0/Ekm * fermiDirac(Ekm/z)
	end function dme2_e_i1

	pure function dme2_e_i2(x, y, z, k)
		!doi:10.1016/S0370-2693(02)01622-2 eq.12 second integral
		real(dl) :: dme2_e_i2
		real(dl), intent(in) :: x, y, z, k
		real(dl) :: m, T, p, Ekm

		if (abs(y - k).gt.1d-6) then
			Ekm = E_k_m(k, x)
			dme2_e_i2 = k/Ekm * fermiDirac(Ekm/z)*log(abs((y+k)/(y-k)))
		else
			dme2_e_i2 = 0.d0
		end if
	end function dme2_e_i2

	pure function dme2_electronFull(x, y, z, logt)
	!doi:10.1016/S0370-2693(02)01622-2 eq.12
		real(dl) :: dme2_electronFull, integr_1, integr_2
		real(dl), intent(in) :: x, y, z
		logical, intent(in), optional :: logt
		integer :: i
		logical :: uselog

		if (present(logt)) then
			uselog = ftqed_log_term .and. logt
		else
			uselog = ftqed_log_term
		end if
		if (ftqed_temperature_corr) then
			integr_1 = 0.d0
			do i=1, N_opt_y
				integr_1 = integr_1 + opt_y_w(i)*dme2_e_i1(x, z, opt_y(i))
			end do
			dme2_electronFull = 2. * alpha_fine * (z*z * PID3 + integr_1/PID2)
			if (uselog) then
				integr_2 = 0.d0
				do i=1, N_opt_y
					integr_2 = integr_2 + opt_y_w(i)*dme2_e_i2(x, y, z, opt_y(i))/(opt_y(i)**2)
				end do
				dme2_electronFull = dme2_electronFull &
					- x*x* alpha_fine / y / PID2 * integr_2
			end if
		else
			dme2_electronFull = 0.d0
		end if
	end function dme2_electronFull

#ifndef NOINTERPOLATION
	function dme2_electron(x, y, z)
		real(dl) :: dme2_electron
		real(dl), intent(in) :: x, y, z
		integer :: iflag

		call dmeCorr%evaluate(x, y, z, dme2_electron)
	end function dme2_electron
#endif

#ifndef NOINTERPOLATION
	function dme2_nolog(x, z)
		real(dl) :: dme2_nolog
		real(dl), intent(in) :: x, z
		integer :: iflag

		call dmeNLCorr%evaluate(x, z, dme2_nolog)
	end function dme2_nolog
#endif

end module ftqed
