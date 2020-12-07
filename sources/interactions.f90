module fpInteractions
	use precision
	use constants
	use variables
	use utilities
	use ftqed
	use fpErrors
	use fpMatrices
	use linear_interpolation_module
	implicit none

	contains

	!fitting formula for the above function
	elemental function dy_damping_fit(y)
		real(dl) :: dy_damping_fit
		real(dl), intent(in) :: y
		real(dl), parameter :: d0 = 129.875
		real(dl), parameter :: dinf = 100.999
		real(dl), parameter :: a0 = 90.7332
		real(dl), parameter :: a1 = -48.4473
		real(dl), parameter :: a2 = 20.1219
		real(dl), parameter :: b1 = -0.529157
		real(dl), parameter :: b2 = 0.20649
		real(dl) :: e101y, e001y, ly
		e101y = exp(-1.01d0*y)
		e001y = exp(-0.01d0*y)
		ly = log(y)
		dy_damping_fit = d0*e101y &
		    + dinf*(1.d0-e001y) &
		    + (e001y - e101y) &
		    * ( &
			(a0 + a1*ly + a2*ly*ly) &
			/ (1.d0 + b1*ly + b2*ly*ly) &
		    )
	end function dy_damping_fit

	!damping factor coefficients are defined here
	subroutine setDampingFactorCoeffs
		integer :: ix, iy
		real(dl) :: nue_nue_nux, nue_nue_nus, nue_numu_nutau, nue_nux_nus
		real(dl) :: nunu_nue_nux, nunu_nue_nus, nunu_numu_nutau, nunu_nux_nus

		dampTermMatrixCoeffNue = 0.d0
		dampTermMatrixCoeffNunu = 0.d0

#ifdef NO_NUE_ANNIHILATION
		if (collint_offdiag_damping) &
			call criticalError("Nu-e annihilation channel can only be disabled in full integral calculations at the moment")
#else
		!numbers from McKellar:1992ja
		!terms for scattering, annihilation with electrons
		!nu_e - nu_X
		nue_nue_nux = &
			8.d0 &!e+nu -> e+nu
			+ (8.d0*sin2thW**2 + 1.d0)!nu+bnu -> e+e-
		!nu_mu - nu_tau
		nue_numu_nutau = &
			(8.d0*sin2thW**2 - 4.d0*sin2thW + 1.d0)!nu+bnu -> e+e-
		!nu_e - nu_s
		nue_nue_nus = 2.d0*(&
			(8.d0*sin2thW**2 + 4.d0*sin2thW + 1.d0) &!e+nu -> e+nu
			+ (4.d0*sin2thW**2 + 2.d0*sin2thW + 0.5d0) &!nu+bnu -> e+e-
		)
		!nu_X - nu_s
		nue_nux_nus = 2.d0*(&
			(8.d0*sin2thW**2 - 4.d0*sin2thW + 1.d0) &!e+nu -> e+nu
			+ (4.d0*sin2thW**2 - 2.d0*sin2thW + 0.5d0) &!nu+bnu -> e+e-
		)
#endif
		!terms for nunu scattering
		!nu_e - nu_X
		nunu_nue_nux = 6.d0 !nu+(b)nu -> nu+(b)nu
		!nu_mu - nu_tau
		nunu_numu_nutau = 6.d0 !nu+(b)nu -> nu+(b)nu
		!nu_e - nu_s
		nunu_nue_nus = 2.d0 * 13.d0 !nu+(b)nu -> nu+(b)nu
		!nu_X - nu_s
		nunu_nux_nus = 2.d0 * 13.d0 !nu+(b)nu -> nu+(b)nu

		if (collint_offdiag_damping .and. collint_damping_type.eq.1) then
			!formulas from Bennett:2020zkv
			nunu_nue_nux = 1.d0
			nunu_numu_nutau = 1.d0
			nue_nue_nux = 2.d0*sin2thW**2 + 0.25d0
			nue_numu_nutau = 2.d0*sin2thW**2 - sin2thW + 0.25d0
			nunu_nue_nus = 0.d0!check
			nunu_nux_nus = 0.d0!check
			nue_nue_nus = 3.d0*sin2thW**2 + 1.d0*sin2thW + 0.25d0!check
			nue_nux_nus = 3.d0*sin2thW**2 - 1.d0*sin2thW + 0.25d0!check
			if (any(sterile)) &
				call criticalError("Error: damping terms not yet implemented with sterile neutrinos.")
		end if
		if (flavorNumber .ge. 2) then
			if (sterile(2)) then
				dampTermMatrixCoeffNue(1, 2) = nue_nue_nus
				dampTermMatrixCoeffNunu(1, 2) = nunu_nue_nus
			else
				dampTermMatrixCoeffNue(1, 2) = nue_nue_nux
				dampTermMatrixCoeffNunu(1, 2) = nunu_nue_nux
			end if
		end if
		if (flavorNumber .ge. 3) then
			if (sterile(3)) then
				dampTermMatrixCoeffNue(1, 3) = nue_nue_nus
				dampTermMatrixCoeffNunu(1, 3) = nunu_nue_nus
				if (.not.sterile(2)) then
					dampTermMatrixCoeffNue(2, 3) = nue_nux_nus
					dampTermMatrixCoeffNunu(2, 3) = nunu_nux_nus
				end if
			else
				dampTermMatrixCoeffNue(1, 3) = nue_nue_nux
				dampTermMatrixCoeffNunu(1, 3) = nunu_nue_nux
				dampTermMatrixCoeffNue(2, 3) = nue_numu_nutau
				dampTermMatrixCoeffNunu(2, 3) = nunu_numu_nutau
			end if
		end if
		do ix=4, flavorNumber
			dampTermMatrixCoeffNue(1, ix) = nue_nue_nus
			dampTermMatrixCoeffNunu(1, ix) = nunu_nue_nus
			do iy=2, 3
				if (.not.sterile(iy)) then
					dampTermMatrixCoeffNue(iy, ix) = nue_nux_nus
					dampTermMatrixCoeffNunu(iy, ix) = nunu_nux_nus
				end if
			end do
		end do

		if (collint_offdiag_damping .and. collint_damping_type.eq.1) then
			!formulas from Bennett:2020zkv
			write(*,*) "[collint] Example d(y) for damping factors a la Bennett:2020zkv..."
			write(*,"(3A14)") "y", "f_eq(y)", "d(y)"
			do ix=1, Ny
				write(*,"(3E14.6)") y_arr(ix), feq_arr(ix), dy_damping_fit(y_arr(ix)) * y_arr(ix)**3
			end do
		end if
	end subroutine setDampingFactorCoeffs

	!phase space
	pure function F_ab_ann_re(n1, n2, f3, f4, a, b, i, j)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.5
		real(dl) :: F_ab_ann_re
		type(cmplxMatNN), intent(in) :: n1, n2
		real(dl), intent(in) :: f3, f4
		integer, intent(in) :: a, b, i, j
		integer :: k
		real(dl) :: t1a, t1b, t2a, t2b
#ifdef FULL_F_AB
		integer :: l, m
#endif

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_ann_re] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
#ifdef FULL_F_AB
		do k=1, flavorNumber
		do l=1, flavorNumber
		do m=1, flavorNumber
			!Ga (1-rho2) Gb (1-rho1)
			t1a = t1a + GLR_vec(a,i,l) * GLR_vec(b,k,m) * ( &
				(idMat(l,k)-n2%re(l,k))*(idMat(m,j)-n1%re(m,j)) - &
				(-n2%im(l,k))*(-n1%im(m,j))&	!idMat has im=0
				)
			!(1-rho1) Gb (1-rho2) Ga
			t1b = t1b + GLR_vec(a,m,j) * GLR_vec(b,l,k) * ( &
				(idMat(i,l)-n1%re(i,l))*(idMat(k,m)-n2%re(k,m)) - &
				(-n1%im(i,l))*(-n2%im(k,m))&
				)
			!rho1 Gb rho2 Ga
			t2a = t2a + GLR_vec(a,m,j) * GLR_vec(b,l,k) * ( &
				(n1%re(i,l))*(n2%re(k,m)) - &
				(n1%im(i,l))*(n2%im(k,m))&
				)
			!Ga Rho2 Gb rho1
			t2b = t2b + GLR_vec(a,i,l) * GLR_vec(b,k,m) * ( &
				(n2%re(l,k))*(n1%re(m,j)) - &
				(n2%im(l,k))*(n1%im(m,j))&
				)
		end do
		end do
		end do
#else
		do k=1, flavorNumber
			!Ga (1-rho2) Gb (1-rho1)
			t1a = t1a + GLR_vec(b,k,k) * ( &
				(idMat(i,k)-n2%re(i,k))*(idMat(k,j)-n1%re(k,j)) - &
				(-n2%im(i,k))*(-n1%im(k,j))&	!idMat ha im=0
				)
			!(1-rho1) Gb (1-rho2) Ga
			t1b = t1b + GLR_vec(b,k,k) * ( &
				(idMat(i,k)-n1%re(i,k))*(idMat(k,j)-n2%re(k,j)) - &
				(-n1%im(i,k))*(-n2%im(k,j))&
				)
			!rho1 Gb rho2 Ga
			t2a = t2a + GLR_vec(b,k,k) * ( &
				(n1%re(i,k))*(n2%re(k,j)) - &
				(n1%im(i,k))*(n2%im(k,j))&
				)
			!Ga Rho2 Gb rho1
			t2b = t2b + GLR_vec(b,k,k) * ( &
				(n2%re(i,k))*(n1%re(k,j)) - &
				(n2%im(i,k))*(n1%im(k,j))&
				)
		end do
		t1a = t1a * GLR_vec(a,i,i)
		t1b = t1b * GLR_vec(a,j,j)
		t2a = t2a * GLR_vec(a,j,j)
		t2b = t2b * GLR_vec(a,i,i)
#endif
		F_ab_ann_re = f3 * f4 * (t1a + t1b)&
			- (1.d0 - f3) * (1.d0 - f4) * (t2a + t2b)
	end function F_ab_ann_re

	pure function F_ab_ann_im(n1, n2, f3, f4, a, b, i, j)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.5
		real(dl) :: F_ab_ann_im
		type(cmplxMatNN), intent(in) :: n1, n2
		real(dl), intent(in) :: f3, f4
		integer, intent(in) :: a, b, i, j
		integer :: k
		real(dl) :: t1a, t1b, t2a, t2b
#ifdef FULL_F_AB
		integer :: l, m
#endif

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_ann_im] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
#ifdef FULL_F_AB
		do k=1, flavorNumber
		do l=1, flavorNumber
		do m=1, flavorNumber
			!Ga (1-rho2) Gb (1-rho1)
			t1a = t1a + GLR_vec(a,i,l) * GLR_vec(b,k,m) * ( &
				(-n2%im(l,k))*(idMat(m,j)-n1%re(m,j)) + &
				(idMat(l,k)-n2%re(l,k))*(-n1%im(m,j))&
				)
			!(1-rho1) Gb (1-rho2) Ga
			t1b = t1b + GLR_vec(a,m,j) * GLR_vec(b,l,k) * ( &
				(-n1%im(i,l))*(idMat(k,m)-n2%re(k,m)) + &
				(idMat(i,l)-n1%re(i,l))*(-n2%im(k,m))&
				)
			!rho1 Gb rho2 Ga
			t2a = t2a + GLR_vec(a,m,j) * GLR_vec(b,l,k) * ( &
				(n1%im(i,l))*(n2%re(k,m)) + &
				(n1%re(i,l))*(n2%im(k,m))&
				)
			!Ga Rho2 Gb rho1
			t2b = t2b + GLR_vec(a,i,l) * GLR_vec(b,k,m) * ( &
				(n2%im(l,k))*(n1%re(m,j)) + &
				(n2%re(l,k))*(n1%im(m,j))&
				)
		end do
		end do
		end do
#else
		do k=1, flavorNumber
			!Ga (1-rho2) Gb (1-rho1)
			t1a = t1a + GLR_vec(b,k,k) * ( &
				(-n2%im(i,k))*(idMat(k,j)-n1%re(k,j)) + &
				(idMat(i,k)-n2%re(i,k))*(-n1%im(k,j))&
				)
			!(1-rho1) Gb (1-rho2) Ga
			t1b = t1b + GLR_vec(b,k,k) * ( &
				(-n1%im(i,k))*(idMat(k,j)-n2%re(k,j)) + &
				(idMat(i,k)-n1%re(i,k))*(-n2%im(k,j))&
				)
			!rho1 Gb rho2 Ga
			t2a = t2a + GLR_vec(b,k,k) * ( &
				(n1%im(i,k))*(n2%re(k,j)) + &
				(n1%re(i,k))*(n2%im(k,j))&
				)
			!Ga Rho2 Gb rho1
			t2b = t2b + GLR_vec(b,k,k) * ( &
				(n2%im(i,k))*(n1%re(k,j)) + &
				(n2%re(i,k))*(n1%im(k,j))&
				)
		end do
		t1a = t1a * GLR_vec(a,i,i)
		t1b = t1b * GLR_vec(a,j,j)
		t2a = t2a * GLR_vec(a,j,j)
		t2b = t2b * GLR_vec(a,i,i)
#endif
		F_ab_ann_im = f3 * f4 * (t1a + t1b)&
			- (1.d0 - f3) * (1.d0 - f4) * (t2a + t2b)
	end function F_ab_ann_im

	pure function F_ab_sc_re (n1, n3, f2, f4, a, b, i, j)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.10
		real(dl) :: F_ab_sc_re
		type(cmplxMatNN), intent(in) :: n1, n3
		real(dl), intent(in) :: f2, f4
		integer, intent(in) :: a, b, i, j
		integer :: k
		real(dl) :: t1a, t1b, t2a, t2b
#ifdef FULL_F_AB
		integer :: l, m
#endif

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_sc_re] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
#ifdef FULL_F_AB
		do k=1, flavorNumber
		do l=1, flavorNumber
		do m=1, flavorNumber
			!Ga rho3 Gb (1-rho1)
			t1a = t1a + GLR_vec(a,i,l) * GLR_vec(b,k,m) * ( &
				(n3%re(l,k))*(idMat(m,j)-n1%re(m,j)) - &
				(n3%im(l,k))*(-n1%im(m,j))&	!idMat has im=0
				)
			!(1-rho1) Gb rho3 Ga
			t1b = t1b + GLR_vec(a,m,j) * GLR_vec(b,l,k) * ( &
				(idMat(i,l)-n1%re(i,l))*(n3%re(k,m)) - &
				(-n1%im(i,l))*(n3%im(k,m))&
				)
			!rho1 Gb (1-rho3) Ga
			t2a = t2a + GLR_vec(a,m,j) * GLR_vec(b,l,k) * ( &
				(n1%re(i,l))*(idMat(k,m)-n3%re(k,m)) - &
				(n1%im(i,l))*(-n3%im(k,m))&
				)
			!Ga (1-rho3) Gb rho1
			t2b = t2b + GLR_vec(a,i,l) * GLR_vec(b,k,m) * ( &
				(idMat(l,k)-n3%re(l,k))*(n1%re(m,j)) - &
				(-n3%im(l,k))*(n1%im(m,j))&
				)
		end do
		end do
		end do
#else
		do k=1, flavorNumber
			!Ga rho3 Gb (1-rho1)
			t1a = t1a + GLR_vec(b,k,k) * ( &
				(n3%re(i,k))*(idMat(k,j)-n1%re(k,j)) - &
				(n3%im(i,k))*(-n1%im(k,j))&	!idMat ha im=0
				)
			!(1-rho1) Gb rho3 Ga
			t1b = t1b + GLR_vec(b,k,k) * ( &
				(idMat(i,k)-n1%re(i,k))*(n3%re(k,j)) - &
				(-n1%im(i,k))*(n3%im(k,j))&
				)
			!rho1 Gb (1-rho3) Ga
			t2a = t2a + GLR_vec(b,k,k) * ( &
				(n1%re(i,k))*(idMat(k,j)-n3%re(k,j)) - &
				(n1%im(i,k))*(-n3%im(k,j))&
				)
			!Ga (1-rho3) Gb rho1
			t2b = t2b + GLR_vec(b,k,k) * ( &
				(idMat(i,k)-n3%re(i,k))*(n1%re(k,j)) - &
				(-n3%im(i,k))*(n1%im(k,j))&
				)
		end do
		t1a = t1a * GLR_vec(a,i,i)
		t1b = t1b * GLR_vec(a,j,j)
		t2a = t2a * GLR_vec(a,j,j)
		t2b = t2b * GLR_vec(a,i,i)
#endif
		F_ab_sc_re = (1.d0 - f2) * f4 * (t1a + t1b)&
			- f2 * (1.d0 - f4) * (t2a + t2b)
	end function F_ab_sc_re

	pure function F_ab_sc_im (n1, n3, f2, f4, a, b, i, j)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.10
		real(dl) :: F_ab_sc_im
		type(cmplxMatNN), intent(in) :: n1, n3
		real(dl), intent(in) :: f2, f4
		integer, intent(in) :: a, b, i, j
		integer :: k
		real(dl) :: t1a, t1b, t2a, t2b
#ifdef FULL_F_AB
		integer :: l, m
#endif

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_sc_im] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
#ifdef FULL_F_AB
		do k=1, flavorNumber
		do l=1, flavorNumber
		do m=1, flavorNumber
			!Ga rho3 Gb (1-rho1)
			t1a = t1a + GLR_vec(a,i,l) * GLR_vec(b,k,m) * ( &
				(n3%im(l,k))*(idMat(m,j)-n1%re(m,j)) + &
				(n3%re(l,k))*(-n1%im(m,j))&
				)
			!(1-rho1) Gb rho3 Ga
			t1b = t1b + GLR_vec(a,m,j) * GLR_vec(b,l,k) * ( &
				(-n1%im(i,l))*(n3%re(k,m)) + &
				(idMat(i,l)-n1%re(i,l))*(n3%im(k,m))&
				)
			!rho1 Gb (1-rho3) Ga
			t2a = t2a + GLR_vec(a,m,j) * GLR_vec(b,l,k) * ( &
				(n1%im(i,l))*(idMat(k,m)-n3%re(k,m)) + &
				(n1%re(i,l))*(-n3%im(k,m))&
				)
			!Ga (1-rho3) Gb rho1
			t2b = t2b + GLR_vec(a,i,l) * GLR_vec(b,k,m) * ( &
				(-n3%im(l,k))*(n1%re(m,j)) + &
				(idMat(l,k)-n3%re(l,k))*(n1%im(m,j))&
				)
		end do
		end do
		end do
#else
		do k=1, flavorNumber
			!Ga rho3 Gb (1-rho1)
			t1a = t1a + GLR_vec(b,k,k) * ( &
				(n3%im(i,k))*(idMat(k,j)-n1%re(k,j)) + &
				(n3%re(i,k))*(-n1%im(k,j))&
				)
			!(1-rho1) Gb rho3 Ga
			t1b = t1b + GLR_vec(b,k,k) * ( &
				(-n1%im(i,k))*(n3%re(k,j)) + &
				(idMat(i,k)-n1%re(i,k))*(n3%im(k,j))&
				)
			!rho1 Gb (1-rho3) Ga
			t2a = t2a + GLR_vec(b,k,k) * ( &
				(n1%im(i,k))*(idMat(k,j)-n3%re(k,j)) + &
				(n1%re(i,k))*(-n3%im(k,j))&
				)
			!Ga (1-rho3) Gb rho1
			t2b = t2b + GLR_vec(b,k,k) * ( &
				(-n3%im(i,k))*(n1%re(k,j)) + &
				(idMat(i,k)-n3%re(i,k))*(n1%im(k,j))&
				)
		end do
		t1a = t1a * GLR_vec(a,i,i)
		t1b = t1b * GLR_vec(a,j,j)
		t2a = t2a * GLR_vec(a,j,j)
		t2b = t2b * GLR_vec(a,i,i)
#endif
		F_ab_sc_im = (1.d0 - f2) * f4 * (t1a + t1b)&
			- f2 * (1.d0 - f4) * (t2a + t2b)
	end function F_ab_sc_im

	!phase space functions for nunu interactions
	!scattering
#ifdef FULL_F_NU
	pure subroutine F_nu_sc_1324(n1, n2, n3, n4, i, j, t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i, tr2, tr4)
		!...
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		real(dl), dimension(:), allocatable, intent(out) :: t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i
		real(dl), intent(out) :: tr2, tr4
		integer :: k, m
		real(dl) :: s13r, s13i, s24r, s24i, d24r

		allocate(t1r(flavorNumber), t3r(flavorNumber)) !t?r(m) -> t?r(i,m)
		allocate(t1i(flavorNumber), t3i(flavorNumber))
		allocate(t2r(flavorNumber), t4r(flavorNumber)) !t?r(m) -> t?r(m,j)
		allocate(t2i(flavorNumber), t4i(flavorNumber))
		tr2 = 0.d0
		tr4 = 0.d0
		do m=1, flavorNumber
			s13r = 0.d0
			s13i = 0.d0
			s24r = 0.d0
			s24i = 0.d0
			d24r = 0.d0
			!compute some terms that are repeated in terms 1/3 and 2/4, to save time
			do k=1, flavorNumber
				s13r = s13r + n1%im(i,k)*n3%im(k,m) - n1%re(i,k)*n3%re(k,m)
				s13i = s13i - n1%re(i,k)*n3%im(k,m) - n1%im(i,k)*n3%re(k,m)
				s24r = s24r + n2%im(m,k)*n4%im(k,j) - n2%re(m,k)*n4%re(k,j)
				s24i = s24i - n2%re(m,k)*n4%im(k,j) - n2%im(m,k)*n4%re(k,j)
				d24r = d24r + n2%im(m,k)*n4%im(k,m) - n2%re(m,k)*n4%re(k,m)
			end do
			!compute the products of two matrices
			t1r(m) = n3%re(i,m) + s13r
			t1i(m) = n3%im(i,m) + s13i
			t3r(m) = n1%re(i,m) + s13r
			t3i(m) = n1%im(i,m) + s13i
			t2r(m) = n4%re(m,j) + s24r
			t2i(m) = n4%im(m,j) + s24i
			t4r(m) = n2%re(m,j) + s24r
			t4i(m) = n2%im(m,j) + s24i
			!compute Tr()
			tr2 = tr2 + n4%re(m,m) + d24r
			tr4 = tr4 + n2%re(m,m) + d24r
		end do
	end subroutine F_nu_sc_1324

	pure function F_nu_sc_re_prod(n1, n2, n3, n4, i, j)
		!...
		real(dl) :: F_nu_sc_re_prod
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		integer :: k, m
		real(dl), dimension(:), allocatable :: t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i
		real(dl) :: s13r, s13i, s24r, s24i, d24r
		real(dl) :: tr2, tr4

		call F_nu_sc_1324(n1, n2, n3, n4, i, j, t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i, tr2, tr4)
		F_nu_sc_re_prod = &
			t1r(j) * tr2 + sum(t1r(:)*t2r(:) - t1i(:)*t2i(:)) &
			- (t3r(j) * tr4 + sum(t3r(:)*t4r(:) - t3i(:)*t4i(:)))
		deallocate(t1r, t2r, t3r, t4r)
		deallocate(t1i, t2i, t3i, t4i)
	end function F_nu_sc_re_prod

	pure function F_nu_sc_im_prod(n1, n2, n3, n4, i, j)
		!...
		real(dl) :: F_nu_sc_im_prod
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		integer :: k, m
		real(dl), dimension(:), allocatable :: t1r, t2r, t3r, t4r
		real(dl), dimension(:), allocatable :: t1i, t2i, t3i, t4i
		real(dl) :: s13r, s13i, s24r, s24i, d24r
		real(dl) :: tr2, tr4

		F_nu_sc_im_prod = 0.d0
		!there can be no imaginary part in the diagonal
		if (i.eq.j) &
			return

		call F_nu_sc_1324(n1, n2, n3, n4, i, j, t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i, tr2, tr4)
		F_nu_sc_im_prod = &
			t1i(j) * tr2 + sum(t1r(:)*t2i(:) + t1i(:)*t2r(:)) &
			- (t3i(j) * tr4 + sum(t3r(:)*t4i(:) + t3i(:)*t4r(:)))
		deallocate(t1r, t2r, t3r, t4r)
		deallocate(t1i, t2i, t3i, t4i)
	end function F_nu_sc_im_prod
#endif

	pure function F_nu_sc_re(n1, n2, n3, n4, i, j)
		!...
		real(dl) :: F_nu_sc_re
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		integer :: k
#ifndef FULL_F_NU
		real(dl), dimension(:), allocatable :: t2r, t4r
		real(dl) :: s13r, s24r
#endif

		F_nu_sc_re = 0.d0

#ifdef FULL_F_NU
		F_nu_sc_re = F_nu_sc_re_prod(n1, n2, n3, n4, i, j)
		if (i.eq.j) then ! +h.c. is the same for diagonal entries
			F_nu_sc_re = 2.d0 * F_nu_sc_re
			return
		end if
		!now repeat to compute the h.c. (off-diagonal entries)
		F_nu_sc_re = F_nu_sc_re + F_nu_sc_re_prod(n1, n2, n3, n4, j, i)
#else
		!only consider diagonal rho, it's easier because all the products involve diagonal matrices only
		!offdiagonal elements vanish here
		if (i.ne.j) &
			return
		allocate(t2r(flavorNumber), t4r(flavorNumber))
		do k=1, flavorNumber
			s24r = n2%re(k, k)*n4%re(k, k)
			t2r(k) = n4%re(k, k) - s24r
			t4r(k) = n2%re(k, k) - s24r
		end do
		s13r = n1%re(i, i) * n3%re(i, i)
		F_nu_sc_re = &
			(n3%re(i,i) - s13r) * (t2r(i) + sum(t2r)) &
			- (n1%re(i,i) - s13r) * (t4r(i) + sum(t4r))
		deallocate(t2r, t4r)
		F_nu_sc_re = 2.d0 * F_nu_sc_re ! +h.c.
#endif
	end function F_nu_sc_re

	pure function F_nu_sc_im(n1, n2, n3, n4, i, j)
		!...
		real(dl) :: F_nu_sc_im
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j

		F_nu_sc_im = 0.d0

#ifdef FULL_F_NU
		!there can be no imaginary part in the diagonal
		if (i.eq.j) &
			return

		F_nu_sc_im = F_nu_sc_im_prod(n1, n2, n3, n4, i, j) - F_nu_sc_im_prod(n1, n2, n3, n4, j, i)
#endif
	end function F_nu_sc_im

	!pair
#ifdef FULL_F_NU
	pure subroutine F_nu_pa_1234(n1, n2, n3, n4, i, j, t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i, tr2, tr4)
		!...
		!first line: sAr -> 12, sB->34
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		real(dl), dimension(:), allocatable, intent(inout) :: t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i
		real(dl), intent(out) :: tr2, tr4
		integer :: k, m
		real(dl) :: sAr, sAi, sBr, sBi, sBd

		tr2 = 0.d0
		tr4 = 0.d0
		do m=1, flavorNumber
			sAr = 0.d0
			sAi = 0.d0
			sBr = 0.d0
			sBi = 0.d0
			sBd = 0.d0
			!compute some terms that are repeated in terms 1/3 and 2/4, to save time
			do k=1, flavorNumber
				sAr = sAr + n1%re(i,k)*n2%re(k,m) - n1%im(i,k)*n2%im(k,m)
				sAi = sAi + n1%im(i,k)*n2%re(k,m) + n1%re(i,k)*n2%im(k,m)
				sBr = sBr + n4%re(m,k)*n3%re(k,j) - n4%im(m,k)*n3%im(k,j)
				sBi = sBi + n4%im(m,k)*n3%re(k,j) + n4%re(m,k)*n3%im(k,j)
				sBd = sBd + n4%re(m,k)*n3%re(k,m) - n4%im(m,k)*n3%im(k,m)
			end do
			!compute the products of two matrices
			t1r(m) = sAr + idMat(i,m) - n1%re(i,m) - n2%re(i,m)
			t1i(m) = sAi - n1%im(i,m) - n2%im(i,m)
			t3r(m) = sAr
			t3i(m) = sAi
			t2r(m) = sBr
			t2i(m) = sBi
			t4r(m) = sBr + idMat(m,j) - n4%re(m,j) - n3%re(m,j)
			t4i(m) = sBi - n4%im(m,j) - n3%im(m,j)
			!compute Tr()
			tr2 = tr2 + sBd
			tr4 = tr4 + sBd + idMat(m,m) - n4%re(m,m) - n3%re(m,m)
		end do
	end subroutine F_nu_pa_1234

	pure subroutine F_nu_pa_1342(n1, n2, n3, n4, i, j, t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i, tr2, tr4)
		!...
		!second line: sAr -> 13, sB->24
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		real(dl), dimension(:), allocatable, intent(inout) :: t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i
		real(dl), intent(out) :: tr2, tr4
		integer :: k, m
		real(dl) :: sAr, sAi, sBr, sBi, sBd

		tr2 = 0.d0
		tr4 = 0.d0
		do m=1, flavorNumber
			sAr = 0.d0
			sAi = 0.d0
			sBr = 0.d0
			sBi = 0.d0
			sBd = 0.d0
			!compute some terms that are repeated in terms 1/3 and 2/4, to save time
			do k=1, flavorNumber
				sAr = sAr + n1%im(i,k)*n3%im(k,m) - n1%re(i,k)*n3%re(k,m)
				sAi = sAi - n1%re(i,k)*n3%im(k,m) - n1%im(i,k)*n3%re(k,m)
				sBr = sBr + n4%im(m,k)*n2%im(k,j) - n4%re(m,k)*n2%re(k,j)
				sBi = sBi - n4%re(m,k)*n2%im(k,j) - n4%im(m,k)*n2%re(k,j)
				sBd = sBd + n4%im(m,k)*n2%im(k,m) - n4%re(m,k)*n2%re(k,m)
			end do
			!compute the products of two matrices
			t1r(m) = n3%re(i,m) + sAr
			t1i(m) = n3%im(i,m) + sAi
			t3r(m) = n1%re(i,m) + sAr
			t3i(m) = n1%im(i,m) + sAi
			t2r(m) = n4%re(m,j) + sBr
			t2i(m) = n4%im(m,j) + sBi
			t4r(m) = n2%re(m,j) + sBr
			t4i(m) = n2%im(m,j) + sBi
			!compute Tr()
			tr2 = tr2 + n4%re(m,m) + sBd
			tr4 = tr4 + n2%re(m,m) + sBd
		end do
	end subroutine F_nu_pa_1342

	pure function F_nu_pa_re_prod(n1, n2, n3, n4, i, j)
		!...
		real(dl) :: F_nu_pa_re_prod
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		integer :: k, m
		real(dl), dimension(:), allocatable :: t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i
		real(dl) :: sAr, sAi, sBr, sBi, sBd
		real(dl) :: tr2, tr4

		F_nu_pa_re_prod = 0.d0

		allocate(t1r(flavorNumber), t3r(flavorNumber)) !t?r(m) -> t?r(i,m)
		allocate(t1i(flavorNumber), t3i(flavorNumber)) !t?r(m) -> t?r(i,m)
		allocate(t2r(flavorNumber), t4r(flavorNumber)) !t?r(m) -> t?r(m,j)
		allocate(t2i(flavorNumber), t4i(flavorNumber)) !t?r(m) -> t?r(m,j)

		!first line: sAr -> 12, sB->34
		call F_nu_pa_1234(n1, n2, n3, n4, i, j, t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i, tr2, tr4)
		F_nu_pa_re_prod = F_nu_pa_re_prod &
			+ t1r(j) * tr2 + sum(t1r(:)*t2r(:) - t1i(:)*t2i(:)) &
			- (t3r(j) * tr4 + sum(t3r(:)*t4r(:) - t3i(:)*t4i(:)))

		!second line: sAr -> 13, sB->24
		call F_nu_pa_1342(n1, n2, n3, n4, i, j, t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i, tr2, tr4)
		F_nu_pa_re_prod = F_nu_pa_re_prod &
			+ t1r(j) * tr2 + sum(t1r(:)*t2r(:) - t1i(:)*t2i(:)) &
			- (t3r(j) * tr4 + sum(t3r(:)*t4r(:) - t3i(:)*t4i(:)))

		deallocate(t1r, t2r, t3r, t4r)
		deallocate(t1i, t2i, t3i, t4i)
	end function F_nu_pa_re_prod

	pure function F_nu_pa_im_prod(n1, n2, n3, n4, i, j)
		!...
		real(dl) :: F_nu_pa_im_prod
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		integer :: k, m
		real(dl), dimension(:), allocatable :: t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i
		real(dl) :: sAr, sAi, sBr, sBi, sBd
		real(dl) :: tr2, tr4

		F_nu_pa_im_prod = 0.d0

		!there can be no imaginary part in the diagonal
		if (i.eq.j) &
			return

		allocate(t1r(flavorNumber), t3r(flavorNumber)) !t?r(m) -> t?r(i,m)
		allocate(t1i(flavorNumber), t3i(flavorNumber)) !t?r(m) -> t?r(i,m)
		allocate(t2r(flavorNumber), t4r(flavorNumber)) !t?r(m) -> t?r(m,j)
		allocate(t2i(flavorNumber), t4i(flavorNumber)) !t?r(m) -> t?r(m,j)

		!first line: sAr -> 12, sB->34
		call F_nu_pa_1234(n1, n2, n3, n4, i, j, t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i, tr2, tr4)
		F_nu_pa_im_prod = F_nu_pa_im_prod &
			+ t1i(j) * tr2 + sum(t1r(:)*t2i(:) + t1i(:)*t2r(:)) &
			- (t3i(j) * tr4 + sum(t3r(:)*t4i(:) + t3i(:)*t4r(:)))

		!second line: sAr -> 13, sB->24
		call F_nu_pa_1342(n1, n2, n3, n4, i, j, t1r, t2r, t3r, t4r, t1i, t2i, t3i, t4i, tr2, tr4)
		F_nu_pa_im_prod = F_nu_pa_im_prod &
			+ t1i(j) * tr2 + sum(t1r(:)*t2i(:) + t1i(:)*t2r(:)) &
			- (t3i(j) * tr4 + sum(t3r(:)*t4i(:) + t3i(:)*t4r(:)))

		deallocate(t1r, t2r, t3r, t4r)
		deallocate(t1i, t2i, t3i, t4i)
	end function F_nu_pa_im_prod
#endif

	pure function F_nu_pa_re(n1, n2, n3, n4, i, j)
		!...
		real(dl) :: F_nu_pa_re
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j
		integer :: k
#ifndef FULL_F_NU
		real(dl), dimension(:), allocatable :: t2r, t4r
		real(dl) :: sAr, sBr
#endif

		F_nu_pa_re = 0.d0

#ifdef FULL_F_NU
		F_nu_pa_re = F_nu_pa_re_prod(n1, n2, n3, n4, i, j)
		if (i.eq.j) then ! +h.c. is the same for diagonal entries
			F_nu_pa_re = 2.d0 * F_nu_pa_re
			return
		end if
		!now repeat to compute the h.c. (off-diagonal entries)
		F_nu_pa_re = F_nu_pa_re + F_nu_pa_re_prod(n1, n2, n3, n4, j, i)
#else
		!only consider diagonal rho, it's easier because all the products involve diagonal matrices only
		!offdiagonal elements vanish here
		if (i.ne.j) &
			return

		allocate(t2r(flavorNumber), t4r(flavorNumber))
		!first line: sAr -> 12, sBr->34
		do k=1, flavorNumber
			sBr = n4%re(k, k)*n3%re(k, k)
			t2r(k) = sBr
			t4r(k) = sBr + 1.d0 - n4%re(k, k) - n3%re(k, k)
		end do
		sAr = n1%re(i, i) * n2%re(i, i)
		F_nu_pa_re = F_nu_pa_re &
			+ (sAr + 1.d0 - n1%re(i,i) - n2%re(i,i)) * (t2r(i) + sum(t2r)) &
			- sAr * (t4r(i) + sum(t4r))
		!second line: sAr -> 13, sBr->24
		do k=1, flavorNumber
			sBr = n4%re(k, k)*n2%re(k, k)
			t2r(k) = n4%re(k, k) - sBr
			t4r(k) = n2%re(k, k) - sBr
		end do
		sAr = n1%re(i, i) * n3%re(i, i)
		F_nu_pa_re = F_nu_pa_re &
			+ (n3%re(i,i) - sAr) * (t2r(i) + sum(t2r)) &
			- (n1%re(i,i) - sAr) * (t4r(i) + sum(t4r))
		deallocate(t2r, t4r)
		F_nu_pa_re = 2.d0 * F_nu_pa_re ! +h.c.
#endif
	end function F_nu_pa_re

	pure function F_nu_pa_im(n1, n2, n3, n4, i, j)
		!...
		real(dl) :: F_nu_pa_im
		type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
		integer, intent(in) :: i, j

		F_nu_pa_im = 0.d0

#ifdef FULL_F_NU
		!there can be no imaginary part in the diagonal
		if (i.eq.j) &
			return

		F_nu_pa_im = F_nu_pa_im_prod(n1, n2, n3, n4, i, j) - F_nu_pa_im_prod(n1, n2, n3, n4, j, i)
#endif
	end function F_nu_pa_im

	!D functions
	elemental function D1_full(y1, y2, y3, y4)
	!10.1103/PhysRevD.94.033009 eq.D1
		implicit none
		real(dl) :: D1_full
		real(dl), intent(in) :: y1, y2, y3, y4

		D1_full = &
			abs( y1+y2+y3-y4) + &
			abs( y1+y2-y3+y4) + &
			abs( y1-y2+y3+y4) + &
			abs(-y1+y2+y3+y4) - &
			abs( y1+y2-y3-y4) - &
			abs( y1-y2-y3+y4) - &
			abs( y1-y2+y3-y4) - &
			abs( y1+y2+y3+y4)
	end function D1_full

	elemental function D2_full(y1, y2, y3, y4)
	!10.1103/PhysRevD.94.033009 eq.D2
		implicit none
		real(dl) :: D2_full
		real(dl), intent(in) :: y1, y2, y3, y4
		real(dl) :: &
			p1p2p3p4, p1p2m3m4, p1m2m3p4, p1m2p3m4, &
			m1p2p3p4, p1m2p3p4, p1p2m3p4, p1p2p3m4
		real(dl) :: a, b
		real(dl) :: &
			Ap1p2m3m4, Ap1m2p3m4, Ap1p2p3m4, Ap1m2m3p4, &
			Ap1p2m3p4, Ap1m2p3p4, Am1p2p3p4, Ap1p2p3p4

		p1p2m3m4 =  y1 +y2 -y3 -y4 
		p1m2m3p4 =  y1 -y2 -y3 +y4 
		p1m2p3m4 =  y1 -y2 +y3 -y4 

		m1p2p3p4 = -y1 +y2 +y3 +y4
		p1m2p3p4 =  y1 -y2 +y3 +y4
		p1p2m3p4 =  y1 +y2 -y3 +y4
		p1p2p3m4 =  y1 +y2 +y3 -y4

		p1p2p3p4 =  y1 +y2 +y3 +y4 

		Ap1p2m3m4 = Abs(p1p2m3m4)
		Ap1m2p3m4 = Abs(p1m2p3m4)
		Ap1p2p3m4 = Abs(p1p2p3m4)
		Ap1m2m3p4 = Abs(p1m2m3p4)
		Ap1p2m3p4 = Abs(p1p2m3p4)
		Ap1m2p3p4 = Abs(p1m2p3p4)
		Am1p2p3p4 = Abs(m1p2p3p4)
		Ap1p2p3p4 = Abs(p1p2p3p4)

		a =   Sign(p1p2m3m4**2, p1p2m3m4) - Sign(p1p2p3m4**2, p1p2p3m4) - &
			  Sign(p1p2m3p4**2, p1p2m3p4) + Sign(p1p2p3p4**2, p1p2p3p4)
		b = - Sign(m1p2p3p4**2,-m1p2p3p4) + Sign(p1m2p3m4**2, p1m2p3m4) &
			+ Sign(p1m2m3p4**2, p1m2m3p4) - Sign(p1m2p3p4**2, p1m2p3p4)

		D2_full = - ( &
			Ap1p2m3m4**3 + Ap1m2p3m4**3 - &
			Ap1p2p3m4**3 + Ap1m2m3p4**3 - &
			Ap1p2m3p4**3 - Ap1m2p3p4**3 - &
			Am1p2p3p4**3 + Ap1p2p3p4**3 + &
			6.d0 * y1 * y2 * ( &
				Ap1p2m3m4 - Ap1m2p3m4 - &
				Ap1p2p3m4 - Ap1m2m3p4 - &
				Ap1p2m3p4 + Ap1m2p3p4 + &
				Am1p2p3p4 + Ap1p2p3p4 &
			) - &
			3.d0 * (y1 * ( a + b ) + y2 * ( a - b )) &
			) / 6.d0
	end function D2_full

	elemental function D3_full(y1, y2, y3, y4)
	!10.1103/PhysRevD.94.033009 eq.D3
		implicit none
		real(dl) :: D3_full
		real(dl), intent(in) :: y1, y2, y3, y4
		real(dl) :: &
			p1p2m3m4, p1m2m3p4, p1m2p3m4, &
			m1p2p3p4, p1m2p3p4, p1p2m3p4, p1p2p3m4
		real(dl) :: p123, p234, p12, p23, p34, m34
		real(dl) :: y1234, y12, y23, y34, y13, y24
		real(dl) :: &
			Ap1p2m3m4, Ap1m2p3m4, Ap1p2p3m4, Ap1m2m3p4, &
			Ap1p2m3p4, Ap1m2p3p4, Am1p2p3p4

		p1p2m3m4 =  y1 +y2 -y3 -y4
		p1m2m3p4 =  y1 -y2 -y3 +y4
		p1m2p3m4 =  y1 -y2 +y3 -y4

		m1p2p3p4 = -y1 +y2 +y3 +y4
		p1m2p3p4 =  y1 -y2 +y3 +y4
		p1p2m3p4 =  y1 +y2 -y3 +y4
		p1p2p3m4 =  y1 +y2 +y3 -y4

		p123 = y1+y2+y3
		p234 = y2+y3+y4
		p12 = y1+y2
		p23 = y2+y3
		p34 = y3+y4
		m34 = y3-y4

		y1234 = y1*y2*y3*y4
		y12 = y1*y2
		y13 = y1*y3
		y23 = y2*y3
		y24 = y2*y4
		y34 = y3*y4

		Ap1p2m3m4 = Abs(p1p2m3m4)
		Ap1m2p3m4 = Abs(p1m2p3m4)
		Ap1p2p3m4 = Abs(p1p2p3m4)
		Ap1m2m3p4 = Abs(p1m2m3p4)
		Ap1p2m3p4 = Abs(p1p2m3p4)
		Ap1m2p3p4 = Abs(p1m2p3p4)
		Am1p2p3p4 = Abs(m1p2p3p4)

		D3_full = &
			(     y1**5 + y2**5 &
				- 5*( &
					y1**2 * ( &
						(y2**2 + y3**2 + y4**2) * y1 + &
						(y2**3 + y3**3 + y4**3) ) + &
					y2**2 * ( &
						(y3**2 + y4**2) * y2 + &
						(y3**3 + y4**3) ) )&
				+  p34**3 * (y3**2 - 3*y34 + y4**2)&
			)/30.d0 + &
			( 	  am1p2p3p4**5 - ap1p2m3m4**5 - ap1m2p3m4**5 + ap1p2p3m4**5 &
				- ap1m2m3p4**5 + ap1p2m3p4**5 + ap1m2p3p4**5 &
			)/120.d0 + &
			(	- 6 * y1234 * ( &
					am1p2p3p4 + ap1m2m3p4 + ap1p2p3m4 + ap1m2p3m4 + &
					ap1p2m3p4 + ap1m2p3p4 + ap1p2m3m4 &
				) + &
				am1p2p3p4**3 * (y34 + y2*p34 - y1*p234) + &
				ap1m2m3p4**3 * (y1*(p23-y4) - y23 + y4*p23) + &
				ap1p2p3m4**3 * (y12 + y23 + y13 - y4*p123) + &
				ap1m2p3m4**3 * ( y2*(m34) + y34 + y1*(y2-m34)) + &
				ap1p2m3p4**3 * (-y2*(m34) - y34 + y1*(y2-m34)) + &
				ap1m2p3p4**3 * (y34 - y2*p34 + y1*(p34-y2)) + &
				ap1p2m3m4**3 * (y2*p34 - y34 + y1*(p34-y2)) &
			)/6.d0 + &
			(	sign(p1p2m3m4**2, p1p2m3m4) * ( &
					y1**3+y2**3+3*y1*(y2**2+y3**2+y4**2+6*(y34-p34*y2)) - &
					p34**3 + 3*(y2*(y3**2+y4**2+6*y34)-y2**2*p34+y1**2*(y2-p34))&
				) + &
				sign(p1p2m3p4**2, p1p2m3p4) * ( &
					- p1p2m3p4**3 + 12*(y12*m34+p12*y34) &
				) + &
				sign(p1m2p3p4**2, p1m2p3p4) * ( &
					- p1m2p3p4**3 + 12*(y12*p34+(y2-y1)*y34) &
				) + &
				sign(p1p2p3m4**2, p1p2p3m4) * ( &
					3*((y2**2+6*y23+y3**2)*y4 - (p23-y4)*y1**2 - p23*y4**2 - &
						(y2**2+6*y23+y3**2+y4**2-6*p23*y4)*y1) +&
					y4**3-y1**3-p23**3 &
				) + &
				sign(p1m2m3p4**2, p1m2m3p4) * ( &
					3*((y2**2+6*y23+y3**2)*y4 - (p23-y4)*y1**2 - p23*y4**2 + &
						(y2**2+6*y23+y3**2+y4**2-6*p23*y4)*y1) +&
					y1**3+y4**3-p23**3 &
				) + &
				sign(p1m2p3m4**2, p1m2p3m4) * ( &
					3*(-(y3**2-6*y34+y4**2)*y2 - (y2-m34)*y1**2 + m34*y2**2 + &
						(y2**2-6*(m34)*y2+y3**2+y4**2-6*y34)*y1) +&
					y1**3-y2**3+(m34)**3 &
				) + &
				sign(m1p2p3p4**2,-m1p2p3p4) * ( &
					3*((y3**2+6*y34+y4**2)*y2 + (p234)*y1**2 + p34*y2**2 - &
						(y2**2+6*p34*y2+y3**2+y4**2+6*y34)*y1) +&
					y2**3-y1**3+p34**3 &
				) &
			)/24.d0
	end function D3_full

	!Pi functions
	elemental function PI1_12_full(y1, y2, y3, y4) !(y1, y2)
		real(dl) :: PI1_12_full
		real(dl), intent(in) :: y1, y2, y3, y4

		!pi_1(y1, y2): nu+nu -> e+e
		PI1_12_full = y1 * y2 * D1_full(y1, y2, y3, y4) - D2_full(y1, y2, y3, y4)
	end function PI1_12_full

	elemental function PI1_13_full(y1, y2, y3, y4) !(y1, y3)
		real(dl) :: PI1_13_full
		real(dl), intent(in) :: y1, y2, y3, y4

		!pi_1(y1, y3): nu+e -> nu+e
		PI1_13_full = y1 * y3 * D1_full(y1, y2, y3, y4) + D2_full(y1, y3, y2, y4)
	end function PI1_13_full

	pure function PI2_nn_f(y1, y2, y3, y4, e3, e4) !1: (y1, y4),	2: (y1, y3)
		real(dl), dimension(2) :: PI2_nn_f
		real(dl), intent(in) :: y1, y2, y3, y4, e3, e4
		real(dl) :: e1, e2, t

		e1 = y1
		e2 = y2
		t = e1 * e2 * e3 * e4 * D1_full(y1, y2, y3, y4) + D3_full(y1, y2, y3, y4)

		!pi_2(y1, y4)
		PI2_nn_f(1) = 2.d0 * (t &
			+ e2 * e3 * D2_full(y1, y4, y2, y3) &
			+ e1 * e4 * D2_full(y2, y3, y1, y4))

		!pi_2(y1, y3)
		PI2_nn_f(2) = 2.d0 * (t &
			+ e1 * e3 * D2_full(y2, y4, y1, y3) &
			+ e2 * e4 * D2_full(y1, y3, y2, y4))
	end function PI2_nn_f

	pure function PI2_ne_f(y1, y2, y3, y4, e2, e4) !1: (y1, y4),	2: (y1, y2)
		real(dl), dimension(2) :: PI2_ne_f
		real(dl), intent(in) :: y1, y2, y3, y4, e2, e4
		real(dl) :: e1, e3, t

		e1 = y1
		e3 = y3
		t = e1 * e2 * e3 * e4 * D1_full(y1, y2, y3, y4) + D3_full(y1, y2, y3, y4)

		!pi_2(y1, y4)
		PI2_ne_f(1) = 2.d0 * (t &
			+ e2 * e3 * D2_full(y1, y4, y2, y3) &
			+ e1 * e4 * D2_full(y2, y3, y1, y4))

		!pi_2(y1, y2)
		PI2_ne_f(2) = 2.d0 * (t &
			- e1 * e2 * D2_full(y3, y4, y1, y2) &
			- e3 * e4 * D2_full(y1, y2, y3, y4))
	end function PI2_ne_f

	!integrands of collision terms
	pure function coll_nue_ann_int(iy2, y4, obj, F_ab)
		!annihilation
		use fpInterfaces1
		procedure (F_annihilation) :: F_ab
		integer, intent(in) :: iy2
		real(dl), intent(in) :: y4
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_ann_int
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y2, y3, f3, f4, E3, E4, dme2, t1, t2

		coll_nue_ann_int = 0.d0

		dme2 = obj%dme2
		y2 = y_arr(iy2)
		E4 = Ebare_i_dme(obj%x, y4, dme2)
		E3 = obj%y1 + y2 - E4
		t1 = E3*E3
		t2 = obj%x*obj%x + dme2
		if (t1<t2) then
			y3 = -1.d0
		else
			y3 = sqrt(t1 - t2)
		endif
		if (.not.(y3.lt.0.d0 &
				.or. obj%y1.gt.y3+y2+y4 &
				.or. y2.gt.obj%y1+y3+y4 &
				.or. y3.gt.obj%y1+y2+y4 &
				.or. y4.gt.obj%y1+y2+y3)) then
			f3 = fermiDirac(E3 / obj%z)
			f4 = fermiDirac(E4 / obj%z)
			pi2_vec = PI2_nn_f(obj%y1, y2, y3, y4, E3, E4)
			coll_nue_ann_int = coll_nue_ann_int + &
				y3/E3 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 1, 1, obj%ix1,obj%ix2) &
					+ pi2_vec(2) * F_ab(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 2, 2, obj%ix1,obj%ix2) &
					+ (obj%x*obj%x + dme2) * PI1_12_full(obj%y1, y2, y3, y4) * ( &
						F_ab(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 2, 1, obj%ix1,obj%ix2) &
						+ F_ab(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_ann_int

	pure function coll_nue_sc_int(iy3, y2, obj, F_ab)
		!scattering, summing positron and electrons
		use fpInterfaces1
		procedure (F_scattering) :: F_ab
		integer, intent(in) :: iy3
		real(dl), intent(in) :: y2
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_sc_int
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y3, y4, f2, f4, E2, E4, dme2, t1, t2

		coll_nue_sc_int = 0.d0

		dme2 = obj%dme2
		E2 = Ebare_i_dme(obj%x,y2, dme2)
		y3 = y_arr(iy3)
		E4 = obj%y1 + E2 - y3
		t1 = E4*E4
		t2 = obj%x*obj%x + dme2
		if (t1<t2) then
			y4 = -1.d0
		else
			y4 = sqrt(t1 - t2)
		endif
		if (.not.(y4.lt.0.d0 &
				.or. obj%y1.gt.y2+y3+y4 &
				.or. y2.gt.obj%y1+y3+y4 &
				.or. y3.gt.obj%y1+y2+y4 &
				.or. y4.gt.obj%y1+y2+y3)) then
			f2 = fermiDirac(E2 / obj%z)
			f4 = fermiDirac(E4 / obj%z)
			pi2_vec = PI2_ne_f (obj%y1, y2, y3, y4, E2, E4)
			coll_nue_sc_int = coll_nue_sc_int + &
				y2/E2 * &
				y4/E4 * &
				( &
					( pi2_vec(1) + pi2_vec(2) ) * ( & !F_sc^LL + F_sc^RR
						F_ab(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy3),f2,f4, 1, 1, obj%ix1,obj%ix2) &
						+ F_ab(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy3),f2,f4, 2, 2, obj%ix1,obj%ix2) &
					) &
					- 2.d0 * (obj%x*obj%x + dme2) * PI1_13_full(obj%y1, y2, y3, y4) * ( & !F_sc^RL and F_sc^LR
						F_ab(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy3),f2,f4, 2, 1, obj%ix1,obj%ix2) &
						+ F_ab(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy3),f2,f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_sc_int

	pure function coll_nue_int(iy, yx, obj, F_ab_ann, F_ab_sc)
		use fpInterfaces1
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		integer, intent(in) :: iy
		real(dl), intent(in) :: yx
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_int
#ifdef NO_NUE_ANNIHILATION
		coll_nue_int = &
			coll_nue_sc_int(iy, yx, obj, F_ab_sc)
#else
		coll_nue_int = &
			coll_nue_sc_int(iy, yx, obj, F_ab_sc) &
			+ coll_nue_ann_int(iy, yx, obj, F_ab_ann)
#endif
	end function coll_nue_int

	pure function coll_nunu_int(iy2, iy3, obj, F_nu_sc, F_nu_pa)
		!nunu interaction
		use fpInterfaces1
		real(dl) :: coll_nunu_int
		integer, intent(in) :: iy2, iy3
		type(coll_args), intent(in) :: obj
		procedure (Fnunu) :: F_nu_sc, F_nu_pa
		type(cmplxMatNN) :: n4
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y2, y3, y4

		coll_nunu_int = 0.d0

		y2 = y_arr(iy2)
		y3 = y_arr(iy3)
		y4 = obj%y1 + y2 - y3
		if (.not.(y4.lt.y_arr(1) &
				.or. y4.gt.y_arr(Ny) & !assume that if y4 is within the y_max range or there is no contribution...otherwise I should extrapolate! maybe later...
				.or. obj%y1.gt.y2+y3+y4 &
				.or. y2.gt.obj%y1+y3+y4 &
				.or. y3.gt.obj%y1+y2+y4 &
				.or. y4.gt.obj%y1+y2+y3) &
		) then
			pi2_vec = PI2_ne_f(obj%y1, y2, y3, y4, y2, y4)
			n4 = get_interpolated_nudens(nuDensMatVecFD, y4, flavorNumber, Ny)
			coll_nunu_int = coll_nunu_int &
				+ pi2_vec(1) * F_nu_pa(nuDensMatVecFD(obj%iy), nuDensMatVecFD(iy2), nuDensMatVecFD(iy3), n4, obj%ix1,obj%ix2) &
				+ pi2_vec(2) * F_nu_sc(nuDensMatVecFD(obj%iy), nuDensMatVecFD(iy2), nuDensMatVecFD(iy3), n4, obj%ix1,obj%ix2)
		end if
	end function coll_nunu_int

	!integrate collision terms
	pure function integrate_collint_nue_NC(f, obj, F_ab_ann, F_ab_sc)
		use fpInterfaces1
		use fpInterfaces2
		implicit None
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		procedure (collision_integrand_nue) :: f
		real(dl) :: integrate_collint_nue_NC
		type(coll_args), intent(in) :: obj
		integer :: ia, ib
		real(dl), dimension(:,:), allocatable :: fy2_arr
		real(dl), dimension(:), allocatable :: tmpfy2

		allocate(fy2_arr(Ny, Ny))
		fy2_arr = 0.d0
		do ia=1, Ny
			do ib=1, Ny
				fy2_arr(ia, ib) = f(ia, y_arr(ib), obj, F_ab_ann, F_ab_sc)
			end do
		end do
		integrate_collint_nue_NC = integral_NC_2d(Ny, Ny, dy_arr, dy_arr, fy2_arr)
		deallocate(fy2_arr)
	end function integrate_collint_nue_NC

	pure function integrate_collint_nue_GL(f, obj, F_ab_ann, F_ab_sc)
		use fpInterfaces1
		use fpInterfaces2
		implicit None
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		procedure (collision_integrand_nue) :: f
		real(dl) :: integrate_collint_nue_GL
		type(coll_args), intent(in) :: obj
		integer :: ia, ib
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(fy2_arr(Ny, Ny))
		fy2_arr = 0.d0
		do ia=1, Ny
			do ib=1, Ny
				fy2_arr(ia, ib) = f(ia, y_arr(ib), obj, F_ab_ann, F_ab_sc)
			end do
		end do
		integrate_collint_nue_GL = integral_GL_2d(Ny, w_gl_arr2, w_gl_arr2, fy2_arr)
		deallocate(fy2_arr)
	end function integrate_collint_nue_GL

	pure function integrate_collint_nunu_NC(f, obj, F_nu_sc, F_nu_pa)
		use fpInterfaces1
		use fpInterfaces2
		implicit None
		procedure (Fnunu) :: F_nu_sc, F_nu_pa
		procedure (collision_integrand_nunu) :: f
		real(dl) :: integrate_collint_nunu_NC
		type(coll_args), intent(in) :: obj
		integer :: ia, ib
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(fy2_arr(Ny, Ny))
		fy2_arr = 0.d0
		do ia=1, Ny
			do ib=1, Ny
				fy2_arr(ia, ib) = f(ia, ib, obj, F_nu_sc, F_nu_pa)
			end do
		end do
		integrate_collint_nunu_NC = integral_NC_2d(Ny, Ny, dy_arr, dy_arr, fy2_arr)
		deallocate(fy2_arr)
	end function integrate_collint_nunu_NC

	pure function integrate_collint_nunu_GL(f, obj, F_nu_sc, F_nu_pa)
		use fpInterfaces1
		use fpInterfaces2
		implicit None
		procedure (Fnunu) :: F_nu_sc, F_nu_pa
		procedure (collision_integrand_nunu) :: f
		real(dl) :: integrate_collint_nunu_GL
		type(coll_args), intent(in) :: obj
		integer :: ia, ib
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(fy2_arr(Ny, Ny))
		fy2_arr = 0.d0
		do ia=1, Ny
			do ib=1, Ny
				fy2_arr(ia, ib) = f(ia, ib, obj, F_nu_sc, F_nu_pa)
			end do
		end do
		integrate_collint_nunu_GL = integral_GL_2d(Ny, w_gl_arr2, w_gl_arr2, fy2_arr)
		deallocate(fy2_arr)
	end function integrate_collint_nunu_GL

	pure function get_collision_terms(collArgsIn, Fint_nue, Fint_nunu)
		use fpInterfaces2
		use fpInterfaces3
		implicit None
		procedure (collision_integrand_nue) :: Fint_nue
		procedure (collision_integrand_nunu) :: Fint_nunu
		type(cmplxMatNN) :: get_collision_terms
		procedure (collision_integrator_nue), pointer :: integrator_nue
		procedure (collision_integrator_nunu), pointer :: integrator_nunu
		real(dl) :: x, w, z, y1, cf, dampfact
		integer :: iy1
		type(coll_args), intent(in) :: collArgsIn
		type(coll_args) :: collArgs
		integer :: i, j

		call allocateCmplxMat(get_collision_terms)

		collArgs=collArgsIn

		x = collArgs%x
		w = collArgs%w
		z = collArgs%z
		y1 = collArgs%y1
		iy1 = collArgs%iy

		get_collision_terms%y = y1
		get_collision_terms%x = x
		get_collision_terms%z = z

		get_collision_terms%re = 0.d0
		get_collision_terms%im = 0.d0

		if (.not. has_offdiagonal() .and. collint_diagonal_zero) then
			return
		end if

		if (use_gauss_laguerre) then
			integrator_nue => integrate_collint_nue_GL
			integrator_nunu => integrate_collint_nunu_GL
		else
			integrator_nue => integrate_collint_nue_NC
			integrator_nunu => integrate_collint_nunu_NC
		end if
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			!diagonal elements:
			if (.not.sterile(i)) then
				if (.not. collint_diagonal_zero) then
					if (.not.collint_d_no_nue) &
						get_collision_terms%re(i,i) = get_collision_terms%re(i,i) &
							+ integrator_nue(Fint_nue, collArgs, F_ab_ann_re, F_ab_sc_re)
					if (.not.collint_d_no_nunu) &
						get_collision_terms%re(i,i) = get_collision_terms%re(i,i) &
							+ integrator_nunu(Fint_nunu, collArgs, F_nu_sc_re, F_nu_pa_re)/4.d0
				end if
			end if
			!coefficient for damping terms
			if (collint_damping_type.eq.2) then !dampings from McKellar:1992ja
				dampfact = z*z*z*z * y1*y1*y1 * dampTermFactor
			else if (collint_damping_type.eq.1) then !dampings from Bennett:2020zkv
				dampfact = z*z*z*z * y1*y1*y1 * 2.d0 * dy_damping_fit(y1/z)
			end if
			!off-diagonal elements:
			if (.not. collint_offdiag_damping) then !full integration for nue, for nunu integration only if FULL_F_NU is defined, damping if not (disabled for tests)
				do j=i+1, flavorNumber
					collArgs%ix2 = j
					if (.not.collint_od_no_nue) then
						get_collision_terms%re(i,j) = get_collision_terms%re(i,j) &
							+ integrator_nue(Fint_nue, collArgs, F_ab_ann_re, F_ab_sc_re)
						get_collision_terms%im(i,j) = get_collision_terms%im(i,j) &
							+ integrator_nue(Fint_nue, collArgs, F_ab_ann_im, F_ab_sc_im)
					end if
#ifndef DO_TESTS
#ifdef FULL_F_NU
					if (.not.collint_od_no_nunu) then
						get_collision_terms%re(i,j) = get_collision_terms%re(i,j) &
							+ integrator_nunu(Fint_nunu, collArgs, F_nu_sc_re, F_nu_pa_re)/4.d0
						get_collision_terms%im(i,j) = get_collision_terms%im(i,j) &
							+ integrator_nunu(Fint_nunu, collArgs, F_nu_sc_im, F_nu_pa_im)/4.d0
					end if
#endif
#endif
				end do
#ifndef DO_TESTS
#ifndef FULL_F_NU
				do j=i+1, flavorNumber
					get_collision_terms%re(i,j) = &
						get_collision_terms%re(i,j) &
						- dampTermMatrixCoeffNunu(i,j) * dampfact * nuDensMatVecFD(iy1)%re(i,j)
					get_collision_terms%im(i,j) = &
						get_collision_terms%im(i,j) &
						- dampTermMatrixCoeffNunu(i,j) * dampfact * nuDensMatVecFD(iy1)%im(i,j)
				end do
#endif
#endif
			else if (collint_damping_type.eq.2 .or. collint_damping_type.eq.1) then
				do j=i+1, flavorNumber
					get_collision_terms%re(i,j) = &
						- (dampTermMatrixCoeffNue(i,j)+dampTermMatrixCoeffNunu(i,j)) &
						* dampfact * nuDensMatVecFD(iy1)%re(i,j)
					get_collision_terms%im(i,j) = &
						- (dampTermMatrixCoeffNue(i,j)+dampTermMatrixCoeffNunu(i,j)) &
						* dampfact * nuDensMatVecFD(iy1)%im(i,j)
				end do
			end if
			!fill other half of the matrix:
			do j=i+1, flavorNumber
				get_collision_terms%re(j,i) = get_collision_terms%re(i,j)
				get_collision_terms%im(j,i) = - get_collision_terms%im(i,j)
			end do
		end do
		!fix coefficients:
		cf = (y1*y1*x**4)
		get_collision_terms%re(:,:) = get_collision_terms%re(:,:) * collTermFactor / cf
		get_collision_terms%im(:,:) = get_collision_terms%im(:,:) * collTermFactor / cf
	end function get_collision_terms
end module
