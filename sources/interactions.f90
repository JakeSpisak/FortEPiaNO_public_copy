module ndInteractions
	use precision
	use constants
	use variables
	use utilities
	use ndErrors
	use ndMatrices
	use linear_interpolation_module
	implicit none

	type(linear_interp_2d) :: dmgCorr, dmeNLCorr
	type(linear_interp_3d) :: dmeCorr

	contains

	subroutine init_interp_dme2_e
		real(dl), dimension(:,:), allocatable :: dmg_vec
		real(dl), dimension(:,:,:), allocatable :: dme_vec
		integer :: ix, iy, iz, iflag
		real(dl) :: x, z, t1, t2
		real(8) :: timer1
		logical :: initial

		call addToLog("[interactions] Initializing interpolation for photon mass corrections...")
		allocate(dmg_vec(interp_nx, interp_nz))
		!$omp parallel do default(shared) private(ix, iz) schedule(dynamic)
		do ix=1, interp_nx
			do iz=1, interp_nz
				dmg_vec(ix,iz) = dmg2_full(interp_xvec(ix),interp_zvec(iz))
			end do
		end do
		!$omp end parallel do
		call dmgCorr%initialize(interp_xvec, interp_zvec, dmg_vec, iflag)!linear

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
		initial = dme2_log_term
		dme2_log_term = .false.
		!$omp parallel do default(shared) private(ix, iz) schedule(dynamic)
		do ix=1, interp_nx
			do iz=1, interp_nz
				dmg_vec(ix,iz) = dme2_electronFull(interp_xvec(ix), 0.d0, interp_zvec(iz))
			end do
		end do
		!$omp end parallel do
		call dmeNLCorr%initialize(interp_xvec, interp_zvec, dmg_vec, iflag)!linear
		dme2_log_term = initial

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
				t1 = dme2_nolog(x,z)
			end do
			call toc(timer1, "<full>")
		end if
		call random_number(x)
		call random_number(z)
		x=(x_fin-x_in)*x + x_in
		z=0.4d0*z + z_in
		write(*,"(' [interactions] test dmg2_interp in ',*(E12.5))") x,z
		t1 = dmg2_full(x,z)
		t2 = dmg2_interp(x,z)
		write(*,"(' [interactions] comparison (true vs interp): ',*(E17.10))") t1,t2

		write(*,"(' [interactions] test dme2_electronInterp in ',*(E12.5))") x,z
		t1 = dme2_electronFull(x,0.d0,z)
		t2 = dme2_electron(x,0.d0,z)
		write(*,"(' [interactions] comparison (true vs interp): ',*(E17.10))") t1,t2

		deallocate(dmg_vec)
		deallocate(dme_vec)
		call addToLog("[interactions] ...done!")
	end subroutine init_interp_dme2_e

	elemental function E_k_m(k,m)
		real(dl), intent(in) :: k, m
		real(dl) :: E_k_m
		E_k_m = sqrt(k*k+m*m)
	end function E_k_m

	elemental function boseEinstein(x)
		real(dl) :: boseEinstein
		real(dl), intent(in) :: x
		boseEinstein = 1.d0/(exp(x) - 1.d0)
	end function boseEinstein

	elemental function fermiDirac(x)
		real(dl) :: fermiDirac
		real(dl), intent(in) :: x
		fermiDirac = 1.d0/(exp(x) + 1.d0)
	end function fermiDirac

	pure function dme2_e_i1(x, z, k)
	!doi:10.1016/S0370-2693(02)01622-2 eq.12 first integral
	!equal to integral in eq.13
		real(dl) :: dme2_e_i1
		real(dl), intent(in) :: x, z, k
		real(dl) :: Ekm
		Ekm = E_k_m(k, x)
		dme2_e_i1 = 1.d0/Ekm * fermiDirac(Ekm/z)
	end function dme2_e_i1

	pure function dme2_e_i2(x, y, z, k)
		!doi:10.1016/S0370-2693(02)01622-2 eq.12 second integral
		real(dl) :: dme2_e_i2
		real(dl), intent(in) :: x, y, z, k
		real(dl) :: m, T, p, Ekm
		Ekm = E_k_m(k, x)
		dme2_e_i2 = k/Ekm * fermiDirac(Ekm/z)*log(abs((y+k)/(y-k)))
	end function dme2_e_i2

	pure function dmg2_full(x, z)
	!doi:10.1016/S0370-2693(02)01622-2 eq.13
		real(dl) :: dmg2_full, integr_1
		real(dl), intent(in) :: x, z
		integer :: i

		if (dme2_temperature_corr) then
			integr_1 = 0.d0
			do i=1, N_opt_y
				integr_1 = integr_1 + opt_y_w(i)*dme2_e_i1(x, z, opt_y(i))
			end do
			dmg2_full = 8. * alpha_fine * integr_1 / PI
		else
			dmg2_full = 0.d0
		end if
	end function dmg2_full

	pure function dme2_electronFull(x, y, z)
	!doi:10.1016/S0370-2693(02)01622-2 eq.12
		real(dl) :: dme2_electronFull, integr_1, integr_2
		real(dl), intent(in) :: x, y, z
		integer :: i

		if (dme2_temperature_corr) then
			integr_1 = 0.d0
			do i=1, N_opt_y
				integr_1 = integr_1 + opt_y_w(i)*dme2_e_i1(x, z, opt_y(i))
			end do
			dme2_electronFull = 2. * alpha_fine * (z*z * PID3 + integr_1/PID2)
			if (dme2_log_term) then
				integr_2 = 0.d0
				do i=1, N_opt_y
					if (y .ne. opt_y(i)) &
					integr_2 = integr_2 + opt_y_w(i)*dme2_e_i2(x, y, z, opt_y(i))/(opt_y(i)**2)
				end do
				dme2_electronFull = dme2_electronFull &
					- x*x* alpha_fine / y / PID2 * integr_2
			end if
		else
			dme2_electronFull = 0.d0
		end if
	end function dme2_electronFull

	function dmg2_interp(x, z)
		real(dl) :: dmg2_interp
		real(dl), intent(in) :: x, z
		integer :: iflag

		call dmgCorr%evaluate(x, z, dmg2_interp)
	end function dmg2_interp

	function dme2_electron(x, y, z)
		real(dl) :: dme2_electron
		real(dl), intent(in) :: x, y, z
		integer :: iflag

		call dmeCorr%evaluate(x, y, z, dme2_electron)
	end function dme2_electron

	function dme2_nolog(x, z)
		real(dl) :: dme2_nolog
		real(dl), intent(in) :: x, z
		integer :: iflag

		call dmeNLCorr%evaluate(x, z, dme2_nolog)
	end function dme2_nolog

	elemental function Ebare_i_dme(x, y, dme2)!for electrons
		real(dl) :: Ebare_i_dme
		real(dl), intent(in) :: x, y, dme2

		Ebare_i_dme = sqrt(x*x+y*y+dme2)
	end function Ebare_i_dme

	!phase space
	pure function F_ab_ann_re(n1, n2, f3, f4, a, b, i, j)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.5
		real(dl) :: F_ab_ann_re
		type(cmplxMatNN), intent(in) :: n1, n2
		real(dl), intent(in) :: f3, f4
		integer, intent(in) :: a, b, i, j
		integer :: k
		real(dl) :: t1a, t1b, t2a, t2b
#ifdef FULLFAB
		integer :: l, m
#endif

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_ann_re] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
#ifdef FULLFAB
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
#ifdef FULLFAB
		integer :: l, m
#endif

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_ann_im] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
#ifdef FULLFAB
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
#ifdef FULLFAB
		integer :: l, m
#endif

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_sc_re] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
#ifdef FULLFAB
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
#ifdef FULLFAB
		integer :: l, m
#endif

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_sc_im] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
#ifdef FULLFAB
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
	pure function coll_nue_3_ann_int(iy2, y4, obj, F_ab)
		!annihilation
		use ndInterfaces1
		procedure (F_annihilation) :: F_ab
		integer, intent(in) :: iy2
		real(dl), intent(in) :: y4
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_3_ann_int
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y2, y3, f3, f4, E3, E4, dme2, t1, t2

		coll_nue_3_ann_int = 0.d0

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
			coll_nue_3_ann_int = coll_nue_3_ann_int + &
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
	end function coll_nue_3_ann_int

	pure function coll_nue_3_sc_int(iy3, y2, obj, F_ab)
		!scattering, summing positron and electrons
		use ndInterfaces1
		procedure (F_scattering) :: F_ab
		integer, intent(in) :: iy3
		real(dl), intent(in) :: y2
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_3_sc_int
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y3, y4, f2, f4, E2, E4, dme2, t1, t2

		coll_nue_3_sc_int = 0.d0

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
			coll_nue_3_sc_int = coll_nue_3_sc_int + &
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
	end function coll_nue_3_sc_int

	pure function coll_nue_3_int(iy, yx, obj, F_ab_ann, F_ab_sc)
		use ndInterfaces1
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		integer, intent(in) :: iy
		real(dl), intent(in) :: yx
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_3_int
		coll_nue_3_int = &
			coll_nue_3_sc_int(iy, yx, obj, F_ab_sc) &
			+ coll_nue_3_ann_int(iy, yx, obj, F_ab_ann)
	end function coll_nue_3_int

	!integrate collision terms
	pure function integrate_coll_int_NC(f, obj, F_ab_ann, F_ab_sc)
		use ndInterfaces1
		use ndInterfaces2
		implicit None
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		procedure (collision_integrand) :: f
		real(dl) :: integrate_coll_int_NC
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
		integrate_coll_int_NC = integral_NC_2d(Ny, Ny, dy_arr, dy_arr, fy2_arr)
		deallocate(fy2_arr)
	end function integrate_coll_int_NC

	pure function integrate_coll_int_GL(f, obj, F_ab_ann, F_ab_sc)
		use ndInterfaces1
		use ndInterfaces2
		implicit None
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		procedure (collision_integrand) :: f
		real(dl) :: integrate_coll_int_GL
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
		integrate_coll_int_GL = integral_GL_2d(Ny, w_gl_arr2, w_gl_arr2, fy2_arr)
		deallocate(fy2_arr)
	end function integrate_coll_int_GL

	pure function get_collision_terms(collArgsIn, Fint)
		use ndInterfaces2
		use ndInterfaces3
		implicit None
		procedure (collision_integrand) :: Fint
		type(cmplxMatNN) :: get_collision_terms
		procedure (collision_integrator), pointer :: integrator
		real(dl) :: x, z, y1, cf
		integer :: iy1
		type(coll_args), intent(in) :: collArgsIn
		type(coll_args) :: collArgs
		integer :: i, j

		call allocateCmplxMat(get_collision_terms)

		collArgs=collArgsIn

		x = collArgs%x
		z = collArgs%z
		y1 = collArgs%y1
		iy1 = collArgs%iy

		get_collision_terms%y = y1
		get_collision_terms%x = x
		get_collision_terms%z = z

		get_collision_terms%re = 0.d0
		get_collision_terms%im = 0.d0

		if (collision_offdiag.eq.0) then
			return
		end if

		if (use_gauss_laguerre) then
			integrator => integrate_coll_int_GL
		else
			integrator => integrate_coll_int_NC
		end if
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			if (.not.sterile(i)) then
				get_collision_terms%re(i,i) = integrator(Fint, collArgs, F_ab_ann_re, F_ab_sc_re)
			end if
			if (collision_offdiag.eq.1) then
				do j=i+1, flavorNumber
					collArgs%ix2 = j
					get_collision_terms%re(i,j) = integrator(Fint, collArgs, F_ab_ann_re, F_ab_sc_re)
					get_collision_terms%im(i,j) = integrator(Fint, collArgs, F_ab_ann_im, F_ab_sc_im)
				end do
			else if (collision_offdiag.eq.2) then
				do j=i+1, flavorNumber
					collArgs%ix2 = j
					get_collision_terms%re(i,j) = dampTermFactor * dampTermMatrixCoeff(i,j) * y1*y1*y1 * nuDensMatVecFD(iy1)%re(i,j)
					get_collision_terms%im(i,j) = dampTermFactor * dampTermMatrixCoeff(i,j) * y1*y1*y1 * nuDensMatVecFD(iy1)%im(i,j)
				end do
			end if
			do j=i+1, flavorNumber
				get_collision_terms%re(j,i) = get_collision_terms%re(i,j)
				get_collision_terms%im(j,i) = - get_collision_terms%im(i,j)
			end do
		end do
		cf = (y1*y1*x**4)
		get_collision_terms%re(:,:) = get_collision_terms%re(:,:) * collTermFactor / cf
		get_collision_terms%im(:,:) = get_collision_terms%im(:,:) * collTermFactor / cf
	end function get_collision_terms
end module
