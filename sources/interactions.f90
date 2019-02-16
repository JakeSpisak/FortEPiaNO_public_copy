module ndInteractions
	use precision
	use constants
	use variables
	use utilities
	use ndErrors
	use ndMatrices
	use linear_interpolation_module
	implicit none

	type interpNuDens_obj
		type(spline_class), dimension(:,:), allocatable :: re, im
	end type interpNuDens_obj
	type(interpNuDens_obj) :: interpNuDens

	type(spline_class) :: FD_interp
	type(linear_interp_2d) :: dmeCorr
	
	type(linear_interp_4d) :: D2_interp, D3_interp, PI1_12_interp, PI1_13_interp

	contains

	subroutine allocate_interpNuDens
		real(dl), dimension(:), allocatable :: ndmv_re, ndmv_im
		integer :: i, j, iy
		allocate(interpNuDens%re(flavorNumber, flavorNumber), &
			interpNuDens%im(flavorNumber, flavorNumber))
		allocate(ndmv_re(Ny), ndmv_im(Ny))
		do i=1, flavorNumber
			do iy=1, Ny
				ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, i)
			end do
			call interpNuDens%re(i, i)%initialize(Ny, y_arr, ndmv_re)
			do j=i+1, flavorNumber
				do iy=1, Ny
					ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, j)
					ndmv_im(iy) = nuDensMatVecFD(iy)%im(i, j)
				end do
				call interpNuDens%re(i, j)%initialize(Ny, y_arr, ndmv_re)
				call interpNuDens%im(i, j)%initialize(Ny, y_arr, ndmv_im)
			end do
		end do
	end subroutine allocate_interpNuDens

	subroutine init_interpNuDens
		real(dl), dimension(:), allocatable :: ndmv_re, ndmv_im
		integer :: i, j, iy
		allocate(ndmv_re(Ny), ndmv_im(Ny))
		do i=1, flavorNumber
			do iy=1, Ny
				ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, i)
			end do
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
			do j=i+1, flavorNumber
				do iy=1, Ny
					ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, j)
					ndmv_im(iy) = nuDensMatVecFD(iy)%im(i, j)
				end do
				call interpNuDens%re(i, j)%replace(Ny, y_arr, ndmv_re)
				call interpNuDens%im(i, j)%replace(Ny, y_arr, ndmv_im)
			end do
		end do
	end subroutine init_interpNuDens

	subroutine init_interp_dme2_e
		real(dl), dimension(:,:), allocatable :: dme_vec
		integer :: ix, iz, iflag
		real(dl) :: x,z,t1,t2
		real(8) :: timer1
		
		call addToLog("[interactions] Initializing interpolation for electron mass corrections...")
		allocate(dme_vec(interp_nx,interp_nz))
		do ix=1, interp_nx
			do iz=1, interp_nz
				dme_vec(ix,iz) = dme2_electronFull(interp_xvec(ix),0.d0,interp_zvec(iz))
			end do
		end do
		call dmeCorr%initialize(interp_xvec,interp_zvec,dme_vec,iflag)!linear
		
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
				t1 = dme2_electron(x,0.d0,z)
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
		call random_number(x)
		call random_number(z)
		x=(x_fin-x_in)*x + x_in
		z=0.4d0*z + z_in
		write(*,"(' [interactions] test dme2_electronInterp in ',*(E12.5))") x,z
		t1 = dme2_electronFull(x,0.d0,z)
		t2 = dme2_electron(x,0.d0,z)
		write(*,"(' [interactions] comparison (true vs interp): ',*(E17.10))") t1,t2

		deallocate(dme_vec)
		call addToLog("[interactions] ...done!")
	end subroutine init_interp_dme2_e
	
	elemental function E_k_m(k,m)
		real(dl), intent(in) :: k,m
		real(dl) :: E_k_m
		E_k_m = sqrt(k*k+m*m)
	end function E_k_m
	
	elemental function fermiDirac_full(x)
		real(dl) :: fermiDirac_full
		real(dl), intent(in) :: x

		fermiDirac_full = 1.d0/(exp(x) + 1.d0)
	end function fermiDirac_full

	elemental function fermiDirac(x)
		real(dl) :: fermiDirac
		real(dl), intent(in) :: x

		fermiDirac = fermiDirac_full(x)
!		fermiDirac = fermiDirac_i(x)
	end function fermiDirac

	pure function fermiDirac_i(x)
		real(dl) :: fermiDirac_i
		real(dl), intent(in) :: x

		fermiDirac_i = FD_interp%evaluate(x)
	end function fermiDirac_i

	subroutine init_interp_FD
		real(dl), dimension(:), allocatable :: fd_vec, fd_x
		integer :: ix, iflag, interp_nfd
		real(dl) :: x, t1,t2
		real(8) :: timer1

		call addToLog("[interactions] Initializing interpolation for FermiDirac...")
		interp_nfd=100
		allocate(fd_x(interp_nfd), fd_vec(interp_nfd))
		fd_x=logspace(-3.d0, 2.d0, interp_nfd)
		do ix=1, interp_nfd
			fd_vec(ix) = fermiDirac_full(fd_x(ix))
		end do
		call FD_interp%initialize(interp_nfd, fd_x, fd_vec)

		call random_seed()
		call random_number(x)
		x=10.d0**(x*(2.d0-(-3.d0)) - 3.d0)
		write(*,"(' [interactions] test fermiDirac in ',*(E12.5))") x
		t1 = fermiDirac_full(x)
		t2 = fermiDirac_i(x)
		write(*,"(' [interactions] comparison (true vs interp): ',*(E17.10))") t1,t2

		if (timing_tests) then
			write (*,*) "[interactions] now doing some timing..."
			call tic(timer1)
			do ix=1, 100000000
				call random_number(x)
				x=10.d0**(x*(2.d0-(-3.d0)) - 3.d0)
				t1 = fermiDirac_i(x)
			end do
			call toc(timer1, "<reset>")

			call tic(timer1)
			do ix=1, 100000000
				call random_number(x)
				x=10.d0**(x*(2.d0-(-3.d0)) - 3.d0)
				t1 = fermiDirac_i(x)
			end do
			call toc(timer1, "<interpolated>")

			call tic(timer1)
			do ix=1, 100000000
				call random_number(x)
				x=10.d0**(x*(2.d0-(-3.d0)) - 3.d0)
				t1 = fermiDirac_full(x)
			end do
			call toc(timer1, "<full>")
		end if

		call addToLog("[interactions] ...done!")
	end subroutine init_interp_FD

	pure function dme2_e_i1(vec, k)
	!doi:10.1016/S0370-2693(02)01622-2 eq.12 first integral
		real(dl) :: dme2_e_i1
		real(dl), intent(in) :: k
		real(dl), dimension(3), intent(in) :: vec
		real(dl) :: m,T,p, Ekm
		m=vec(1)
		T=vec(2)
		Ekm = E_k_m(k,m)
		dme2_e_i1 = k*k/Ekm * fermiDirac(Ekm/T)
	end function dme2_e_i1
	
	pure function dme2_e_i2(vec, k) !not used
	!doi:10.1016/S0370-2693(02)01622-2 eq.12 second integral
		real(dl) :: dme2_e_i2
		real(dl), intent(in) :: k
		real(dl), dimension(3), intent(in) :: vec
		real(dl) :: m,T,p, Ekm
		m=vec(1)
		T=vec(2)
		p=vec(3)
		Ekm = E_k_m(k,m)
		dme2_e_i2 = k/Ekm * fermiDirac(Ekm/T)*log(abs((p+k)/(p-k)))
	end function dme2_e_i2
	
	function dme2_electronFull(x,y,z)
	!doi:10.1016/S0370-2693(02)01622-2 eq.12
		real(dl) :: dme2_electronFull, tmp
		real(dl), intent(in) :: x,y,z
		real(dl), dimension(3) :: vec
		
		if (dme2_temperature_corr) then
			vec(1) = x
			vec(2) = z
			vec(3) = y !not used
			
			tmp = rombint_vec(vec, dme2_e_i1, fe_l, fe_u, toler_dme2, maxiter)
			dme2_electronFull = 2. * alpha_fine * z*z * (PID3 + tmp/PID2)
		else
			dme2_electronFull = 0.d0
		end if
	end function dme2_electronFull
	
	function dme2_electron(x,y,z)
		real(dl) :: dme2_electron
		real(dl), intent(in) :: x,y,z
		integer :: iflag
		
		call dmeCorr%evaluate(x,z,dme2_electron)!linear
	end function dme2_electron
	
	elemental function Ebare_i_dme(x,y,dme2)!for electrons
		real(dl) :: Ebare_i_dme
		real(dl), intent(in) :: x,y,dme2
		
		Ebare_i_dme = sqrt(x*x+y*y+dme2)
	end function Ebare_i_dme
	
	pure function F_ab_ann_re(n1,n2,f3,f4, a,b, i, j)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.5
		real(dl) :: F_ab_ann_re
		type(cmplxMatNN), intent(in) :: n1,n2
		real(dl), intent(in) :: f3, f4
		integer, intent(in) :: a,b, i, j
		integer :: k
		real(dl) :: t1a, t1b, t2a, t2b
		
!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_ann_re] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
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
		F_ab_ann_re = f3 * f4 * &
				(t1a * GLR_vec(a,i,i) + t1b * GLR_vec(a,j,j))&
			- (1.d0 - f3) * (1.d0 - f4) * &
				(t2a * GLR_vec(a,j,j) + t2b * GLR_vec(a,i,i))
	end function F_ab_ann_re
	
	pure function F_ab_ann_im(n1,n2,f3,f4, a,b, i, j)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.5
		real(dl) :: F_ab_ann_im
		type(cmplxMatNN), intent(in) :: n1,n2
		real(dl), intent(in) :: f3,f4
		integer, intent(in) :: a,b, i, j
		integer :: k
		real(dl) :: t1a, t1b, t2a, t2b

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_ann_im] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
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
		F_ab_ann_im = f3 * f4 * &
				(t1a * GLR_vec(a,i,i) + t1b * GLR_vec(a,j,j))&
			- (1.d0 - f3) * (1.d0 - f4) * &
				(t2a * GLR_vec(a,j,j) + t2b * GLR_vec(a,i,i))
	end function F_ab_ann_im

	pure function F_ab_sc_re (n1,f2,n3,f4, a,b, i, j)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.10
		real(dl) :: F_ab_sc_re
		type(cmplxMatNN), intent(in) :: n1,n3
		real(dl), intent(in) :: f2,f4
		integer, intent(in) :: a,b, i, j
		integer :: k
		real(dl) :: t1a, t1b, t2a, t2b

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_sc_re] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
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
		F_ab_sc_re = (1.d0 - f2) * f4 * &
				((t1a) * GLR_vec(a,i,i) + (t1b) * GLR_vec(a,j,j))&
			- f2 * (1.d0 - f4) * &
				((t2a) * GLR_vec(a,j,j) + (t2b) * GLR_vec(a,i,i))
	end function F_ab_sc_re

	pure function F_ab_sc_im (n1,f2,n3,f4, a,b, i, j)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.10
		real(dl) :: F_ab_sc_im
		type(cmplxMatNN), intent(in) :: n1,n3
		real(dl), intent(in) :: f2,f4
		integer, intent(in) :: a,b, i, j
		integer :: k
		real(dl) :: t1a, t1b, t2a, t2b

!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_sc_im] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
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
		F_ab_sc_im = (1.d0 - f2) * f4 * &
				((t1a) * GLR_vec(a,i,i) + (t1b) * GLR_vec(a,j,j))&
			- f2 * (1.d0 - f4) * &
				((t2a) * GLR_vec(a,j,j) + (t2b) * GLR_vec(a,i,i))
	end function F_ab_sc_im

	subroutine init_interp_d123
		real(dl), dimension(:,:,:,:), allocatable :: d2, d3, pi1_12, pi1_13
		real(dl) :: y1,  y2,  y3,  y4, d1a, d1b, d2a, d2b, d3a, d3b
		integer  :: ix1, ix2, ix3, ix4, iflag
		real(dl), dimension(2) :: PI1
		real(dl), dimension(:), allocatable :: yarr
		real(8) :: timer1
		
		call addToLog("[interactions] Initializing interpolation for D1/2/3 functions...")
		allocate(yarr(interp_ny))
		yarr = logspace(log10(y_min), log10(y_max), interp_ny)
		allocate(d2(interp_ny,interp_ny,interp_ny,interp_ny), &
			d3(interp_ny,interp_ny,interp_ny,interp_ny), &
			pi1_12(interp_ny,interp_ny,interp_ny,interp_ny), &
			pi1_13(interp_ny,interp_ny,interp_ny,interp_ny))
		do ix1=1, interp_ny
			do ix2=1, interp_ny
				do ix3=1, interp_ny
					do ix4=1, interp_ny
						d2(ix1,ix2,ix3,ix4) = d2_full(yarr(ix1), yarr(ix2), yarr(ix3), yarr(ix4))
						d3(ix1,ix2,ix3,ix4) = d3_full(yarr(ix1), yarr(ix2), yarr(ix3), yarr(ix4))

						pi1_12(ix1,ix2,ix3,ix4) = PI1_12_full(yarr(ix1), yarr(ix2), yarr(ix3), yarr(ix4))
						pi1_13(ix1,ix2,ix3,ix4) = PI1_13_full(yarr(ix1), yarr(ix2), yarr(ix3), yarr(ix4))
					end do
				end do
			end do
		end do
		call PI1_12_interp%initialize(yarr, yarr, yarr, yarr, pi1_12, iflag)
		deallocate(pi1_12)
		call PI1_13_interp%initialize(yarr, yarr, yarr, yarr, pi1_13, iflag)
		deallocate(pi1_13)

		call D2_interp%initialize(yarr, yarr, yarr, yarr, d2, iflag)
		deallocate(d2)
		call D3_interp%initialize(yarr, yarr, yarr, yarr, d3, iflag)
		deallocate(d3)

		call random_seed()
		call random_number(y1)
		call random_number(y2)
		call random_number(y3)
		call random_number(y4)
		y1=(y_max-y_min)*y1 + y_min
		y2=(y_max-y_min)*y2 + y_min
		y3=(y_max-y_min)*y3 + y_min
		y4=(y_max-y_min)*y4 + y_min

		write(*,"(' [interactions] test PI_1, D2_f, D3_f in ',*(E12.5))") y1,y2,y3,y4
		d2a = D2_f(y1,y2,y3,y4)
		d3a = D3_f(y1,y2,y3,y4)
		d2b = D2_full(y1,y2,y3,y4)
		d3b = D3_full(y1,y2,y3,y4)
		pi1 = PI1_full(y1,y2,y3,y4)
		write(*,"(' [interactions] comparison pi1_12_i (true vs interp): ',*(E17.10))") PI1_12_full(y1,y2,y3,y4), pi1_12_i(y1,y2,y3,y4)
		write(*,"(' [interactions] comparison pi1_13_i (true vs interp): ',*(E17.10))") PI1_13_full(y1,y2,y3,y4), pi1_13_i(y1,y2,y3,y4)
		write(*,"(' [interactions] comparison D2_f (true vs interp): ',*(E17.10))") d2a,d2b
		write(*,"(' [interactions] comparison D3_f (true vs interp): ',*(E17.10))") d3a,d3b

		if (timing_tests) then
			call tic(timer1)
			write (*,*) "[interactions] now doing some timing..."
			call tic(timer1)
			do ix1=1, 10000000
				call random_number(y1)
				call random_number(y2)
				call random_number(y3)
				call random_number(y4)
				y1=(y_max-y_min)*y1 + y_min
				y2=(y_max-y_min)*y2 + y_min
				y3=(y_max-y_min)*y3 + y_min
				y4=(y_max-y_min)*y4 + y_min
!				d1a = d1_cases(y1,y2,y3,y4)
				d1a = d2_cases(y1,y2,y3,y4)
				d1b = d3_cases(y1,y2,y3,y4)
!				d2a = PI1_12_i(y1,y2,y3,y4)
!				d3a = PI1_13_i(y1,y2,y3,y4)
			end do
			call toc(timer1, "<reset>")

			call tic(timer1)
			do ix1=1, 100000000
				call random_number(y1)
				call random_number(y2)
				call random_number(y3)
				call random_number(y4)
				y1=(y_max-y_min)*y1 + y_min
				y2=(y_max-y_min)*y2 + y_min
				y3=(y_max-y_min)*y3 + y_min
				y4=(y_max-y_min)*y4 + y_min
!				d1a = d1_cases(y1,y2,y3,y4)
				d1a = d2_cases(y1,y2,y3,y4)
				d1b = d3_cases(y1,y2,y3,y4)
!				d2a = PI1_12_i(y1,y2,y3,y4)
!				d3a = PI1_13_i(y1,y2,y3,y4)
			end do
			call toc(timer1, "<cases>")

			call tic(timer1)
			do ix1=1, 100000000
				call random_number(y1)
				call random_number(y2)
				call random_number(y3)
				call random_number(y4)
				y1=(y_max-y_min)*y1 + y_min
				y2=(y_max-y_min)*y2 + y_min
				y3=(y_max-y_min)*y3 + y_min
				y4=(y_max-y_min)*y4 + y_min
				d1a = d1_full(y1,y2,y3,y4)
				d1a = d2_full(y1,y2,y3,y4)
				d1b = d3_full(y1,y2,y3,y4)
!				d2a = PI1_12_full(y1,y2,y3,y4)
!				d3a = PI1_13_full(y1,y2,y3,y4)
			end do
			call toc(timer1, "<full>")

			call tic(timer1)
			do ix1=1, 100000000
				call random_number(y1)
				call random_number(y2)
				call random_number(y3)
				call random_number(y4)
				y1=(y_max-y_min)*y1 + y_min
				y2=(y_max-y_min)*y2 + y_min
				y3=(y_max-y_min)*y3 + y_min
				y4=(y_max-y_min)*y4 + y_min
				d1a = d1_pablo(y1,y2,y3,y4)
				d1a = d2_pablo(y1,y2,y3,y4)
				d1b = d3_pablo(y1,y2,y3,y4)
!				d2a = PI1_12_full(y1,y2,y3,y4)
!				d3a = PI1_13_full(y1,y2,y3,y4)
			end do
			call toc(timer1, "<pablo>")
		end if

		call addToLog("[interactions] ...done!")
	end subroutine init_interp_d123

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

	elemental function D1_pablo(y1, y2, y3, y4)
		implicit none
		real(dl) :: D1_pablo
		real(dl), intent(in) :: y1, y2, y3, y4
          real(dl) a12,a34
          a12=y1+y2
          a34=y3+y4
          D1_pablo=a12+a34-dabs(a12-a34)-dabs(y1-y2+y3-y4)-dabs(y1-y2-y3+y4)
	end function D1_pablo

	function D2_f(y1, y2, y3, y4)
		implicit none
		real(dl) :: D2_f
		real(dl), intent(in) :: y1, y2, y3, y4
		integer :: iflag

		D2_f=0.d0
		call D2_interp%evaluate(y1,y2,y3,y4, D2_f)
	end function D2_f
	
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

	elemental function D2_cases(q1, q2, q3, q4)
	!10.1103/PhysRevD.94.033009 eq.D2
		implicit none
		real(dl) :: D2_cases
		real(dl), intent(in) :: q1, q2, q3, q4
		real(dl) :: y1, y2, y3, y4
		real(dl) :: a

		if (y1.le.y2) then
			y1=q2
			y2=q1
		else
			y1=q1
			y2=q2
		end if
		if (y3.le.y4) then
			y3=q4
			y4=q3
		else
			y3=q3
			y4=q4
		end if
		if (y1+y2 .ge. y3+y4) then
			if (y1+y4 .ge. y2+y3) then
				if (y1 .le. y2+y3+y4) then
					a = (y1-y2)
					D2_cases = ( a**3 &
						+ 2.d0 * (y3**3 + y4**3) &
						- 3.d0 * a * (y3**2 + y4**2) &
						) / 12.d0
				else
					D2_cases = 0.d0
				end if
			else
				D2_cases = y4**3 / 3.d0
			end if
		else
			if (y1+y4 .le. y2+y3) then
				if (y3 .le. y2+y1+y4) then
					a = (y1+y2)
					D2_cases = ( - a**3 &
						- 2.d0 * (y3**3 - y4**3) &
						+ 3.d0 * a * (y3**2 + y4**2) &
						) / 12.d0
				else
					D2_cases = 0.d0
				end if
			else
				D2_cases = y2 / 6.d0 * ( &
					3.d0 * (y3**2 + y4**2 - y1**2) &
					- y2**2 )
			end if
		end if

	end function D2_cases

    elemental real(dl) function D2_pablo(y1,y2,y3,y4)
          implicit none
          real(dl),intent(in)::y1,y2,y3,y4
          real(dl) a12,a13,a14,a34,a123,a124,a134,a234,a1234,s12,s13,s14
          real(dl) a12_2,sum1,sum2

          a12_2=y1+y2
          a34=y3+y4
          
          a123=a12_2+y3-y4
          a124=a12_2-y3+y4
          a134=y1-y2+a34
          a234=-y1+y2+a34
          a1234=a12_2+a34
          a12=a12_2-a34
          a13=y1-y2+y3-y4
          a14=y1-y2-y3+y4
          s12=a12*a12*dsign(1.d0,a12)
          s13=a13*a13*dsign(1.d0,a13)
          s14=a14*a14*dsign(1.d0,a14)

          sum1=-a123*a123 + a1234*a1234 - a124*a124 + s12
          sum2=- a134*a134 + a234*a234 + s13 + s14

          D2_pablo=(a123**3 - a1234**3 + a124**3 + a134**3 + a234**3 + &
           3.d0*(sum1 + sum2)*y1 + &
           3.d0*(sum1 - sum2)*y2 - dabs(a12)**3 - dabs(a13)**3 - &
           dabs(a14)**3 + 6.d0*y1*y2*(a12_2 - 3.d0*a34 - dabs(a12) + &
           dabs(a13) + dabs(a14)))/6.d0

        end function D2_pablo

	function D3_f(y1, y2, y3, y4)
		implicit none
		real(dl) :: D3_f
		real(dl), intent(in) :: y1, y2, y3, y4
		integer :: iflag

		D3_f=0.d0
		call D3_interp%evaluate(y1,y2,y3,y4, D3_f)
	end function D3_f
	
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

	elemental function D3_cases(q1, q2, q3, q4)
	!10.1103/PhysRevD.94.033009
		implicit none
		real(dl) :: D3_cases
		real(dl), intent(in) :: q1, q2, q3, q4
		real(dl) :: y1, y2, y3, y4
		real(dl) :: a

		if (y1.le.y2) then
			y1=q2
			y2=q1
		else
			y1=q1
			y2=q2
		end if
		if (y3.le.y4) then
			y3=q4
			y4=q3
		else
			y3=q3
			y4=q4
		end if
		if (y1+y2 .ge. y3+y4) then
			if (y1+y4 .ge. y2+y3) then
				if (y1 .le. y2+y3+y4) then
					D3_cases = 1.d0 / 60.d0 * ( &
						y1**5 - y2**5 - y3**5 - y4**5 + 5.d0 * &
						( y1**2 * (y2**3 + y3**3 + y4**3) &
						- y1**3 * (y2**2 + y3**2 + y4**2) &
						+ y2**2 * y3**2 * (y2 + y3) &
						+ y2**2 * y4**2 * (y2 + y4) &
						+ y4**2 * y3**2 * (y4 + y3) ) )
				else
					D3_cases = 0.d0
				end if
			else
				D3_cases = y4**3 / 30.d0 * (5.d0 * (y1**2 + y2**2 + y3**2) -y4**2 )
			end if
		else
			if (y1+y4 .le. y2+y3) then
				if (y3 .le. y2+y1+y4) then
					D3_cases = 1.d0 / 60.d0 * ( &
						y3**5 - y4**5 - y1**5 - y2**5 + 5.d0 * &
						( y3**2 * (y4**3 + y1**3 + y2**3) &
						- y3**3 * (y4**2 + y1**2 + y2**2) &
						+ y4**2 * y1**2 * (y4 + y1) &
						+ y4**2 * y2**2 * (y4 + y2) &
						+ y2**2 * y1**2 * (y2 + y1) ) )
					D3_cases = 1.d0
				else
					D3_cases = 0.d0
				end if
			else
				D3_cases = y2**3 / 30.d0 * (5.d0 * (y1**2 + y4**2 + y3**2) -y2**2 )
			end if
		end if

	end function D3_cases

    elemental real(dl) function D3_pablo(y1,y2,y3,y4)
          implicit none
          real(dl),intent(in)::y1,y2,y3,y4
          real(dl) a12,a13,a14,a123,a124,a134,a234,a1234,P1234,s12,s13,s14
          real(dl) absa12,absa13,absa14,absa12e3,absa13e3,absa14e3
          real(dl) ampp,appm,apmp,a123e2,a124e2,a134e2,a234e2,a1234e2
          real(dl) a123e3,a124e3,a134e3,a234e3,a1234e3
          real(dl) a123e5,a124e5,a134e5,a234e5,a1234e5

          a12=y1+y2-y3-y4
          a13=y1-y2+y3-y4
          a14=y1-y2-y3+y4
          a123=y1+y2+y3-y4
          a124=y1+y2-y3+y4
          a134=y1-y2+y3+y4
          a234=-y1+y2+y3+y4
          a1234=y1+y2+y3+y4
          s12=a12**2*dsign(1.0d0,a12)
          s13=a13**2*dsign(1.0d0,a13)
          s14=a14**2*dsign(1.0d0,a14)
          P1234=-120.d0*y1*y2*y3*y4

          absa12=dabs(a12)
          absa13=dabs(a13)
          absa14=dabs(a14)
          absa12e3=absa12**3
          absa13e3=absa13**3
          absa14e3=absa14**3
          ampp=-absa12e3 + absa13e3 + absa14e3
          appm=absa12e3 + absa13e3 - absa14e3
          apmp=absa12e3 - absa13e3 + absa14e3
          a123e2=a123*a123
          a124e2=a124*a124
          a134e2=a134*a134
          a234e2=a234*a234
          a1234e2=a1234*a1234
          a123e3=a123e2*a123
          a124e3=a124e2*a124
          a134e3=a134e2*a134
          a234e3=a234e2*a234
          a1234e3=a1234e2*a1234

          a123e5=a123e3*a123e2
          a124e5=a124e3*a124e2
          a134e5=a134e3*a134e2
          a234e5=a234e3*a234e2
          a1234e5=a1234e3*a1234e2

          D3_pablo= (4.d0*(a1234e5 - a123e5 - a124e5 - a134e5 - a234e5) - &
           absa12**5 - absa13**5 - absa14**5 + &
           (a123 + a1234 + a124 + a134 + a234 + absa12 + absa13 + absa14)*P1234 + &
           5.d0*(a12**3*s12+a13**3*s13 +a14**3*s14) + &
           20.d0*((-a1234e3 + a123e3 + a124e3 - a134e3 - a234e3 + ampp)*y1*y2 + &
            (-a1234e3 + a123e3 - a124e3 + a134e3 - a234e3 + apmp)*y1*y3 + &
            (-a1234e3 + a123e3 - a124e3 - a134e3 + a234e3 + appm)*y2*y3 + &
            (-a1234e3 - a123e3 + a124e3 + a134e3 - a234e3 + appm)*y1*y4 + &
            (-a1234e3 - a123e3 + a124e3 - a134e3 + a234e3 + apmp)*y2*y4 + &
            (-a1234e3 - a123e3 - a124e3 + a134e3 + a234e3 + ampp)*y3*y4) + &
           60.d0*((a1234e2 + a123e2 - a124e2 + a134e2 + a234e2 - s12 + s13 - s14)*y1*y2*y4 + &
            (a1234e2 - a123e2 + a124e2 + a134e2 + a234e2 - s12 - s13 + s14)*y1*y2*y3 + &
            (a1234e2 + a123e2 + a124e2 - a134e2 + a234e2 + s12 - s13 - s14)*y1*y3*y4 + &
            (a1234e2 + a123e2 + a124e2 + a134e2 - a234e2 + s12 + s13 + s14)*y2*y3*y4))/120.d0

        end function D3_pablo

	function PI1_12_i(y1, y2, y3, y4)
		implicit none
		real(dl) :: PI1_12_i
		real(dl), intent(in) :: y1, y2, y3, y4
		integer :: iflag

		PI1_12_i=0.d0
		call PI1_12_interp%evaluate(y1,y2,y3,y4, PI1_12_i)
	end function PI1_12_i

	function PI1_13_i(y1, y2, y3, y4)
		implicit none
		real(dl) :: PI1_13_i
		real(dl), intent(in) :: y1, y2, y3, y4
		integer :: iflag

		PI1_13_i=0.d0
		call PI1_13_interp%evaluate(y1,y2,y3,y4, PI1_13_i)
	end function PI1_13_i

	elemental function PI1_12_full (y1, y2, y3, y4) !(y1, y2)
		real(dl) :: PI1_12_full
		real(dl), intent(in) :: y1, y2, y3, y4

		!pi_1(y1, y2): nu+nu -> e+e
		PI1_12_full = y1 * y2 * D1_full(y1, y2, y3, y4) - D2_full(y1, y2, y3, y4)
	end function PI1_12_full

	elemental function PI1_13_full (y1, y2, y3, y4) !(y1,y3)
		real(dl) :: PI1_13_full
		real(dl), intent(in) :: y1, y2, y3, y4

		!pi_1(y1, y3): nu+e -> nu+e
		PI1_13_full = y1 * y3 * D1_full(y1, y2, y3, y4) + D2_full(y1, y3, y2, y4)
	end function PI1_13_full

	pure function PI1_full (y1, y2, y3, y4) !1: (y1,y3), 	2: (y1, y2)
		real(dl), dimension(2) :: PI1_full
		real(dl), intent(in) :: y1, y2, y3, y4
		real(dl) :: d1

		d1 = D1_full(y1, y2, y3, y4)

		!pi_1(y1, y3): nu+e -> nu+e
		PI1_full(1) = y1 * y3 * d1 + D2_full(y1, y3, y2, y4)

		!pi_1(y1, y2): nu+nu -> e+e
		PI1_full(2) = y1 * y2 * d1 - D2_full(y1, y2, y3, y4)
	end function PI1_full
	
	pure function PI2_nn_f (y1, y2, y3, y4, e3, e4) !1: (y1, y4),	2: (y1, y3)
		real(dl), dimension(2) :: PI2_nn_f
		real(dl), intent(in) :: y1, y2, y3, y4, e3, e4
		real(dl) :: e1, e2, t

		e1 = y1
		e2 = y2
		t = e1 * e2 * e3 * e4 * D1_full(y1, y2, y3, y4) + D3_full(y1, y2, y3, y4)

		!pi_2(y1, y4)
		PI2_nn_f(1) = 2.d0 * (t &
			+ e2 * e3 * D2_full(y1, y4, y2, y3) &
			+ e1 * e4 * D2_full(y2, y3, y1, y4) )

		!pi_2(y1, y3)
		PI2_nn_f(2) = 2.d0 * (t &
			+ e1 * e3 * D2_full(y2, y4, y1, y3) &
			+ e2 * e4 * D2_full(y1, y3, y2, y4) )
	end function PI2_nn_f
	
	pure function PI2_ne_f (y1, y2, y3, y4, e2, e4) !1: (y1, y4),	2: (y1, y2)
		real(dl), dimension(2) :: PI2_ne_f
		real(dl), intent(in) :: y1, y2, y3, y4, e2, e4
		real(dl) :: e1, e3, t
		
		e1 = y1
		e3 = y3
		t = e1 * e2 * e3 * e4 * D1_full(y1, y2, y3, y4) + D3_full(y1, y2, y3, y4)
		
		!pi_2(y1, y4)
		PI2_ne_f(1) = 2.d0 * (t &
			+ e2 * e3 * D2_full(y1, y4, y2, y3) &
			+ e1 * e4 * D2_full(y2, y3, y1, y4) )
		
		!pi_2(y1, y2)
		PI2_ne_f(2) = 2.d0 * (t &
			- e1 * e2 * D2_full(y3, y4, y1, y2) &
			- e3 * e4 * D2_full(y1, y2, y3, y4) )
		
	end function PI2_ne_f
	
	function coll_nue_4_ann_int_re(ndim, ve, obj)
		!annihilation
		integer, intent(in) :: ndim
		real(dl), dimension(2), intent(in) :: ve
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_4_ann_int_re
		integer :: ix, iy
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y2,y3,y4, f3,f4, E3, E4, dme2
		type(cmplxMatNN) :: nB

		coll_nue_4_ann_int_re = 0.d0

		dme2 = obj%dme2
		y3 = ve(1)
		E3 = Ebare_i_dme(obj%x,y3, dme2)
		f3 = fermiDirac(E3 / obj%z)
		y4 = ve(2)
		E4 = Ebare_i_dme(obj%x,y4, dme2)
		f4 = fermiDirac(E4 / obj%z)

		y2 = E3 + E4 - obj%y1
		if (.not. (y2.lt.0.d0 &
				.or. obj%y1.gt.y3+y2+y4 &
				.or. y2.gt.obj%y1+y3+y4 &
				.or. y3.gt.obj%y1+y2+y4 &
				.or. y4.gt.obj%y1+y2+y3)) then
			call allocateCmplxMat(nB)
			do ix=1, flavorNumber
				nB%re(ix,ix) = interpNuDens%re(ix,ix)%evaluate(y2)
				nB%im(ix,ix) = 0.d0
				do iy=ix+1, flavorNumber
					nB%re(ix,iy) = interpNuDens%re(ix,iy)%evaluate(y2)
					nB%im(ix,iy) = interpNuDens%im(ix,iy)%evaluate(y2)
					nB%re(iy,ix) = nB%re(ix,iy)
					nB%im(iy,ix) = -nB%im(ix,iy)
				end do
			end do
			pi2_vec = PI2_nn_f (obj%y1, y2, y3, y4, E3, E4)
			coll_nue_4_ann_int_re = coll_nue_4_ann_int_re + &
				y3/E3 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab_ann_re(nuDensMatVecFD(obj%iy),nB,f3,f4, 1, 1, obj%ix1,obj%ix2) + &
					pi2_vec(2) * F_ab_ann_re(nuDensMatVecFD(obj%iy),nB,f3,f4, 2, 2, obj%ix1,obj%ix2) + &
					(obj%x*obj%x + dme2) * PI1_12_full(obj%y1, y2, y3, y4) * ( &
						F_ab_ann_re(nuDensMatVecFD(obj%iy),nB,f3,f4, 2, 1, obj%ix1,obj%ix2) + &
						F_ab_ann_re(nuDensMatVecFD(obj%iy),nB,f3,f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_4_ann_int_re

	function coll_nue_4_sc_int_re(ndim, ve, obj)
		!scattering, summing positron and electrons
		integer, intent(in) :: ndim
		real(dl), dimension(2), intent(in) :: ve
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_4_sc_int_re
		integer :: ix, iy
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y2,y3,y4, f2,f4, E2, E4, dme2
		type(cmplxMatNN) :: nB

		coll_nue_4_sc_int_re = 0.d0

		dme2 = obj%dme2
		y2 = ve(1)
		E2 = Ebare_i_dme(obj%x,y2, dme2)
		f2 = fermiDirac(E2 / obj%z)
		y4 = ve(2)
		E4 = Ebare_i_dme(obj%x,y4, dme2)
		f4 = fermiDirac(E4 / obj%z)

		y3 = obj%y1 + E2 - E4
		if (.not. (y3.lt.0.d0 &
				.or. obj%y1.gt.y2+y3+y4 &
				.or. y2.gt.obj%y1+y3+y4 &
				.or. y3.gt.obj%y1+y2+y4 &
				.or. y4.gt.obj%y1+y2+y3)) then
			call allocateCmplxMat(nB)
			do ix=1, flavorNumber
				nB%re(ix,ix) = interpNuDens%re(ix,ix)%evaluate(y3)
				nB%im(ix,ix) = 0.d0
				do iy=ix+1, flavorNumber
					nB%re(ix,iy) = interpNuDens%re(ix,iy)%evaluate(y3)
					nB%im(ix,iy) = interpNuDens%im(ix,iy)%evaluate(y3)
					nB%re(iy,ix) = nB%re(ix,iy)
					nB%im(iy,ix) = -nB%im(ix,iy)
				end do
			end do
			pi2_vec = PI2_ne_f (obj%y1, y2, y3, y4, E2, E4)
			coll_nue_4_sc_int_re = coll_nue_4_sc_int_re + &
				y2/E2 * &
				y4/E4 * &
				( &
					( pi2_vec(1) + pi2_vec(2) ) * ( & !F_sc^LL + F_sc^RR
						F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nB,f4, 1, 1, obj%ix1,obj%ix2) + &
						F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nB,f4, 2, 2, obj%ix1,obj%ix2) &
					) - &
					2.d0 * (obj%x*obj%x + dme2) * PI1_13_full(obj%y1, y2, y3, y4) * ( & !F_sc^RL and F_sc^LR
						F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nB,f4, 2, 1, obj%ix1,obj%ix2) + &
						F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nB,f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_4_sc_int_re

	function coll_nue_4_int_re(ndim, ve, obj)
		integer, intent(in) :: ndim
		real(dl), dimension(2), intent(in) :: ve
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_4_int_re
		coll_nue_4_int_re = &
			coll_nue_4_sc_int_re(ndim, ve, obj) &
			+ coll_nue_4_ann_int_re(ndim, ve, obj)
	end function coll_nue_4_int_re

	pure function coll_nue_3_ann_int_re(iy2, y4, obj)
		!annihilation
		integer, intent(in) :: iy2
		real(dl), intent(in) :: y4
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_3_ann_int_re
		integer :: ix, iy
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y2,y3, f3,f4, E3, E4, dme2, t1, t2

		coll_nue_3_ann_int_re = 0.d0

		dme2 = obj%dme2
		y2 = y_arr(iy2)
		E4 = Ebare_i_dme(obj%x,y4, dme2)
		E3 = obj%y1 + y2 - E4
		t1 = E3*E3
		t2 = obj%x*obj%x + dme2
		if (t1<t2) then
			y3 = -1.d0
		else
			y3 = sqrt(t1 - t2)
		endif
!		write(*,multidblfmt) obj%y1, y2, y3, y4
!		print*, obj%iy, iy2
!		write(*,multidblfmt) nuDensMatVecFD(obj%iy)
!		write(*,multidblfmt) nuDensMatVecFD(iy2)
		if (.not.(y3.lt.0.d0 &
				.or. obj%y1.gt.y3+y2+y4 &
				.or. y2.gt.obj%y1+y3+y4 &
				.or. y3.gt.obj%y1+y2+y4 &
				.or. y4.gt.obj%y1+y2+y3)) then
			f3 = fermiDirac(E3 / obj%z)
			f4 = fermiDirac(E4 / obj%z)
			pi2_vec = PI2_nn_f(obj%y1, y2, y3, y4, E3, E4)
!			write(*,multidblfmt) pi2_vec, PI1_12_full(obj%y1, y2, y3, y4)
!			write(*,multidblfmt)F_ab_ann_re(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 1, 1, obj%ix1,obj%ix2), F_ab_ann_re(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 2, 2, obj%ix1,obj%ix2),F_ab_ann_re(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 2, 1, obj%ix1,obj%ix2), F_ab_ann_re(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 1, 2, obj%ix1,obj%ix2)
			coll_nue_3_ann_int_re = coll_nue_3_ann_int_re + &
				y3/E3 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab_ann_re(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 1, 1, obj%ix1,obj%ix2) &
					+ pi2_vec(2) * F_ab_ann_re(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 2, 2, obj%ix1,obj%ix2) &
					+ (obj%x*obj%x + dme2) * PI1_12_full(obj%y1, y2, y3, y4) * ( &
						F_ab_ann_re(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 2, 1, obj%ix1,obj%ix2) &
						+ F_ab_ann_re(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_3_ann_int_re

	pure function coll_nue_3_sc_int_re(iy3, y2, obj)
		!scattering, summing positron and electrons
		integer, intent(in) :: iy3
		real(dl), intent(in) :: y2
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_3_sc_int_re
		integer :: ix, iy
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y3,y4, f2,f4, E2, E4, dme2, t1, t2

		coll_nue_3_sc_int_re = 0.d0

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
!		write(*,multidblfmt) obj%y1, y2, y3, y4, dme2, E2, E4
!		print*, obj%iy, iy3
!		write(*,multidblfmt) nuDensMatVecFD(obj%iy)
!		write(*,multidblfmt) nuDensMatVecFD(iy3)
		if (.not.(y4.lt.0.d0 &
				.or. obj%y1.gt.y2+y3+y4 &
				.or. y2.gt.obj%y1+y3+y4 &
				.or. y3.gt.obj%y1+y2+y4 &
				.or. y4.gt.obj%y1+y2+y3)) then
			f2 = fermiDirac(E2 / obj%z)
			f4 = fermiDirac(E4 / obj%z)
			pi2_vec = PI2_ne_f (obj%y1, y2, y3, y4, E2, E4)
!			write(*,multidblfmt) pi2_vec, PI1_13_full(obj%y1, y2, y3, y4)
!			write(*,multidblfmt)F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 1, 1, obj%ix1,obj%ix2),F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 2, 2, obj%ix1,obj%ix2),F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 2, 1, obj%ix1,obj%ix2),F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 1, 2, obj%ix1,obj%ix2)
			coll_nue_3_sc_int_re = coll_nue_3_sc_int_re + &
				y2/E2 * &
				y4/E4 * &
				( &
					( pi2_vec(1) + pi2_vec(2) ) * ( & !F_sc^LL + F_sc^RR
						F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 1, 1, obj%ix1,obj%ix2) &
						+ F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 2, 2, obj%ix1,obj%ix2) &
					) &
					- 2.d0 * (obj%x*obj%x + dme2) * PI1_13_full(obj%y1, y2, y3, y4) * ( & !F_sc^RL and F_sc^LR
						F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 2, 1, obj%ix1,obj%ix2) &
						+ F_ab_sc_re(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_3_sc_int_re

	pure function coll_nue_3_int_re(iy, yx, obj)
		integer, intent(in) :: iy
		real(dl), intent(in) :: yx
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_3_int_re
		coll_nue_3_int_re = &
			coll_nue_3_sc_int_re(iy, yx, obj) &
				+ coll_nue_3_ann_int_re(iy, yx, obj)
	end function coll_nue_3_int_re

	pure function coll_nue_3_ann_int_im(iy2, y4, obj)
		!annihilation
		integer, intent(in) :: iy2
		real(dl), intent(in) :: y4
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_3_ann_int_im
		integer :: ix, iy
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y2,y3, f3,f4, E3, E4, dme2, t1, t2

		coll_nue_3_ann_int_im = 0.d0

		dme2 = obj%dme2
		y2 = y_arr(iy2)
		E4 = Ebare_i_dme(obj%x,y4, dme2)

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
			coll_nue_3_ann_int_im = coll_nue_3_ann_int_im + &
				y3/E3 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab_ann_im(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 1, 1, obj%ix1,obj%ix2) &
					+ pi2_vec(2) * F_ab_ann_im(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 2, 2, obj%ix1,obj%ix2) &
					+ (obj%x*obj%x + dme2) * PI1_12_full(obj%y1, y2, y3, y4) * ( &
						F_ab_ann_im(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 2, 1, obj%ix1,obj%ix2) &
						+ F_ab_ann_im(nuDensMatVecFD(obj%iy),nuDensMatVecFD(iy2),f3,f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_3_ann_int_im

	pure function coll_nue_3_sc_int_im(iy3, y2, obj)
		!scattering, summing positron and electrons
		integer, intent(in) :: iy3
		real(dl), intent(in) :: y2
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_3_sc_int_im
		integer :: ix, iy
		real(dl), dimension(2) :: pi2_vec
		real(dl) :: y3,y4, f2,f4, E2, E4, dme2, t1, t2

		coll_nue_3_sc_int_im = 0.d0

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
			coll_nue_3_sc_int_im = coll_nue_3_sc_int_im + &
				y2/E2 * &
				y4/E4 * &
				( &
					( pi2_vec(1) + pi2_vec(2) ) * ( & !F_sc^LL + F_sc^RR
						F_ab_sc_im(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 1, 1, obj%ix1,obj%ix2) &
						+ F_ab_sc_im(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 2, 2, obj%ix1,obj%ix2) &
					) &
					- 2.d0 * (obj%x*obj%x + dme2) * PI1_13_full(obj%y1, y2, y3, y4) * ( & !F_sc^RL and F_sc^LR
						F_ab_sc_im(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 2, 1, obj%ix1,obj%ix2) &
						+ F_ab_sc_im(nuDensMatVecFD(obj%iy),f2,nuDensMatVecFD(iy3),f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_3_sc_int_im

	pure function coll_nue_3_int_im(iy, yx, obj)
		integer, intent(in) :: iy
		real(dl), intent(in) :: yx
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_3_int_im
		coll_nue_3_int_im = &
			coll_nue_3_sc_int_im(iy, yx, obj) &
				+ coll_nue_3_ann_int_im(iy, yx, obj)
	end function coll_nue_3_int_im

	pure function integrate_coll_int_3(f, obj)
		interface
			pure real(dl) function f(a, b, o)
				use variables
				integer, intent(in) :: a
				real(dl), intent(in) :: b
				type(coll_args), intent(in) :: o
			end function
		end interface
		real(dl) :: integrate_coll_int_3
		type(coll_args), intent(in) :: obj
		integer :: ia, ib
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(fy2_arr(Ny, Ny))

		do ia=1, Ny
			do ib=1, Ny
				fy2_arr(ia,ib) = f(ia, y_arr(ib), obj)
			end do
		end do
		integrate_coll_int_3 = integral_linearized_2d(Ny, Ny, dy_arr, dy_arr, fy2_arr)
		deallocate(fy2_arr)
	end function integrate_coll_int_3

	function coll_nue_int_im(ndim, ve, obj)
		integer, intent(in) :: ndim
		real(dl), dimension(2), intent(in) :: ve
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_int_im
		integer :: ix, iy
		real(dl), dimension(2) :: pi2_vec
		type(cmplxMatNN) :: nX
		real(dl) :: yA, y2,y3,y4, f2,f3,f4, E2, E3, E4, dme2

		coll_nue_int_im = 0.0
		call allocateCmplxMat(nX)

		yA = ve(1)
		y4 = ve(2)
		dme2 = obj%dme2
		E4 = Ebare_i_dme(obj%x,y4, dme2)

		!scattering, summing positrons and electrons
		y2 = yA
		E2 = Ebare_i_dme(obj%x,y2, dme2)
		y3 = obj%y1 + E2 - E4
		f2 = fermiDirac(E2 / obj%z)
		f4 = fermiDirac(E4 / obj%z)
		if (y3.lt.0.d0 &
			.or. obj%y1.gt.y2+y3+y4 &
			.or. y2.gt.obj%y1+y3+y4 &
			.or. y3.gt.obj%y1+y2+y4 &
			.or. y4.gt.obj%y1+y2+y3) then
			coll_nue_int_im = coll_nue_int_im + 0.d0
		else
!			fd = fermiDirac(y3 / obj%z)
			do ix=1, flavorNumber
				nX%re(ix,ix) = interpNuDens%re(ix,ix)%evaluate(y3)! * fd
				nX%im(ix,ix) = 0.d0
				do iy=ix+1, flavorNumber
					nX%re(ix,iy) = interpNuDens%re(ix,iy)%evaluate(y3)
					nX%im(ix,iy) = interpNuDens%im(ix,iy)%evaluate(y3)
					nX%re(iy,ix) = nX%re(ix,iy)
					nX%im(iy,ix) = -nX%im(ix,iy)
				end do
			end do

			pi2_vec = PI2_ne_f (obj%y1, y2, y3, y4, E2, E4)

			coll_nue_int_im = coll_nue_int_im + &
				y2/E2 * &
				y4/E4 * &
				( &
					( pi2_vec(1) * pi2_vec(2) ) * ( & !F_sc^LL + F_sc^RR
						F_ab_sc_im(nuDensMatVecFD(obj%iy),f2,nX,f4, 1, 1, obj%ix1,obj%ix2) + &
						F_ab_sc_im(nuDensMatVecFD(obj%iy),f2,nX,f4, 2, 2, obj%ix1,obj%ix2) &
					) - &
					2.d0 * (obj%x*obj%x + dme2) * PI1_13_full(obj%y1, y2, y3, y4) * ( &!F_sc^RL and F_sc^LR
						F_ab_sc_im(nuDensMatVecFD(obj%iy),f2,nX,f4, 2,1, obj%ix1,obj%ix2) + &
						F_ab_sc_im(nuDensMatVecFD(obj%iy),f2,nX,f4, 1,2, obj%ix1,obj%ix2) ) &
				)
		end if

		!annihilation
		y3=yA
		E3 = Ebare_i_dme(obj%x,y3, dme2)
		y2 = E3 + E4 - obj%y1
		f3 = fermiDirac(E3 / obj%z)
		f4 = fermiDirac(E4 / obj%z)
		if (y2.lt.0.d0 &
			.or. obj%y1.gt.y3+y2+y4 &
			.or. y2.gt.obj%y1+y3+y4 &
			.or. y3.gt.obj%y1+y2+y4 &
			.or. y4.gt.obj%y1+y2+y3) then
			coll_nue_int_im = coll_nue_int_im + 0.d0
		else
!			fd = fermiDirac(y2 / obj%z)
			do ix=1, flavorNumber
				nX%re(ix,ix) = interpNuDens%re(ix,ix)%evaluate(y2)! * fd
				nX%im(ix,ix) = 0.d0
				do iy=ix+1, flavorNumber
					nX%re(ix,iy) = interpNuDens%re(ix,iy)%evaluate(y2)
					nX%im(ix,iy) = interpNuDens%im(ix,iy)%evaluate(y2)
					nX%re(iy,ix) = nX%re(ix,iy)
					nX%im(iy,ix) = -nX%im(ix,iy)
				end do
			end do

			pi2_vec = PI2_nn_f (obj%y1, y2, y3, y4, E3, E4)

			coll_nue_int_im = coll_nue_int_im + &
				y3/E3 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab_ann_im(nuDensMatVecFD(obj%iy),nX,f3,f4, 1, 1, obj%ix1,obj%ix2) + &
					pi2_vec(2) * F_ab_ann_im(nuDensMatVecFD(obj%iy),nX,f3,f4, 2, 2, obj%ix1,obj%ix2) + &
					(obj%x*obj%x + dme2) * PI1_12_full(obj%y1, y2, y3, y4) * ( &
						F_ab_ann_im(nuDensMatVecFD(obj%iy),nX,f3,f4, 2, 1, obj%ix1,obj%ix2) + &
						F_ab_ann_im(nuDensMatVecFD(obj%iy),nX,f3,f4, 1, 2, obj%ix1,obj%ix2) ) &
				)
		end if
	end function coll_nue_int_im

	pure SUBROUTINE region(ndim,x,j,c,d)
		use precision
		IMPLICIT NONE
		REAL (dl), INTENT (OUT) :: c, d
		INTEGER, INTENT (IN)    :: j, ndim
		REAL (dl), INTENT (IN)  :: x(ndim)
		c = y_min
		d = y_max
		RETURN
	END SUBROUTINE region

	pure function get_collision_terms(collArgsIn, Fre, Fim)
		interface
			pure real(dl) function Fre(a, b, o)
				use variables
				integer, intent(in) :: a
				real(dl), intent(in) :: b
				type(coll_args), intent(in) :: o
			end function
			pure real(dl) function Fim(a, b, o)
				use variables
				integer, intent(in) :: a
				real(dl), intent(in) :: b
				type(coll_args), intent(in) :: o
			end function
		end interface
		type(cmplxMatNN) :: get_collision_terms
		real(dl) :: x,z, y1
		integer :: iy1
		type(coll_args), intent(in) :: collArgsIn
		type(coll_args) :: collArgs
		integer :: i,j, k
		real(dl) :: cf, res1, res2

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

		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			if (.not.sterile(i)) then
				get_collision_terms%re(i,i) = integrate_coll_int_3(Fre, collArgs)
			end if
			if (collision_offdiag.eq.1) then
				do j=i+1, flavorNumber
					collArgs%ix2 = j
					get_collision_terms%re(i,j) = integrate_coll_int_3(Fre, collArgs)

					get_collision_terms%im(i,j) = integrate_coll_int_3(Fim, collArgs)
				end do
			else if (collision_offdiag.eq.2) then
				!damping terms from Dolgov:2002ab, eq A.11
				do j=i+1, flavorNumber
					collArgs%ix2 = j
					get_collision_terms%re(i,j) = 0.d0
					get_collision_terms%im(i,j) = 0.d0
!						if (i .eq. 1 .and. (j.eq.2 .or. j.eq.3)) then
!							get_collision_terms%re(i,j) = dampTermFactor * Damp_ex * y1*y1*y1 * nuDensMatVecFD(iy1)%re(i,j)
!						elseif (i .eq. 2 .and. j.eq.3) then
!							get_collision_terms%re(i,j) = dampTermFactor * Damp_mt * y1*y1*y1 * nuDensMatVecFD(iy1)%re(i,j)
!						end if
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
