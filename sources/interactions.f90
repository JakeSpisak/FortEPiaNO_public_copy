module ndInteractions
	use precision
	use constants
	use variables
	use utilities
	use ndErrors
	use ndMatrices
!	use ndconfig
	use linear_interpolation_module
!	use bspline_module
	implicit none
	
	type(spline_class) :: FD_interp
	type(linear_interp_2d) :: dmeCorr
!	type(bspline_2d) :: dmeCorr
	
	type(linear_interp_4d) :: D2_interp, D3_interp, PI1_12_interp, PI1_13_interp

	contains
	
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
!		call dmeCorr%initialize(interp_xvec,interp_zvec,dme_vec,4,4,iflag)!bspline
		
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
			call tic(timer1)
			write (*,*) "[interactions] now doing some timing..."
			call tic(timer1)
			do ix=1, 1000000
				call random_number(x)
				x=10.d0**(x*(2.d0-(-3.d0)) - 3.d0)
				t1 = fermiDirac_i(x)
			end do
			call toc(timer1, "<interpolated>")

			call tic(timer1)
			do ix=1, 1000000
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
			
			tmp = rombint_vec(vec, dme2_e_i1, fe_l, fe_u, 1d-3, maxiter)
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
!		call dmeCorr%evaluate(x,z,0,0,dme2_electron,iflag)!bspline
	end function dme2_electron
	
	elemental function Ebare_i_dme(x,y,dme2)!for electrons
		real(dl) :: Ebare_i_dme
		real(dl), intent(in) :: x,y,dme2
		
		Ebare_i_dme = sqrt(x*x+y*y+dme2)
	end function Ebare_i_dme
	
	pure function F_ab_ann(x,z, n1,n2,e3,e4, a,b, i, j, rI)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.5
		real(dl) :: F_ab_ann
		type(cmplxMatNN), intent(in) :: n1,n2
		real(dl), intent(in) :: e3,e4, x,z
		integer, intent(in) :: a,b, i, j, rI
		integer :: k
		real(dl) :: f3, f4
		real(dl) :: t1a, t1b, t2a, t2b
		
!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_ann] a and b must be either 1(=L) or 2(=R)")

		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
		if (rI.eq.1) then !take real part
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
		else if (rI .eq.2) then!take imaginary part
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
!		else
!			call criticalError("[F_ab_ann] invalid rI")
		end if
		f3 = fermiDirac(E_k_m(e3, x) / z)
		f4 = fermiDirac(E_k_m(e4, x) / z)
		F_ab_ann = f3 * f4 * &
				(t1a * GLR_vec(a,i,i) + t1b * GLR_vec(a,j,j))&
			- (1.d0 - f3) * (1.d0 - f4) * &
				(t2a * GLR_vec(a,j,j) + t2b * GLR_vec(a,i,i))
	end function F_ab_ann
	
	pure function F_ab_sc (x,z, n1,e2,n3,e4, a,b, i, j, rI)!a, b must be either 1(=L) or 2(=R)
	!doi:10.1088/1475-7516/2016/07/051 eq. 2.10
		real(dl) :: F_ab_sc
		type(cmplxMatNN), intent(in) :: n1,n3
		real(dl), intent(in) :: e2,e4, x,z
		integer, intent(in) :: a,b, i, j, rI
		integer :: k
		real(dl) :: f2, f4
		real(dl) :: t1a, t1b, t2a, t2b
		
!		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
!			call criticalError("[F_ab_sc] a and b must be either 1(=L) or 2(=R)")
		
		t1a=0.d0
		t1b=0.d0
		t2a=0.d0
		t2b=0.d0
		if (rI.eq.1) then !take real part
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
		else if (rI .eq.2) then!take imaginary part
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
!		else
!			call criticalError("[F_ab_sc] invalid rI")
		end if
		f2 = fermiDirac(E_k_m(e2, x) / z)
		f4 = fermiDirac(E_k_m(e4, x) / z)
		F_ab_sc = (1.d0 - f2) * f4 * &
				((t1a) * GLR_vec(a,i,i) + (t1b) * GLR_vec(a,j,j))&
			- f2 * (1.d0 - f4) * &
				((t2a) * GLR_vec(a,j,j) + (t2b) * GLR_vec(a,i,i))
	end function F_ab_sc
	
	!same as F_ab_ann, but with full matrix products (useful for nonstandard gL, gR)
	function F_ab_ann_mp (x,z, n1,n2,e3,e4, a,b, i, j, rI)!a, b must be either 1(=L) or 2(=R)
		real(dl) :: F_ab_ann_mp
		type(cmplxMatNN), intent(in) :: n1,n2
		real(dl), dimension(:,:), allocatable, save :: t1,t2,t3,t4,p1,p2,idm1,idm2,su!
		real(dl), intent(in) :: e3,e4, x,z
		integer, intent(in) :: a,b, i, j, rI
		integer :: k,f
		real(dl) :: f3, f4
		
		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
			call criticalError("[F_ab_ann_mp] a and b must be either 1(=L) or 2(=R)")
		
		f=flavorNumber
		if (.not.allocated(t1)) &
			allocate(t1(f,f))
		if (.not.allocated(t2)) &
			allocate(t2(f,f))
		if (.not.allocated(t3)) &
			allocate(t3(f,f))
		if (.not.allocated(t4)) &
			allocate(t4(f,f))
		if (.not.allocated(p1)) &
			allocate(p1(f,f))
		if (.not.allocated(p2)) &
			allocate(p2(f,f))
		if (.not.allocated(idm1)) &
			allocate(idm1(f,f))
		if (.not.allocated(idm2)) &
			allocate(idm2(f,f))
		if (.not.allocated(su)) &
			allocate(su(f,f))
			
		su=0.d0
		idm1=idMat-n1%re
		idm2=idMat-n2%re
		if (rI.eq.1) then !take real part
			!t1 = Ga (1-rho2) Gb (1-rho1)
			call quadrupleProdMat(GLR_vec(a,:,:), idm2 , GLR_vec(b,:,:), idm1 , p1)
			call quadrupleProdMat(GLR_vec(a,:,:), n2%im, GLR_vec(b,:,:), n1%im, p2)
			t1 = p1 - p2
			!t2 = (1-rho1) Gb (1-rho2) Ga
			call quadrupleProdMat(idm1 , GLR_vec(b,:,:), idm2 , GLR_vec(a,:,:), p1)
			call quadrupleProdMat(n1%im, GLR_vec(b,:,:), n2%im, GLR_vec(a,:,:), p2)
			t2 = p1 - p2
			!t3 = rho1 Gb rho2 Ga
			call quadrupleProdMat(n1%re, GLR_vec(b,:,:), n2%re, GLR_vec(a,:,:), p1)
			call quadrupleProdMat(n1%im, GLR_vec(b,:,:), n2%im, GLR_vec(a,:,:), p2)
			t3 = p1 - p2
			!t4 = Ga rho2 Gb rho1
			call quadrupleProdMat(GLR_vec(a,:,:), n2%re, GLR_vec(b,:,:), n1%re, p1)
			call quadrupleProdMat(GLR_vec(a,:,:), n2%im, GLR_vec(b,:,:), n1%im, p2)
			t4 = p1 - p2
		else if (rI .eq.2) then!take imaginary part
			!t1 = Ga (1-rho2) Gb (1-rho1)
			call quadrupleProdMat(GLR_vec(a,:,:), idm2 , GLR_vec(b,:,:), n1%im, p1)
			call quadrupleProdMat(GLR_vec(a,:,:), n2%im, GLR_vec(b,:,:), idm1 , p2)
			t1 = - p1 - p2
			!t2 = (1-rho1) Gb (1-rho2) Ga
			call quadrupleProdMat(idm1 , GLR_vec(b,:,:), n2%im, GLR_vec(a,:,:), p1)
			call quadrupleProdMat(n1%im, GLR_vec(b,:,:), idm2 , GLR_vec(a,:,:), p2)
			t2 = - p1 - p2
			!t3 = rho1 Gb rho2 Ga
			call quadrupleProdMat(n1%re, GLR_vec(b,:,:), n2%im, GLR_vec(a,:,:), p1)
			call quadrupleProdMat(n1%im, GLR_vec(b,:,:), n2%re, GLR_vec(a,:,:), p2)
			t3 = p1 + p2
			!t4 = Ga rho2 Gb rho1
			call quadrupleProdMat(GLR_vec(a,:,:), n2%re, GLR_vec(b,:,:), n1%im, p1)
			call quadrupleProdMat(GLR_vec(a,:,:), n2%im, GLR_vec(b,:,:), n1%re, p2)
			t4 = p1 + p2
		else
			call criticalError("[F_ab_ann_mp] invalid rI")
		end if
		f3 = fermiDirac(E_k_m(e3, x) / z)
		f4 = fermiDirac(E_k_m(e4, x) / z)
		su =  f3 * f4 * (t1 + t2) &
			- (1.d0 - f3) * (1.d0 - f4) * (t3 + t4)
		F_ab_ann_mp = su(i,j)
	end function F_ab_ann_mp
	
	!same as F_ab_sc, but with full matrix products (useful for nonstandard gL, gR)
	function F_ab_sc_mp (x,z, n1,e2,n3,e4, a,b, i, j, rI)!a, b must be either 1(=L) or 2(=R)
		real(dl) :: F_ab_sc_mp
		type(cmplxMatNN), intent(in) :: n1,n3
		real(dl), dimension(:,:), allocatable, save :: t1,t2,t3,t4,p1,p2,idm1,idm3,su!
		real(dl), intent(in) :: e2,e4, x,z
		integer, intent(in) :: a,b, i, j, rI
		integer :: k,f
		real(dl) :: f2, f4
		
		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
			call criticalError("[F_ab_sc_mp] a and b must be either 1(=L) or 2(=R)")
		
		f=flavorNumber
		if (.not.allocated(t1)) &
			allocate(t1(f,f))
		if (.not.allocated(t2)) &
			allocate(t2(f,f))
		if (.not.allocated(t3)) &
			allocate(t3(f,f))
		if (.not.allocated(t4)) &
			allocate(t4(f,f))
		if (.not.allocated(p1)) &
			allocate(p1(f,f))
		if (.not.allocated(p2)) &
			allocate(p2(f,f))
		if (.not.allocated(idm1)) &
			allocate(idm1(f,f))
		if (.not.allocated(idm3)) &
			allocate(idm3(f,f))
		if (.not.allocated(su)) &
			allocate(su(f,f))
			
		su=0.d0
		idm1=idMat-n1%re
		idm3=idMat-n3%re
		if (rI.eq.1) then !take real part
			!t1 = Ga rho3 Gb (1-rho1)
			call quadrupleProdMat(GLR_vec(a,:,:), n3%re, GLR_vec(b,:,:), idm1 , p1)
			call quadrupleProdMat(GLR_vec(a,:,:), n3%im, GLR_vec(b,:,:), n1%im, p2)
			t1 = p1 + p2
			!t2 = (1-rho1) Gb rho3 Ga
			call quadrupleProdMat(idm1 , GLR_vec(b,:,:), n3%re, GLR_vec(a,:,:), p1)
			call quadrupleProdMat(n1%im, GLR_vec(b,:,:), n3%im, GLR_vec(a,:,:), p2)
			t2 = p1 + p2
			!t3 = rho1 Gb (1-rho3) Ga
			call quadrupleProdMat(n1%re, GLR_vec(b,:,:), idm3 , GLR_vec(a,:,:), p1)
			call quadrupleProdMat(n1%im, GLR_vec(b,:,:), n3%im, GLR_vec(a,:,:), p2)
			t3 = p1 + p2
			!t4 = Ga (1-rho3) Gb rho1
			call quadrupleProdMat(GLR_vec(a,:,:), idm3 , GLR_vec(b,:,:), n1%re, p1)
			call quadrupleProdMat(GLR_vec(a,:,:), n3%im, GLR_vec(b,:,:), n1%im, p2)
			t4 = p1 + p2
		else if (rI .eq.2) then!take imaginary part
			!t1 = Ga rho3 Gb (1-rho1)
			call quadrupleProdMat(GLR_vec(a,:,:), n3%re, GLR_vec(b,:,:), n1%im, p1)
			call quadrupleProdMat(GLR_vec(a,:,:), n3%im, GLR_vec(b,:,:), idm1 , p2)
			t1 = - p1 + p2
			!t2 = (1-rho1) Gb rho3 Ga
			call quadrupleProdMat(idm1 , GLR_vec(b,:,:), n3%im, GLR_vec(a,:,:), p1)
			call quadrupleProdMat(n1%im, GLR_vec(b,:,:), n3%re, GLR_vec(a,:,:), p2)
			t2 = p1 - p2
			!t3 = rho1 Gb (1-rho3) Ga
			call quadrupleProdMat(n1%re, GLR_vec(b,:,:), n3%im, GLR_vec(a,:,:), p1)
			call quadrupleProdMat(n1%im, GLR_vec(b,:,:), idm3 , GLR_vec(a,:,:), p2)
			t3 = - p1 + p2
			!t4 = Ga (1-rho3) Gb rho1
			call quadrupleProdMat(GLR_vec(a,:,:), idm3 , GLR_vec(b,:,:), n1%im, p1)
			call quadrupleProdMat(GLR_vec(a,:,:), n3%im, GLR_vec(b,:,:), n1%re, p2)
			t4 = p1 - p2
		else
			call criticalError("[F_ab_sc_mp] invalid rI")
		end if
		f2 = fermiDirac(E_k_m(e2, x) / z)
		f4 = fermiDirac(E_k_m(e4, x) / z)
		su =  (1.d0 - f2) * f4 * (t1 + t2) &
			- f2 * (1.d0 - f4) * (t3 + t4)
		F_ab_sc_mp = su(i,j)
	end function F_ab_sc_mp

	subroutine init_interp_d123
		real(dl), dimension(:,:,:,:), allocatable :: d2, d3, pi1_12, pi1_13
		real(dl) :: y1,  y2,  y3,  y4, d1a, d1b, d2a, d2b, d3a, d3b
		integer  :: ix1, ix2, ix3, ix4, iflag
		real(dl), dimension(2) :: PI1
		real(dl), dimension(:), allocatable :: yarr
		real(8) :: timer1
		
		call addToLog("[interactions] Initializing interpolation for D1/2/3 functions...")
		allocate(yarr(interp_ny))
		yarr = logspace(logy_min, logy_max, interp_ny)
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
				d2a = d2_f(y1,y2,y3,y4)
				d2a = d3_f(y1,y2,y3,y4)
				d2a = PI1_12_i(y1,y2,y3,y4)
				d3a = PI1_13_i(y1,y2,y3,y4)
			end do
			call toc(timer1, "<interpolated>")

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
				d2a = d2_full(y1,y2,y3,y4)
				d2a = d3_full(y1,y2,y3,y4)
				d2a = PI1_12_full(y1,y2,y3,y4)
				d3a = PI1_13_full(y1,y2,y3,y4)
			end do
			call toc(timer1, "<full>")
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

	elemental function PI1_12_full (y1, y2, y3, y4) !1: (y1,y3), 	2: (y1, y2)
		real(dl) :: PI1_12_full
		real(dl), intent(in) :: y1, y2, y3, y4

		!pi_1(y1, y2): nu+nu -> e+e
		PI1_12_full = y1 * y2 * D1_full(y1, y2, y3, y4) - D2_full(y1, y2, y3, y4)
	end function PI1_12_full

	elemental function PI1_13_full (y1, y2, y3, y4) !1: (y1,y3)
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
		PI2_nn_f(1) = 2 * (t &
			+ e2 * e3 * D2_full(y1, y4, y2, y3) &
			+ e1 * e4 * D2_full(y2, y3, y1, y4) )

		!pi_2(y1, y3)
		PI2_nn_f(2) = 2 * (t &
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
		PI2_ne_f(1) = 2 * (t &
			+ e2 * e3 * D2_full(y1, y4, y2, y3) &
			+ e1 * e4 * D2_full(y2, y3, y1, y4) )
		
		!pi_2(y1, y2)
		PI2_ne_f(2) = 2 * (t &
			- e1 * e2 * D2_full(y3, y4, y1, y2) &
			- e3 * e4 * D2_full(y1, y2, y3, y4) )
		
	end function PI2_ne_f

	function interp_nuDens(y)
		type(cmplxMatNN) :: interp_nuDens
		real(dl), intent(in) :: y
		integer :: ix,i,j,k
		real(dl) :: ixr, yr, logy
		character(len=100) :: vals
		
		call allocateCmplxMat(interp_nuDens)
		if (y.gt.nuDensMatVec(Ny)%y .or. y .lt. nuDensMatVec(1)%y) then
			interp_nuDens%re =0.d0
			interp_nuDens%im =0.d0
		elseif (abs(y-nuDensMatVec(Ny)%y) .lt. 1.d-6) then
			interp_nuDens = nuDensMatVec(Ny)
		else
			logy=log10(y)
			ixr=(Ny-1)*(logy-logy_min)/(logy_max-logy_min)+1
			ix=int(ixr)
			if (ix .gt. 0 .and. ix .lt. Ny) then
				if (abs(ixr-ix) .lt. 1d-5) then
					interp_nuDens = nuDensMatVec(ix)
				else
					yr = (logy - nuDensMatVec(ix)%logy) / (nuDensMatVec(ix+1)%logy - nuDensMatVec(ix)%logy)
					interp_nuDens%y = y
					interp_nuDens%re(:,:) = (1.d0 - yr) * nuDensMatVec(ix)%re(:,:) + yr * nuDensMatVec(ix+1)%re(:,:)
					interp_nuDens%im(:,:) = (1.d0 - yr) * nuDensMatVec(ix)%im(:,:) + yr * nuDensMatVec(ix+1)%im(:,:)
				end if
			else
				write(vals,"('y=',"//dblfmt//",'    ix=',I3)") y, ix
				call Error("k out of range or errors occurred in interpolation"//NEW_LINE('A')//vals)
			end if
		end if
	end function interp_nuDens
	
	function interp_nuDensIJre(y, i, j)
		real(dl) :: interp_nuDensIJre
		real(dl), intent(in) :: y
		integer :: ix,i,j,k
		real(dl) :: ixr, yr, logy
		character(len=100) :: vals
		
		if (y.gt.nuDensMatVec(Ny)%y .or. y .lt. nuDensMatVec(1)%y) then
			interp_nuDensIJre =0.d0
		elseif (abs(y - nuDensMatVec(Ny)%y).lt.1d-5) then
			interp_nuDensIJre = nuDensMatVec(Ny)%re(i,j)
		else
			logy=log10(y)
			ixr=(Ny-1)*(logy-logy_min)/(logy_max-logy_min)+1
			ix=int(ixr)
!			if (ix .gt. 0 .and. ix .lt. Ny) then
!				if (abs(ixr-ix) .lt. 1d-6) then
!					interp_nuDensIJre = nuDensMatVec(ix)%re(i,j)
!				else
					yr = (logy - nuDensMatVec(ix)%logy) / (nuDensMatVec(ix+1)%logy - nuDensMatVec(ix)%logy)
					interp_nuDensIJre = nuDensMatVec(ix)%re(i,j) + &
						yr * (nuDensMatVec(ix+1)%re(i,j) - nuDensMatVec(ix)%re(i,j))
!				end if
!			else
!				write(vals,"('y=',"//dblfmt//",'    ix=',I3)") y, ix
!				call Error("k out of range or errors occurred in interpolation"//NEW_LINE('A')//vals)
!			end if
		end if
!      dataIndex = floor((xVal(inputIndex)-minXData)/xRange);

!      weight = (xVal - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex));
!      yVal(inputIndex) = (1.0-weight)*yData(dataIndex) + ...
!                         weight*yData(dataIndex+1);
	end function interp_nuDensIJre
	
	function coll_nue_sc_int(ndim, ve, obj)
		integer, intent(in) :: ndim
		real(dl), dimension(2), intent(in) :: ve
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_sc_int
		integer :: ix
		real(dl), dimension(2) :: pi2_vec
		type(cmplxMatNN) :: n3
		real(dl) :: y2,y3,y4, E2, E4, dme2, fd
		
		y2 = ve(1)
		y4 = ve(2)
		dme2 = obj%dme2
		E2 = Ebare_i_dme(obj%x,y2, dme2)
		E4 = Ebare_i_dme(obj%x,y4, dme2)
		y3 = obj%y1 + E2 - E4
		if (y3.lt.0.d0 &
			.or. obj%y1.gt.y2+y3+y4 &
			.or. y2.gt.obj%y1+y3+y4 &
			.or. y3.gt.obj%y1+y2+y4 &
			.or. y4.gt.obj%y1+y2+y3) then
			coll_nue_sc_int = 0.0
		else
			n3 = interp_nuDens(y3)
			fd = fermiDirac(y3 / obj%z)
			do ix=1, flavorNumber
				n3%re(ix,ix) = n3%re(ix,ix) * fd
				n3%im(ix,ix) = n3%im(ix,ix) * fd
			end do

			pi2_vec = PI2_ne_f (obj%y1, y2, y3, y4, E2, E4)

			coll_nue_sc_int = &
				y2/E2 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab_sc(obj%x,obj%z,obj%n,y2,n3,y4, obj%a,obj%a, obj%ix1,obj%ix2, obj%rI) + & !F_sc^LL or F_sc^RR
					pi2_vec(2) * F_ab_sc(obj%x,obj%z,obj%n,y2,n3,y4, obj%b,obj%b, obj%ix1,obj%ix2, obj%rI) - & !F_sc^RR or F_sc^LL
					(obj%x*obj%x + dme2) * PI1_13_full(obj%y1, y2, y3, y4) * ( &
						F_ab_sc(obj%x,obj%z,obj%n,y2,n3,y4, 2,1, obj%ix1,obj%ix2, obj%rI) + &!F_sc^RL and F_sc^LR
						F_ab_sc(obj%x,obj%z,obj%n,y2,n3,y4, 1,2, obj%ix1,obj%ix2, obj%rI) ) &
				)
		end if
	end function coll_nue_sc_int
	
	function coll_nue_ann_int(ndim, ve, obj)
		integer, intent(in) :: ndim
		real(dl), dimension(2), intent(in) :: ve
		type(coll_args), intent(in) :: obj
		real(dl) :: coll_nue_ann_int
		integer :: ix
		real(dl), dimension(2) :: pi2_vec
		type(cmplxMatNN) :: n2
		real(dl) :: y2,y3,y4, E3, E4, dme2, fd
		
		y3=ve(1)
		y4=ve(2)
		dme2 = obj%dme2
		E3 = Ebare_i_dme(obj%x,y3, dme2)
		E4 = Ebare_i_dme(obj%x,y4, dme2)
		y2 = E3 + E4 - obj%y1
		if (y2.lt.0.d0 &
			.or. obj%y1.gt.y3+y2+y4 &
			.or. y2.gt.obj%y1+y3+y4 &
			.or. y3.gt.obj%y1+y2+y4 &
			.or. y4.gt.obj%y1+y2+y3) then
			coll_nue_ann_int = 0.0
		else
			n2 = interp_nuDens(y2)
			fd = fermiDirac(y2 / obj%z)
			do ix=1, flavorNumber
				n2%re(ix,ix) = n2%re(ix,ix) * fd
				n2%im(ix,ix) = n2%im(ix,ix) * fd
			end do

			pi2_vec = PI2_nn_f (obj%y1, y2, y3, y4, E3, E4)
			
			coll_nue_ann_int = &
				y3/E3 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab_ann(obj%x,obj%z,obj%n,n2,y3,y4, obj%a,obj%a, obj%ix1,obj%ix2, obj%rI) + &
					pi2_vec(2) * F_ab_ann(obj%x,obj%z,obj%n,n2,y3,y4, obj%b,obj%b, obj%ix1,obj%ix2, obj%rI) + &
					(obj%x*obj%x + dme2) * PI1_12_full(obj%y1, y2, y3, y4) * ( &
						F_ab_ann(obj%x,obj%z,obj%n,n2,y3,y4, 2,1, obj%ix1,obj%ix2, obj%rI) + &
						F_ab_ann(obj%x,obj%z,obj%n,n2,y3,y4, 1,2, obj%ix1,obj%ix2, obj%rI) ) &
				)
		end if
	end function coll_nue_ann_int
	
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
       
	function collision_terms (x, z, y1, n1)
		type(cmplxMatNN) :: collision_terms
		real(dl), intent(in) :: x,z, y1
		type(cmplxMatNN) :: n1, tmp
		type(coll_args), save :: collArgs
		integer :: i,j, k
		real(dl) :: errr1,errr2, res1,res2, cf
		INTEGER :: IFAIL, ITRANS, N, NPTS, NRAND
		real(dl) ::  VK(2)
		real(dl) ::  VKa(3)

		npts=1
		nrand=1
		n=2
		
		call allocateCmplxMat(collision_terms)
		call allocateCmplxMat(collArgs%n)
		call allocateCmplxMat(tmp)
		
		collision_terms%y = y1
		collision_terms%x = x
		collision_terms%z = z
		collision_terms%re = 0.d0
		collision_terms%im = 0.d0
		
!		!$omp parallel &
!		!$omp default(shared) &
!		!$omp private(collArgs,tmp,ifail,itrans,i,j,res1,res2,errr1,errr2)
		collArgs%x = x
		collArgs%z = z
		collArgs%y1 = y1
		collArgs%dme2 = dme2_electron(x, 0.d0, z)
		collArgs%n = n1
		tmp%re = 0.d0
		tmp%im = 0.d0
!		!$omp sections
!		!$omp section
		!scattering of neutrinos vs electrons:
		if (coll_scatt_em) then
			collArgs%a = 2
			collArgs%b = 1
			do i=1, flavorNumber
				collArgs%ix1 = i
				do j=1, flavorNumber
					if (i.eq.j .or. collision_offdiag.eq.1) then
						collArgs%ix2 = j
						collArgs%rI = 1
						ifail=0
						itrans=0
						call D01GCF(n,coll_nue_sc_int, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
						tmp%re(i,j) = tmp%re(i,j) + res1

						collArgs%rI = 2
						ifail=0
						itrans=0
						call D01GCF(n,coll_nue_sc_int, region, npts, vk, nrand,itrans,res2,ERRr2,ifail, collArgs)
						tmp%im(i,j) = tmp%im(i,j) + res2
					else if (i.ne.j .and. collision_offdiag.eq.2) then
						print *,"Warning: should use damping factor...NYI"
					end if
				end do
			end do
			collision_terms%re(:,:) = collision_terms%re(:,:) + tmp%re(:,:)
			collision_terms%im(:,:) = collision_terms%im(:,:) + tmp%im(:,:)
		end if

!		!$omp section
!		!scattering of neutrinos vs positrons:
		if (coll_scatt_ep) then
			collArgs%a = 1
			collArgs%b = 2
!			!$omp do
			do i=1, flavorNumber
				collArgs%ix1 = i
				do j=1, flavorNumber
					if (i.eq.j .or. collision_offdiag.eq.1) then
						collArgs%ix2 = j
						collArgs%rI = 1
						ifail=0
						itrans=0
						call D01GCF(n,coll_nue_sc_int, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
						tmp%re(i,j) = tmp%re(i,j) + res1

						collArgs%rI = 2
						ifail=0
						itrans=0
						call D01GCF(n,coll_nue_sc_int, region, npts, vk, nrand,itrans,res2,ERRr2,ifail, collArgs)
						tmp%im(i,j) = tmp%im(i,j) + res2
					else if (i.ne.j .and. collision_offdiag.eq.2) then
						print *,"Warning: should use damping factor...NYI"
					end if
				end do
			end do
			collision_terms%re(:,:) = collision_terms%re(:,:) + tmp%re(:,:)
			collision_terms%im(:,:) = collision_terms%im(:,:) + tmp%im(:,:)
!			!$omp end do
		end if
		
!		!$omp section
		!annihilation in e+e-
		if (coll_annih_epem) then
			collArgs%a = 1
			collArgs%b = 2
!			!$omp do
			do i=1, flavorNumber
				collArgs%ix1 = i
				do j=1, flavorNumber
					if (i.eq.j .or. collision_offdiag.eq.1) then
						collArgs%ix2 = j
						collArgs%rI = 1
						ifail=0
						itrans=0
						call D01GCF(n,coll_nue_ann_int, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
						tmp%re(i,j) = tmp%re(i,j) + res1

						collArgs%rI = 2
						ifail=0
						itrans=0
						call D01GCF(n,coll_nue_ann_int, region, npts, vk, nrand,itrans,res2,ERRr2,ifail, collArgs)
						tmp%im(i,j) = tmp%im(i,j) + res2
					else if (i.ne.j .and. collision_offdiag.eq.2) then
						print *,"Warning: should use damping factor...NYI"
					end if
				end do
			end do
			collision_terms%re(:,:) = collision_terms%re(:,:) + tmp%re(:,:)
			collision_terms%im(:,:) = collision_terms%im(:,:) + tmp%im(:,:)
!			!$omp end do
		end if
!		!$omp end sections
!		!$omp end parallel
		
		cf = (y1*y1*x**4)
		collision_terms%re(:,:) = collision_terms%re(:,:) * collTermFactor / cf
		collision_terms%im(:,:) = collision_terms%im(:,:) * collTermFactor / cf
		
	end function collision_terms
end module
