module ndInteractions
	use precision
	use constants
	use variables
	use utilities
	use ndErrors
	use ndMatrices
!	use ndconfig
	use bspline_module
	implicit none
	
	logical :: dme2_e_loaded = .false.
	type(bspline_2d) :: dmeCorr
	procedure (funcXYZ), pointer :: dme2_electron => null ()
	contains
	
	subroutine dme2_e_load()
		real(dl), dimension(:,:), allocatable :: dme_vec
		integer :: ix, iz, iflag
		real(dl) :: x,z
		
		call addToLog("[interactions] Initializing interpolation for electron mass corrections...")
		allocate(dme_vec(interp_nx,interp_nz))
		do ix=1, interp_nx
			do iz=1, interp_nz
				dme_vec(ix,iz) = dme2_electronFull(interp_xvec(ix),0.d0,interp_zvec(iz))
			end do
		end do
		call dmeCorr%initialize(interp_xvec,interp_zvec,dme_vec,4,4,iflag)
		dme2_electron => dme2_electronInterp
		
		call random_seed()
		call random_number(x)!=0.13151
		call random_number(z)!=0.13151
		x=(x_fin-x_in)*x + x_in
		z=0.4d0*z + z_in
		write(*,"(' [interactions] test dme2_electronInterp in ',*(E12.5))") x,z
		write(*,"(' [interactions] comparison (true vs interp): ',*(E17.10))") dme2_electronFull(x,0.d0,z), dme2_electron(x,0.d0,z)
		
		dme2_e_loaded = .true.
		deallocate(dme_vec)
		call addToLog("[interactions] ...done!")
	end subroutine 
	
	function dme2_e_interp(x,z)
		real(dl) :: dme2_e_interp
		real(dl) :: x,z
		
		dme2_e_interp = 0.
	end function dme2_e_interp
	
	function E_k_m(k,m)
		real(dl), intent(in) :: k,m
		real(dl) :: E_k_m
		E_k_m = sqrt(k*k+m*m)
	end function E_k_m
	
	function fermiDirac(k,m,T)!it's equivalent to fermiDirac(y,x,z)
		real(dl) :: fermiDirac, tmp
		real(dl), intent(in) :: k,m,T
		
		tmp = exp(E_k_m(k,m)/T) + 1.d0
		fermiDirac = 1.d0/tmp
	end function fermiDirac
	
	function fermiDirac_massless(k,T)!it's equivalent to fermiDirac(y,z)
		real(dl) :: fermiDirac_massless, tmp
		real(dl), intent(in) :: k,T
		
		fermiDirac_massless = 1.d0/(exp(k/T) + 1.d0)
	end function fermiDirac_massless
	
	function dme2_e_i1(vec, k)
		real(dl) :: dme2_e_i1
		real(dl), intent(in) :: k
		real(dl), dimension(3) :: vec
		real(dl) :: m,T,p
		m=vec(1)
		T=vec(2)
		dme2_e_i1 = k*k/E_k_m(k,m)*fermiDirac(k,m,T)
	end function dme2_e_i1
	
	function dme2_e_i2(vec, k) !not used
		real(dl) :: dme2_e_i2
		real(dl), intent(in) :: k
		real(dl), dimension(3) :: vec
		real(dl) :: m,T,p
		m=vec(1)
		T=vec(2)
		p=vec(3)
		dme2_e_i2 = k/E_k_m(k,m)*fermiDirac(k,m,T)*log(abs((p+k)/(p-k)))
	end function dme2_e_i2
	
	function dme2_electronFull(x,y,z)
		real(dl) :: dme2_electronFull, tmp, rombint_obj
		real(dl), intent(in) :: x,y,z
		real(dl), dimension(3) :: vec
		external rombint_obj
		
!		if (dme2_e_loaded) then
!			dme2_electronFull = dme2_e_interp(x,z)
!		else
			vec(1) = x
			vec(2) = z
			vec(3) = y !not used
		
			tmp = rombint_obj(vec, dme2_e_i1, 0.d0, 60.d0, 1d-3, 200)
			dme2_electronFull = 2. * alpha_fine * z*z * (PID3 + tmp/PID2)
!		end if
	end function dme2_electronFull
	
	function dme2_electronInterp(x,y,z)
		real(dl) :: dme2_electronInterp
		real(dl), intent(in) :: x,y,z
		integer :: iflag
		
		call dmeCorr%evaluate(x,z,0,0,dme2_electronInterp,iflag)
	end function dme2_electronInterp
	
	function Ebare_i(x,y,z)!for electrons
		real(dl) :: Ebare_i
		real(dl), intent(in) :: x,y,z
		
		Ebare_i = sqrt(x*x+y*y+dme2_electron(x,y,z))
	end function Ebare_i
	
	function Ebare_i_dme(x,y,z,dme2)!for electrons
		real(dl) :: Ebare_i_dme
		real(dl), intent(in) :: x,y,z,dme2
		
		Ebare_i_dme = sqrt(x*x+y*y+dme2)
	end function Ebare_i_dme
	
	function Ebarn_i(x,y,z)!for neutrinos
		real(dl) :: Ebarn_i
		real(dl), intent(in) :: x,y,z
		
		Ebarn_i = y
	end function Ebarn_i
	
	function Ebar_i(s,x,y,z)!s: True for e+/e-, False for nu
		real(dl) :: Ebar_i
		real(dl), intent(in) :: x,y,z
		logical, intent(in) :: s !True for e+/e-, False for nu
		
		if (s) then
			Ebar_i = Ebare_i(x,y,z)
		else
			Ebar_i = Ebarn_i(x,y,z)
		end if
	end function Ebar_i
	
	function F_ab_ann(x,z, n1,n2,e3,e4, a,b, i, j, rI)!a, b must be either 1(=L) or 2(=R)
		real(dl) :: F_ab_ann
		type(cmplxMatNN), intent(in) :: n1,n2
		real(dl), intent(in) :: e3,e4, x,z
		integer, intent(in) :: a,b, i, j, rI
		integer :: k
		real(dl) :: term1a, term1b, term2a, term2b
		
		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
			call criticalError("a and b must be either 1(=L) or 2(=R)")
		
		term1a=0.d0
		term1b=0.d0
		term2a=0.d0
		term2b=0.d0
		if (rI.eq.1) then !take real part
			do k=1, flavorNumber
				!Ga (1-rho2) Gb (1-rho1)
				term1a = term1a + GLR_vec(b,k,k) * ( &
					(idMat(i,k)-n2%re(i,k))*(idMat(k,j)-n1%re(k,j)) - &
					(-n2%im(i,k))*(-n1%im(k,j))&	!idMat ha im=0
					)
				!(1-rho1) Gb (1-rho2) Ga
				term1b = term1b + GLR_vec(b,k,k) * ( &
					(idMat(i,k)-n1%re(i,k))*(idMat(k,j)-n2%re(k,j)) - &
					(-n1%im(i,k))*(-n2%im(k,j))&
					)
				!rho1 Gb rho2 Ga
				term2a = term2a + GLR_vec(b,k,k) * ( &
					(n1%re(i,k))*(n2%re(k,j)) - &
					(n1%im(i,k))*(n2%im(k,j))&
					)
				!Ga Rho2 Gb rho1
				term2b = term2b + GLR_vec(b,k,k) * ( &
					(n2%re(i,k))*(n1%re(k,j)) - &
					(n2%im(i,k))*(n1%im(k,j))&
					)
			end do
		else if (rI .eq.2) then!take imaginary part
			do k=1, flavorNumber
				!Ga (1-rho2) Gb (1-rho1)
				term1a = term1a + GLR_vec(b,k,k) * ( &
					(-n2%im(i,k))*(idMat(k,j)-n1%re(k,j)) - &
					(idMat(i,k)-n2%re(i,k))*(-n1%im(k,j))&
					)
				!(1-rho1) Gb (1-rho2) Ga
				term1b = term1b + GLR_vec(b,k,k) * ( &
					(-n1%im(i,k))*(idMat(k,j)-n2%re(k,j)) - &
					(idMat(i,k)-n1%re(i,k))*(-n2%im(k,j))&
					)
				!rho1 Gb rho2 Ga
				term2a = term2a + GLR_vec(b,k,k) * ( &
					(n1%im(i,k))*(n2%re(k,j)) - &
					(n1%re(i,k))*(n2%im(k,j))&
					)
				!Ga Rho2 Gb rho1
				term2b = term2b + GLR_vec(b,k,k) * ( &
					(n2%im(i,k))*(n1%re(k,j)) - &
					(n2%re(i,k))*(n1%im(k,j))&
					)
			end do
		end if
		F_ab_ann = fermiDirac(e3,x,z)*fermiDirac(e4,x,z) * &
				(term1a * GLR_vec(a,i,i) + term1b * GLR_vec(a,j,j))&
			- (1-fermiDirac(e3,x,z))*(1-fermiDirac(e4,x,z)) * &
				(term2a * GLR_vec(a,j,j) + term2b * GLR_vec(a,i,i))
	end function F_ab_ann
	
	function F_ab_sc (x,z, n1,e2,n3,e4, a,b, i, j, rI)!a, b must be either 1(=L) or 2(=R)
		real(dl) :: F_ab_sc
		type(cmplxMatNN), intent(in) :: n1,n3
		real(dl), intent(in) :: e2,e4, x,z
		integer, intent(in) :: a,b, i, j, rI
		integer :: k
		real(dl) :: term1a, term1b, term2a, term2b
		
		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
			call criticalError("a and b must be either 1(=L) or 2(=R)")
		
		term1a=0.d0
		term1b=0.d0
		term2a=0.d0
		term2b=0.d0
		if (rI.eq.1) then !take real part
			do k=1, flavorNumber
				!Ga rho3 Gb (1-rho1)
				term1a = term1a + GLR_vec(b,k,k) * ( &
					(n3%re(i,k))*(idMat(k,j)-n1%re(k,j)) - &
					(n3%im(i,k))*(-n1%im(k,j))&	!idMat ha im=0
					)
				!(1-rho1) Gb rho3 Ga
				term1b = term1b + GLR_vec(b,k,k) * ( &
					(idMat(i,k)-n1%re(i,k))*(n3%re(k,j)) - &
					(-n1%im(i,k))*(n3%im(k,j))&
					)
				!rho1 Gb (1-rho3) Ga
				term2a = term2a + GLR_vec(b,k,k) * ( &
					(n1%re(i,k))*(idMat(k,j)-n3%re(k,j)) - &
					(n1%im(i,k))*(-n3%im(k,j))&
					)
				!Ga (1-rho3) Gb rho1
				term2b = term2b + GLR_vec(b,k,k) * ( &
					(idMat(i,k)-n3%re(i,k))*(n1%re(k,j)) - &
					(-n3%im(i,k))*(n1%im(k,j))&
					)
			end do
		else if (rI .eq.2) then!take imaginary part
			do k=1, flavorNumber
				!Ga rho3 Gb (1-rho1)
				term1a = term1a + GLR_vec(b,k,k) * ( &
					(n3%im(i,k))*(idMat(k,j)-n1%re(k,j)) - &
					(n3%re(i,k))*(-n1%im(k,j))&
					)
				!(1-rho1) Gb rho3 Ga
				term1b = term1b + GLR_vec(b,k,k) * ( &
					(-n1%im(i,k))*(n3%re(k,j)) - &
					(idMat(i,k)-n1%re(i,k))*(n3%im(k,j))&
					)
				!rho1 Gb (1-rho3) Ga
				term2a = term2a + GLR_vec(b,k,k) * ( &
					(n1%im(i,k))*(idMat(k,j)-n3%re(k,j)) - &
					(n1%re(i,k))*(-n3%im(k,j))&
					)
				!Ga (1-rho3) Gb rho1
				term2b = term2b + GLR_vec(b,k,k) * ( &
					(-n3%im(i,k))*(n1%re(k,j)) - &
					(idMat(i,k)-n3%re(i,k))*(n1%im(k,j))&
					)
			end do
		end if
		F_ab_sc = (1-fermiDirac(e2,x,z))*fermiDirac(e4,x,z) * &
				(term1a * GLR_vec(a,i,i) + term1b * GLR_vec(a,j,j))&
			- fermiDirac(e2,x,z)*(1-fermiDirac(e4,x,z)) * &
				(term2a * GLR_vec(a,j,j) + term2b * GLR_vec(a,i,i))
	end function F_ab_sc
	
	function D1_f(y1, y2, y3, y4)
		real(dl) :: D1_f
		real(dl), intent(in) :: y1, y2, y3, y4
		
		D1_f = &
			abs( y1+y2+y3-y4) + &
			abs( y1+y2-y3+y4) + &
			abs( y1-y2+y3+y4) + &
			abs(-y1+y2+y3+y4) - &
			abs( y1+y2-y3-y4) - &
			abs( y1-y2-y3+y4) - &
			abs( y1-y2+y3-y4) - &
			abs( y1+y2+y3+y4)
	end function D1_f
	
	function D2_f(y1, y2, y3, y4)
		implicit none
		real(dl) :: D2_f
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

		a = Sign(p1p2m3m4, (p1p2m3m4)**2) - Sign(p1p2p3m4, (p1p2p3m4)**2) - &
			Sign(p1p2m3p4, (p1p2m3p4)**2) + Sign(p1p2p3p4, (p1p2p3p4)**2)
		b = - Sign(m1p2p3p4, (m1p2p3p4)**2) + Sign(p1m2p3m4, (p1m2p3m4)**2) &
			+ Sign(p1m2m3p4, (p1m2m3p4)**2) - Sign(p1m2p3p4, (p1m2p3p4)**2)

		D2_f = - ( &
			Ap1p2m3m4**3 + Ap1m2p3m4**3 - &
			Ap1p2p3m4**3 + Ap1m2m3p4**3 - &
			Ap1p2m3p4**3 - Ap1m2p3p4**3 - &
			Am1p2p3p4**3 + Ap1p2p3p4**3 + &
			6 * y1 * y2 * ( &
				Ap1p2m3m4 - Ap1m2p3m4 - &
				Ap1p2p3m4 - Ap1m2m3p4 - &
				Ap1p2m3p4 + Ap1m2p3p4 + &
				Am1p2p3p4 + Ap1p2p3p4 &
			) - &
			3 * y1 * ( a + b ) &
			- 3 * y2 * ( a - b )) / 6.d0
	end function D2_f
	
	function D3_f(y1, y2, y3, y4)
	!from 1605.09383
		real(dl) :: D3_f
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

		D3_f = &
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
				ap1m2m3p4**3 * (y1 + (p23-y4) - y23 + y4*p23) + &
				ap1p2p3m4**3 * (y12 + y23 + y13 - y4*p123) + &
				ap1m2p3m4**3 * (y2*(m34) + y34 + y1*(y24-y3)) + &
				ap1p2m3p4**3 * (-y2*(m34) - y34 + y1*(y24-y3)) + &
				ap1m2p3p4**3 * (y34 - y2*p34 + y1*(p34-y2)) + &
				ap1p2m3m4**3 * (y2*p34 - y34 + y1*(p34-y2)) &
			)/6.d0 + &
			(	sign(p1p2m3m4,p1p2m3m4**2) * ( &
					y1**3+y2**3+3*y1*(y2**2+y3**2+y4**2+6*(y34-p34*y2)) + &
					p34**3 + 3*(y2*(y3**2+y4**2+6*y34)-y2**2*p34+y1**2*(y2-p34))&
				) + &
				sign(p1p2m3p4,p1p2m3p4**2) * ( &
					- p1p2m3p4**3 + 12*(y12*m34+p12*y34) &
				) + &
				sign(p1m2p3p4,p1m2p3p4**2) * ( &
					- p1m2p3p4**3 + 12*(y12*p34+(y2-y1)*y34) &
				) + &
				sign(p1p2p3m4,p1p2p3m4**2) * ( &
					3*((y2**2+6*y23+y3**2)*y4 - (p23-y4)*y1**2 - p23*y4**2 - &
						(y2**2+6*y23+y3**2+y4**2-6*p23*y4)*y1) +&
					y4**3-y1**3-p23**3 &
				) + &
				sign(p1m2m3p4,p1m2m3p4**2) * ( &
					3*((y2**2+6*y23+y3**2)*y4 - (p23-y4)*y1**2 - p23*y4**2 + &
						(y2**2+6*y23+y3**2+y4**2-6*p23*y4)*y1) +&
					y1**3+y4**3-p23**3 &
				) + &
				sign(p1m2p3m4,p1m2p3m4**2) * ( &
					3*(-(y3**2-6*y34+y4**2)*y4 - (p23-y3)*y1**2 + m34*y2**2 + &
						(y2**2-6*(m34)*y2+y3**2+y4**2-6*y23)*y1) +&
					y1**3-y2**3-(m34)**3 &
				) + &
				sign(-m1p2p3p4,m1p2p3p4**2) * ( &
					3*((y3**2+6*y34+y4**2)*y2 + (p234)*y1**2 + p34*y2**2 - &
						(y2**2+6*p34*y2+y3**2+y4**2+6*y34)*y1) +&
					y2**3-y1**3+p34**3 &
				) &
			)/24.d0

	end function D3_f
	
	function PI1_f (x, z, y1, y2, y3, y4, s1, s2, s3) !1: (y1,y3), 	2: (y1, y2)
		real(dl), dimension(2) :: PI1_f
		real(dl), intent(in) :: y1, y2, y3, y4, x, z
		logical, intent(in) :: s1, s2, s3
		real(dl) :: e1, d1
		
		e1 = Ebar_i(s1,x,y1,z)
		d1 = D1_f(y1, y2, y3, y4)
		
		!pi_1(y1, y3)
		PI1_f(1) = e1 * Ebar_i(s3,x,y3,z) * d1 &
			+ D2_f(y1, y3, y2, y4)
			
		!pi_1(y1, y2)
		PI1_f(2) = e1 * Ebar_i(s2,x,y2,z) * d1 &
			- D2_f(y1, y2, y3, y4)
	end function PI1_f
	
	function PI2_f (x, z, y1, y2, y3, y4, s1, s2, s3, s4) !1: (y1, y4),	2: (y1, y2), 3:	(y1, y3)
		real(dl), dimension(3) :: PI2_f
		real(dl), intent(in) :: y1, y2, y3, y4, x, z
		logical, intent(in) :: s1, s2, s3, s4
		real(dl) :: e1, e2, e3, e4, t
		
		e1 = Ebar_i(s1,x,y1,z)
		e2 = Ebar_i(s2,x,y2,z)
		e3 = Ebar_i(s3,x,y3,z)
		e4 = Ebar_i(s4,x,y4,z)
		t = e1 * e2 * e3 * e4 * D1_f(y1, y2, y3, y4) + D3_f(y1, y2, y3, y4)
		
		!pi_2(y1, y4)
		PI2_f(1) = 2 * (t &
			+ e2 * e3 * D2_f(y1, y4, y2, y3) &
			+ e1 * e4 * D2_f(y2, y3, y1, y4) )
		
		!pi_2(y1, y2)
		PI2_f(2) = 2 * (t &
			- e1 * e2 * D2_f(y3, y4, y1, y2) &
			- e3 * e4 * D2_f(y1, y2, y3, y4) )
		
		!pi_1(y1, y3)
		PI2_f(3) = 2 * (t &
			+ e1 * e3 * D2_f(y2, y4, y1, y3) &
			+ e2 * e4 * D2_f(y1, y3, y2, y4) )
	end function PI2_f

	function interp_nuDens(y)
		type(cmplxMatNN) :: interp_nuDens, mu, ml
		real(dl), intent(in) :: y
		integer :: ix,i,j,k
		real(dl) :: ixr, yr
		logical :: exact
		character(len=100) :: vals
		
		call allocateCmplxMat(interp_nuDens)
		call allocateCmplxMat(mu)
		call allocateCmplxMat(ml)
		if (y.gt.nuDensMatVec(Ny)%y .or. y .lt. nuDensMatVec(1)%y) then
			interp_nuDens%re =0.d0
			interp_nuDens%im =0.d0
		elseif (abs(y-nuDensMatVec(Ny)%y) .lt. 1.d-6) then
			interp_nuDens = nuDensMatVec(Ny)
		else
			ixr=(Ny-1)*(log10(y)-logy_min)/(logy_max-logy_min)+1
			ix=int(ixr)
			exact = abs(ixr-ix) .lt. 1d-5
			if (ix .gt. 0 .and. ix .lt. Ny) then
				if (exact) then
					interp_nuDens = nuDensMatVec(ix)
				else
					ml = nuDensMatVec(ix)
					mu = nuDensMatVec(ix+1)
					yr = (y - nuDensMatVec(ix)%y) / (nuDensMatVec(ix+1)%y - nuDensMatVec(ix)%y)
					interp_nuDens%y = y
					interp_nuDens%re(:,:) = ml%re(:,:) + yr * (mu%re(:,:) - ml%re(:,:))
					interp_nuDens%im(:,:) = ml%im(:,:) + yr * (mu%im(:,:) - ml%im(:,:))
				end if
			else
				write(vals,"('y=',E14.7,'    ix=',I3)") y, ix
				call Error("k out of range or errors occurred in interpolation"//NEW_LINE('A')//vals)
			end if
		end if
	end function interp_nuDens
	
	function interp_nuDensIJ(y, i, j)
		type(cmplxMatNN) :: interp_nuDensIJ, mu, ml
		real(dl), intent(in) :: y
		integer :: ix,i,j,k
		real(dl) :: ixr, yr
		logical :: exact
		character(len=100) :: vals
		
		call allocateCmplxMat(interp_nuDensIJ)
		call allocateCmplxMat(mu)
		call allocateCmplxMat(ml)
		if (y.gt.nuDensMatVec(Ny)%y .or. y .lt. nuDensMatVec(1)%y) then
			interp_nuDensIJ%re =0.d0
			interp_nuDensIJ%im =0.d0
		elseif (abs(y - nuDensMatVec(Ny)%y).lt.1d-6) then
			interp_nuDensIJ = nuDensMatVec(Ny)
		else
			ixr=(Ny-1)*(log10(y)-logy_min)/(logy_max-logy_min)+1
			ix=int(ixr)
			exact = abs(ixr-ix) .lt. 1d-6
			if (ix .gt. 0 .and. ix .lt. Ny) then
				if (exact) then
					interp_nuDensIJ = nuDensMatVec(ix)
				else
					ml = nuDensMatVec(ix)
					mu = nuDensMatVec(ix+1)
					yr = (y - nuDensMatVec(ix)%y) / (nuDensMatVec(ix+1)%y - nuDensMatVec(ix)%y)
					interp_nuDensIJ%y = y
					interp_nuDensIJ%re(i,j) = ml%re(i,j) + yr * (mu%re(i,j) - ml%re(i,j))
					interp_nuDensIJ%im(i,j) = ml%im(i,j) + yr * (mu%im(i,j) - ml%im(i,j))
				end if
			else
				write(vals,"('y=',E14.7,'    ix=',I3)") y, ix
				call Error("k out of range or errors occurred in interpolation"//NEW_LINE('A')//vals)
			end if
		end if
	end function interp_nuDensIJ
	
	function coll_nue_sc_int(ndim, ve, obj)
		integer :: ndim
		type(coll_args) :: obj
		real(dl), dimension(20) :: ve
		real(dl) :: coll_nue_sc_int
		real(dl), dimension(2) :: pi1_vec
		real(dl), dimension(3) :: pi2_vec
		type(cmplxMatNN) :: n3
		real(dl) :: y2,y3,y4, E2, E4, dme2
		
		call allocateCmplxMat(n3)
		
		y2 = ve(1)
		y4 = ve(2)
		dme2 = obj%dme2
		E2 = Ebare_i_dme(obj%x,y2,obj%z, dme2)
		E4 = Ebare_i_dme(obj%x,y4,obj%z, dme2)
		y3 = obj%y1 + E2 - E4
		if (y3.lt.0.d0 &
			.or. obj%y1.gt.y3+y2+y4 &
			.or. y2.gt.obj%y1+y3+y4 &
			.or. y3.gt.obj%y1+y2+y4 &
			.or. y4.gt.obj%y1+y2+y3) then
			coll_nue_sc_int = 0.0
		else
			n3 = interp_nuDens(y3)
			
			pi1_vec = PI1_f (obj%x, obj%z, obj%y1, y2, y3, y4, obj%s1, obj%s2, obj%s3)
			pi2_vec = PI2_f (obj%x, obj%z, obj%y1, y2, y3, y4, obj%s1, obj%s2, obj%s3, obj%s4)
			
			coll_nue_sc_int = &
				y2/E2 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab_sc(obj%x,obj%z,obj%na,y2,n3,y4, obj%a,obj%a, obj%ix1,obj%ix2, obj%rI) + &
					pi2_vec(2) * F_ab_sc(obj%x,obj%z,obj%na,y2,n3,y4, obj%b,obj%b, obj%ix1,obj%ix2, obj%rI) - &
					(obj%x*obj%x + dme2) * pi1_vec(1)* (&
						F_ab_sc(obj%x,obj%z,obj%na,y2,n3,y4, 2,1, obj%ix1,obj%ix2, obj%rI) + &
						F_ab_sc(obj%x,obj%z,obj%na,y2,n3,y4, 1,2, obj%ix1,obj%ix2, obj%rI)) &
				)
		end if
	end function coll_nue_sc_int
	
	function coll_nue_ann_int(n, ve, obj)
		real(dl) :: coll_nue_ann_int
		real(dl), dimension(2) :: pi1_vec
		real(dl), dimension(3) :: pi2_vec
		real(dl), intent(in), dimension(2) :: ve
		type(cmplxMatNN) :: n2
		real(dl) :: y2,y3,y4, E3, E4, dme2
		integer :: n
		type(coll_args) :: obj
		
		call allocateCmplxMat(n2)
		
		y3=ve(1)
		y4=ve(2)
		dme2 = obj%dme2
		E3 = Ebare_i_dme(obj%x,y4,obj%z,dme2)
		E4 = Ebare_i_dme(obj%x,y4,obj%z,dme2)
		y2 = E3 + E4 - obj%y1
		if (y2.lt.0.d0 &
			.or. obj%y1.gt.y3+y2+y4 &
			.or. y2.gt.obj%y1+y3+y4 &
			.or. y3.gt.obj%y1+y2+y4 &
			.or. y4.gt.obj%y1+y2+y3) then
			coll_nue_ann_int = 0.0
		else
			n2 = interp_nuDens(y2)
			
			pi1_vec = PI1_f (obj%x, obj%z, obj%y1, y2, y3, y4, obj%s1, obj%s2, obj%s3)
			pi2_vec = PI2_f (obj%x, obj%z, obj%y1, y2, y3, y4, obj%s1, obj%s2, obj%s3, obj%s4)
			
			coll_nue_ann_int = &
				y3/E3 * &
				y4/E4 * &
				( &
					pi2_vec(1) * F_ab_ann(obj%x,obj%z,obj%na,n2,y3,y4, obj%a,obj%a, obj%ix1,obj%ix2, obj%rI) + &
					pi2_vec(2) * F_ab_ann(obj%x,obj%z,obj%na,n2,y3,y4, obj%b,obj%b, obj%ix1,obj%ix2, obj%rI) + &
					(obj%x*obj%x + dme2) * pi1_vec(1)* (&
						F_ab_ann(obj%x,obj%z,obj%na,n2,y3,y4, 2,1, obj%ix1,obj%ix2, obj%rI) + &
						F_ab_ann(obj%x,obj%z,obj%na,n2,y3,y4, 1,2, obj%ix1,obj%ix2, obj%rI)) &
				)
		end if
	end function coll_nue_ann_int
	
	SUBROUTINE region(ndim,x,j,c,d)
		use precision
		IMPLICIT NONE
		REAL (dl), INTENT (OUT) :: c, d
		INTEGER, INTENT (IN)    :: j, ndim
		REAL (dl), INTENT (IN)  :: x(ndim)
		c = y_min
		d = y_max
		RETURN
	END SUBROUTINE region
       
	function collision_terms (x, z, y1)
		type(cmplxMatNN) :: collision_terms
		real(dl), intent(in) :: x,z, y1
		type(cmplxMatNN) :: n1
		type(coll_args) :: collArgs
		integer :: i,j, k
		real(dl) :: rombint_obj
		external rombint_obj
		real(dl) :: ERRr, RES
		INTEGER :: IFAIL, ITRANS, N, NPTS, NRAND
		real(dl) ::  VK(2)

		npts=5
		n=2
		
		call allocateCmplxMat(n1)
		call allocateCmplxMat(collision_terms)
		call allocateCmplxMat(collArgs%na)
		call allocateCmplxMat(collArgs%nb)
		
		collArgs%s1 = .false.
		collArgs%s2 = .true.
		collArgs%s3 = .false.
		collArgs%s4 = .true.
		collArgs%x = x
		collArgs%z = z
		collArgs%y1 = y1
		collArgs%dme2 = dme2_electron(x, 0.d0, z)
		do k=1, Ny
			if (nuDensMatVec(k)%y .eq. y1) &
				collArgs%na = nuDensMatVec(k)
		end do
		collision_terms%y = y1
		collision_terms%x = x
		collision_terms%z = z
		
		!scattering of neutrinos vs electrons:
		collArgs%a = 2
		collArgs%b = 1
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				collArgs%rI = 1
				nrand=1
				ifail=0
				itrans=0
				call D01GCF(n,coll_nue_sc_int, region, npts, vk, nrand,itrans,res,ERRr,ifail, collArgs)
				collision_terms%re(i,j) = res
				
				collArgs%rI = 2
				nrand=1
				ifail=0
				itrans=0
				call D01GCF(n,coll_nue_sc_int, region, npts, vk, nrand,itrans,res,ERRr,ifail, collArgs)
				collision_terms%im(i,j) = res
			end do
		end do

		!scattering of neutrinos vs positrons:
		collArgs%a = 1
		collArgs%b = 2
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				collArgs%rI = 1
				nrand=1
				ifail=0
				itrans=0
				call D01GCF(n,coll_nue_sc_int, region, npts, vk, nrand,itrans,res,ERRr,ifail, collArgs)
				collision_terms%re(i,j) = collision_terms%re(i,j) + res
				
				collArgs%rI = 2
				nrand=1
				ifail=0
				itrans=0
				call D01GCF(n,coll_nue_sc_int, region, npts, vk, nrand,itrans,res,ERRr,ifail, collArgs)
				collision_terms%im(i,j) = collision_terms%im(i,j) + res
			end do
		end do
		
		collArgs%s2 = .false.
		collArgs%s3 = .true.
		!annihilation in e+e-
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				collArgs%rI = 1
				nrand=1
				ifail=0
				itrans=0
				call D01GCF(n,coll_nue_ann_int, region, npts, vk, nrand,itrans,res,ERRr,ifail, collArgs)
				collision_terms%re(i,j) = collision_terms%re(i,j) + res
				
				collArgs%rI = 2
				nrand=1
				ifail=0
				itrans=0
				call D01GCF(n,coll_nue_ann_int, region, npts, vk, nrand,itrans,res,ERRr,ifail, collArgs)
				collision_terms%im(i,j) = collision_terms%im(i,j) + res
			end do
		end do
		
		collision_terms%re(:,:) = collision_terms%re(:,:) * G_Fsq/(y1*y1*8.d0*PICub) * (m_e**3) / (x**4)
		collision_terms%im(:,:) = collision_terms%im(:,:) * G_Fsq/(y1*y1*8.d0*PICub) * (m_e**3) / (x**4)
		
	end function collision_terms
end module
