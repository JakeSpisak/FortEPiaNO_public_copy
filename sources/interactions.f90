module ndInteractions
	use precision
	use constants
	use variables
	use ndErrors
	use ndMatrices
!	use ndconfig
	implicit none
	
	type coll_args
		type(cmplxMatNN) :: na, nb
		real(dl) :: y1, y2, y3, y4, x, z
		logical :: s1, s2, s3, s4
		integer :: ix1, ix2, a, b
	end type coll_args
	
	contains
	
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
	
	function dme2_electron(x,y,z)
		real(dl) :: dme2_electron, tmp, rombint_obj
		real(dl), intent(in) :: x,y,z
		real(dl), dimension(3) :: vec
		external rombint_obj
		
		vec(1) = x
		vec(2) = z
		vec(3) = y !not used
		
		tmp = rombint_obj(vec, dme2_e_i1, 0.d0, 1d8, 1d-3, 200, 5)
		dme2_electron = 2. * alpha_fine * z*z * (PID3 + tmp/PID2)
	end function dme2_electron
	
	function Ebare_i(x,y,z)!for electrons
		real(dl) :: Ebare_i
		real(dl), intent(in) :: x,y,z
		
		Ebare_i = sqrt(x*x+y*y+dme2_electron(x,y,z)**2)
	end function Ebare_i
	
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
	
	function F_ab_ann(x,z, n1,n2,e3,e4, a,b, i, j)!a, b must be either 1(=L) or 2(=R)
		real(dl), dimension(2) :: F_ab_ann
		type(cmplxMatNN), intent(in) :: n1,n2
		real(dl), intent(in) :: e3,e4, x,z
		integer, intent(in) :: a,b, i, j
		integer :: k
		real(dl), dimension(2) :: term1a, term1b, term2a, term2b
		
		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
			call criticalError("a and b must be either 1(=L) or 2(=R)")
		
		term1a=0.d0
		term1b=0.d0
		term2a=0.d0
		term2b=0.d0
		do k=1, flavorNumber
			!Ga (1-rho2) Gb (1-rho1)
			term1a(1) = term1a(1) + GLR_vec(b,k,k) * ( &
				(idMat(i,k)-n2%re(i,k))*(idMat(k,j)-n1%re(k,j)) - &
				(-n2%im(i,k))*(-n1%im(k,j))&	!idMat ha im=0
				)
			term1a(2) = term1a(2) + GLR_vec(b,k,k) * ( &
				(-n2%im(i,k))*(idMat(k,j)-n1%re(k,j)) - &
				(idMat(i,k)-n2%re(i,k))*(-n1%im(k,j))&
				)
			!(1-rho1) Gb (1-rho2) Ga
			term1b(1) = term1b(1) + GLR_vec(b,k,k) * ( &
				(idMat(i,k)-n1%re(i,k))*(idMat(k,j)-n2%re(k,j)) - &
				(-n1%im(i,k))*(-n2%im(k,j))&
				)
			term1b(2) = term1b(2) + GLR_vec(b,k,k) * ( &
				(-n1%im(i,k))*(idMat(k,j)-n2%re(k,j)) - &
				(idMat(i,k)-n1%re(i,k))*(-n2%im(k,j))&
				)
			!rho1 Gb rho2 Ga
			term2a(1) = term2a(1) + GLR_vec(b,k,k) * ( &
				(n1%re(i,k))*(n2%re(k,j)) - &
				(n1%im(i,k))*(n2%im(k,j))&
				)
			term2a(2) = term2a(2) + GLR_vec(b,k,k) * ( &
				(n1%im(i,k))*(n2%re(k,j)) - &
				(n1%re(i,k))*(n2%im(k,j))&
				)
			!Ga Rho2 Gb rho1
			term2b(1) = term2b(1) + GLR_vec(b,k,k) * ( &
				(n2%re(i,k))*(n1%re(k,j)) - &
				(n2%im(i,k))*(n1%im(k,j))&
				)
			term2b(2) = term2b(2) + GLR_vec(b,k,k) * ( &
				(n2%im(i,k))*(n1%re(k,j)) - &
				(n2%re(i,k))*(n1%im(k,j))&
				)
		end do
		F_ab_ann(:) = fermiDirac(e3,x,z)*fermiDirac(e4,x,z) * &
				(term1a(:) * GLR_vec(a,i,i) + term1b(:) * GLR_vec(a,j,j))&
			- (1-fermiDirac(e3,x,z))*(1-fermiDirac(e4,x,z)) * &
				(term2a(:) * GLR_vec(a,j,j) + term2b(:) * GLR_vec(a,i,i))
	end function F_ab_ann
	
	function F_ab_sc (x,z, n1,e2,n3,e4, a,b, i, j)!a, b must be either 1(=L) or 2(=R)
		real(dl), dimension(2) :: F_ab_sc
		type(cmplxMatNN), intent(in) :: n1,n3
		real(dl), intent(in) :: e2,e4, x,z
		integer, intent(in) :: a,b, i, j
		integer :: k
		real(dl), dimension(2) :: term1a, term1b, term2a, term2b
		
		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
			call criticalError("a and b must be either 1(=L) or 2(=R)")
		
		term1a=0.d0
		term1b=0.d0
		term2a=0.d0
		term2b=0.d0
		do k=1, flavorNumber
			!Ga rho3 Gb (1-rho1)
			term1a(1) = term1a(1) + GLR_vec(b,k,k) * ( &
				(n3%re(i,k))*(idMat(k,j)-n1%re(k,j)) - &
				(n3%im(i,k))*(-n1%im(k,j))&	!idMat ha im=0
				)
			term1a(2) = term1a(2) + GLR_vec(b,k,k) * ( &
				(n3%im(i,k))*(idMat(k,j)-n1%re(k,j)) - &
				(n3%re(i,k))*(-n1%im(k,j))&
				)
			!(1-rho1) Gb rho3 Ga
			term1b(1) = term1b(1) + GLR_vec(b,k,k) * ( &
				(idMat(i,k)-n1%re(i,k))*(n3%re(k,j)) - &
				(-n1%im(i,k))*(n3%im(k,j))&
				)
			term1b(2) = term1b(2) + GLR_vec(b,k,k) * ( &
				(-n1%im(i,k))*(n3%re(k,j)) - &
				(idMat(i,k)-n1%re(i,k))*(n3%im(k,j))&
				)
			!rho1 Gb (1-rho3) Ga
			term2a(1) = term2a(1) + GLR_vec(b,k,k) * ( &
				(n1%re(i,k))*(idMat(k,j)-n3%re(k,j)) - &
				(n1%im(i,k))*(-n3%im(k,j))&
				)
			term2a(2) = term2a(2) + GLR_vec(b,k,k) * ( &
				(n1%im(i,k))*(idMat(k,j)-n3%re(k,j)) - &
				(n1%re(i,k))*(-n3%im(k,j))&
				)
			!Ga (1-rho3) Gb rho1
			term2b(1) = term2b(1) + GLR_vec(b,k,k) * ( &
				(idMat(i,k)-n3%re(i,k))*(n1%re(k,j)) - &
				(-n3%im(i,k))*(n1%im(k,j))&
				)
			term2b(2) = term2b(2) + GLR_vec(b,k,k) * ( &
				(-n3%im(i,k))*(n1%re(k,j)) - &
				(idMat(i,k)-n3%re(i,k))*(n1%im(k,j))&
				)
		end do
		F_ab_sc(:) = (1-fermiDirac(e2,x,z))*fermiDirac(e4,x,z) * &
				(term1a(:) * GLR_vec(a,i,i) + term1b(:) * GLR_vec(a,j,j))&
			- fermiDirac(e2,x,z)*(1-fermiDirac(e4,x,z)) * &
				(term2a(:) * GLR_vec(a,j,j) + term2b(:) * GLR_vec(a,i,i))
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
		real(dl) :: D2_f
		real(dl), intent(in) :: y1, y2, y3, y4
		real(dl) :: &
			p1m2m3m4, m1p2m3m4, m1m2p3m4, m1m2m3p4, &
			p1p2m3m4, m1p2p3m4, m1m2p3p4, p1m2m3p4, &
			p1m2p3m4, m1p2m3p4, m1p2p3p4, p1m2p3p4, &
			p1p2m3p4, p1p2p3m4, p1p2p3p4, m1m2m3m4 
		
		p1m2m3m4 =  y1 -y2 -y3 -y4 
		m1p2m3m4 = -y1 +y2 -y3 -y4 
		m1m2p3m4 = -y1 -y2 +y3 -y4 
		m1m2m3p4 = -y1 -y2 -y3 +y4 
		
		p1p2m3m4 =  y1 +y2 -y3 -y4 
		m1p2p3m4 = -y1 +y2 +y3 -y4 
		m1m2p3p4 = -y1 -y2 +y3 +y4 
		p1m2m3p4 =  y1 -y2 -y3 +y4 
		p1m2p3m4 =  y1 -y2 +y3 -y4 
		m1p2m3p4 = -y1 +y2 -y3 +y4 
		
		m1p2p3p4 = -p1m2m3m4
		p1m2p3p4 = -m1p2m3m4
		p1p2m3p4 = -m1m2p3m4
		p1p2p3m4 = -m1m2m3p4
		
		p1p2p3p4 =  y1 +y2 +y3 +y4 
		m1m2m3m4 = - p1p2p3p4
				
		D2_f = ( &
			Abs(p1p2m3m4)**3 + &
			Abs(p1m2p3m4)**3 - &
			Abs(p1p2p3m4)**3 + &
			Abs(p1m2m3p4)**3 - &
			Abs(p1p2m3p4)**3 - &
			Abs(p1m2p3p4)**3 - &
			Abs(m1p2p3p4)**3 + &
			Abs(p1p2p3p4)**3 + &
			6 * y1 * y2 * ( &
				Abs(p1p2m3m4) - &
				Abs(p1m2p3m4) - &
				Abs(p1p2p3m4) - &
				Abs(p1m2m3p4) - &
				Abs(p1p2m3p4) + &
				Abs(p1m2p3p4) + &
				Abs(m1p2p3p4) + &
				Abs(p1p2p3p4) &
			) - &
			3 * y1 * ( &
				- Sign(p1m2m3m4, (m1p2p3p4)**2) &
				+ Sign(p1p2m3m4, (p1p2m3m4)**2) &
				+ Sign(p1m2p3m4, (p1m2p3m4)**2) &
				- Sign(p1p2p3m4, (p1p2p3m4)**2) &
				+ Sign(p1m2m3p4, (p1m2m3p4)**2) &
				- Sign(p1p2m3p4, (p1p2m3p4)**2) &
				- Sign(p1m2p3p4, (p1m2p3p4)**2) &
				+ Sign(p1p2p3p4, (p1p2p3p4)**2) &
			) &
			- 3 * y2 * ( &
				  Sign(p1m2m3m4, (m1p2p3p4)**2) &
				+ Sign(p1p2m3m4, (p1p2m3m4)**2) &
				- Sign(p1m2p3m4, (p1m2p3m4)**2) &
				- Sign(p1p2p3m4, (p1p2p3m4)**2) &
				- Sign(p1m2m3p4, (p1m2m3p4)**2) &
				- Sign(p1p2m3p4, (p1p2m3p4)**2) &
				+ Sign(p1m2p3p4, (p1m2p3p4)**2) &
				+ Sign(p1p2p3p4, (p1p2p3p4)**2) &
			)) / 6.d0
	end function D2_f
	
	function D3_f(y1, y2, y3, y4)
	!copied from Pablo
		real(dl) :: D3_f
		real(dl), intent(in) :: y1, y2, y3, y4
		
		real(dl) :: a12,a13,a14,a123,a124,a134,a234,a1234,P1234,s12,s13,s14
		real(dl) :: absa12,absa13,absa14,absa12e3,absa13e3,absa14e3
		real(dl) :: ampp,appm,apmp,a123e2,a124e2,a134e2,a234e2,a1234e2
		real(dl) :: a123e3,a124e3,a134e3,a234e3,a1234e3
		real(dl) :: a123e5,a124e5,a134e5,a234e5,a1234e5

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
          
		D3_f = &
			(4.d0*(a1234e5 - a123e5 - a124e5 - a134e5 - a234e5) - &
			absa12**5 - absa13**5 - absa14**5 + &
			(a123 + a1234 + a124 + a134 + a234 + absa12 + absa13 + absa14)*P1234 + &
			5.d0*(a12**3*s12+a13**3*s13 +a14**3*s14) + &
			20.d0*(&
				(-a1234e3 + a123e3 + a124e3 - a134e3 - a234e3 + ampp)*y1*y2 + &
				(-a1234e3 + a123e3 - a124e3 + a134e3 - a234e3 + apmp)*y1*y3 + &
				(-a1234e3 + a123e3 - a124e3 - a134e3 + a234e3 + appm)*y2*y3 + &
				(-a1234e3 - a123e3 + a124e3 + a134e3 - a234e3 + appm)*y1*y4 + &
				(-a1234e3 - a123e3 + a124e3 - a134e3 + a234e3 + apmp)*y2*y4 + &
				(-a1234e3 - a123e3 - a124e3 + a134e3 + a234e3 + ampp)*y3*y4) + &
			60.d0*(&
				(a1234e2 + a123e2 - a124e2 + a134e2 + a234e2 - s12 + s13 - s14)*y1*y2*y4 + &
				(a1234e2 - a123e2 + a124e2 + a134e2 + a234e2 - s12 - s13 + s14)*y1*y2*y3 + &
				(a1234e2 + a123e2 + a124e2 - a134e2 + a234e2 + s12 - s13 - s14)*y1*y3*y4 + &
				(a1234e2 + a123e2 + a124e2 + a134e2 - a234e2 + s12 + s13 + s14)*y2*y3*y4)&
			)/120.d0

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
		type(cmplxMatNN) :: interp_nuDens
		real(dl), intent(in) :: y
	
	end function interp_nuDens
	
	function coll_nue_sc_int(n, ve, obj)
		real(dl), dimension(2) :: coll_nue_sc_int
		real(dl), dimension(2) :: pi1_vec
		real(dl), dimension(3) :: pi2_vec
		real(dl), intent(in), dimension(2) :: ve
		type(cmplxMatNN) :: n3
		real(dl) :: y2,y3,y4
		integer :: n
		type(coll_args) :: obj
		
		y2=ve(1)
		y4=ve(2)
		y3 = obj%y1 + y2 - y4
		n3 = interp_nuDens(y3)
		
		pi1_vec = PI1_f (obj%x, obj%z, obj%y1, y2, y3, y4, obj%s1, obj%s2, obj%s3)
		pi2_vec = PI2_f (obj%x, obj%z, obj%y1, y2, y3, y4, obj%s1, obj%s2, obj%s3, obj%s4)
		
		coll_nue_sc_int = &
			obj%y2/Ebare_i(obj%x, y2, obj%z) * &
			    y4/Ebare_i(obj%x, y4, obj%z) * &
			( &
				pi2_vec(1) * F_ab_sc(obj%x,obj%z,obj%na,y2,n3,y4, obj%a,obj%a, obj%ix1,obj%ix2) + &
				pi2_vec(2) * F_ab_sc(obj%x,obj%z,obj%na,y2,n3,y4, obj%b,obj%b, obj%ix1,obj%ix2) - &
				(obj%x*obj%x + dme2_electron(obj%x,0.d0, obj%z)) * pi1_vec(1)* (&
					F_ab_sc(obj%x,obj%z,obj%na,y2,n3,y4, 2,1, obj%ix1,obj%ix2) + F_ab_sc(obj%x,obj%z,obj%na,y2,n3,y4, 1,2, obj%ix1,obj%ix2)) &
			)
	end function coll_nue_sc_int
	
	function coll_nue_ann_int(n, ve, obj)
		real(dl), dimension(2) :: coll_nue_ann_int
		real(dl), dimension(2) :: pi1_vec
		real(dl), dimension(3) :: pi2_vec
		real(dl), intent(in), dimension(2) :: ve
		type(cmplxMatNN) :: n2
		real(dl) :: y2,y3,y4
		integer :: n
		type(coll_args) :: obj
		
		y3=ve(1)
		y4=ve(2)
		y2 = obj%y1 - y3 - y4
		n2 = interp_nuDens(y2)
		
		pi1_vec = PI1_f (obj%x, obj%z, obj%y1, y2, y3, y4, obj%s1, obj%s2, obj%s3)
		pi2_vec = PI2_f (obj%x, obj%z, obj%y1, y2, y3, y4, obj%s1, obj%s2, obj%s3, obj%s4)
		
		coll_nue_ann_int = &
			obj%y2/Ebare_i(obj%x, y2, obj%z) * &
			    y4/Ebare_i(obj%x, y4, obj%z) * &
			( &
				pi2_vec(1) * F_ab_ann(obj%x,obj%z,obj%na,n2,y3,y4, obj%a,obj%a, obj%ix1,obj%ix2) + &
				pi2_vec(2) * F_ab_ann(obj%x,obj%z,obj%na,n2,y3,y4, obj%b,obj%b, obj%ix1,obj%ix2) - &
				(obj%x*obj%x + dme2_electron(obj%x,0.d0, obj%z)) * pi1_vec(1)* (&
					F_ab_ann(obj%x,obj%z,obj%na,n2,y3,y4, 2,1, obj%ix1,obj%ix2) + F_ab_ann(obj%x,obj%z,obj%na,n2,y3,y4, 1,2, obj%ix1,obj%ix2)) &
			)
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
       
	function collision_terms (x, z, y1, n1, n2, n3)
		real(dl), dimension(:,:), allocatable :: collision_terms
		real(dl), intent(in) :: x,z, y1
		type(cmplxMatNN), intent(in) :: n1, n2, n3
		type(coll_args) :: collArgs
		integer :: i,j
		real(dl) :: rombint_obj
		external rombint_obj
		real(dl) :: ERRr, RES
		INTEGER :: IFAIL, ITRANS, N, NPTS, NRAND
		real(dl) ::  VK(2)
		
		collArgs%s1 = .false.
		collArgs%s2 = .true.
		collArgs%s3 = .false.
		collArgs%s4 = .true.
		collArgs%x = x
		collArgs%z = z
		collArgs%y1 = y1
		
		collArgs%na = n1
		collArgs%nb = n2
		
		!scattering of neutrinos vs electrons:
		collArgs%a = 2
		collArgs%b = 1
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				nrand=1
				ifail=0
				itrans=0
				call D01GCF(2,coll_nue_sc_int, region, 6, vk, nrand,itrans,res,ERRr,ifail)
				collision_terms(i,j) = res
			end do
		end do

		!scattering of neutrinos vs positrons:
		collArgs%a = 1
		collArgs%b = 2
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				nrand=1
				ifail=0
				itrans=0
				call D01GCF(2,coll_nue_sc_int, region, 6, vk, nrand,itrans,res,ERRr,ifail)
				collision_terms(i,j) = collision_terms(i,j) + res
			end do
		end do
		
		!annihilation in e+e-
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				nrand=1
				ifail=0
				itrans=0
				call D01GCF(2,coll_nue_ann_int, region, 6, vk, nrand,itrans,res,ERRr,ifail)
				collision_terms(i,j) = collision_terms(i,j) + res
			end do
		end do
		
		collision_terms(:,:) = collision_terms(:,:) * G_Fsq/(y1*y1*8.d0*PICub)
		
	end function collision_terms
end module
