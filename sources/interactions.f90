module ndInteractions
	use precision
	use constants
	use ndErrors
	use ndconfig
	implicit none
	
	
	contains
	
	function E_k_m(k,m)
		real(dl), intent(in) :: k,m
		real(dl) :: E_k_m
		E_k_m = sqrt(k*k+m*m)
	end function E_k_m
	
	function fermiDirac(k,m,T)
		real(dl) :: fermiDirac, tmp
		real(dl), intent(in) :: k,m,T
		
		tmp = exp(E_k_m(k,m)/T) + 1.d0
		fermiDirac = 1.d0/tmp
	end function fermiDirac
	
	function dme2_e_i1(vec, k)
		real(dl) :: dme2_e_i1
		real(dl), intent(in) :: k
		real(dl), dimension(3) :: vec
		real(dl) :: m,T,p
		m=vec(1)
		T=vec(2)
		dme2_e_i1 = k*k/E_k_m(k,m)*fermiDirac(k,m,T)
	end function dme2_e_i1
	
!	function dme2_e_i2(vec, k)
!		real(dl) :: dme2_e_i2
!		real(dl), intent(in) :: k
!		real(dl), dimension(3) :: vec
!		real(dl) :: m,T,p
!		m=vec(1)
!		T=vec(2)
!		p=vec(3)
!		dme2_e_i2 = k/E_k_m(k,m)*fermiDirac(k,m,T)*log(abs((p+k)/(p-k)))
!	end function dme2_e_i2
	
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
	
	function Ebar_i(s,x,y,z)
		real(dl) :: Ebar_i
		real(dl), intent(in) :: x,y,z
		logical, intent(in) :: s !True for e+/e-, False for nu
		
		if (s) then
			Ebar_i = sqrt(x*x+y*y+dme2_electron(x,y,z)**2)
		else
			Ebar_i = y
		end if
	
	end function Ebar_i
	
	function F_ab_ann(n1,n2,e3,e4,a,b)
		real(dl), dimension(:,:), allocatable :: F_ab_ann
		real(dl), dimension(:,:), intent(in) :: n1,n2
		real(dl), intent(in) :: e3,e4
		integer, intent(in) :: a,b
		integer :: nf
		real(dl), dimension(:,:), allocatable :: tm1, tm2, tm3, tm4, tm5, tm6, tm7
		
		nf = flavorNumber
		allocate(tm1(nf,nf), tm2(nf,nf), tm3(nf,nf), tm4(nf,nf), tm5(nf,nf), tm6(nf,nf), tm7(nf,nf))
		allocate( F_ab_ann(nf,nf) )
		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
			call criticalError("a and b must be either 1(=L) or 2(=R)")
		
		tm4 = GLR_vec(a,:,:)
		tm3 = GLR_vec(b,:,:)
		tm2(:,:) = 1.d0 - n2(:,:)
		tm1(:,:) = 1.d0 - n1(:,:)
		call quadrupleProdMat(tm4, tm2, tm3, tm1, tm5)
		call quadrupleProdMat(tm1, tm3, tm2, tm4, tm6)
		tm7=tm5+tm6

		call quadrupleProdMat(tm4, n2, tm3, n1, tm5)
		call quadrupleProdMat(n1, tm3, n2, tm4, tm6)
		tm1=tm5+tm6
		
		F_ab_ann = e3*e4 * tm7 - (1-e3)*(1-e4)*tm1
	end function F_ab_ann
	
	function F_ab_sc(n1,e2,n3,e4,a,b)
		real(dl), dimension(:,:), allocatable :: F_ab_sc
		real(dl), dimension(:,:), intent(in) :: n1,n3
		real(dl), intent(in) :: e2,e4
		integer, intent(in) :: a,b
		integer :: nf
		real(dl), dimension(:,:), allocatable :: tm1, tm2, tm3, tm4, tm5, tm6, tm7
		
		nf = flavorNumber
		allocate(tm1(nf,nf), tm2(nf,nf), tm3(nf,nf), tm4(nf,nf), tm5(nf,nf), tm6(nf,nf), tm7(nf,nf))
		allocate( F_ab_sc(nf,nf) )
		if ((a.ne.1 .and. a.ne.2) .or. (b.ne.1 .and. b.ne.2)) &
			call criticalError("a and b must be either 1(=L) or 2(=R)")
		
		tm4 = GLR_vec(a,:,:)
		tm3(:,:) = 1.d0 - n3(:,:)
		tm2 = GLR_vec(b,:,:)
		tm1(:,:) = 1.d0 - n1(:,:)
		call quadrupleProdMat(tm4, n3, tm2, tm1, tm5)
		call quadrupleProdMat(tm1, tm2, n3, tm4, tm6)
		tm7=tm5+tm6

		call quadrupleProdMat(n1, tm2, tm3, tm4, tm5)
		call quadrupleProdMat(tm4, tm3, tm2, n1, tm6)
		tm1=tm5+tm6
		
		F_ab_sc = (1-e2)*e4 * tm7 - e2*(1-e4)*tm1
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
				
		D2_f = 0. + &
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
			)
	end function D2_f
	
	function D3_f(y1, y2, y3, y4)
		real(dl) :: D3_f
		real(dl), intent(in) :: y1, y2, y3, y4
		
		D3_f = 0.
	end function D3_f
	
	function PI1_f (y1, y2, y3, y4)
		real(dl) :: PI1_f
		real(dl), intent(in) :: y1, y2, y3, y4
		
		PI1_f=0
	end function PI1_f

end module
