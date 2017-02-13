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
	
	function dme2_e_i2(vec, k)
		real(dl) :: dme2_e_i2
		real(dl), intent(in) :: k
		real(dl), dimension(3) :: vec
		real(dl) :: m,T,p
		m=vec(1)
		T=vec(2)
		p=vec(3)
		dme2_e_i2 = k/E_k_m(k,m)*fermiDirac(k,m,T)*log(abs((p+k)/(p-k)))
	end function dme2_e_i2
	
	function dme2_electron(p,T)
		real(dl) :: dme2_electron, tmp, rombint_obj
		real(dl), intent(in) :: p,T
		real(dl), dimension(3) :: vec
		external rombint_obj
		
		vec(1) = m_e
		vec(2) = T
		vec(3) = p
		
		dme2_electron = PIx2 * alpha_fine * T*T/3.d0
		tmp = rombint_obj(vec, dme2_e_i1, 0.d0, 1d8, 1d-3, 200, 5)
		dme2_electron = dme2_electron + tmp*alpha_fine/PId4
		tmp = rombint_obj(vec, dme2_e_i2, 0.d0, 1d8, 1d-3, 200, 5)
		dme2_electron = dme2_electron + tmp*alpha_fine*m_e*m_e/(PId2*p)
	end function dme2_electron
	
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

end module
