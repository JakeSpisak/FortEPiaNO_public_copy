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
