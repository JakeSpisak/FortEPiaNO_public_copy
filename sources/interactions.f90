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

	!kappa function for damping terms (YYYW)
	elemental function kappa_damp(a, b, c)
		real(dl) :: kappa_damp
		real(dl), intent(in) :: a, b, c
		kappa_damp = a*a*(c-b) - 2.d0/3.d0*a* (c**3-b**3) + (c**5-b**5)/5.d0
	end function kappa_damp

	!function that returns the integrand of the nunu damping coefficient in YYYW terms
	elemental function nunu_damp_integrand(p, p1, q)
		real(dl) :: nunu_damp_integrand
		real(dl), intent(in) :: p, q, p1
		real(dl) :: q1
		real(dl) :: as, bs, cs, ap, bp, cp
		real(dl) :: fq, fp1, fq1
		q1 = p + q - p1

		cs = p + q
		as = cs**2
		bs = max(abs(p-q), abs(p1-q1))

		ap = (q - p1)**2
		bp = abs(p1 - q)
		cp = min(p+q1, q+p1)

		fq = fermiDirac(q)
		fp1 = fermiDirac(p1)
		fq1 = fermiDirac(q1)

		nunu_damp_integrand = &
			(kappa_damp(as, bs, cs) + 2.d0 * kappa_damp(ap, bp, cp)) &
			* ((1.d0-fq)*fp1*fq1 + fq*(1.d0-fp1)*(1.d0-fq1))
	end function nunu_damp_integrand

	!integral of the above function
	function dy_damping(y)
		real(dl) :: dy_damping
		real(dl), intent(in) :: y
		integer :: ia, ib
		integer, parameter :: nydamp = 2**13
		real(dl), dimension(:), allocatable :: ya, dya
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(ya(nydamp), dya(nydamp))
		ya = linspace(0.d0, 100.d0, nydamp)
		do ia=1, nydamp-1
			dya(ia) = ya(ia+1) - ya(ia)
		end do
		allocate(fy2_arr(nydamp, nydamp))
		fy2_arr = 0.d0
		do ia=1, nydamp
			do ib=1, nydamp
				if (ya(ib) .gt. ya(ia)-y) &
					fy2_arr(ia, ib) = nunu_damp_integrand(y, ya(ia), ya(ib))
			end do
		end do
		dy_damping = integral_NC_2d(nydamp, nydamp, dya, dya, fy2_arr) / y**3
		deallocate(ya, dya, fy2_arr)
	end function dy_damping

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
		use fpInterfaces1
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
		use fpInterfaces1
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
		use fpInterfaces1
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
		use fpInterfaces1
		use fpInterfaces2
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
		use fpInterfaces1
		use fpInterfaces2
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
		use fpInterfaces2
		use fpInterfaces3
		implicit None
		procedure (collision_integrand) :: Fint
		type(cmplxMatNN) :: get_collision_terms
		procedure (collision_integrator), pointer :: integrator
		real(dl) :: x, w, z, y1, cf, w4, z4
		integer :: iy1
		type(coll_args), intent(in) :: collArgsIn
		type(coll_args) :: collArgs
		integer :: i, j

		call allocateCmplxMat(get_collision_terms)

		collArgs=collArgsIn

		x = collArgs%x
		w = collArgs%w
		z = collArgs%z
		w4 = w**4
		z4 = z**4
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
			!diagonal elements:
			if (.not.sterile(i)) then
				if (collision_offdiag.eq.5) then !only damping factors (YYYW)
					get_collision_terms%re(i,i) = &
						- ( &
							dampTermMatrixCoeffNue(i,i) * (nuDensMatVecFD(iy1)%re(i,i) - fermiDirac(y_arr(iy1)/z))&
							+ dampTermMatrixCoeffNunu(i,i) * (nuDensMatVecFD(iy1)%re(i,i) - fermiDirac(y_arr(iy1)/w))&
						) &
						* dampTermYYYWdy(iy1)
				else !full integration for nue
					get_collision_terms%re(i,i) = &
						integrator(Fint, collArgs, F_ab_ann_re, F_ab_sc_re)
					if (collision_offdiag.eq.4) then !add nunu diagonal damping (YYYW)
						get_collision_terms%re(i,i) = get_collision_terms%re(i,i) &
							- (dampTermMatrixCoeffNunu(i,i)) &
							* dampTermYYYWdy(iy1) &
							* (nuDensMatVecFD(iy1)%re(i,i) - fermiDirac(y_arr(iy1)/w))
					end if
				end if
			end if
			!off-diagonal elements:
			if (collision_offdiag.eq.1) then !full integration for nue, dampings for nunu (disabled for tests)
				do j=i+1, flavorNumber
					collArgs%ix2 = j
					get_collision_terms%re(i,j) = integrator(Fint, collArgs, F_ab_ann_re, F_ab_sc_re)
					get_collision_terms%im(i,j) = integrator(Fint, collArgs, F_ab_ann_im, F_ab_sc_im)
#ifndef DO_TESTS
					get_collision_terms%re(i,j) = &
						get_collision_terms%re(i,j) &
						+ dampTermFactor * dampTermMatrixCoeffNunu(i,j) &
						* z4 * y1*y1*y1 * nuDensMatVecFD(iy1)%re(i,j)
					get_collision_terms%im(i,j) = &
						get_collision_terms%im(i,j) &
						+ dampTermFactor * dampTermMatrixCoeffNunu(i,j) &
						* z4 * y1*y1*y1 * nuDensMatVecFD(iy1)%im(i,j)
#endif
				end do
			else if (collision_offdiag.eq.2) then !nue and nunu dampings from McKellar:1992ja
				do j=i+1, flavorNumber
					get_collision_terms%re(i,j) = &
						dampTermFactor &
						* (dampTermMatrixCoeffNue(i,j)+dampTermMatrixCoeffNunu(i,j)) &
						* z4 * y1*y1*y1 * nuDensMatVecFD(iy1)%re(i,j)
					get_collision_terms%im(i,j) = &
						dampTermFactor &
						* (dampTermMatrixCoeffNue(i,j)+dampTermMatrixCoeffNunu(i,j)) &
						* z4 * y1*y1*y1 * nuDensMatVecFD(iy1)%im(i,j)
				end do
			else if (collision_offdiag.eq.4 .or. collision_offdiag.eq.5) then !nue and nunu dampings from YYYW
				do j=i+1, flavorNumber
					get_collision_terms%re(i,j) = &
						- (dampTermMatrixCoeffNue(i,j)+dampTermMatrixCoeffNunu(i,j)) &
						* dampTermYYYWdy(iy1) * nuDensMatVecFD(iy1)%re(i,j)
					get_collision_terms%im(i,j) = &
						- (dampTermMatrixCoeffNue(i,j)+dampTermMatrixCoeffNunu(i,j)) &
						* dampTermYYYWdy(iy1) * nuDensMatVecFD(iy1)%im(i,j)
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
