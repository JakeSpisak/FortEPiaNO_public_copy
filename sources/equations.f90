module ndEquations
	use precision
	use constants
	use ndErrors
	use ndinteractions
	use ndcosmology
	implicit none
	
	type collTerm
		type(cmplxMatNN) :: mat
		real(dl) :: x,y,z
	end type collTerm
	type(collTerm) :: lastColl
	
	type integrRhoNuPar
		real(dl) :: x,z
		integer :: flavor
	end type integrRhoNuPar
	
	real(dl), parameter :: upper = 1.d2

	contains 

	subroutine densMat_2_vec(vec)
		real(dL), dimension(:), intent(out) :: vec
		integer :: i,j,k,m
		
		k=1
		do m=1, Ny
			do i=1, flavorNumber
				do j=i, flavorNumber
					vec(k) = nuDensMatVec(m)%re(i,j)
					k=k+1
				end do
				if (i.lt.flavorNumber) then
					do j=i+1, flavorNumber
						vec(k) = nuDensMatVec(m)%im(i,j)
						k=k+1
					end do
				end if
			end do
		end do
	end subroutine densMat_2_vec

	subroutine vec_2_densMat(vec)
		real(dL), dimension(:), intent(in) :: vec
		integer :: i,j,k,m
	
		k=1
		do m=1, Ny
			do i=1, flavorNumber
				do j=i, flavorNumber
					nuDensMatVec(m)%re(i,j) = vec(k)
					k=k+1
				end do
				if (i.lt.flavorNumber) then
					do j=i+1, flavorNumber
						nuDensMatVec(m)%im(i,j) = vec(k)
						k=k+1
					end do
				end if
			end do
		end do
	end subroutine vec_2_densMat
	
	function J_int(o, u)
		real(dl) :: J_int, esuo
		real(dl), intent(in) :: o, u
		
		esuo=exp(sqrt(u*u+o*o))
		J_int = u*u * esuo / ((1.d0+esuo)**2)
	end function J_int
	
	function J_func(o)
		real(dl) :: J_func
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
		
		J_func = rombint_obj(o, J_int,0,upper,toler)/PISQ
	end function J_func
	
	function Jprime_int(o, u)
		real(dl) :: Jprime_int
		real(dl), intent(in) :: o, u
		real(dl) :: uuoo, sqrtuuoo, expsqrtuuoo
		
		uuoo=u*u+o*o
		sqrtuuoo = sqrt(uuoo)
		expsqrtuuoo = exp(sqrtuuoo)
		
		Jprime_int = u*u * expsqrtuuoo * (1.d0-expsqrtuuoo) /(sqrtuuoo*(expsqrtuuoo+1)**3)
	end function Jprime_int
	
	function Jprime(o)
		real(dl) :: Jprime
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
	
		Jprime = rombint_obj(o, Jprime_int,0,upper,toler) * o / PISQ
	end function Jprime
	
	function K_int(o, u)
		real(dl) :: K_int,suo
		real(dl), intent(in) :: o, u
		
		suo=sqrt(u*u+o*o)
		K_int = u*u / (suo * (1.d0+exp(suo)))
	end function K_int
	
	function K_func(o)
		real(dl) :: K_func
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
		
		k_func = rombint_obj(o, k_int,0,upper,toler)/PISQ
	end function K_func
	
	function Kprime_int(o, u)
		real(dl) :: Kprime_int
		real(dl), intent(in) :: o, u
		real(dl) :: uuoo, sqrtuuoo, expsqrtuuoo
		
		uuoo=u*u+o*o
		sqrtuuoo = sqrt(uuoo)
		expsqrtuuoo = exp(sqrtuuoo)
		
		Kprime_int = - u*u / (uuoo*sqrtuuoo*(expsqrtuuoo+1)**2) * &
			(1.d0 + expsqrtuuoo*(sqrtuuoo+1.d0))
	end function Kprime_int
	
	function Kprime(o)
		real(dl) :: Kprime
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
	
		Kprime = rombint_obj(o, Kprime_int,0,upper,toler) * o / PISQ
	end function Kprime
	
	function Y_int(o, u)
		real(dl) :: Y_int, esuo
		real(dl), intent(in) :: o, u
		
		esuo=exp(sqrt(u*u+o*o))
		Y_int = u**4 * esuo / ((1.d0+esuo)**2)
	end function Y_int
	
	function Y_func(o)
		real(dl) :: Y_func
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
		
		Y_func = rombint_obj(o, Y_int,0,upper,toler)/PISQ
	end function Y_func
	
	function G12_func(o)
		real(dl), dimension(2) :: G12_func
		real(dl), intent(in) :: o
		real(dl) :: ko, jo, kp, jp, tmp
		ko=k_func(o)
		jo=j_func(o)
		kp=kprime(o)
		jp=jprime(o)
		tmp = (kp/6.d0 - ko*kp + jp/6.d0 +jp*ko + jo*kp)
		
		G12_func(1) = PIx2*alpha_fine *(&
			(ko/3.d0 + 2*ko*ko - jo/6.d0 -ko*jo)/o + &
			tmp )
		G12_func(2) = PIx2*alpha_fine*( &
			o * tmp &
			- 4.d0*( (ko+jo)/6.d0 + ko*jo -ko*ko/2.d0) )
	end function G12_func
	
	subroutine init_interp_jkyg12
		!must be written
	end subroutine init_interp_jkyg12
	
	function integrate_dRhoNu(params, y)
		real(dl) :: integrate_dRhoNu, y
		type(integrRhoNuPar) :: params
	
		integrate_dRhoNu = y**3 * drhoy_dx_ij_rI(params%x,y,params%z, params%flavor, params%flavor, 1)
	end function integrate_dRhoNu
	
	function dz_o_dx(x,z)
		real(dl) :: dz_o_dx, rombint_obj
		real(dl), intent(in) :: x,z
		real(dl) :: x_o_z, jxoz, tmp
		real(dl), dimension(2) :: g12
		type(integrRhoNuPar) :: params
		integer :: ix
		external rombint_obj
		
		x_o_z = x/z
		jxoz = J_func(x_o_z)
		g12 = g12_func(x_o_z)
		
		params%x=x
		params%z=z
		tmp=0.d0

		do ix=1, flavorNumber
			params%flavor = ix
			tmp = tmp + rombint_obj(params, integrate_dRhoNu, y_min, y_max, toler)
		end do
		
		dz_o_dx = (x_o_z * jxoz + g12(1) &
			- 1.d0/(2 * z**3 *PISQ)* tmp) &
			/ (x_o_z**2*jxoz +Y_func(x_o_z) + PISQ/7.5d0 +g12(2))
		
	end function dz_o_dx
	
	function drhoy_dx_ij_rI (x,y,z, i, j, rI)
		type(cmplxMatNN) :: matrix
		real(dl) :: drhoy_dx_ij_rI
		real(dl) :: x,y,z, overallNorm
		integer :: i,j,k, rI
		type(cmplxMatNN) :: n1, term1, comm
		
		drhoy_dx_ij_rI = 0.d0
		overallNorm = planck_mass / (sqrt(radDensity(x,y,z)*PIx8D3))
		
		n1 = interp_nuDens(y)
		leptonDensities(1,1) = leptDensFactor * y / x**6 * 2 * electronDensity(x,z)
		term1%im = 0
		term1%re(:,:) = nuMassesMat/(2*y) + leptonDensities
		!switch imaginary and real parts because of the "-i" factor
		call Commutator(term1%re, n1%re, comm%im)
		call Commutator(term1%re, n1%im, comm%re)
		
		comm%im = - comm%im * x**2/m_e_cub
		comm%re = comm%re * x**2/m_e_cub
		
		matrix%re = overallNorm * comm%re
		matrix%im = overallNorm * comm%im 
		
		if (x .ne. lastColl%x .or. y .ne. lastColl%y .or. z .ne. lastColl%z) then
			lastColl%mat = collision_terms(x, z, y)
			lastColl%x = x
			lastColl%y = y
			lastColl%z = z
		end if
		matrix%re = matrix%re + lastColl%mat%re * overallNorm
		matrix%im = matrix%im + lastColl%mat%im * overallNorm

		if (rI.eq.1) then
			drhoy_dx_ij_rI = matrix%re(i,j)
		else if (rI.eq.2) then
			drhoy_dx_ij_rI = matrix%im(i,j)
		end if

	end function drhoy_dx_ij_rI
	
	function drhoy_dx_fullMat (x,z,iy)
		type(cmplxMatNN) :: drhoy_dx_fullMat
		real(dl) :: drhoy_dx_ij_rI
		real(dl) :: x,y,z, overallNorm
		integer :: iy,k
		type(cmplxMatNN) :: n1, term1, comm
		
		n1 = nuDensMatVec(iy)
		y = nuDensMatVec(iy)%y
		
		overallNorm = planck_mass / (sqrt(radDensity(x,y,z)*PIx8D3))
		
		leptonDensities(1,1) = leptDensFactor * y / x**6 * 2 * electronDensity(x,z)
		term1%im = 0
		term1%re(:,:) = nuMassesMat/(2*y) + leptonDensities
		!switch imaginary and real parts because of the "-i" factor
		call Commutator(term1%re, n1%re, comm%im)
		call Commutator(term1%re, n1%im, comm%re)
		
		comm%im = - comm%im * x**2/m_e_cub
		comm%re = comm%re * x**2/m_e_cub
		
		drhoy_dx_fullMat%re = overallNorm * comm%re
		drhoy_dx_fullMat%im = overallNorm * comm%im 
		
		if (x .ne. lastColl%x .or. y .ne. lastColl%y .or. z .ne. lastColl%z) then
			lastColl%mat = collision_terms(x, z, y)
			lastColl%x = x
			lastColl%y = y
			lastColl%z = z
		end if
		drhoy_dx_fullMat%re = drhoy_dx_fullMat%re + lastColl%mat%re * overallNorm
		drhoy_dx_fullMat%im = drhoy_dx_fullMat%im + lastColl%mat%im * overallNorm

	end function drhoy_dx_fullMat
	
	subroutine saveRelevantInfo(x, vec)
		real(dl) :: x
		real(dl), dimension(:) :: vec
	
	end subroutine saveRelevantInfo
	
	subroutine solver
		real(dl) :: xstart, xend
		integer :: ix
		!for dlsoda:
		real(dl) :: rtol
		integer :: itol, itask, istate, iopt, lrw, liw, jt
		real(dl), dimension(:), allocatable :: rwork, atol
		integer, dimension(:), allocatable :: iwork
		external derivatives, jdum
		
		itol=2
		rtol=1.d-3
		itask=1
		istate=1
		iopt=0
		
		lrw=22+ntot*(ntot+9)
		liw=20+ntot
		allocate(atol(ntot), rwork(lrw), iwork(liw))
		atol=1.d-6
		rwork=0.
		iwork=0
		jt=2
		
		call densMat_2_vec(nuDensVec)
		nuDensVec(ntot)=z_in
		print *,nuDensVec
		
		xstart=x_arr(1)
		do ix=1, Nx/printEveryNIter
			xend   = x_arr((ix)*printEveryNIter)
			call dlsoda(derivatives,ntot,nuDensVec,xstart,xend,&
						itol,rtol,atol,itask,istate, &
						iopt,rwork,lrw,iwork,liw,jdum,jt)
			if (istate.lt.0) &
				call criticalError('istate=')!,istate
				
			call saveRelevantInfo(xend, nuDensVec)
			xstart=xend
		end do
	end subroutine solver
end module ndEquations
