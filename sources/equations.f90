module ndEquations
	use precision
	use variables
	use constants
	use ndErrors
	use ndinteractions
	use ndcosmology
	implicit none
	
	type(bspline_1d) :: Jfunc_interp, JPrime_interp, Kfunc_interp, Kprime_interp, yfunc_interp
	type(bspline_1d) :: g1func_interp, g2func_interp
	
	procedure (funcX), pointer :: J_func => null ()
	procedure (funcX), pointer :: K_func => null ()
	procedure (funcX), pointer :: Jprime => null ()
	procedure (funcX), pointer :: Kprime => null ()
	procedure (funcX), pointer :: Y_func => null ()
	procedure (func2_X), pointer :: G12_func => null ()
	
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
	
	function J_funcFull(o)
		real(dl) :: J_funcFull
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj

		J_funcFull = rombint_obj(o, J_int,0.d0,upper,toler,maxiter)/PISQ
	end function J_funcFull
	
	function Jprime_int(o, u)
		real(dl) :: Jprime_int
		real(dl), intent(in) :: o, u
		real(dl) :: uuoo, sqrtuuoo, expsqrtuuoo
		
		uuoo=u*u+o*o
		sqrtuuoo = sqrt(uuoo)
		expsqrtuuoo = exp(sqrtuuoo)
		
		Jprime_int = u*u * expsqrtuuoo * (1.d0-expsqrtuuoo) /(sqrtuuoo*(expsqrtuuoo+1)**3)
	end function Jprime_int
	
	function JprimeFull(o)
		real(dl) :: JprimeFull
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
	
		JprimeFull = rombint_obj(o, Jprime_int,0.d0,upper,toler,maxiter) * o / PISQ
	end function JprimeFull
	
	function K_int(o, u)
		real(dl) :: K_int,suo
		real(dl), intent(in) :: o, u
		
		suo=sqrt(u*u+o*o)
		K_int = u*u / (suo * (1.d0+exp(suo)))
	end function K_int
	
	function K_funcFull(o)
		real(dl) :: K_funcFull
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
		
		k_funcFull = rombint_obj(o, k_int,0.d0,upper,toler,maxiter)/PISQ
	end function K_funcFull
	
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
	
	function KprimeFull(o)
		real(dl) :: KprimeFull
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
	
		KprimeFull = rombint_obj(o, Kprime_int,0.d0,upper,toler, maxiter) * o / PISQ
	end function KprimeFull
	
	function Y_int(o, u)
		real(dl) :: Y_int, esuo
		real(dl), intent(in) :: o, u
		
		esuo=exp(sqrt(u*u+o*o))
		Y_int = u**4 * esuo / ((1.d0+esuo)**2)
	end function Y_int
	
	function Y_funcFull(o)
		real(dl) :: Y_funcFull
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
		
		Y_funcFull = rombint_obj(o, Y_int,0.d0,upper,toler, maxiter)/PISQ
	end function Y_funcFull
	
	function G12_funcFull(o)
		real(dl), dimension(2) :: G12_funcFull
		real(dl), intent(in) :: o
		real(dl) :: ko, jo, kp, jp, tmp
		ko=k_func(o)
		jo=j_func(o)
		kp=kprime(o)
		jp=jprime(o)
		tmp = (kp/6.d0 - ko*kp + jp/6.d0 +jp*ko + jo*kp)
		
		G12_funcFull(1) = PIx2*alpha_fine *(&
			(ko/3.d0 + 2*ko*ko - jo/6.d0 -ko*jo)/o + &
			tmp )
		G12_funcFull(2) = PIx2*alpha_fine*( &
			o * tmp &
			- 4.d0*( (ko+jo)/6.d0 + ko*jo -ko*ko/2.d0) )
	end function G12_funcFull
	
	subroutine init_interp_jkyg12
		real(dl), dimension(:),allocatable :: j, k, jp, kp, y, g1, g2
		real(dl), dimension(2) :: g12
		integer::ix, nx, iflag
		
		call addToLog("[equations] Initializing interpolation for J, J', K, K', Y, G_1, G_2...")
		nx=interp_nx
		allocate(j(nx),jp(nx),k(nx),kp(nx),y(nx),g1(nx),g2(nx))
		do ix=1, nx
			j(ix)  = J_funcFull(interp_xozvec(ix))
			jp(ix) = JprimeFull(interp_xozvec(ix))
			k(ix)  = K_funcFull(interp_xozvec(ix))
			kp(ix) = KprimeFull(interp_xozvec(ix))
			y(ix)  = Y_funcFull(interp_xozvec(ix))
			g12    = G12_funcFull(interp_xozvec(ix))
			g1(ix) = g12(1)
			g2(ix) = g12(2)
		end do
		call  Jfunc_interp%initialize(interp_xozvec, j,4,iflag)
		call JPrime_interp%initialize(interp_xozvec,jp,4,iflag)
		call  Kfunc_interp%initialize(interp_xozvec, k,4,iflag)
		call Kprime_interp%initialize(interp_xozvec,kp,4,iflag)
		call  yfunc_interp%initialize(interp_xozvec, y,4,iflag)
		call g1func_interp%initialize(interp_xozvec,g1,4,iflag)
		call g2func_interp%initialize(interp_xozvec,g2,4,iflag)

		J_func   => J_funcInterp
		K_func   => K_funcInterp
		Jprime   => JprimeInterp
		Kprime   => KprimeInterp
		Y_func   => Y_funcInterp
		G12_func => G12_funcInterp
		deallocate(j, k, jp, kp, y, g1, g2)
		call addToLog("[equations] ...done!")
	end subroutine init_interp_jkyg12
	
	function J_funcInterp(o)
		real(dl) :: J_funcInterp
		real(dl), intent(in) :: o
		integer :: iflag
		call Jfunc_interp%evaluate(o,0,J_funcInterp,iflag)
	end function J_funcInterp
	
	function K_funcInterp(o)
		real(dl) :: K_funcInterp
		real(dl), intent(in) :: o
		integer :: iflag
		call Kfunc_interp%evaluate(o,0,K_funcInterp,iflag)
	end function K_funcInterp
	
	function jprimeInterp(o)
		real(dl) :: jprimeInterp
		real(dl), intent(in) :: o
		integer :: iflag
		call JPrime_interp%evaluate(o,0,jprimeInterp,iflag)
	end function jprimeInterp
	
	function kprimeInterp(o)
		real(dl) :: kprimeInterp
		real(dl), intent(in) :: o
		integer :: iflag
		call Kprime_interp%evaluate(o,0,kprimeInterp,iflag)
	end function kprimeInterp
	
	function Y_funcInterp(o)
		real(dl) :: Y_funcInterp
		real(dl), intent(in) :: o
		integer :: iflag
		call yfunc_interp%evaluate(o,0,Y_funcInterp,iflag)
	end function Y_funcInterp
	
	function G12_funcInterp(o)
		real(dl), dimension(2) :: G12_funcInterp
		real(dl), intent(in) :: o
		integer :: iflag
		call g1func_interp%evaluate(o,0,G12_funcInterp(1),iflag)
		call g2func_interp%evaluate(o,0,G12_funcInterp(2),iflag)
	end function G12_funcInterp
	
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
		g12 = G12_func(x_o_z)
		
		params%x=x
		params%z=z
		tmp=0.d0

		do ix=1, flavorNumber
			params%flavor = ix
			tmp = tmp + rombint_obj(params, integrate_dRhoNu, y_min, y_max, toler, maxiter)*nuFactor(ix)
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
		
		call allocateCmplxMat(n1)
		call allocateCmplxMat(term1)
		call allocateCmplxMat(comm)
		call allocateCmplxMat(matrix)
		
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
		real(dl) :: x,y,z, overallNorm
		integer :: iy,k
		type(cmplxMatNN) :: n1, term1, comm
		
		call allocateCmplxMat(n1)
		call allocateCmplxMat(term1)
		call allocateCmplxMat(comm)
		call allocateCmplxMat(drhoy_dx_fullMat)
		
		n1 = nuDensMatVec(iy)
		y = nuDensMatVec(iy)%y
		
		overallNorm = planck_mass / (sqrt(radDensity(x,y,z)*PIx8D3))
		
		leptonDensities(1,1) = leptDensFactor * y / x**6 * 2 * electronDensity(x,z)
		term1%im = 0
		term1%re(:,:) = nuMassesMat/(2*y) + leptonDensities
!		print *,leptDensFactor, y , x**6 , electronDensity(x,z)
!		print *,'comm',leptonDensities, nuMassesMat, term1%re, n1%re
		
		!switch imaginary and real parts because of the "-i" factor
		call Commutator(term1%re, n1%re, comm%im)
		call Commutator(term1%re, n1%im, comm%re)
		
		comm%im = - comm%im * x**2/m_e_cub
		comm%re = comm%re * x**2/m_e_cub
		
!		print *,'comm2',comm%re, comm%im
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

!		print *,'fullMat',lastColl%mat%re,lastColl%mat%im
	end function drhoy_dx_fullMat
	
	subroutine saveRelevantInfo(x, vec)
		real(dl) :: x
		real(dl), dimension(:) :: vec
		real(dl), dimension(:), allocatable :: tmpvec
		integer :: k, i, j, iy, totFiles
		integer, dimension(:), allocatable :: units
		character(len=200) :: fname
		
		totFiles=flavorNumber**2+1
		allocate(units(totFiles),tmpvec(Ny))
		do i=1, totFiles
			units(i) = 8972 + i
		end do
		
		write(fname, '(A,E14.7)') '[output] Saving info at x=',x
		call addToLog(trim(fname))!not a filename but the above string
		
		do k=1, flavorNumber
			write(fname, '(A,I1,A)') trim(outputFolder)//'/nuDens_diag',k,'.dat'
			call openFile(units(k), trim(fname),firstWrite)
			do iy=1, nY
				tmpvec(iy)=nuDensMatVec(iy)%y**2*nuDensMatVec(iy)%re(k,k)
			end do
			write(units(k), '(*(E14.7))') x, tmpvec
		end do
		do i=1, flavorNumber
			do j=i,flavorNumber
				if (j.gt.i) then
					write(fname, '(A,I1,I1,A)') trim(outputFolder)//'/nuDens_nd_',i,j,'_re.dat'
					call openFile(units(k), trim(fname),firstWrite)
					tmpvec=0
					do iy=1, nY
						tmpvec(iy)=nuDensMatVec(iy)%re(i,j)
					end do
					write(units(k), '(*(E14.7))') x, tmpvec
					
					k=k+1
					write(fname, '(A,I1,I1,A)') trim(outputFolder)//'/nuDens_nd_',i,j,'_im.dat'
					call openFile(units(k), trim(fname),firstWrite)
					tmpvec=0
					do iy=1, nY
						tmpvec(iy)=nuDensMatVec(iy)%im(i,j)
					end do
					write(units(k), '(*(E14.7))') x, tmpvec
					k=k+1
				end if
			end do
		end do
		call openFile(units(k), trim(outputFolder)//'/z.dat', firstWrite)
		write(units(k), '(*(E14.7))') x, vec(ntot)
		
		do i=1, totFiles
			close(units(i))
		end do
		deallocate(units,tmpvec)
		
		firstWrite=.false.
	end subroutine saveRelevantInfo
	
	subroutine solver
		real(dl) :: xstart, xend, xchk
		real(dl), dimension(:), allocatable :: y, ydot
		integer :: ix, nchk, ix_in
		character(len=3) :: istchar
		character(len=100) :: tmpstring
		real(dl), dimension(:), allocatable :: ychk
		logical :: chk
		!for dlsoda:
		real(dl) :: rtol
		integer :: itol, itask, istate, iopt, lrw, liw, jt
		real(dl), dimension(:), allocatable :: rwork, atol
		integer, dimension(:), allocatable :: iwork
!		external jdum!derivatives, 
		
		itol=2
		rtol=dlsoda_rtol
		itask=1
		istate=1
		iopt=1
		
		lrw=22+ntot*(ntot+9)
		liw=20+ntot
		allocate(atol(ntot), rwork(lrw), iwork(liw))
		atol=dlsoda_atol
		rwork=0.
		iwork=0
		jt=2
		
		call densMat_2_vec(nuDensVec)
		nuDensVec(ntot)=z_in
		
		call readCheckpoints(nchk, xchk, ychk, chk)
		
		if (chk .and. &
			nchk.eq.ntot) then
			xstart=xchk
			nuDensVec=ychk
			firstWrite=.false.
			firstPoint=.true.
			ix_in=1 + int((log10(xchk)-logx_in)/(logx_fin-logx_in)*(Nx-1))
			write(tmpstring,"('ntot =',I4,' - x =',E14.7,' (i=',I4,') - z =',E14.7)"), nchk, xchk, ix_in, ychk(ntot)
			call addToLog("[ckpt] ##### Checkpoint file found. Will start from there. #####")
			call addToLog(trim(tmpstring))
		else
			xstart=x_arr(1)
			ix_in=1
		end if
		
		call addToLog("[solver] Starting DLSODA...")
		call saveRelevantInfo(xstart, nuDensVec)
		do ix=ix_in, Nx
			xend   = x_arr(ix)
			call dlsoda(derivatives,ntot,nuDensVec,xstart,xend,&
						itol,rtol,atol,itask,istate, &
						iopt,rwork,lrw,iwork,liw,jdum,jt)
			
			call writeCheckpoints(ntot, xend, nuDensVec)
			if (istate.lt.0) then
				write(istchar, "(I3)") istate
				call criticalError('istate='//istchar)
			end if
			
			call vec_2_densMat(nuDensVec)
			call saveRelevantInfo(xend, nuDensVec)
			xstart=xend
		end do
	end subroutine solver
	
	subroutine derivatives(n, x, vars, ydot)
		use omp_lib
		
		type(cmplxMatNN) :: mat
		integer :: n, i, j, k, m
		real(dl) :: x, z
		real(dl), dimension(n), intent(in) :: vars
		real(dl), dimension(n), intent(out) :: ydot
		integer :: flsq
		real(dl), dimension(:,:), allocatable :: tmpvec
		character(len=15) :: tmpstr
		
		flsq=flavorNumber**2
		allocate(tmpvec(Ny,flsq))
		write (tmpstr, '(E14.7)') x
		call addToLog("[eq] Calling derivatives...x="//trim(tmpstr))
		call openFile(6587, "tmpfile.txt",.false.)
		write(6587, '(*(E14.7))') x, vars
		close(6587)

		call allocateCmplxMat(mat)
		z = vars(n)
		call vec_2_densMat(vars(1:n-1))
		tmpvec=0
		!$omp parallel &
		!$omp default(shared) &
		!$omp private(mat,i,j)
		do m=1, Ny
!			print *,OMP_GET_THREAD_NUM(),"/",OMP_GET_NUM_THREADS()
			k=1
			mat = drhoy_dx_fullMat(x,z,m)
			do i=1, flavorNumber
				do j=i, flavorNumber
					tmpvec(m,k) = mat%re(i,j)
					k=k+1
				end do
				if (i.lt.flavorNumber) then
					do j=i+1, flavorNumber
						tmpvec(m,k) = mat%im(i,j)
						k=k+1
					end do
				end if
			end do
		end do
		!$omp end parallel
		do m=1, Ny
			do k=1,flsq
				ydot(k+(m-1)*flsq)=tmpvec(m,k)
			end do
		end do
		deallocate(tmpvec)
		
		ydot(ntot) = dz_o_dx(x,z)
		call sleep(1)
		call densMat_2_vec(nuDensVec)
		call deallocateCmplxMat(mat)
	end subroutine derivatives
	
	subroutine jdum

	end subroutine jdum
end module ndEquations
