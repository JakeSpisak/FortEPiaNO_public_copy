module ndEquations
	use precision
	use variables
	use utilities
	use constants
	use ndErrors
	use ndinteractions
	use ndcosmology
	implicit none
	
	type(bspline_1d) :: Jfunc_interp, JPrime_interp, Kfunc_interp, Kprime_interp, yfunc_interp
	type(bspline_1d) :: g1func_interp, g2func_interp
	
!	procedure (funcX), pointer :: J_func => null ()
!	procedure (funcX), pointer :: K_func => null ()
!	procedure (funcX), pointer :: Jprime => null ()
!	procedure (funcX), pointer :: Kprime => null ()
!	procedure (funcX), pointer :: Y_func => null ()
!	procedure (func2_X), pointer :: G12_func => null ()
	
	type collTerm
		type(cmplxMatNN) :: mat
		real(dl) :: x,y,z
	end type collTerm
	type(collTerm) :: lastColl
	
	real(dl), parameter :: upper = 1.d2
	
	real(dl) :: deriv_counter

	contains 

	subroutine densMat_2_vec(vec)
		real(dL), dimension(:), intent(out) :: vec
		integer :: i,j,k,m
		
		k=1
		do m=1, Ny
			do i=1, flavorNumber
				vec(k+i-1) = nuDensMatVec(m)%re(i,i)
			end do
			k=k+flavorNumber
			if (collision_offdiag.ne.3) then
				do i=1, flavorNumber-1
					do j=i+1, flavorNumber
						vec(k) = nuDensMatVec(m)%re(i,j)
						vec(k+1) = nuDensMatVec(m)%im(i,j)
						k=k+2
					end do
				end do
			end if
		end do
	end subroutine densMat_2_vec

	subroutine vec_2_densMat(vec)
		real(dL), dimension(:), intent(in) :: vec
		integer :: i,j,k,m
	
		k=1
		do m=1, Ny
			do i=1, flavorNumber
				nuDensMatVec(m)%re(i,i) = vec(k+i-1)
			end do
			k=k+flavorNumber
			if (collision_offdiag.ne.3) then
				do i=1, flavorNumber-1
					do j=i+1, flavorNumber
						nuDensMatVec(m)%re(i,j) = vec(k)
						nuDensMatVec(m)%im(i,j) = vec(k+1)
						nuDensMatVec(m)%re(j,i) = vec(k)
						nuDensMatVec(m)%im(j,i) = -vec(k+1)
						k=k+2
					end do
				end do
			end if
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

		J_funcFull = rombint_re(o, J_int,0.d0,upper,toler,maxiter)/PISQ
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
	
		JprimeFull = rombint_re(o, Jprime_int,0.d0,upper,toler,maxiter) * o / PISQ
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
		
		k_funcFull = rombint_re(o, k_int,0.d0,upper,toler,maxiter)/PISQ
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
	
		KprimeFull = rombint_re(o, Kprime_int,0.d0,upper,toler, maxiter) * o / PISQ
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
		
		Y_funcFull = rombint_re(o, Y_int,0.d0,upper,toler, maxiter)/PISQ
	end function Y_funcFull
	
	function G12_funcFull(o)
		real(dl), dimension(2) :: G12_funcFull
		real(dl), intent(in) :: o
		real(dl) :: ko, jo, kp, jp, tmp
		
		if (dme2_temperature_corr) then
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
		else
			G12_funcFull = 0.d0
		end if
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

!		J_func   => J_funcInterp
!		K_func   => K_funcInterp
!		Jprime   => JprimeInterp
!		Kprime   => KprimeInterp
!		Y_func   => Y_funcInterp
!		G12_func => G12_funcInterp
		deallocate(j, k, jp, kp, y, g1, g2)
		call addToLog("[equations] ...done!")
	end subroutine init_interp_jkyg12
	
	function J_func(o)
		real(dl) :: J_func
		real(dl), intent(in) :: o
		integer :: iflag
		call Jfunc_interp%evaluate(o,0,J_func,iflag)
	end function J_func
	
	function K_func(o)
		real(dl) :: K_func
		real(dl), intent(in) :: o
		integer :: iflag
		call Kfunc_interp%evaluate(o,0,K_func,iflag)
	end function K_func
	
	function jprime(o)
		real(dl) :: jprime
		real(dl), intent(in) :: o
		integer :: iflag
		call JPrime_interp%evaluate(o,0,jprime,iflag)
	end function jprime
	
	function kprime(o)
		real(dl) :: kprime
		real(dl), intent(in) :: o
		integer :: iflag
		call Kprime_interp%evaluate(o,0,kprime,iflag)
	end function kprime
	
	function Y_func(o)
		real(dl) :: Y_func
		real(dl), intent(in) :: o
		integer :: iflag
		call yfunc_interp%evaluate(o,0,Y_func,iflag)
	end function Y_func
	
	function G12_func(o)
		real(dl), dimension(2) :: G12_func
		real(dl), intent(in) :: o
		integer :: iflag
		call g1func_interp%evaluate(o,0,G12_func(1),iflag)
		call g2func_interp%evaluate(o,0,G12_func(2),iflag)
	end function G12_func
	
	function integrate_dRhoNu(params, y)
		real(dl) :: integrate_dRhoNu, y, dfnudx
		type(integrRhoNuPar) :: params
	
        call splint(y_arr, params%der, params%y2, params%ntot, y, dfnudx) !x, y(x), derivatives, N, xout, yout)
		integrate_dRhoNu = y**3 * dfnudx
	end function integrate_dRhoNu
	
	subroutine dz_o_dx(x,z, ydot, n)
		real(dl) :: dzodx
		real(dl), intent(in) :: x,z
		real(dl) :: x_o_z, jxoz, tmp
		real(dl), dimension(2) :: g12
		type(integrRhoNuPar) :: params
		real(dl), dimension(n) :: ydot
		real(dl), dimension(:), allocatable :: dertemp
		integer :: ix, n, m
		
		x_o_z = x/z
		jxoz = J_func(x_o_z)
		g12 = G12_func(x_o_z)
		
		params%x=x
		params%z=z
		params%ntot=Ny
		if (.not. allocated(params%der)) &
			allocate(params%der(Ny))
		if (.not. allocated(params%y2)) &
			allocate(params%y2(Ny))
!		params%der=ydot 
		tmp=0.d0
		
		do ix=1, flavorNumber
			params%flavor = ix
			params%der = 0.d0
			do m=1, Ny
!			print *,(m-1)*flavNumSqu + ix,m,ydot((m-1)*flavNumSqu + ix)
				params%der(m) = ydot((m-1)*flavNumSqu + ix) * fermiDirac_massless(y_arr(m), z)
			end do
			call spline(y_arr, params%der, Ny, 1.d30, 1.d30, params%y2) !x, y(x), N, ?, ?, derivatives)
			tmp = tmp + rombint_dRN(params, integrate_dRhoNu, y_min, y_max, toler, maxiter)*nuFactor(ix)
		end do
		
		write(*,"(*(E15.7))"), x_o_z * jxoz, g12(1), -1.d0/(2 * z**3 *PISQ)* tmp, &
			x_o_z**2*jxoz ,Y_func(x_o_z) , PISQ/7.5d0 ,g12(2)
		call sleep(1)
		dzodx = (x_o_z * jxoz + g12(1) &
			- 1.d0/(2 * z**3 *PISQ)* tmp) &
			/ (x_o_z**2*jxoz +Y_func(x_o_z) + PISQ/7.5d0 +g12(2))

		ydot(n)=dzodx
		
	end subroutine dz_o_dx
	
	subroutine drhoy_dx_fullMat (matrix,x,z,iy)
		type(cmplxMatNN),intent(out) :: matrix
		real(dl) :: x,y,z, overallNorm, a
		integer :: iy,k, ix
		type(cmplxMatNN), save :: n1, term1, comm
		
		call allocateCmplxMat(n1)
		call allocateCmplxMat(term1)
		call allocateCmplxMat(comm)
		call allocateCmplxMat(matrix)
		
		n1 = nuDensMatVec(iy)
		y = nuDensMatVec(iy)%y
!		print *,'a', n1%re!, n1%re
		do ix=1, flavorNumber
			n1%re(ix,ix) = n1%re(ix,ix) * fermiDirac_massless(y,z)
			n1%im(ix,ix) = n1%im(ix,ix) * fermiDirac_massless(y,z)
		end do
		
		overallNorm = planck_mass / (sqrt(radDensity(x,y,z)*PIx8D3))
		
		leptonDensities=0.d0
		leptonDensities(1,1) = leptDensFactor * y / x**6 * electronDensity(x,z)
		term1%im = 0
		term1%re(:,:) = nuMassesMat(:,:)/(2*y) + leptonDensities(:,:)
!		print *,leptDensFactor, y , x**6 , electronDensity(x,z)
!		print *,'comm',leptonDensities, nuMassesMat/(2*y)
!		print *,'b', n1%re!, term1%re
		
		!switch imaginary and real parts because of the "-i" factor
		call Commutator(term1%re, n1%re, comm%im)
		call Commutator(term1%re, n1%im, comm%re)
!		print *,'comm1',comm%re,comm%im
		
		matrix%im = - comm%im * x**2/m_e_cub
		matrix%re = comm%re * x**2/m_e_cub
		
!		print *,'comm2',matrix%re, matrix%im
		if (x .ne. lastColl%x .or. y .ne. lastColl%y .or. z .ne. lastColl%z) then
!			print *,'xyz',x,y,z
			lastColl%mat = collision_terms(x, z, y, n1)
			lastColl%x = x
			lastColl%y = y
			lastColl%z = z
		end if
!		print *,'coll',lastColl%mat%re
		matrix%re = matrix%re + lastColl%mat%re
		matrix%im = matrix%im + lastColl%mat%im
!		print *,'full',matrix%re

		matrix%re = matrix%re * overallNorm
		matrix%im = matrix%im * overallNorm
		do ix=1, flavorNumber
			matrix%re(ix,ix) = matrix%re(ix,ix) / fermiDirac_massless(y,z)
			matrix%im(ix,ix) = matrix%im(ix,ix) / fermiDirac_massless(y,z)
		end do
!		print *,'x,y:',x,y,z
!		print *,'fullMat',matrix%re!,matrix%im
!		call sleep(2)
	end subroutine drhoy_dx_fullMat
	
	subroutine saveRelevantInfo(x, vec)
		real(dl) :: x
		real(dl), dimension(:) :: vec
		real(dl), dimension(:), allocatable :: tmpvec
		integer :: k, i, j, iy, totFiles
		integer, dimension(:), allocatable :: units
		character(len=200) :: fname
		
		call vec_2_densMat(vec)
		totFiles=flavNumSqu+1
		allocate(units(totFiles),tmpvec(Ny))
		do i=1, totFiles
			units(i) = 8972 + i
		end do
		
		write(fname, '(A,'//dblfmt//')') '[output] Saving info at x=',x
		call addToLog(trim(fname))!not a filename but the above string
		
		do k=1, flavorNumber
			write(fname, '(A,I1,A)') trim(outputFolder)//'/nuDens_diag',k,'.dat'
			call openFile(units(k), trim(fname),firstWrite)
			do iy=1, nY
				tmpvec(iy)=nuDensMatVec(iy)%y**2*nuDensMatVec(iy)%re(k,k)! * &
!					fermiDirac_massless(nuDensMatVec(iy)%y, vec(ntot)) 
			end do
			write(units(k), '(*('//dblfmt//'))') x, tmpvec
		end do
		if (collision_offdiag.ne.3) then
			do i=1, flavorNumber-1
				do j=i+1,flavorNumber
					write(fname, '(A,I1,I1,A)') trim(outputFolder)//'/nuDens_nd_',i,j,'_re.dat'
					call openFile(units(k), trim(fname),firstWrite)
					tmpvec=0
					do iy=1, nY
						tmpvec(iy)=nuDensMatVec(iy)%re(i,j)
					end do
					write(units(k), '(*('//dblfmt//'))') x, tmpvec
					k=k+1
					
					write(fname, '(A,I1,I1,A)') trim(outputFolder)//'/nuDens_nd_',i,j,'_im.dat'
					call openFile(units(k), trim(fname),firstWrite)
					tmpvec=0
					do iy=1, nY
						tmpvec(iy)=nuDensMatVec(iy)%im(i,j)
					end do
					write(units(k), '(*('//dblfmt//'))') x, tmpvec
					k=k+1
				end do
			end do
		end if
		call openFile(units(k), trim(outputFolder)//'/z.dat', firstWrite)
		write(units(k), '(*('//dblfmt//'))') x, vec(ntot)
		
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
		
		deriv_counter = 0
		
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
		iwork(6)=99999999
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
			write(tmpstring,"('ntot =',I4,' - x =',"//dblfmt//",' (i=',I4,') - z =',"//dblfmt//")") &
				nchk, xchk, ix_in, ychk(ntot)
			call addToLog("[ckpt] ##### Checkpoint file found. Will start from there. #####")
			call addToLog(trim(tmpstring))
		else
			xstart=x_arr(1)
			ix_in=1
		end if
		
		call saveRelevantInfo(xstart, nuDensVec)
		do ix=ix_in+1, Nx
			xend   = x_arr(ix)
			write(tmpstring,"('x_start =',"//dblfmt//",' - x_end =',"//dblfmt//")"), xstart, xend
			call addToLog("[solver] Starting DLSODA..."//trim(tmpstring))
			call dlsoda(derivatives,ntot,nuDensVec,xstart,xend,&
						itol,rtol,atol,itask,istate, &
						iopt,rwork,lrw,iwork,liw,jdum,jt)
			
			if (istate.lt.0) then
				write(istchar, "(I3)") istate
				call criticalError('istate='//istchar)
			end if
			call writeCheckpoints(ntot, xend, nuDensVec)
			
			call saveRelevantInfo(xend, nuDensVec)
			xstart=xend
		end do
		write(tmpstring,"('x_end =',"//dblfmt//",' - z_end =',"//dblfmt//")"), xend, nuDensVec(ntot)
		call addToLog("[solver] Solver ended. "//trim(tmpstring))
	end subroutine solver
	
	subroutine derivatives(n, x, vars, ydot)
		use omp_lib
		
		type(cmplxMatNN), save :: mat
		integer :: n, i, j, k, m, s
		real(dl) :: x, z
		real(dl), dimension(n), intent(in) :: vars
		real(dl), dimension(n), intent(out) :: ydot
!		real(dl), dimension(:,:), allocatable, save :: tmpvec
		real(dl), dimension(:), allocatable, save :: tmpv
		character(len=100) :: tmpstr
!		!$omp threadprivate(leptonDensities,lastColl)
		
		deriv_counter = deriv_counter+1
!		if (.not. allocated(tmpvec)) &
!			allocate(tmpvec(Ny,flavNumSqu))
		if (.not. allocated(tmpv)) &
			allocate(tmpv(flavNumSqu))
		write (tmpstr, "('[eq] Calling derivatives (',"//dblfmt//",' x=',"//dblfmt//",')')") deriv_counter,x
		call printVerbose(trim(tmpstr),2)

		call allocateCmplxMat(mat)
		z = vars(n)
		call vec_2_densMat(vars(1:n-1))
!		tmpvec=0
!		!$omp parallel do &
!		!$omp default(shared) &
!		!$omp private(i,j,k,mat,tmpv)
		do m=1, Ny
			call drhoy_dx_fullMat(mat, x,z,m)
			do i=1, flavorNumber
				tmpv(i) = mat%re(i,i)
			end do
			k=flavorNumber+1
			if (collision_offdiag.ne.3) then
				do i=1, flavorNumber-1
					do j=i+1, flavorNumber
						tmpv(k) = mat%re(i,j)
						tmpv(k+1) = mat%im(i,j)
						k=k+2
					end do
				end do
			end if
!			!$omp critical
!			tmpvec(m,k) = tmpv(k)
			s=(m-1)*flavNumSqu
			do k=1,flavNumSqu
!				ydot(s+k)=tmpvec(m,k)
				ydot(s+k)=tmpv(k)
			end do
!			!$omp end critical
		end do
!		!$omp end parallel do
		
!		call openFile(6587, "tmpfile.txt",.false.)
!        write(6587, '(*(E14.7))') x, ydot
!        close(6587)
        
		call dz_o_dx(x,z, ydot, ntot)
		
		call densMat_2_vec(nuDensVec)
	end subroutine derivatives
	
	subroutine jdum

	end subroutine jdum
end module ndEquations
