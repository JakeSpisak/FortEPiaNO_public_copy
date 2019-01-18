module ndEquations
	use precision
	use variables
	use utilities
	use constants
	use ndErrors
	use ndinteractions
	use ndcosmology
	use bspline_module
	use linear_interpolation_module
	use sg_interpolate
	implicit none

	type(bspline_1d) :: dzodx_A_interp, dzodx_B_interp
	
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
	
	elemental function J_int(o, u)
		real(dl) :: J_int, esuo
		real(dl), intent(in) :: o, u

		esuo=exp(sqrt(u*u+o*o))
		J_int = u*u * esuo / ((1.d0+esuo)**2)
	end function J_int
	
	pure function J_funcFull(o)
		real(dl) :: J_funcFull
		real(dl), intent(in) :: o

		J_funcFull = rombint_re(o, J_int,0.d0,upper,toler_jkyg,maxiter)/PISQ
	end function J_funcFull
	
	elemental function Jprime_int(o, u)
		real(dl) :: Jprime_int
		real(dl), intent(in) :: o, u
		real(dl) :: uuoo, sqrtuuoo, expsqrtuuoo
		
		uuoo=u*u+o*o
		sqrtuuoo = sqrt(uuoo)
		expsqrtuuoo = exp(sqrtuuoo)
		
		Jprime_int = u*u * expsqrtuuoo * (1.d0-expsqrtuuoo) /(sqrtuuoo*(expsqrtuuoo+1)**3)
	end function Jprime_int
	
	pure function JprimeFull(o)
		real(dl) :: JprimeFull
		real(dl), intent(in) :: o
	
		JprimeFull = rombint_re(o, Jprime_int,0.d0,upper,toler_jkyg,maxiter) * o / PISQ
	end function JprimeFull
	
	elemental function K_int(o, u)
		real(dl) :: K_int,suo
		real(dl), intent(in) :: o, u
		
		suo=sqrt(u*u+o*o)
		K_int = u*u / (suo * (1.d0+exp(suo)))
	end function K_int
	
	pure function K_funcFull(o)
		real(dl) :: K_funcFull
		real(dl), intent(in) :: o
		
		k_funcFull = rombint_re(o, k_int,0.d0,upper,toler_jkyg,maxiter)/PISQ
	end function K_funcFull
	
	elemental function Kprime_int(o, u)
		real(dl) :: Kprime_int
		real(dl), intent(in) :: o, u
		real(dl) :: uuoo, sqrtuuoo, expsqrtuuoo
		
		uuoo=u*u+o*o
		sqrtuuoo = sqrt(uuoo)
		expsqrtuuoo = exp(sqrtuuoo)
		
		Kprime_int = u*u / (uuoo*sqrtuuoo*(expsqrtuuoo+1)**2) * &
			(1.d0 + expsqrtuuoo*(sqrtuuoo+1.d0))
	end function Kprime_int
	
	pure function KprimeFull(o)
		real(dl) :: KprimeFull
		real(dl), intent(in) :: o
	
		KprimeFull = - rombint_re(o, Kprime_int,0.d0,upper,toler_jkyg, maxiter) * o / PISQ
	end function KprimeFull
	
	elemental function Y_int(o, u)
		real(dl) :: Y_int, esuo
		real(dl), intent(in) :: o, u
		
		esuo=exp(sqrt(u*u+o*o))
		Y_int = u**4 * esuo / ((1.d0+esuo)**2)
	end function Y_int
	
	pure function Y_funcFull(o)
		real(dl) :: Y_funcFull
		real(dl), intent(in) :: o
		
		Y_funcFull = rombint_re(o, Y_int,0.d0,upper,toler_jkyg, maxiter)/PISQ
	end function Y_funcFull
	
	pure function G12_funcFull(o)
		real(dl), dimension(2) :: G12_funcFull
		real(dl), intent(in) :: o
		real(dl) :: ko, jo, kp, jp, tmp
		
		if (dme2_temperature_corr) then
			ko=k_funcFull(o)
			jo=j_funcFull(o)
			kp=kprimeFull(o)
			jp=jprimeFull(o)
			tmp = (kp/6.d0 - ko*kp + jp/6.d0 +jp*ko + jo*kp)
			
			G12_funcFull(1) = PIx2*alpha_fine *(&
				(ko/3.d0 + 2.d0*ko*ko - jo/6.d0 -ko*jo)/o + &
				tmp )
			G12_funcFull(2) = PIx2*alpha_fine*( &
				o * tmp &
				- 4.d0*( (ko+jo)/6.d0 + ko*jo -ko*ko/2.d0) )
		else
			G12_funcFull = 0.d0
		end if
	end function G12_funcFull
	
	subroutine init_interp_jkyg12
		real(dl) :: num, den, j, y
		real(dl), dimension(:),allocatable :: A, B
		real(dl), dimension(2) :: g12
		integer::ix, nx, iflag
		
		call addToLog("[equations] Initializing interpolation for coefficients in dz/dx...")
		nx=interp_nx
		allocate(A(nx),B(nx))
		do ix=1, nx
			j      = J_funcFull(interp_xozvec(ix))
			y      = Y_funcFull(interp_xozvec(ix))
			g12    = G12_funcFull(interp_xozvec(ix))
			num= interp_xozvec(ix) * j + g12(1)
			den= interp_xozvec(ix)**2 * j + y + PISQ/7.5d0 + g12(2)
			A(ix) = num / den
			B(ix) = 1./(2.d0*PISQ*den)
		end do
		call dzodx_A_interp%initialize(interp_xozvec,A,4,iflag)
		call dzodx_B_interp%initialize(interp_xozvec,B,4,iflag)

		deallocate(A, B)
		call addToLog("[equations] ...done!")
	end subroutine init_interp_jkyg12
	
	function dzodxcoef_interp_func(o)
		real(dl), dimension(2) :: dzodxcoef_interp_func
		real(dl), intent(in) :: o
		integer :: iflag
		call dzodx_A_interp%evaluate(o,0,dzodxcoef_interp_func(1),iflag)
		call dzodx_B_interp%evaluate(o,0,dzodxcoef_interp_func(2),iflag)
	end function dzodxcoef_interp_func
	
	pure function integrate_dRhoNu(params, y)
		real(dl) :: integrate_dRhoNu
		real(dl), intent(in) :: y
		type(spline_class), intent(in) :: params
	
		integrate_dRhoNu = params%evaluate(y)
	end function integrate_dRhoNu
	
	subroutine dz_o_dx(x,z, ydot, n)!eq 17 from doi:10.1016/S0370-2693(02)01622-2
		real(dl), intent(in) :: x,z
		integer, intent(in) :: n
		real(dl), dimension(n), intent(inout) :: ydot
		real(dl) :: dzodx, tmp
		real(dl), dimension(2) :: coeffs
		type(spline_class) :: params
		real(dl), dimension(:), allocatable :: der
		integer :: ix, m

		coeffs = dzodxcoef_interp_func(x/z)

		allocate(der(Ny))

		der = 0.d0
		do m=1, Ny
			do ix=1, flavorNumber
				der(m) = der(m) + ydot((m-1)*flavNumSqu + ix) * nuFactor(ix)
			end do
			der(m) = y_arr(m)**3 * der(m) * fermiDirac(y_arr(m) / z)
		end do
		call params%initialize(Ny, y_arr, der)
		tmp = rombint_spli(params, integrate_dRhoNu, y_min, y_max, toler, maxiter)

		dzodx = coeffs(1) - coeffs(2)/ z**3 * tmp

		ydot(n)=dzodx

	end subroutine dz_o_dx

	subroutine dz_o_dx_lin(x,z, ydot, n)
		!eq 17 from doi:10.1016/S0370-2693(02)01622-2
		!integral without need of interpolation, just use linear approx in y bins
		real(dl), intent(in) :: x,z
		integer, intent(in) :: n
		real(dl), dimension(n), intent(inout) :: ydot
		real(dl) :: tmp
		real(dl), dimension(2) :: coeffs
		integer :: ix, m

		do m=1, Ny
			tmp = 0.d0
			do ix=1, flavorNumber
				tmp = tmp + ydot((m-1)*flavNumSqu + ix) * nuFactor(ix)
			end do
			fy_arr(m) = y_arr(m)**3 * tmp * fermiDirac(y_arr(m) / z)
		end do
		tmp = integral_linearized_1d(Ny, dy_arr, fy_arr)
		coeffs = dzodxcoef_interp_func(x/z)
		ydot(n) = coeffs(1) - coeffs(2)/ z**3 * tmp
	end subroutine dz_o_dx_lin

	subroutine test_dzodx_speed
		real(dl), dimension(:), allocatable :: fd_vec, fd_x
		integer :: ix, iflag, interp_nfd
		real(dl) :: x,z
		real(8) :: timer1
		integer, parameter :: n=901
		real(dl), dimension(n) :: ydot
		integer :: m

		ydot = 0.d0
		do m=1, Ny
			ydot((m-1)*flavNumSqu + 1) = cos(0.02d0*y_arr(m))
			ydot((m-1)*flavNumSqu + 2) = y_arr(m)/20.d0
			ydot((m-1)*flavNumSqu + 3) = 1.d0
		end do

		call addToLog("[equations] Testing dz_o_dx speed...")
		call random_seed()
		call random_number(x)
		call random_number(z)
		call tic(timer1)
		do ix=1, 1000000
			call random_number(x)
			call random_number(z)
			x=(x_fin-x_in)*x + x_in
			z=0.4d0*z + z_in
			call dz_o_dx(x, z, ydot, n)
			call dz_o_dx_lin(x, z, ydot, n)
		end do
		call toc(timer1, "<reset>")

		call tic(timer1)
		do ix=1, 1000000
			call random_number(x)
			call random_number(z)
			x=(x_fin-x_in)*x + x_in
			z=0.4d0*z + z_in
			call dz_o_dx(x, z, ydot, n)
		end do
		call toc(timer1, "<rombint>")

		call tic(timer1)
		do ix=1, 1000000
			call random_number(x)
			call random_number(z)
			x=(x_fin-x_in)*x + x_in
			z=0.4d0*z + z_in
			call dz_o_dx_lin(x, z, ydot, n)
		end do
		call toc(timer1, "<linearized>")

		call addToLog("[equations] ...done!")
	end subroutine test_dzodx_speed

	subroutine drhoy_dx_fullMat (matrix,x,z,iy)
		type(cmplxMatNN), intent(out) :: matrix
		real(dl) :: x,y,z, overallNorm, a, fd, cf
		integer :: iy,k, ix
		type(cmplxMatNN), save :: n1, comm
		
		call allocateCmplxMat(n1)
		call allocateCmplxMat(comm)
		call allocateCmplxMat(matrix)
		
		n1 = nuDensMatVec(iy)
		y = nuDensMatVec(iy)%y
		fd = fermiDirac(y / z)
		do ix=1, flavorNumber
			n1%re(ix,ix) = n1%re(ix,ix) * fd
			n1%im(ix,ix) = n1%im(ix,ix) * fd
		end do
		
		overallNorm = planck_mass / (sqrt(radDensity(x,z)*PIx8D3))
		
		leptonDensities=0.d0
		leptonDensities(1,1) = leptDensFactor * y / x**6 * electronDensity(x,z)
!		matrix%im = 0.d0
		matrix%re(:,:) = nuMassesMat(:,:)/(2*y) + leptonDensities(:,:)
		
		!switch imaginary and real parts because of the "-i" factor
		call Commutator(matrix%re, n1%re, comm%im)
		call Commutator(matrix%re, n1%im, comm%re)
		
		!matrix is now redefined
		cf = x**2/m_e_cub
		matrix%im = - comm%im * cf
		matrix%re = comm%re * cf
		
		if (x .ne. lastColl%x .or. y .ne. lastColl%y .or. z .ne. lastColl%z) then
			lastColl%mat = collision_terms(x, z, y, n1)
			lastColl%x = x
			lastColl%y = y
			lastColl%z = z
		end if
		matrix%re = matrix%re + lastColl%mat%re
		matrix%im = matrix%im + lastColl%mat%im

		matrix%re = matrix%re * overallNorm
		matrix%im = matrix%im * overallNorm
		do ix=1, flavorNumber
			matrix%re(ix,ix) = matrix%re(ix,ix) / fd
			matrix%im(ix,ix) = matrix%im(ix,ix) / fd
		end do
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
!					fermiDirac(nuDensMatVec(iy)%y / vec(ntot))
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
		real(dl) :: rtol
		integer :: itol, itask, istate, iopt, lrw, liw, jt
		real(dl), dimension(:), allocatable :: rwork, atol
		integer, dimension(:), allocatable :: iwork
		
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
			call saveRelevantInfo(xstart, nuDensVec)
		end if
		
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

	subroutine derivative (x, z, m, mat, n, ydot)
!		compute rho derivatives for a given momentum y_arr(m), save to ydot
		real(dl), intent(in) :: x, z
		integer, intent(in) :: m, n
		real(dl), dimension(n), intent(out) :: ydot
		integer :: i, j, k, s
		real(dl), dimension(:), allocatable :: tmpv
		type(cmplxMatNN) :: mat
	
		allocate(tmpv(flavNumSqu))
		call drhoy_dx_fullMat(mat, x,z,m)
		tmpv=0.d0
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
		s=(m-1)*flavNumSqu
		do i=1,flavNumSqu
			ydot(s+i)=tmpv(i)
		end do
	end subroutine derivative

	subroutine derivatives(n, x, vars, ydot)
!		compute all the rho derivatives (drho/dx for all y, dz/dx)
!		needs allocation and interpolation of density matrix
		use omp_lib
		
		type(cmplxMatNN) :: mat
		integer :: n, i, j, k, m, s
		real(dl) :: x, z
		real(dl), dimension(n), intent(in) :: vars
		real(dl), dimension(n), intent(out) :: ydot
!		real(dl), dimension(:), allocatable, save :: tmpv
		character(len=100) :: tmpstr
		logical :: file_exist
		integer :: stat
		real(8) :: timer1

		inquire(file="terminate", exist=file_exist)
		if (file_exist) then
			open(unit=1234, iostat=stat, file="terminate", status='old')
			if (stat == 0) close(1234, status='delete')
			call criticalError("Termination request received")
		end if

#ifdef TESTSPEED
		if (deriv_counter.eq.0) call tic(timer1)
		if (deriv_counter.eq.1000) call toc(timer1, "1000 derivatives")
#endif
		deriv_counter = deriv_counter+1
		write (tmpstr, "('[eq] Calling derivatives (',"//dblfmt//",' x=',"//dblfmt//",')')") deriv_counter,x
		call printVerbose(trim(tmpstr),2)

		call allocateCmplxMat(mat)
		call init_interpNuDens
		z = vars(n)
		call vec_2_densMat(vars(1:n-1))
		do m=1, Ny
			call derivative(x, z, m, mat, n, ydot)
		end do
		
		call dz_o_dx(x,z, ydot, ntot)
		
		call densMat_2_vec(nuDensVec)
	end subroutine derivatives
	
	subroutine jdum
		!necessary but not needed
	end subroutine jdum

	function Neff_from_rho_z(z)
		real(dl) :: Neff_from_rho_z
		real(dl), intent(in) :: z
		real(dl) :: rhototnu, ndeq, totf
		integer :: ix

		totf=0.d0
		do ix=1, flavorNumber
			totf = totf + nuFactor(ix)
		end do
		ndeq=nuDensityLinEq(z)
		rhototnu = (allNuDensity(z) - totf*ndeq)/ndeq
		Neff_from_rho_z = (zid/z)**4 * &
			(3.d0 + rhototnu)
	end function Neff_from_rho_z
end module ndEquations
