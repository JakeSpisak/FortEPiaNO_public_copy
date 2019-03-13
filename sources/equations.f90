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
		logical :: f5
	end type collTerm
	type(collTerm) :: lastColl

	real(dl), parameter :: upper = 1.d2

	real(dl) :: deriv_counter

	contains

	subroutine updateLeptonDensities(x,y,z)
		real(dl), intent(in) :: x,y,z
		real(dl) :: ldf

		leptonDensities=0.d0
		ldf = leptDensFactor * y / x**6
		leptonDensities(1,1) = ldf * electronDensity(x,z)
		leptonDensities(2,2) = ldf * muonDensity(x,z)
	end subroutine updateLeptonDensities

	subroutine densMat_2_vec(vec)
		real(dL), dimension(:), intent(out) :: vec
		integer :: i,j,k,m

		k=1
		do m=1, Ny
			do i=1, flavorNumber
				vec(k+i-1) = nuDensMatVec(m)%re(i,i)
			end do
			k=k+flavorNumber
			if (collision_offdiag.ne.0 .and. collision_offdiag.ne.3) then
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
		real(dL) :: fd

		k=1
		do m=1, Ny
			do i=1, flavorNumber
				nuDensMatVec(m)%re(i,i) = vec(k+i-1)
				nuDensMatVec(m)%im(i,i) = 0.d0
			end do
			k=k+flavorNumber
			if (collision_offdiag.ne.0 .and. collision_offdiag.ne.3) then
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
			nuDensMatVecFD(m)%re = nuDensMatVec(m)%re
			nuDensMatVecFD(m)%im = nuDensMatVec(m)%im
			fd = fermiDirac(y_arr(m))!/z_in)
			do i=1, flavorNumber
				nuDensMatVecFD(m)%re(i,i) = (1.d0 + nuDensMatVec(m)%re(i,i)) * fd
			end do
		end do
	end subroutine vec_2_densMat

	pure function J_int(o, u)
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

	pure function Jprime_int(o, u)
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

	pure function K_int(o, u)
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

	pure function Kprime_int(o, u)
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

	pure function Y_int(o, u)
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
		nx=interp_nxz
		allocate(A(nx),B(nx))
		!$omp parallel do default(shared) private(ix, j, y, g12, num, den) schedule(dynamic)
		do ix=1, nx
			j      = J_funcFull(interp_xozvec(ix))
			y      = Y_funcFull(interp_xozvec(ix))
			g12    = G12_funcFull(interp_xozvec(ix))
			num= interp_xozvec(ix) * j + g12(1)
			den= interp_xozvec(ix)**2 * j + y + PISQ/7.5d0 + g12(2)
			A(ix) = num / den
			B(ix) = 1./(2.d0*PISQ*den)
		end do
		!$omp end parallel do
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
			fy_arr(m) = y_arr(m)**3 * tmp * fermiDirac(y_arr(m))!/z_in)
		end do
		tmp = integral_linearized_1d(Ny, dy_arr, fy_arr)
		coeffs = dzodxcoef_interp_func(x/z)
		ydot(n) = coeffs(1) - coeffs(2) * tmp/ z**3
	end subroutine dz_o_dx_lin

	subroutine dz_o_dx_eq_lin(n, x, vars, ydot)
		integer, intent(in) :: n
		real(dl), intent(in) :: x
		real(dl), dimension(n), intent(in) :: vars
		real(dl), dimension(n), intent(out) :: ydot
		!eq 17 from doi:10.1016/S0370-2693(02)01622-2
		!integral without need of interpolation, just use linear approx in y bins
		real(dl) :: z, j, y, xoz, num, den
		real(dl), dimension(2) :: coeffs
		integer :: ix, m
		real(dl), dimension(2) :: g12

		z = vars(n)
		xoz=x/z

		j      = J_funcFull(xoz)
		y      = Y_funcFull(xoz)
		g12    = G12_funcFull(xoz)
		num= xoz * j + g12(1)
		den= xoz**2 * j + y + PISQ/7.5d0 + g12(2)

		ydot(n) = num/den
	end subroutine dz_o_dx_eq_lin

	subroutine zin_solver
		integer :: n
		character(len=3) :: istchar
		character(len=100) :: tmpstring
		real(dl) :: rtol, xvh
		integer :: itol, itask, istate, iopt, lrw, liw, jt
		real(dl), dimension(:), allocatable :: rwork, atol, cvec
		integer, dimension(:), allocatable :: iwork

		call addToLog("[zin] Computation of z_in started. ")
		n=1

		itol=2
		rtol=dlsoda_rtol
		itask=1
		istate=1
		iopt=1

		lrw=60
		liw=60
		allocate(atol(n), cvec(n), rwork(lrw), iwork(liw))
		atol=1d-7
		rwork=0.
		iwork=0
		iwork(6)=99999999
		jt=2

		cvec(1) = 1.d0
		xvh = 0.0001d0

		call dlsoda(dz_o_dx_eq_lin,n,cvec,xvh,x_in,&
					itol,rtol,atol,itask,istate, &
					iopt,rwork,lrw,iwork,liw,jdum,jt)

		if (istate.lt.0) then
			write(istchar, "(I3)") istate
			call criticalError('[zin] istate='//istchar)
		end if

		write(tmpstring,"('[zin] ended with z_in =',E16.9,'.')") cvec(n)
		z_in = cvec(n)
		call addToLog(trim(tmpstring))
	end subroutine zin_solver

	subroutine drhoy_dx_fullMat(matrix,x,z,iy, Fre, Fim)
		interface
			pure real(dl) function Fre(a, b, o)
				use variables
				integer, intent(in) :: a
				real(dl), intent(in) :: b
				type(coll_args), intent(in) :: o
			end function
			pure real(dl) function Fim(a, b, o)
				use variables
				integer, intent(in) :: a
				real(dl), intent(in) :: b
				type(coll_args), intent(in) :: o
			end function
		end interface
		type(cmplxMatNN), intent(out) :: matrix
		real(dl) :: x,y,z, overallNorm, a, fd, cf, ldf
		integer :: iy,k, ix
		type(coll_args) :: collArgs
		type(cmplxMatNN), save :: comm

		call allocateCmplxMat(comm)
		call allocateCmplxMat(matrix)

		y = nuDensMatVecFD(iy)%y
		fd = fermiDirac(y)!/z_in)

		collArgs%x = x
		collArgs%z = z
		collArgs%y1 = y
		collArgs%dme2 = dme2_electron(x, 0.d0, z)
		collArgs%iy = iy

		overallNorm = overallFactor / sqrt(radDensity(x,z))
		call updateLeptonDensities(x,y,z)
		matrix%re(:,:) = nuMassesMat(:,:)/(2*y) + leptonDensities(:,:)

		!switch imaginary and real parts because of the "-i" factor
		call Commutator(matrix%re, nuDensMatVecFD(iy)%re, comm%im)
		call Commutator(matrix%re, nuDensMatVecFD(iy)%im, comm%re)

		!matrix is now redefined
		cf = x**2/m_e_cub
		matrix%im = - comm%im * cf
		matrix%re = comm%re * cf

		if (x .ne. lastColl%x .or. y .ne. lastColl%y .or. z .ne. lastColl%z .or. lastColl%f5) then
			lastColl%mat = get_collision_terms(collArgs, Fre, Fim)
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
			matrix%im(ix,ix) = 0.d0
		end do
	end subroutine drhoy_dx_fullMat

	subroutine saveRelevantInfo(x, vec)
		real(dl) :: x
		real(dl), dimension(:) :: vec
		real(dl), dimension(:), allocatable :: tmpvec
		integer :: k, i, j, iy, totFiles
		integer, dimension(:), allocatable :: units
		character(len=200) :: fname

		totFiles=flavNumSqu+2
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
				tmpvec(iy)=nuDensMatVec(iy)%re(k,k)
			end do
			write(units(k), '(*('//dblfmt//'))') x, tmpvec
		end do
		if (collision_offdiag.ne.0 .and. collision_offdiag.ne.3) then
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
		call openFile(units(k+1), trim(outputFolder)//'/Neff.dat', firstWrite)
		write(units(k+1), '(*('//dblfmt//'))') x, Neff_from_rho_z(vec(ntot))

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
		integer,dimension(8) :: values
		integer, parameter :: timefileu = 8970
		character(len=*), parameter :: timefilen = '/time.log'

		deriv_counter = 0

		call openFile(timefileu, trim(outputFolder)//timefilen,.true.)
		write(timefileu,*) "starting solver"
		close(timefileu)

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

		nuDensVec(ntot)=z_in
		call densMat_2_vec(nuDensVec)

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
			write(tmpstring,"('x_start =',"//dblfmt//",' - x_end =',"//dblfmt//")") xstart, xend
			call addToLog("[solver] Start DLSODA..."//trim(tmpstring))

			call date_and_time(VALUES=values)
			call openFile(timefileu, trim(outputFolder)//timefilen, .false.)
			write(timefileu, &
				'("-- ",I0.2,"/",I0.2,"/",I4," - h",I2,":",I0.2,":",I0.2,'&
				//"' - DLSODA x_start =',"//dblfmt//",' - x_end =',"//dblfmt//")") &
				values(3), values(2), values(1), values(5),values(6),values(7), xstart, xend
			close(timefileu)

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
		write(tmpstring,"('x_end =',"//dblfmt//",' - z_end =',"//dblfmt//")") xend, nuDensVec(ntot)

		call date_and_time(VALUES=values)
		call openFile(timefileu, trim(outputFolder)//timefilen, .false.)
		write(timefileu, &
			'("-- ",I0.2,"/",I0.2,"/",I4," - h",I2,":",I0.2,":",I0.2,' &
			//"' - DLSODA end after ',"//dblfmt//",' derivatives - " &
			//"xend =',"//dblfmt//",' - T =',"//dblfmt//")") &
			values(3), values(2), values(1), values(5),values(6),values(7), deriv_counter, xend, nuDensVec(ntot)
		close(timefileu)

		call addToLog("[solver] Solver ended. "//trim(tmpstring))
	end subroutine solver

	subroutine derivative (x, z, m, mat, n, ydot)
!		compute rho derivatives for a given momentum y_arr(m), save to ydot
		real(dl), intent(in) :: x, z
		integer, intent(in) :: m, n
		real(dl), dimension(n), intent(out) :: ydot
		integer :: i, j, k, s
		type(cmplxMatNN) :: mat

		call drhoy_dx_fullMat(mat, x,z,m, coll_nue_3_int_re, coll_nue_3_int_im)
		s=(m-1)*flavNumSqu
		do i=1, flavorNumber
			ydot(s+i) = mat%re(i,i)
		end do
		k=flavorNumber+1
		if (collision_offdiag.ne.0 .and. collision_offdiag.ne.3) then
			do i=1, flavorNumber-1
				do j=i+1, flavorNumber
					ydot(s+k) = mat%re(i,j)
					ydot(s+k+1) = mat%im(i,j)
					k=k+2
				end do
			end do
		end if
	end subroutine derivative

	subroutine derivatives(n, x, vars, ydot)
!		compute all the rho derivatives (drho/dx for all y, dz/dx)
!		needs allocation and interpolation of density matrix
		use omp_lib

		type(cmplxMatNN) :: mat
		integer :: n, m
		real(dl) :: x, z
		real(dl), dimension(n), intent(in) :: vars
		real(dl), dimension(n), intent(out) :: ydot
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
		call printVerbose(trim(tmpstr), 1+int(mod(deriv_counter, Nprintderivs)))

!		write(*,multidblfmt) nuDensMatVec
!		write(*,multidblfmt) nuDensMatVecFD
		call allocateCmplxMat(mat)
		z = vars(n)
		call vec_2_densMat(vars)
		do m=1, Ny
			call derivative(x, z, m, mat, n, ydot)
		end do

!		write(*,multidblfmt) ydot
		call dz_o_dx_lin(x,z, ydot, ntot)

		call densMat_2_vec(nuDensVec)
	end subroutine derivatives

	subroutine jdum
		!necessary for dlsoda but not needed
	end subroutine jdum

	function Neff_from_rho_z(z)
		real(dl) :: Neff_from_rho_z
		real(dl), intent(in) :: z
		real(dl) :: rhototnu, ndeq, totf
		integer :: ix

		totf=0.d0
		do ix=1, flavorNumber
			if(.not.sterile(ix)) totf = totf + nuFactor(ix)
		end do
		ndeq=nuDensityLinEq(z)
		rhototnu = (allNuDensity(z) - totf*ndeq)/ndeq
		Neff_from_rho_z = (zid/z)**4 * &
			(3.d0 + rhototnu)
	end function Neff_from_rho_z

	subroutine finalresults
		real(dl) :: ndeq, tmp
		integer :: ix

		call openFile(9876, trim(outputFolder)//"/resume.dat", .true.)
		write(*,"('final z = ',F11.8)") nuDensVec(ntot)
		write(9876,"('final z = ',F11.8)") nuDensVec(ntot)
		ndeq=nuDensityLinEq(nuDensVec(ntot))
		do ix=1, flavorNumber
			tmp = (nuDensityLin(nuDensVec(ntot), ix) - ndeq)*nuFactor(ix)/ndeq
			write(*,"('dRho_',I1,'  = ',F9.6)") ix, tmp
			write(9876,"('dRho_',I1,'  = ',F9.6)") ix, tmp
		end do
		tmp = Neff_from_rho_z(nuDensVec(ntot))
		write(*,"('Neff    = ',F9.6)") tmp
		write(9876,"('Neff    = ',F9.6)") tmp
		close(9876)
	end subroutine finalresults
end module ndEquations
