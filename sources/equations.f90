module fpEquations
	use precision
	use variables
	use utilities
	use constants
	use bspline_module
	use linear_interpolation_module
	use fpInterpolate
	use diagonalize
	use fpErrors
	use fpInteractions
	use ftqed
	use fpCosmology
	use fpInterfaces1
	use fpMatter
	use fpOutput
	use sgTestUtils
	implicit none

#ifndef NO_INTERPOLATION
	type(bspline_1d) :: dzodx_A_interp, dzodx_B_interp
	type(bspline_1d) :: dwodx_A_interp, dwodx_B_interp
	type(bspline_1d) :: dzodx_eq_interp
#endif

	contains

#ifndef NO_INTERPOLATION
	pure subroutine do_ABweq(ix, elContr0, num, den, numw, denw)
		integer, intent(in) :: ix
		real(dl), dimension(2), intent(in) :: elContr0
		real(dl), intent(out) :: num, den, numw, denw
		real(dl), dimension(2) :: g12, fContr
		real(dl) :: xoz
		integer :: j

		xoz = interp_xozvec(ix)
		g12 = G12_funcFull(xoz, 1.d0)
		num = g12(1)
		den = PISQ/7.5d0 + g12(2)
		do j=1, fermions_number
			fContr = fermions(j)%dzodx_terms(xoz)
			num = num + fContr(1)
			den = den + fContr(2)
		end do
		numw = g12(1) + elContr0(1)
		denw = PISQ/7.5d0 + g12(2) + elContr0(2)
#ifndef NO_MUONS
		do j=2, fermions_number
			fContr = fermions(j)%dzodx_terms(xoz)
			numw = numw + fContr(1)
			denw = denw + fContr(2)
		end do
#endif
	end subroutine do_ABweq

	subroutine init_interp_jkyg12
		real(dl) :: num, den, numw, denw
		real(dl), dimension(:), allocatable :: A, B
		real(dl), dimension(:), allocatable :: Aw, Bw
		real(dl), dimension(:), allocatable :: eq
		real(dl), dimension(2) :: elContr0
		integer :: ix, nx, iflag
		character(len=300) :: tmpstr
		integer, parameter :: uid = 8324
		logical :: exists

		call addToLog("[equations] Initializing interpolation for coefficients in dz/dx and dw/dx...")
		nx=interp_nxz
		allocate(A(nx), B(nx))
		allocate(Aw(nx), Bw(nx))
		allocate(eq(nx))
		elContr0 = electrons%dzodx_terms(interp_xozvec(1))
		if (ftqed_log_term) then
			call criticalError("log term FTQED corrections not supported when interpolating")
		end if
#ifdef NO_MUONS
		write(tmpstr, "(A,'ftqed_nm_',L,'_',L,'.dat')") trim(get_interpolation_folder()), ftqed_temperature_corr, ftqed_ord3
#else
		write(tmpstr, "(A,'ftqed_',L,'_',L,'.dat')") trim(get_interpolation_folder()), ftqed_temperature_corr, ftqed_ord3
#endif
		inquire(file=trim(tmpstr), exist=exists)
		if (exists) then
			call addToLog("[equations] read values from existing file: "//trim(tmpstr))
			open(file=trim(tmpstr), unit=uid, form="unformatted")
			do ix=1, nx
				read(uid) A(ix), B(ix), Aw(ix), Bw(ix), eq(ix)
			end do
			close(uid)
			call addToLog("[equations] check if few saved values are correct: ")
			ix=123
			call do_ABweq(ix, elContr0, num, den, numw, denw)
			call assert_double_rel_safe("check dz/dx coefficient A 1", A(ix), num / den, 1d-7, 1d-6)
			call assert_double_rel_safe("check dz/dx coefficient B 1", B(ix), 1./(2.d0*PISQ*den), 1d-7, 1d-6)
			call assert_double_rel_safe("check dz/dx coefficient Aw 1", Aw(ix), numw / denw, 1d-7, 1d-6)
			call assert_double_rel_safe("check dz/dx coefficient Bw 1", Bw(ix), 1./(2.d0*PISQ*denw), 1d-7, 1d-6)
			call assert_double_rel_safe("check dz/dx coefficient eq 1", eq(ix), num / (den + PISQ/7.5d0*0.875d0*tot_factor_active_nu), 1d-7, 1d-6)
			ix=1422
			call do_ABweq(ix, elContr0, num, den, numw, denw)
			call assert_double_rel_safe("check dz/dx coefficient A 2", A(ix), num / den, 1d-7, 1d-6)
			call assert_double_rel_safe("check dz/dx coefficient B 2", B(ix), 1./(2.d0*PISQ*den), 1d-7, 1d-6)
			call assert_double_rel_safe("check dz/dx coefficient Aw 2", Aw(ix), numw / denw, 1d-7, 1d-6)
			call assert_double_rel_safe("check dz/dx coefficient Bw 2", Bw(ix), 1./(2.d0*PISQ*denw), 1d-7, 1d-6)
			call assert_double_rel_safe("check dz/dx coefficient eq 2", eq(ix), num / (den + PISQ/7.5d0*0.875d0*tot_factor_active_nu), 1d-7, 1d-6)
			call addToLog("everything works!")
		else
			!$omp parallel do default(private) shared(A, B, Aw, Bw, eq, nx, interp_xozvec, fermions, elContr0, tot_factor_active_nu) schedule(dynamic)
			do ix=1, nx
				call do_ABweq(ix, elContr0, num, den, numw, denw)
				A(ix) = num / den
				B(ix) = 1.d0/(2.d0*PISQ*den)
				Aw(ix) = numw / denw
				Bw(ix) = 1.d0/(2.d0*PISQ*denw)
				eq(ix) = num / (den + PISQ/7.5d0*0.875d0*tot_factor_active_nu)
			end do
			!$omp end parallel do
			open(file=trim(tmpstr), unit=uid, status="unknown", form="unformatted")
			do ix=1, nx
				write(uid) A(ix), B(ix), Aw(ix), Bw(ix), eq(ix)
			end do
			close(uid)
			call addToLog("[equations] values saved to file: "//trim(tmpstr))
		end if
		call dzodx_A_interp%initialize(interp_xozvec, A, 4, iflag)
		call dzodx_B_interp%initialize(interp_xozvec, B, 4, iflag)
		call dwodx_A_interp%initialize(interp_xozvec, Aw, 4, iflag)
		call dwodx_B_interp%initialize(interp_xozvec, Bw, 4, iflag)
		call dzodx_eq_interp%initialize(interp_xozvec, eq, 4, iflag)

		deallocate(A, B)
		deallocate(Aw, Bw)
		deallocate(eq)
		call addToLog("[equations] ...done!")
	end subroutine init_interp_jkyg12

	function dzodxcoef_interp_func(o)
		real(dl), dimension(2) :: dzodxcoef_interp_func
		real(dl), intent(in) :: o
		integer :: iflag
		call dzodx_A_interp%evaluate(o, 0, dzodxcoef_interp_func(1), iflag)
		call dzodx_B_interp%evaluate(o, 0, dzodxcoef_interp_func(2), iflag)
	end function dzodxcoef_interp_func

	function dwodxcoef_interp_func(o)
		real(dl), dimension(2) :: dwodxcoef_interp_func
		real(dl), intent(in) :: o
		integer :: iflag
		call dwodx_A_interp%evaluate(o, 0, dwodxcoef_interp_func(1), iflag)
		call dwodx_B_interp%evaluate(o, 0, dwodxcoef_interp_func(2), iflag)
	end function dwodxcoef_interp_func
#endif

	subroutine dz_o_dx(x, w, z, ydot, n)
		!eq 17 from doi:10.1016/S0370-2693(02)01622-2
		!Newton-Cotes integral without need of interpolation
		real(dl), intent(in) :: w, x, z
		integer, intent(in) :: n
		real(dl), dimension(n), intent(inout) :: ydot
		real(dl) :: nudrho
		real(dl), dimension(2) :: coeffs
		integer :: ix, m
#ifdef NO_INTERPOLATION
		integer :: j
		real(dl) :: num, den, xoz, numw, denw
		real(dl), dimension(2) :: g12, fContr, elContr0

		elContr0 = electrons%dzodx_terms(interp_xozvec(1))
#endif

		if (use_gauss_laguerre) then
			!$omp parallel do shared(fy_arr) private(nudrho, m, ix)
			do m=1, Ny
				nudrho = 0.d0
				do ix=1, flavorNumber
					nudrho = nudrho + ydot((m-1)*flavNumSqu + ix) * nuFactor(ix)
				end do
				fy_arr(m) = nudrho * feq_arr(m)
			end do
			!$omp end parallel do
			nudrho = integral_GL_1d(w_gl_arr, fy_arr)
		else
			!$omp parallel do shared(fy_arr) private(nudrho, m, ix)
			do m=1, Ny
				nudrho = 0.d0
				do ix=1, flavorNumber
					nudrho = nudrho + ydot((m-1)*flavNumSqu + ix) * nuFactor(ix)
				end do
				fy_arr(m) = y_arr(m)**3 * nudrho * feq_arr(m)
			end do
			!$omp end parallel do
			nudrho = integral_NC_1d(Ny, dy_arr, fy_arr)
		end if

#ifdef NO_INTERPOLATION
		xoz = x/z
		g12 = G12_funcFull(x, z)
		num = g12(1)
		den = PISQ/7.5d0 + g12(2)
		do j=1, fermions_number
			fContr = fermions(j)%dzodx_terms(xoz)
			num = num + fContr(1)
			den = den + fContr(2)
		end do
		numw = g12(1) + elContr0(1)
		denw = PISQ/7.5d0 + g12(2) + elContr0(2)
#ifndef NO_MUONS
		do j=2, fermions_number
			fContr = fermions(j)%dzodx_terms(xoz)
			numw = numw + fContr(1)
			denw = denw + fContr(2)
		end do
#endif
		!neutrino effective temperature
		ydot(n-1) = numw / denw - nudrho / z**3 / (2.d0*PISQ*denw)
		!photon temperature
		ydot(n) = num/den - nudrho / z**3 / (2.d0*PISQ*den)
#else
		!neutrino effective temperature
		coeffs = dwodxcoef_interp_func(x/z)
		ydot(n-1) = coeffs(1) - coeffs(2) * nudrho / z**3
		!photon temperature
		coeffs = dzodxcoef_interp_func(x/z)
		ydot(n) = coeffs(1) - coeffs(2) * nudrho / z**3
#endif
	end subroutine dz_o_dx

	subroutine dz_o_dx_eq(n, x, vars, ydot)
		!eq 17 from doi:10.1016/S0370-2693(02)01622-2
		! at equilibrium (i.e. no contribution from neutrino distortions), needed for initial z
		integer, intent(in) :: n
		real(dl), intent(in) :: x
		real(dl), dimension(n), intent(in) :: vars
		real(dl), dimension(n), intent(out) :: ydot
		real(dl) :: z, xoz, coeff
#ifdef NO_INTERPOLATION
		integer :: j
		real(dl) :: num, den
		real(dl), dimension(2) :: g12, fContr
#else
		integer :: iflag
#endif

		z = vars(n)+1.d0
		xoz=x/z

#ifdef NO_INTERPOLATION
		g12 = G12_funcFull(x, z)
		num = g12(1)
		den = PISQ/7.5d0 + g12(2)
		do j=1, fermions_number
			fContr = fermions(j)%dzodx_terms(xoz)
			num = num + fContr(1)
			den = den + fContr(2)
		end do
		coeff = num / (den + PISQ/7.5d0*0.875d0*tot_factor_active_nu)
#else
		call dzodx_eq_interp%evaluate(xoz, 0, coeff, iflag)
#endif
		ydot(n) = coeff
	end subroutine dz_o_dx_eq

	subroutine zin_solver
		integer :: n
		character(len=3) :: istchar
		character(len=100) :: tmpstring
		real(dl) :: rtol, xvh
		integer :: itol, itask, istate, iopt, lrw, liw, jt
		real(dl), dimension(:), allocatable :: rwork, atol, cvec
		integer, dimension(:), allocatable :: iwork

		xvh = very_early_x
		write(tmpstring,"('[zin] Computation of z_in started, x_in=',E10.3,', x_end=',E10.3,'.')") xvh, x_in
		call addToLog(tmpstring)
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
		rwork=0.d0
		iwork=0
		iwork(6)=99999999
		jt=2

		cvec(n) = 0.d0

		call dlsoda(dz_o_dx_eq,n,cvec,xvh,x_in,&
					itol,rtol,atol,itask,istate, &
					iopt,rwork,lrw,iwork,liw,jacobian,jt)

		if (istate.lt.0) then
			write(istchar, "(I3)") istate
			call criticalError('[zin] istate='//istchar)
		end if

		write(tmpstring,"('[zin] ended with z_in-1 =',E16.9,'.')") cvec(n)
		call addToLog(trim(tmpstring))
		z_in = cvec(n) + 1.d0
	end subroutine zin_solver

	pure subroutine drhoy_dx_fullMat(matrix, x, w, z, iy, dme2, sqrtraddens, Fint_nue, Fint_nunu, Heff, comm, collterms)
!		compute rho derivatives for a given momentum y_arr(m), save to a matrix
		use fpInterfaces2
		procedure (collision_integrand_nue) :: Fint_nue
		procedure (collision_integrand_nunu) :: Fint_nunu
		real(dl), intent(in) :: x, w, z, dme2, sqrtraddens
		integer, intent(in) :: iy
		type(cmplxMatNN), intent(out) :: matrix
		type(cmplxMatNN), intent(inout) :: Heff, comm, collterms
		real(dl) :: y, overallNorm, cf
		integer :: ix
		type(coll_args) :: collArgs

		y = nuDensMatVecFD(iy)%y

		collArgs%x = x
		collArgs%w = w
		collArgs%z = z
		collArgs%y1 = y
		collArgs%dme2 = dme2
		collArgs%iy = iy

		overallNorm = overallFactor / sqrtraddens
		Heff = H_eff(y)

		!switch imaginary and real parts because of the "-i" factor
		call Commutator(Heff%re, nuDensMatVecFD(iy)%re, comm%im)
		call Commutator(Heff%re, nuDensMatVecFD(iy)%im, comm%re)

		!matrix is now redefined
		cf = x**2/m_e_cub
		matrix%im = - comm%im * cf
		matrix%re = comm%re * cf

		collterms = get_collision_terms(collArgs, Fint_nue, Fint_nunu)
		matrix%re = matrix%re + collterms%re
		matrix%im = matrix%im + collterms%im

		matrix%re = matrix%re * overallNorm
		matrix%im = matrix%im * overallNorm
		do ix=1, flavorNumber
			matrix%re(ix,ix) = matrix%re(ix,ix) / feq_arr(iy)
			matrix%im(ix,ix) = 0.d0
		end do
	end subroutine drhoy_dx_fullMat

	pure subroutine drho_y_dx(x, w, z, m, dme2, sqrtraddens, n, ydot, Heff, comm, collterms)
!		compute rho derivatives for a given momentum y_arr(m), save to ydot
		real(dl), intent(in) :: x, w, z, dme2, sqrtraddens
		integer, intent(in) :: m, n
		real(dl), dimension(n), intent(out) :: ydot
		type(cmplxMatNN), intent(inout) :: Heff, comm, collterms
		integer :: i, j, k
		type(cmplxMatNN) :: mat

		call drhoy_dx_fullMat(&
			mat, x, w, z, m, dme2, sqrtraddens, &
			coll_nue_int, coll_nunu_int, &
			Heff, comm, collterms &
		)
		do i=1, flavorNumber
			ydot(i) = mat%re(i,i)
		end do
		k=flavorNumber+1
		if (has_offdiagonal()) then
			do i=1, flavorNumber-1
				do j=i+1, flavorNumber
					ydot(k) = mat%re(i,j)
					ydot(k+1) = mat%im(i,j)
					k=k+2
				end do
			end do
		end if
	end subroutine drho_y_dx

	subroutine derivatives(n, x, vars, ydot)
!		compute all the rho derivatives (drho/dx for all y, dz/dx)
!		needs allocation and interpolation of density matrix
		integer, intent(in) :: n
		real(dl), intent(in) :: x
		real(dl), dimension(n), intent(in) :: vars
		real(dl), dimension(n), intent(out) :: ydot
		real(dl) :: w, z, dme2, sqrtraddens
		character(len=100) :: tmpstr
		logical :: file_exist
		integer :: tstat, m, s
		real(8) :: timer1
		real(dl), dimension(:), allocatable :: tmpvec

		inquire(file="terminate", exist=file_exist)
		if (file_exist) then
			open(unit=1234, iostat=tstat, file="terminate", status='old')
			if (tstat == 0) close(1234, status='delete')
			call criticalError("Termination request received")
		end if

#ifdef TESTSPEED
		if (deriv_counter.eq.0) call tic(timer1)
		if (deriv_counter.eq.1000) call toc(timer1, "1000 derivatives")
#endif
		deriv_counter = deriv_counter+1
		write (tmpstr, "('[eq] Calling derivatives (',"//dblfmt//",' x=',"//dblfmt//",')')") deriv_counter,x
		call printVerbose(trim(tmpstr), 1+int(mod(deriv_counter, Nprintderivs)))

		w = vars(n-1) + 1.d0
		z = vars(n) + 1.d0
		call vec_2_densMat(vars)

#ifdef NO_INTERPOLATION
		dme2 = dme2_electronFull(x, 0.d0, z, .false.)
#else
		dme2 = dme2_nolog(x, z)
#endif
		sqrtraddens = sqrt(totalRadiationDensity(x, z))
		call updateMatterDensities(x, z)

		!$omp parallel shared(ydot, x, z, dme2, sqrtraddens, flavNumSqu, intermediateSteps) private(m, s, tmpvec)
		allocate(tmpvec(flavNumSqu))
		tmpvec = 0
		!$omp do schedule(static)
		do m=1, Ny
			call drho_y_dx( &
				x, w, z, m, dme2, sqrtraddens, flavNumSqu, &
				tmpvec, &
				intermediateSteps%Heff(m), intermediateSteps%commutator(m), intermediateSteps%collterms(m) &
			)
			s=(m-1)*flavNumSqu
			ydot(s+1:s+flavNumSqu) = tmpvec(:)
		end do
		!$omp end do
		deallocate(tmpvec)
		!$omp end parallel

		call dz_o_dx(x, w, z, ydot, ntot)

		call densMat_2_vec(nuDensVec)

		!optionally save all output params
		if (intermediateSteps%output) then
			intermediateSteps%x = x
			intermediateSteps%norm = overallFactor / sqrtraddens
			intermediateSteps%yvec = vars
			intermediateSteps%ydot = ydot
			! the other quantities should have been saved previously
			call saveIntermediateSteps
		end if
	end subroutine derivatives

	subroutine solver
		real(dl) :: xstart, xend, xchk, xs
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
		integer :: m, i, j, k

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
		atol = dlsoda_atol_z
		k = 1
		do m = 1, Ny
			do i = 1, flavorNumber
				atol(k+i-1) = dlsoda_atol_d
			end do
			k = k + flavorNumber
			if (has_offdiagonal()) then
				do i = 1, flavorNumber-1
					do j = i+1, flavorNumber
						atol(k) = dlsoda_atol_o
						atol(k+1) = dlsoda_atol_o
						k=k+2
					end do
				end do
			end if
		end do
		atol(ntot-1) = dlsoda_atol_z
		atol(ntot) = dlsoda_atol_z
		rwork=0.d0
		iwork=0
		iwork(6)=99999999
		jt=2

		call densMat_2_vec(nuDensVec)
		call allocateStoreVars

		nuDensVec(ntot-1) = z_in - 1.d0 !neutrino temperature start at same value as photon temperature
		nuDensVec(ntot) = z_in - 1.d0

		call readCheckpoints(nchk, xchk, ychk, chk)

		if (chk .and. &
			nchk.eq.ntot) then
			xstart=xchk
			nuDensVec=ychk
			firstWrite=.false.
			firstPoint=.true.
			ix_in=1 + int((log10(xchk)-logx_in)/(logx_fin-logx_in)*(Nx-1))
			write(tmpstring,"('ntot =',I4,' - x =',"//dblfmt//",' (i=',I4,') - w =',"//dblfmt//",' - z =',"//dblfmt//")") &
				nchk, xchk, ix_in, ychk(ntot-1)+1.d0, ychk(ntot)+1.d0
			call addToLog("[ckpt] ##### Checkpoint file found. Will start from there. #####")
			call addToLog(trim(tmpstring))
		else
			xstart=x_arr(1)
			ix_in=1
			call saveRelevantInfo(xstart, nuDensVec)
		end if

		do ix=ix_in+1, Nx
			xs = xstart
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
						iopt,rwork,lrw,iwork,liw,jacobian,jt)

			if (istate.lt.0) then
				write(istchar, "(I3)") istate
				call criticalError('istate='//istchar)
			end if
			call writeCheckpoints(ntot, xend, nuDensVec)
			call saveRelevantInfo(xend, nuDensVec)
			xstart=xend
		end do
		write(tmpstring,"('x_end =',"//dblfmt//",' - w_end =',"//dblfmt//",' - z_end =',"//dblfmt//")") &
			xend, nuDensVec(ntot-1)+1.d0, nuDensVec(ntot)+1.d0

		call date_and_time(VALUES=values)
		call openFile(timefileu, trim(outputFolder)//timefilen, .false.)
		write(timefileu, &
			'("-- ",I0.2,"/",I0.2,"/",I4," - h",I2,":",I0.2,":",I0.2,' &
			//"' - DLSODA end after ',"//dblfmt//",' derivatives - " &
			//"xend =',"//dblfmt//",' - w =',"//dblfmt//",' - z =',"//dblfmt//")") &
			values(3), values(2), values(1), values(5),values(6),values(7), deriv_counter, xend, &
			nuDensVec(ntot-1)+1.d0, nuDensVec(ntot)+1.d0
		close(timefileu)

		call deleteCheckpoints
		call deallocateStoreVars

		call addToLog("[solver] Solver ended. "//trim(tmpstring))
	end subroutine solver

	subroutine jacobian
		!optional for dlsoda, needed only for time optimization (would gain a factor ~Ny/5)
	end subroutine jacobian

end module fpEquations
