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
	use diagonalize
	use ndInterfaces1
	implicit none

	type(bspline_1d) :: dzodx_A_interp, dzodx_B_interp

	real(dl), parameter :: upper = 1.d2

	real(dl) :: deriv_counter

	contains

	subroutine updateMatterDensities(x, z)
		real(dl), intent(in) :: x, z
		real(dl) :: ldf
		integer :: ix, iy
		procedure (nuDensity_integrator), pointer :: nuDensityInt

		if (use_gauss_laguerre) then
			nuDensityInt => nuDensityGL
		else
			nuDensityInt => nuDensityNC
		end if

		leptonDensities = 0.d0
		ldf = leptDensFactor / x**6
		leptonDensities(1,1) = ldf * electronDensity(x, z)
		leptonDensities(2,2) = ldf * muonDensity(x, z)

		nuDensities%re = 0.d0
		nuDensities%im = 0.d0
		do ix=1, flavorNumber
			if (.not.sterile(ix)) then
				nuDensities%re(ix, ix) = nuDensities%re(ix, ix) + nuDensityInt(ix, ix)
				do iy=ix+1, flavorNumber
					if (.not.sterile(iy)) then
						nuDensities%re(ix, iy) = nuDensities%re(ix, iy) + nuDensityInt(ix, iy)
						nuDensities%im(ix, iy) = nuDensities%im(ix, iy) + nuDensityInt(ix, iy, .false.)
					end if
				end do
			end if
			do iy=ix+1, flavorNumber
				nuDensities%re(iy, ix) = nuDensities%re(ix, iy)
				nuDensities%im(iy, ix) = - nuDensities%im(ix, iy)
			end do
		end do
		nuDensities%re(:,:) = nuDensities%re(:,:) * ldf * (cos2thW)
		nuDensities%im(:,:) = nuDensities%im(:,:) * ldf * (cos2thW)
	end subroutine updateMatterDensities

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
			fd = fermiDirac(y_arr(m))
			do i=1, flavorNumber
				nuDensMatVecFD(m)%re(i,i) = (1.d0 + nuDensMatVec(m)%re(i,i)) * fd
			end do
		end do
	end subroutine vec_2_densMat

	pure function J_int(o, u)
		real(dl) :: J_int, esuo
		real(dl), intent(in) :: o, u

		esuo=exp(sqrt(u*u+o*o))
		J_int = esuo / ((1.d0+esuo)**2)
	end function J_int

	pure function J_funcFull(o)
		real(dl) :: J_funcFull
		real(dl), intent(in) :: o
		integer :: i

		J_funcFull = 0.d0
		do i=1, N_opt_xoz
			J_funcFull = J_funcFull + opt_xoz_w(i)*J_int(o, opt_xoz(i))
		end do
		J_funcFull = J_funcFull / PISQ
	end function J_funcFull

	pure function Jprime_int(o, u)
		real(dl) :: Jprime_int
		real(dl), intent(in) :: o, u
		real(dl) :: uuoo, sqrtuuoo, expsqrtuuoo

		uuoo=u*u+o*o
		sqrtuuoo = sqrt(uuoo)
		expsqrtuuoo = exp(sqrtuuoo)

		Jprime_int = expsqrtuuoo * (1.d0-expsqrtuuoo) /(sqrtuuoo*(expsqrtuuoo+1)**3)
	end function Jprime_int

	pure function JprimeFull(o)
		real(dl) :: JprimeFull
		real(dl), intent(in) :: o
		integer :: i

		JprimeFull = 0.d0
		do i=1, N_opt_xoz
			JprimeFull = JprimeFull + opt_xoz_w(i)*Jprime_int(o, opt_xoz(i))
		end do
		JprimeFull = JprimeFull * o / PISQ
	end function JprimeFull

	pure function K_int(o, u)
		real(dl) :: K_int,suo
		real(dl), intent(in) :: o, u

		suo=sqrt(u*u+o*o)
		K_int = 1.d0 / (suo * (1.d0+exp(suo)))
	end function K_int

	pure function K_funcFull(o)
		real(dl) :: K_funcFull
		real(dl), intent(in) :: o
		integer :: i

		K_funcFull = 0.d0
		do i=1, N_opt_xoz
			K_funcFull = K_funcFull + opt_xoz_w(i)*K_int(o, opt_xoz(i))
		end do
		K_funcFull = K_funcFull / PISQ
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

		KprimeFull = - rombint_re(o, Kprime_int, zero, upper, toler_jkyg, maxiter) * o / PISQ
	end function KprimeFull

	pure function Y_int(o, u)
		real(dl) :: Y_int, esuo
		real(dl), intent(in) :: o, u

		esuo=exp(sqrt(u*u+o*o))
		Y_int = u**2 * esuo / ((1.d0+esuo)**2)
	end function Y_int

	pure function Y_funcFull(o)
		real(dl) :: Y_funcFull
		real(dl), intent(in) :: o
		integer :: i

		Y_funcFull = 0.d0
		do i=1, N_opt_xoz
			Y_funcFull = Y_funcFull + opt_xoz_w(i)*Y_int(o, opt_xoz(i))
		end do
		Y_funcFull = Y_funcFull / PISQ
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
			tmp = (kp/6.d0 - ko*kp + jp/6.d0 + jp*ko + jo*kp)

			G12_funcFull(1) = PIx2*alpha_fine *(&
				(ko/3.d0 + 2.d0*ko*ko - jo/6.d0 - ko*jo)/o + &
				tmp )
			G12_funcFull(2) = PIx2*alpha_fine*( &
				o * tmp &
				- 4.d0*((ko+jo)/6.d0 + ko*jo - ko*ko/2.d0))
		else
			G12_funcFull = 0.d0
		end if
	end function G12_funcFull

	subroutine init_interp_jkyg12
		real(dl) :: num, den, j, y, jmu, ymu
		real(dl), dimension(:), allocatable :: A, B
		real(dl), dimension(2) :: g12
		integer :: ix, nx, iflag

		call addToLog("[equations] Initializing interpolation for coefficients in dz/dx...")
		nx=interp_nxz
		allocate(A(nx), B(nx))
		!$omp parallel do default(shared) private(ix, j, y, g12, num, den) schedule(dynamic)
		do ix=1, nx
			j = J_funcFull(interp_xozvec(ix))
			y = Y_funcFull(interp_xozvec(ix))
			if (interp_xozvec(ix).gt.x_muon_cut) then
				jmu=0.d0
				ymu=0.d0
			else
				jmu = J_funcFull(interp_xozvec(ix)*m_mu_o_m_e)
				ymu = Y_funcFull(interp_xozvec(ix)*m_mu_o_m_e)
			end if
			g12 = G12_funcFull(interp_xozvec(ix))
			num= interp_xozvec(ix)*(j+jmu) + g12(1)
			den= interp_xozvec(ix)**2*(j+jmu) + y+ymu + PISQ/7.5d0 + g12(2)
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

	subroutine dz_o_dx(x, w, z, ydot, n)
		!eq 17 from doi:10.1016/S0370-2693(02)01622-2
		!Newton-Cotes integral without need of interpolation
		real(dl), intent(in) :: w, x, z
		integer, intent(in) :: n
		real(dl), dimension(n), intent(inout) :: ydot
		real(dl) :: tmp
		real(dl), dimension(2) :: coeffs
		integer :: ix, m

		if (use_gauss_laguerre) then
			!$omp parallel do shared(fy_arr) private(tmp, m, ix)
			do m=1, Ny
				tmp = 0.d0
				do ix=1, flavorNumber
					tmp = tmp + ydot((m-1)*flavNumSqu + ix) * nuFactor(ix)
				end do
				fy_arr(m) = tmp * fermiDirac(y_arr(m))
			end do
			!$omp end parallel do
			tmp = integral_GL_1d(w_gl_arr, fy_arr)
		else
			!$omp parallel do shared(fy_arr) private(tmp, m, ix)
			do m=1, Ny
				tmp = 0.d0
				do ix=1, flavorNumber
					tmp = tmp + ydot((m-1)*flavNumSqu + ix) * nuFactor(ix)
				end do
				fy_arr(m) = y_arr(m)**3 * tmp * fermiDirac(y_arr(m))
			end do
			!$omp end parallel do
			tmp = integral_NC_1d(Ny, dy_arr, fy_arr)
		end if
		ydot(n-1) = coeff_dw_dx * tmp / w**3 / tot_factor_active_nu
		coeffs = dzodxcoef_interp_func(x/z)
		ydot(n) = coeffs(1) - coeffs(2) * tmp / z**3
	end subroutine dz_o_dx

	subroutine dz_o_dx_eq(n, x, vars, ydot)
		!eq 17 from doi:10.1016/S0370-2693(02)01622-2
		! at equilibrium (i.e. no contribution from neutrino distortions), needed for initial z
		integer, intent(in) :: n
		real(dl), intent(in) :: x
		real(dl), dimension(n), intent(in) :: vars
		real(dl), dimension(n), intent(out) :: ydot
		real(dl) :: z, xoz
		real(dl), dimension(2) :: coeffs

		z = vars(n)+1.d0
		xoz=x/z
		coeffs = dzodxcoef_interp_func(x/z)
		ydot(n) = coeffs(1)
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
		rwork=0.
		iwork=0
		iwork(6)=99999999
		jt=2

		cvec(n) = 0.d0

		call dlsoda(dz_o_dx_eq,n,cvec,xvh,x_in,&
					itol,rtol,atol,itask,istate, &
					iopt,rwork,lrw,iwork,liw,jdum,jt)

		if (istate.lt.0) then
			write(istchar, "(I3)") istate
			call criticalError('[zin] istate='//istchar)
		end if

		write(tmpstring,"('[zin] ended with z_in-1 =',E16.9,'.')") cvec(n)
		call addToLog(trim(tmpstring))
		z_in = cvec(n) + 1.d0
	end subroutine zin_solver

	pure function H_eff(y)
		real(dl), intent(in) :: y
		type(cmplxMatNN) :: H_eff

		call allocateCmplxMat(H_eff)

		!missing: term for NC!
		H_eff%re = 0.d0 &
			+ nuMassesMat(:,:)/(2*y) &
			+ leptonDensities(:,:) * y &
			+ nuDensities%re(:,:) * y
		H_eff%im = 0.d0 &
			+ nuDensities%im(:,:) * y
	end function H_eff

	pure function H_eff_cmplx(y)
		complex(dl), dimension(maxflavorNumber, maxflavorNumber) :: H_eff_cmplx
		real(dl), intent(in) :: y
		type(cmplxMatNN) :: H
		integer :: i, j

		H = H_eff(y)
		H_eff_cmplx(:,:) = cmplx(0.d0, 0.d0)
		do i=1, flavorNumber
			do j=1, flavorNumber
				H_eff_cmplx(i, j) = cmplx(H%re(i,j), H%im(i,j))
			end do
		end do
	end function H_eff_cmplx

	pure function rho_diag_mass(iy)
		type(cmplxMatNN) :: rho_diag_mass
		integer, intent(in) :: iy
		complex(dl), dimension(maxflavorNumber, maxflavorNumber) :: tmpComplMat, transfMat
		real(dl), dimension(maxflavorNumber) :: tmpvec
		integer :: i, j, k, l

		call allocateCmplxMat(rho_diag_mass)
		rho_diag_mass%re(:,:) = 0.d0
		rho_diag_mass%im(:,:) = 0.d0
		tmpvec = 0.d0

		transfMat(:,:) = cmplx(0.d0, 0.d0)
		tmpComplMat = H_eff_cmplx(y_arr(iy))
		call HEigensystem(flavorNumber, tmpComplMat, flavorNumber, tmpvec, transfMat, flavorNumber, 0)
		do k=1, flavorNumber
			do l=k, flavorNumber
				do i=1, flavorNumber
					do j=1, flavorNumber
						rho_diag_mass%re(k, l) = rho_diag_mass%re(k, l) &
							+ dble(conjg(transfMat(i, k))*transfMat(j, l)&
								*cmplx(nuDensMatVecFD(iy)%re(i, j), nuDensMatVecFD(iy)%im(i, j)))
						rho_diag_mass%im(k, l) = rho_diag_mass%im(k, l) &
							+ dimag(conjg(transfMat(i, k))*transfMat(j, l)&
								*cmplx(nuDensMatVecFD(iy)%re(i, j), nuDensMatVecFD(iy)%im(i, j)))
					end do
				end do
			end do
		end do
	end function rho_diag_mass

	subroutine nuDens_to_file(u, ix, iy, x, mat, reim, fname)
		integer, intent(in) :: u, ix, iy
		real(dl), intent(in) :: x
		logical, intent(in) :: reim!true for real, false for imaginary part
		type(cmplxMatNN), dimension(:), allocatable, intent(in) :: mat
		character(len=*), intent(in) :: fname
		integer :: m
		real(dl), dimension(:), allocatable :: tmpvec

		allocate(tmpvec(Ny))
		call openFile(u, trim(fname), firstWrite)
		if (reim) then
			do m=1, nY
				tmpvec(m)=mat(m)%re(ix, iy)
			end do
		else
			do m=1, nY
				tmpvec(m)=mat(m)%im(ix, iy)
			end do
		end if
		write(u, multidblfmt) x, tmpvec
		deallocate(tmpvec)
		close(u)
	end subroutine nuDens_to_file

	subroutine saveRelevantInfo(x, vec)
		real(dl), intent(in) :: x
		real(dl), dimension(:), intent(in) :: vec
		type(cmplxMatNN), dimension(:), allocatable :: rho_mass
		complex(dl), dimension(:,:), allocatable :: tmpComplMat, transfMat
		integer :: k, i, j, iy
		real(dl) :: neff, z
		integer, parameter :: iu = 8972
		character(len=200) :: fname

		write(fname, '(A,'//dblfmt//')') '[output] Saving info at x=', x
		call addToLog(trim(fname))!not a filename but the above string

		z = vec(ntot)+1.d0
		if (save_nuDens_evolution) then
			!density matrix in flavor space
			do k=1, flavorNumber
				write(fname, '(A,I1,A)') trim(outputFolder)//'/nuDens_diag', k, '.dat'
				call nuDens_to_file(iu, k, k, x, nuDensMatVecFD, .true., trim(fname))
			end do
			if (collision_offdiag.ne.0 .and. collision_offdiag.ne.3) then
				do i=1, flavorNumber-1
					do j=i+1, flavorNumber
						write(fname, '(A,I1,I1,A)') trim(outputFolder)//'/nuDens_nd_', i, j, '_re.dat'
						call nuDens_to_file(iu, i, j, x, nuDensMatVecFD, .true., trim(fname))

						write(fname, '(A,I1,I1,A)') trim(outputFolder)//'/nuDens_nd_', i, j, '_im.dat'
						call nuDens_to_file(iu, i, j, x, nuDensMatVecFD, .false., trim(fname))
					end do
				end do
			end if
			!density matrix in mass space
			allocate(rho_mass(Ny))
			call updateMatterDensities(x, z)
			!$omp parallel do shared(rho_mass) private(iy) schedule(static)
			do iy=1, Ny
				rho_mass(iy) = rho_diag_mass(iy)
			end do
			!$omp end parallel do
			do k=1, flavorNumber
				write(fname, '(A,I1,A)') trim(outputFolder)//'/nuDens_mass', k, '.dat'
				call nuDens_to_file(iu, k, k, x, rho_mass, .true., trim(fname))
			end do
			if (collision_offdiag.ne.0 .and. collision_offdiag.ne.3) then
				do i=1, flavorNumber-1
					do j=i+1, flavorNumber
						write(fname, '(A,I1,I1,A)') trim(outputFolder)//'/nuDens_mass_nd_', i, j, '_re.dat'
						call nuDens_to_file(iu, i, j, x, rho_mass, .true., trim(fname))

						write(fname, '(A,I1,I1,A)') trim(outputFolder)//'/nuDens_mass_nd_', i, j, '_im.dat'
						call nuDens_to_file(iu, i, j, x, rho_mass, .false., trim(fname))
					end do
				end do
			end if
		end if
		if (save_z_evolution) then
			call openFile(iu, trim(outputFolder)//'/z.dat', firstWrite)
			if (save_w_evolution) then
				write(iu, multidblfmt) x, z, vec(ntot-1)+1.d0
			else
				write(iu, multidblfmt) x, z
			end if
			close(iu)
		end if
		if (save_Neff) then
			neff = allNuDensity()/photonDensity(vec(ntot)+1.d0)
			call openFile(iu, trim(outputFolder)//'/Neff.dat', firstWrite)
			write(iu, multidblfmt) &
				x, 8.d0/7.d0*neff, 8.d0/7.d0*zid**4*neff
			close(iu)
		end if
		firstWrite=.false.
	end subroutine saveRelevantInfo

	subroutine solver
		real(dl) :: xstart, xend, xchk
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

		call densMat_2_vec(nuDensVec)
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

		call addToLog("[solver] Solver ended. "//trim(tmpstring))
	end subroutine solver

	pure subroutine drhoy_dx_fullMat(matrix, x, z, iy, dme2, sqrtraddens, Fint)
		use ndInterfaces2
		procedure (collision_integrand) :: Fint
		type(cmplxMatNN), intent(out) :: matrix
		real(dl), intent(in) :: x, z, dme2, sqrtraddens
		integer, intent(in) :: iy
		real(dl) :: y, overallNorm, fd, cf
		integer :: ix
		type(coll_args) :: collArgs
		type(cmplxMatNN) :: tmpmat

		call allocateCmplxMat(tmpmat)

		y = nuDensMatVecFD(iy)%y

		collArgs%x = x
		collArgs%z = z
		collArgs%y1 = y
		collArgs%dme2 = dme2
		collArgs%iy = iy

		overallNorm = overallFactor / sqrtraddens
		matrix = H_eff(y)

		!switch imaginary and real parts because of the "-i" factor
		call Commutator(matrix%re, nuDensMatVecFD(iy)%re, tmpmat%im)
		call Commutator(matrix%re, nuDensMatVecFD(iy)%im, tmpmat%re)

		!matrix is now redefined
		cf = x**2/m_e_cub
		matrix%im = - tmpmat%im * cf
		matrix%re = tmpmat%re * cf

		tmpmat = get_collision_terms(collArgs, Fint)
		matrix%re = matrix%re + tmpmat%re
		matrix%im = matrix%im + tmpmat%im

		matrix%re = matrix%re * overallNorm
		matrix%im = matrix%im * overallNorm
		fd = fermiDirac(y)
		do ix=1, flavorNumber
			matrix%re(ix,ix) = matrix%re(ix,ix) / fd
			matrix%im(ix,ix) = 0.d0
		end do
	end subroutine drhoy_dx_fullMat

	pure subroutine derivative (x, z, m, dme2, sqrtraddens, n, ydot)
!		compute rho derivatives for a given momentum y_arr(m), save to ydot
		real(dl), intent(in) :: x, z, dme2, sqrtraddens
		integer, intent(in) :: m, n
		real(dl), dimension(n), intent(out) :: ydot
		integer :: i, j, k
		type(cmplxMatNN) :: mat

		call drhoy_dx_fullMat(mat, x, z, m, dme2, sqrtraddens, coll_nue_3_int)
		do i=1, flavorNumber
			ydot(i) = mat%re(i,i)
		end do
		k=flavorNumber+1
		if (collision_offdiag.ne.0 .and. collision_offdiag.ne.3) then
			do i=1, flavorNumber-1
				do j=i+1, flavorNumber
					ydot(k) = mat%re(i,j)
					ydot(k+1) = mat%im(i,j)
					k=k+2
				end do
			end do
		end if
	end subroutine derivative

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

		dme2 = dme2_electron(x, 0.d0, z)
		sqrtraddens = sqrt(radDensity(x, z))
		call updateMatterDensities(x, z)

		!$omp parallel shared(ydot, x, z, dme2, sqrtraddens, flavNumSqu) private(m, s, tmpvec)
		allocate(tmpvec(flavNumSqu))
		tmpvec = 0
		!$omp do schedule(static)
		do m=1, Ny
			call derivative(x, z, m, dme2, sqrtraddens, flavNumSqu, tmpvec)
			s=(m-1)*flavNumSqu
			ydot(s+1:s+flavNumSqu) = tmpvec(:)
		end do
		!$omp end do
		deallocate(tmpvec)
		!$omp end parallel

		call dz_o_dx(x, w, z, ydot, ntot)

		call densMat_2_vec(nuDensVec)
	end subroutine derivatives

	subroutine jdum
		!necessary for dlsoda but not needed
	end subroutine jdum

	function Neff_from_rho_z(w, z)
		real(dl) :: Neff_from_rho_z
		real(dl), intent(in) :: w, z
		real(dl) :: rhototnu, ndeq, totf
		integer :: ix

		totf=0.d0
		do ix=1, flavorNumber
			if(.not.sterile(ix)) totf = totf + nuFactor(ix)
		end do
		ndeq=nuDensityEq(w)
		rhototnu = (allNuDensity() - totf*ndeq)/ndeq
		Neff_from_rho_z = (w*zid/z)**4 * &
			(3.d0 + rhototnu)
	end function Neff_from_rho_z

	subroutine finalresults
		use ndInterfaces1
		procedure (nuDensity_integrator), pointer :: nuDensityInt
		real(dl) :: ndeq, tmp, w, z
		real(dl), dimension(:), allocatable :: tmpvec
		integer :: ix, iy

		if (use_gauss_laguerre) then
			nuDensityInt => nuDensityGL
		else
			nuDensityInt => nuDensityNC
		end if

		call openFile(9876, trim(outputFolder)//'/rho_final.dat', .true.)
		allocate(tmpvec(flavorNumber))
		do iy=1, nY
			do ix=1, flavorNumber
				tmpvec(ix)=nuDensMatVecFD(iy)%re(ix, ix)
			end do
			write(9876, multidblfmt) nuDensMatVecFD(iy)%y, tmpvec
		end do
		close(9876)
		deallocate(tmpvec)

		w = nuDensVec(ntot-1) + 1.d0
		z = nuDensVec(ntot) + 1.d0
		call openFile(9876, trim(outputFolder)//'/resume.dat', .true.)
		if (save_w_evolution) then
			write(*,"('final w = ',F11.8)") w
			write(9876,"('final w = ',F11.8)") w
		end if
		write(*,"('final z = ',F11.8)") z
		write(9876,"('final z = ',F11.8)") z
		!since it was never taken into account, w must not be used here to get the delta_rho,
		!otherwise the result is not referring to the same quantity
		ndeq=nuDensityEq(1.d0)
		do ix=1, flavorNumber
			tmp = (nuDensityInt(ix, ix) - ndeq)*nuFactor(ix)/ndeq
			write(*,"('dRho_',I1,'  = ',F9.6)") ix, tmp
			write(9876,"('dRho_',I1,'  = ',F9.6)") ix, tmp
		end do
		tmp = Neff_from_rho_z(1.d0, z)
		write(*,"('Neff    = ',F9.6)") tmp
		write(9876,"('Neff    = ',F9.6)") tmp
		close(9876)
	end subroutine finalresults
end module ndEquations
