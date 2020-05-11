module fpConfig
	use precision
	use constants
	use variables
	use fpErrors
	use FileIni
	use fpMatrices
	use fpInteractions
	use fpEquations
	use utilities
	implicit none

	integer :: num_args, ixa
	character(len=300), dimension(:), allocatable :: args

	contains

	subroutine allocateStuff()
		integer :: nf, ix
		nf = flavorNumber
		allocate(nuMasses(nf), nuFactor(nf), sterile(nf))
		allocate(mixMat(nf,nf), mixMatInv(nf,nf))
		allocate(nuMassesMat(nf,nf), leptonDensities(nf,nf))
		call allocateCmplxMat(nuDensities)
		allocate(dampTermYYYWdy(Ny))
		allocate(dampTermMatrixCoeffNue(nf,nf))
		allocate(dampTermMatrixCoeffNunu(nf,nf))
		allocate(GL_mat(nf,nf), GR_mat(nf,nf), GLR_vec(2, nf,nf))
		allocate(massSplittings(nf))
		allocate(mixingAngles(nf,nf))

		if (.not.allocated(ln_2dint_y)) &
			allocate(ln_2dint_y(ln_2dint_Npts), ln_2dint_dy(ln_2dint_Npts))
		ln_2dint_y = geomspace(ln_2dint_lower, ln_2dint_upper, ln_2dint_Npts)
		do ix = 1, ln_2dint_Npts-1
			ln_2dint_dy(ix) = ln_2dint_y(ix+1) - ln_2dint_y(ix)
		end do
	end subroutine allocateStuff

	subroutine setMassMatrix()
		integer :: i, nf
		real(dl), dimension(:), allocatable :: mv
		real(dl), dimension(:,:), allocatable :: m1

		nf = flavorNumber
		allocate(mv(nf), m1(nf,nf))

		if (nf .ge. 2) then
			mv(1) = 0.d0
			do i=2, nf
				mv(i) = mv(1) + massSplittings(i)
			end do
		end if
		mv = mv - minval(mv)

		call createDiagMat(m1, nf, mv)
		call tripleProdMat(mixMat, m1, mixMatInv, nuMassesMat)
		write(*,*) "masses:"
		call printVec(mv, maxFlavorNumber)
		write(*,*) "rotated mass matrix:"
		call printMat(nuMassesMat)
		deallocate(mv,m1)
	end subroutine setMassMatrix

	subroutine setMixingMatrix()
		integer :: nf, i, j
		real(dl), dimension(:,:), allocatable :: m1, m2
		character(len=1) :: nfstr

		nf = flavorNumber
		write(nfstr,"(I1)") nf

		if (verbose.gt.1) &
			call addToLog("[config] creating mixing matrix, using "//nfstr//" flavors")
		allocate(m1(nf,nf), m2(nf,nf))

		call createIdentityMat(m1, nf)
		do i = flavorNumber, 2, -1
			do j = i-1, 1, -1
				call createRotMat(m2, nf, j, i, mixingAngles(j,i))
				m1 = matmul(m1, m2)
			end do
		end do
		
		mixMat=m1
		call inverseMat(mixMat, mixMatInv)
		write(*,*) "mixing matrix:"
		call printMat(mixMat)
		deallocate(m1,m2)
	end subroutine setMixingMatrix

	subroutine setDampingFactorCoeffs
		integer :: ix, iy
		real(dl) :: nue_nue_nux, nue_nue_nus, nue_numu_nutau, nue_nux_nus
		real(dl) :: nunu_nue_nux, nunu_nue_nus, nunu_numu_nutau, nunu_nux_nus
		real(dl) :: nue_nue, nue_nux, nue_nus
		real(dl) :: nunu_nue, nunu_nux, nunu_nus
		real(dl) :: xW
		character(len=300) :: tmpstr

		write(*,*) "[config] Damping factors settings: ", collision_offdiag
		xW = sin2thW

		dampTermYYYWdy = 0.d0
		dampTermMatrixCoeffNue = 0.d0
		dampTermMatrixCoeffNunu = 0.d0
		nunu_nue = 0.d0
		nunu_nux = 0.d0
		nunu_nus = 0.d0
		nue_nue = 0.d0
		nue_nux = 0.d0
		nue_nus = 0.d0

		!numbers from McKellar:1992ja
		!terms for scattering, annihilation with electrons
		!nu_e - nu_X
		nue_nue_nux = &
			8.d0 &!e+nu -> e+nu
			+ (8.d0*sin2thW**2 + 1.d0)!nu+bnu -> e+e-
		!nu_mu - nu_tau
		nue_numu_nutau = &
			(8.d0*sin2thW**2 - 4.d0*sin2thW + 1.d0)!nu+bnu -> e+e-
		!nu_e - nu_s
		nue_nue_nus = 2.d0*(&
			(8.d0*sin2thW**2 + 4.d0*sin2thW + 1.d0) &!e+nu -> e+nu
			+ (4.d0*sin2thW**2 + 2.d0*sin2thW + 0.5d0) &!nu+bnu -> e+e-
		)
		!nu_X - nu_s
		nue_nux_nus = 2.d0*(&
			(8.d0*sin2thW**2 - 4.d0*sin2thW + 1.d0) &!e+nu -> e+nu
			+ (4.d0*sin2thW**2 - 2.d0*sin2thW + 0.5d0) &!nu+bnu -> e+e-
		)
		!terms for nunu scattering
		!nu_e - nu_X
		nunu_nue_nux = 6.d0 !nu+(b)nu -> nu+(b)nu
		!nu_mu - nu_tau
		nunu_numu_nutau = 6.d0 !nu+(b)nu -> nu+(b)nu
		!nu_e - nu_s
		nunu_nue_nus = 2.d0 * 13.d0 !nu+(b)nu -> nu+(b)nu
		!nu_X - nu_s
		nunu_nux_nus = 2.d0 * 13.d0 !nu+(b)nu -> nu+(b)nu

		if (collision_offdiag.eq.4 .or. collision_offdiag.eq.5) then
			!formulas from YYYW notes
			nunu_nue = 2.d0
			nunu_nux = 2.d0
			nunu_nue_nux = 2.d0
			nunu_numu_nutau = 2.d0
			nue_nue = 4.d0*xW**2 + 2.d0 * xW + 0.5d0
			nue_nux = 4.d0*xW**2 - 2.d0 * xW + 0.5d0
			nue_nue_nux = 4.d0*xW**2 + 0.5d0
			nue_numu_nutau = 4.d0*xW**2 - 2.d0 * xW + 0.5d0
			nunu_nus = 0.d0
			nunu_nue_nus = 0.d0
			nunu_nux_nus = 0.d0
			nue_nus = 0.d0
			nue_nue_nus = 3.d0*xW**2 + 1.d0*xW + 0.25d0
			nue_nux_nus = 3.d0*xW**2 - 1.d0*xW + 0.25d0
		end if
		dampTermMatrixCoeffNue(1, 1) = nue_nue
		dampTermMatrixCoeffNunu(1, 1) = nunu_nue
		if (flavorNumber .ge. 2) then
			if (sterile(2)) then
				dampTermMatrixCoeffNue(1, 2) = nue_nue_nus
				dampTermMatrixCoeffNunu(1, 2) = nunu_nue_nus
				dampTermMatrixCoeffNue(2, 2) = nue_nus
				dampTermMatrixCoeffNunu(2, 2) = nunu_nus
			else
				dampTermMatrixCoeffNue(1, 2) = nue_nue_nux
				dampTermMatrixCoeffNunu(1, 2) = nunu_nue_nux
				dampTermMatrixCoeffNue(2, 2) = nue_nux
				dampTermMatrixCoeffNunu(2, 2) = nunu_nux
			end if
		end if
		if (flavorNumber .ge. 3) then
			if (sterile(3)) then
				dampTermMatrixCoeffNue(1, 3) = nue_nue_nus
				dampTermMatrixCoeffNunu(1, 3) = nunu_nue_nus
				if (.not.sterile(2)) then
					dampTermMatrixCoeffNue(2, 3) = nue_nux_nus
					dampTermMatrixCoeffNunu(2, 3) = nunu_nux_nus
				end if
				dampTermMatrixCoeffNue(3, 3) = nue_nus
				dampTermMatrixCoeffNunu(3, 3) = nunu_nus
			else
				dampTermMatrixCoeffNue(1, 3) = nue_nue_nux
				dampTermMatrixCoeffNunu(1, 3) = nunu_nue_nux
				dampTermMatrixCoeffNue(2, 3) = nue_numu_nutau
				dampTermMatrixCoeffNunu(2, 3) = nunu_numu_nutau
				dampTermMatrixCoeffNue(3, 3) = nue_nux
				dampTermMatrixCoeffNunu(3, 3) = nunu_nux
			end if
		end if
		if (damping_read_zero) then
			do ix=4, flavorNumber
				write(tmpstr,"('damping_',2I1,'_zero')") 1, ix
				if (.not. read_ini_logical(trim(tmpstr), .false.)) &
					dampTermMatrixCoeffNue(1, ix) = nue_nue_nus
					dampTermMatrixCoeffNunu(1, ix) = nunu_nue_nus
				do iy=2, ix-1
					write(tmpstr,"('damping_',2I1,'_zero')") iy, ix
					if (.not. read_ini_logical(trim(tmpstr), .false.)) &
						dampTermMatrixCoeffNue(iy, ix) = nue_nux_nus
						dampTermMatrixCoeffNunu(iy, ix) = nunu_nux_nus
				end do
			end do
		else
			do ix=4, flavorNumber
				dampTermMatrixCoeffNue(1, ix) = nue_nue_nus
				dampTermMatrixCoeffNunu(1, ix) = nunu_nue_nus
				dampTermMatrixCoeffNue(ix, ix) = nue_nus
				dampTermMatrixCoeffNunu(ix, ix) = nunu_nus
				do iy=2, 3
					if (.not.sterile(iy)) then
						dampTermMatrixCoeffNue(iy, ix) = nue_nux_nus
						dampTermMatrixCoeffNunu(iy, ix) = nunu_nux_nus
					end if
				end do
			end do
		end if
		if (damping_no_nue) &
			dampTermMatrixCoeffNue = 0.d0
		if (damping_no_nunu) &
			dampTermMatrixCoeffNunu = 0.d0
		write(*,*)"[config] Damping factors (Nue):"
		call printMat(dampTermMatrixCoeffNue)
		write(*,*)"[config] Damping factors (NuNu):"
		call printMat(dampTermMatrixCoeffNunu)

		if (collision_offdiag.eq.4 .or. collision_offdiag.eq.5) then
			!formulas from YYYW notes
			write(*,*) "[config] Computing d(y) for damping factors a la YYYW..."
			!$omp parallel do default(shared) private(ix) schedule(dynamic)
			do ix=1, Ny
				dampTermYYYWdy(ix) = dy_damping(y_arr(ix)) * y_arr(ix)**3
			end do
			!$omp end parallel do
			write(*,"(3A14)") "y", "f_eq(y)", "d(y)"
			do ix=1, Ny
				write(*,"(3E14.6)") y_arr(ix), feq_vec(ix), dampTermYYYWdy(ix)
			end do
		end if
		write(*,*)"[config] Damping factors done."
	end subroutine setDampingFactorCoeffs

	subroutine init_matrices
		real(dl), dimension(:), allocatable :: diag_el
		integer :: ix, iy
		real(dl) :: fdm, fdme

		allocate(diag_el(flavorNumber))

		!GL
		diag_el = sin2thW - 0.5d0
		diag_el(1) = sin2thW + 0.5d0
		do ix=1, flavorNumber
			if (sterile(ix)) diag_el(ix) = 0.d0
		end do
		call createDiagMat(GL_mat, flavorNumber, diag_el)
		GLR_vec(1,:,:) = GL_mat

		!GR
		diag_el = sin2thW
		do ix=1, flavorNumber
			if (sterile(ix)) diag_el(ix) = 0.d0
		end do
		call createDiagMat(GR_mat, flavorNumber, diag_el)
		GLR_vec(2,:,:) = GR_mat

		call setDampingFactorCoeffs

		!lepton densities
		leptonDensities = 0.d0
		nuDensities%re = 0.d0
		nuDensities%im = 0.d0

		!identity (Matrix complex)
		call createIdentityMat(idMat, flavorNumber)

		!nu density matrix
		allocate(nuDensMatVec(Ny), nuDensMatVecFD(Ny))
		ntot = Ny*(flavNumSqu) + 2 !independent elements in nuDensity(y) + w,z
		allocate(nuDensVec(ntot))
		if(save_fd .and. trim(outputFolder).ne."")&
			call openFile(3154, trim(outputFolder)//"/fd.dat", .true.)
		do ix=1, Ny
			allocate(nuDensMatVec(ix)%re(flavorNumber,flavorNumber), nuDensMatVec(ix)%im(flavorNumber,flavorNumber))
			allocate(nuDensMatVecFD(ix)%re(flavorNumber,flavorNumber), nuDensMatVecFD(ix)%im(flavorNumber,flavorNumber))
			nuDensMatVec(ix)%y = y_arr(ix)
			nuDensMatVecFD(ix)%y = y_arr(ix)
			nuDensMatVec(ix)%re(:,:) = 0.d0
			nuDensMatVec(ix)%im(:,:) = 0.d0
			fdme = fermiDirac(y_arr(ix)/z_in)
			fdm = fermiDirac(y_arr(ix))
			if(save_fd .and. trim(outputFolder).ne."")&
				write(3154, multidblfmt) y_arr(ix), fdm * y_arr(ix)*y_arr(ix)
			do iy=2, flavorNumber
				if (sterile(iy)) &
					nuDensMatVec(ix)%re(iy,iy) = -1.d0
			end do
			nuDensMatVecFD(ix)%re = nuDensMatVec(ix)%re
			nuDensMatVecFD(ix)%im = nuDensMatVec(ix)%im
			do iy=1, flavorNumber
				nuDensMatVecFD(ix)%re(iy,iy) = (1.d0+nuDensMatVec(ix)%re(iy,iy)) * fdme
				nuDensMatVec(ix)%re(iy,iy) = nuDensMatVecFD(ix)%re(iy,iy)/fdm - 1.d0
			end do
		end do
		if(save_fd .and. trim(outputFolder).ne."")&
			close(3154)

		deallocate(diag_el)
	end subroutine init_matrices

	subroutine init_fermions
#ifndef NOINTERPOLATION
		call init_interp_dme2_e
#endif
		call fermions(1)%initialize("electrons", .true., 1.d0, 1d3)
		electrons => fermions(1)
#ifndef NO_MUONS
		call fermions(2)%initialize("muons", .false., m_mu_o_m_e, x_muon_cut)
		muons => fermions(2)
#endif
#ifndef NOINTERPOLATION
		call init_interp_jkyg12
#endif
	end subroutine init_fermions

	subroutine finish_y_arrays
		integer :: ix
		real(dl), dimension(:), allocatable :: fake
		allocate(dy_arr(Ny), fy_arr(Ny), feq_vec(Ny))

		call get_GLq_vectors(N_opt_y, opt_y, opt_y_w, fake, .true., 2, opt_y_cut)
		call get_GLq_vectors(N_opt_xoz, opt_xoz, opt_xoz_w, fake, .true., 2, opt_xoz_cut)
		do ix=1, Ny-1
			dy_arr(ix) = y_arr(ix+1) - y_arr(ix)
		end do
		dy_arr(Ny) = 0.d0

		feq_vec = 0.d0
		do ix=1, Ny
			feq_vec(ix) = fermiDirac(y_arr(ix))
		end do
	end subroutine

	subroutine initConfig()
		use omp_lib
		character(len=300) :: tmparg, tmpstr
		integer :: ix, iy, num_threads
		logical :: file_exist

		if (verbose>0) write(*,*) '[config] init configuration'
		num_args = command_argument_count()

		allocate(args(num_args))
		do ixa = 1, num_args
			call get_command_argument(ixa,args(ixa))
		end do
		if(num_args.gt.1) then
			logFile = trim(args(2))
		end if
		call openLogFile
		if(num_args.eq.0) &
			call criticalError("You are not passing a configuration file, please provide one.")

		call addToLog("[config] reading additional configuration from "//trim(args(1)))
		call ini_file_open(trim(args(1)), trim(args(1))//".log")

		force_replace = read_ini_logical('force_replace', .false.)
		tmparg=trim(read_ini_char('outputFolder'))
		if (trim(tmparg)/="") outputFolder=tmparg
		inquire(file=trim(outputFolder)//"/resume.dat", exist=file_exist)
		if (file_exist .and. .not. force_replace) then
			call criticalError("[config] resume.dat already exists! No need to proceed.")
		end if
		call addToLog("[config] Writing to: "//trim(outputFolder)//"/...")
		call system('mkdir -p '//trim(outputFolder))

		num_threads = read_ini_int("num_threads", 0)
		if (num_threads .gt.0) &
			CALL OMP_SET_NUM_THREADS(num_threads)

		verbose = read_ini_int('verbose', verbose)
		Nprintderivs = read_ini_real('Nprintderivs', Nprintderivs)
		checkpoint = read_ini_logical('checkpoint', .true.)
		maxiter = read_ini_int('maxiter', 100)
		toler_jkyg = read_ini_real('tolerance_jkyg', 1.d-7)
		dlsoda_atol_z = read_ini_real('dlsoda_atol_z', 1.d-6)
		dlsoda_atol_d = read_ini_real('dlsoda_atol_d', 1.d-6)
		dlsoda_atol_o = read_ini_real('dlsoda_atol_o', 1.d-6)
		dlsoda_rtol = read_ini_real('dlsoda_rtol', 1.d-6)

		interp_nx = read_ini_int('interp_nx', interp_nx0)
		interp_nz = read_ini_int('interp_nz', interp_nz0)
		interp_nxz = read_ini_int('interp_nxz', interp_nxz0)
		interp_zmin = read_ini_real('interp_zmin', interp_zmin0)
		interp_zmax = read_ini_real('interp_zmax', interp_zmax0)

		Nx = read_ini_int('Nx',100)
		use_gauss_laguerre = read_ini_logical('use_gauss_laguerre', .true.)
		Ny = read_ini_int('Ny',30)

		Nylog = read_ini_int('Nylog',7)
		allocate(x_arr(Nx))

		x_in    = read_ini_real('x_in', 0.01d0)
		if (x_in.lt.very_early_x) then
			write(tmparg, multidblfmt) very_early_x
			call criticalError("initial x is too early: x_in cannot be smaller than "//tmparg)
		end if
		x_fin   = read_ini_real('x_fin', 40.d0)
		logx_in  = log10(x_in)
		logx_fin = log10(x_fin)
		x_arr = logspace(logx_in, logx_fin, Nx)

		y_min = read_ini_real('y_min', 0.01d0)
		y_max = read_ini_real('y_max', 20.0d0)
		y_cen = read_ini_real('y_cen', 1.0d0)

		if (use_gauss_laguerre) then
			call get_GLq_vectors(Ny, y_arr, w_gl_arr, w_gl_arr2, .true., 3, y_max)
		else
			allocate(y_arr(Ny))
			if (Nylog .lt. 2) then
				y_arr = linspace(y_min, y_max, Ny)
			else
				y_arr = loglinspace(y_min, y_cen, y_max, Ny, Nylog)
			end if
			write(tmpstr,"('[config] Using Ny=',I3,' and Nylog=',I3)") Ny, Nylog
			call addToLog(trim(tmpstr))
		end if

		call finish_y_arrays

		!read mixing parameters and create matrices
		giveSinSq = read_ini_logical('givesinsq', .true.)

		!settings for collisional
		collision_offdiag = read_ini_int("collision_offdiag", 1)
		damping_no_nue = read_ini_logical("damping_no_nue", .false.)
		damping_no_nunu = read_ini_logical("damping_no_nunu", .false.)
		damping_read_zero = .true.
		ftqed_temperature_corr = read_ini_logical("ftqed_temperature_corr",.true.)
		ftqed_log_term = read_ini_logical("ftqed_log_term",.false.)
		ftqed_ord3 = read_ini_logical("ftqed_ord3",.true.)
		ftqed_e_mth_leptondens = read_ini_logical("ftqed_e_mth_leptondens",.true.)

		!settings for saving files
		save_fd = read_ini_logical("save_fd", .true.)
		save_energy_entropy_evolution = read_ini_logical("save_energy_entropy_evolution", .true.)
		save_Neff = read_ini_logical("save_Neff", .true.)
		save_nuDens_evolution = read_ini_logical("save_nuDens_evolution", .true.)
		save_w_evolution = read_ini_logical("save_w_evolution", .true.)
		save_z_evolution = read_ini_logical("save_z_evolution", .true.)

		z_in=1.d0
		allocate(interp_xvec(interp_nx), interp_yvec(interp_ny), interp_zvec(interp_nz), interp_xozvec(interp_nxz))
		interp_xvec = logspace(interp_logx_in, logx_fin, interp_nx)
		interp_yvec = logspace(interp_logy_min, interp_logy_max, interp_ny)
		interp_zvec = linspace(interp_zmin, interp_zmax, interp_nz)
		interp_xozvec = logspace(log10(very_early_x/interp_zmax), logx_fin, interp_nxz)

		flavorNumber = read_ini_int('flavorNumber', i_flavorNumber)
		if (collision_offdiag.ne.0 .and. collision_offdiag.ne.3) then
			flavNumSqu = flavorNumber**2
		else
			flavNumSqu = flavorNumber
		end if
		call allocateStuff

		do ix=1, flavorNumber
			write(tmparg,"('nuFactor',I1)") ix
			nuFactor(ix) = read_ini_real(trim(tmparg), 1.d0)
			write(tmparg,"('sterile',I1)") ix
			if (ix.le.3) then
				sterile(ix) = read_ini_logical(trim(tmparg), .false.)
			else
				sterile(ix) = read_ini_logical(trim(tmparg), .true.)
			end if
		end do
		tot_factor_active_nu = 0.d0
		tot_factor_nu = 0.d0
		do ix=1, flavorNumber
			tot_factor_nu = tot_factor_nu + nuFactor(ix)
			if (.not. sterile(ix)) &
				tot_factor_active_nu = tot_factor_active_nu + nuFactor(ix)
		end do
		write (tmparg, &
			"('[config] using ',I1,' neutrinos, counting as ',*(E10.3))") flavorNumber, nuFactor
		write(tmpstr,"(A,', of which they are steriles:',*(L2))") trim(tmparg), sterile
		call addToLog(trim(tmpstr))

		if (flavorNumber .gt. maxFlavorNumber) then
			write(tmpstr,"('[config] WARNING: only up to ',I1,' neutrino flavors are supported. Using N=',I1)") maxFlavorNumber, maxFlavorNumber
			call error(tmpstr)
			flavorNumber = maxFlavorNumber
		end if

		call init_fermions
		call zin_solver

		massSplittings = 0.d0
		massSplittings(2) = read_ini_real('dm21', i_dm21)
		if (flavorNumber .gt. 2) then
			massSplittings(3) = read_ini_real('dm31', i_dm31)
		end if
		do ix = 4, flavorNumber
			write(tmpstr,"('dm',I1,'1')") ix
			massSplittings(ix) = read_ini_real(trim(tmpstr), zero)
		end do
		if (giveSinSq) then
			mixingAngles(1,2) = asin(sqrt(read_ini_real('theta12', i_theta12)))
			if (flavorNumber .gt. 2) then
				mixingAngles(1,3) = asin(sqrt(read_ini_real('theta13', i_theta13)))
				mixingAngles(2,3) = asin(sqrt(read_ini_real('theta23', i_theta23)))
			end if
			do ix=4, flavorNumber
				do iy=1, ix-1
					write(tmpstr,"('theta',2I1)") iy, ix
					mixingAngles(iy, ix) = asin(sqrt(read_ini_real(trim(tmpstr), zero)))
				end do
			end do
		else
			mixingAngles(1,2) = read_ini_real('theta12', asin(sqrt(i_theta12)))
			if (flavorNumber .gt. 2) then
				mixingAngles(1,3) = read_ini_real('theta13', asin(sqrt(i_theta13)))
				mixingAngles(2,3) = read_ini_real('theta23', asin(sqrt(i_theta23)))
			end if
			do ix=4, flavorNumber
				do iy=1, ix-1
					write(tmpstr,"('theta',2I1)") iy, ix
					mixingAngles(iy, ix) = read_ini_real(trim(tmpstr), zero)
				end do
			end do
		end if

		call setMixingMatrix()
		call setMassMatrix()
		!create other matrices
		call init_matrices

		call ini_file_close()
		call rename(trim(args(1))//".log", trim(outputFolder)//'/ini.log')
		call addToLog("[config] Read configuration from ini file: complete.")
	end subroutine initconfig
end module fpConfig
