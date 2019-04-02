module ndConfig
	use precision
	use constants
	use variables
	use ndErrors
	use FileIni
	use ndMatrices
	use ndInteractions
	use ndEquations
	use utilities
	implicit none

	integer :: num_args, ixa
	character(len=300), dimension(:), allocatable :: args

	contains

	subroutine allocateStuff()
		integer :: nf
		nf = flavorNumber
		allocate(nuMasses(nf), nuFactor(nf), sterile(nf))
		allocate(mixMat(nf,nf), mixMatInv(nf,nf))
		allocate(nuMassesMat(nf,nf), leptonDensities(nf,nf))
		allocate(dampTermMatrixCoeff(nf,nf))
		allocate(GL_mat(nf,nf), GR_mat(nf,nf), GLR_vec(2, nf,nf))
		allocate(massSplittings(nf))
		allocate(mixingAngles(nf,nf))
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
		real(dl) :: nue_nux, nue_nus, numu_nutau, nux_nus
		dampTermMatrixCoeff = 0.d0
		!numbers from McKellar:1992ja
		!nu_e - nu_X
		nue_nux = &
			8.d0 &!e+nu -> e+nu
			+ 6.d0 &!nu+(b)nu -> nu+(b)nu
			+ (8.d0*sin2thW**2 + 1.d0)!nu+bnu -> e+e-
		!nu_mu - nu_tau
		numu_nutau = &
			6.d0 &!nu+(b)nu -> nu+(b)nu
			+ (8.d0*sin2thW**2 - 4.d0*sin2thW + 1.d0)!nu+bnu -> e+e-
		!nu_e - nu_s
		nue_nus = 2.d0*(&
			(8.d0*sin2thW**2 + 4.d0*sin2thW + 1.d0) &!e+nu -> e+nu
			+ 13.d0 &!nu+(b)nu -> nu+(b)nu
			+ (4.d0*sin2thW**2 + 2.d0*sin2thW + 0.5d0))!nu+bnu -> e+e-
		!nu_X - nu_s
		nux_nus = 2.d0*(&
			(8.d0*sin2thW**2 - 4.d0*sin2thW + 1.d0) &!e+nu -> e+nu
			+ 13.d0 &!nu+(b)nu -> nu+(b)nu
			+ (4.d0*sin2thW**2 - 2.d0*sin2thW + 0.5d0))!nu+bnu -> e+e-
		if (flavorNumber .ge. 2) then
			if (sterile(2)) then
				dampTermMatrixCoeff(1, 2) = nue_nus
			else
				dampTermMatrixCoeff(1, 2) = nue_nux
			end if
		end if
		if (flavorNumber .ge. 3) then
			if (sterile(3)) then
				dampTermMatrixCoeff(1, 3) = nue_nus
				dampTermMatrixCoeff(2, 3) = nux_nus
			else
				dampTermMatrixCoeff(1, 3) = nue_nux
				dampTermMatrixCoeff(2, 3) = numu_nutau
			end if
		end if
		if (flavorNumber .ge. 4) then
			dampTermMatrixCoeff(1, 4) = nue_nus
			dampTermMatrixCoeff(2, 4) = nux_nus
			dampTermMatrixCoeff(3, 4) = nux_nus
		end if
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
			nuDensMatVec(ix)%re(1,1) = 0.d0
			if (flavorNumber.gt.1 .and. sterile(2)) &
				nuDensMatVec(ix)%re(2,2) = -1.d0
			if (flavorNumber.gt.2 .and. sterile(3)) &
				nuDensMatVec(ix)%re(3,3) = -1.d0
			if (flavorNumber.gt.3) &
				nuDensMatVec(ix)%re(4,4) = -1.d0
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
		if(num_args.gt.0) then
			call addToLog("[config] reading additional configuration from "//trim(args(1)))
			call ini_file_open(trim(args(1)), trim(args(1))//".log")

			tmparg=trim(read_ini_char('outputFolder'))
			if (trim(tmparg)/="") outputFolder=tmparg
			inquire(file=trim(outputFolder)//"/resume.dat", exist=file_exist)
			if (file_exist) then
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
			toler = read_ini_real('tolerance', 1.d-5)
			toler_dme2 = read_ini_real('tolerance_dme2', 1.d-5)
			toler_ed = read_ini_real('tolerance_ed', 1.d-4)
			toler_jkyg = read_ini_real('tolerance_jkyg', 1.d-7)
			dlsoda_atol = read_ini_real('dlsoda_atol', 1.d-6)
			dlsoda_rtol = read_ini_real('dlsoda_rtol', 1.d-6)

			interp_nx = read_ini_int('interp_nx', interp_nx0)
			interp_nz = read_ini_int('interp_nz', interp_nz0)
			interp_nxz = read_ini_int('interp_nxz', interp_nxz0)
			interp_zmin = read_ini_real('interp_zmin', interp_zmin0)
			interp_zmax = read_ini_real('interp_zmax', interp_zmax0)

			Nx = read_ini_int('Nx',100)
			Ny = read_ini_int('Ny',100)
			Nylog = read_ini_int('Nylog',7)
			allocate(x_arr(Nx), y_arr(Ny))
			allocate(dy_arr(Ny), fy_arr(Ny))

			x_in    = read_ini_real('x_in', 0.01d0)
			x_fin   = read_ini_real('x_fin', 40.d0)
			logx_in  = log10(x_in)
			logx_fin = log10(x_fin)
			x_arr = logspace(logx_in, logx_fin, Nx)

			y_min = read_ini_real('y_min', 0.01d0)
			y_max = read_ini_real('y_max', 20.0d0)
			y_cen = read_ini_real('y_cen', 1.0d0)
			if (Nylog .lt. 2) then
				y_arr = linspace(y_min, y_max, Ny)
			else
				y_arr = loglinspace(y_min, y_cen, y_max, Ny, Nylog)
			end if

			do ix=1, Ny-1
				dy_arr(ix) = y_arr(ix+1) - y_arr(ix)
			end do
			dy_arr(Ny) = 0.d0

			call zin_solver

			!read mixing parameters and create matrices
			giveSinSq = read_ini_logical('givesinsq', .true.)

			!settings for collisional
			collision_offdiag = read_ini_int("collision_offdiag", 1)
			dme2_temperature_corr = read_ini_logical("dme2_temperature_corr",.true.)

			!settings for saving files
			save_fd = read_ini_logical("save_fd",.true.)
			save_nuDens_evolution = read_ini_logical("save_nuDens_evolution",.true.)
			save_z_evolution = read_ini_logical("save_z_evolution",.true.)
			save_w_evolution = read_ini_logical("save_w_evolution",.true.)

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

			if (flavorNumber .gt. maxflavorNumber) then
				write(tmpstr,"('[config] WARNING: only up to ',I1,' neutrino flavors are supported. Using N=',I1)") maxflavorNumber, maxflavorNumber
				call error(tmpstr)
				flavorNumber = maxflavorNumber
			end if

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
			write(tmpstr,"('[config] Using Ny=',I3,' and Nylog=',I3)") Ny, Nylog
			call addToLog(trim(tmpstr))
			call init_matrices

			allocate(interp_xvec(interp_nx), interp_zvec(interp_nz), interp_xozvec(interp_nxz))
			interp_xvec = logspace(logx_in, logx_fin, interp_nx)
			interp_zvec = linspace(interp_zmin, interp_zmax, interp_nz)
			interp_xozvec = logspace(log10(x_in/interp_zmax), logx_fin, interp_nxz)
		else
			call criticalError("You are not passing a configuration file...are you sure?")
		end if
		call ini_file_close()
		call addToLog("[config] Read configuration from ini file: complete.")
	end subroutine initconfig
end module ndConfig
