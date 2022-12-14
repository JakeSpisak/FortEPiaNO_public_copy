module fpConfig
	use precision
	use constants
	use variables
	use fpErrors
	use FileIni
	use fpinput
	use fpMatrices
	use fpInteractions
	use fpEquations
	use utilities
	implicit none

	integer :: num_args, ixa
	character(len=300), dimension(:), allocatable :: args
	character(len=300) :: inilogname

	contains

	subroutine initialOperations(logname)
		character(len=*), intent(in) :: logname

		logFile = logname
		call openLogFile

#ifdef FULL_F_AB
		call addToLog("[precompiler] Compiled to compute the full products of rho, G_L and G_R matrices")
#endif
#ifdef FULL_F_NU
		call addToLog("[precompiler] Compiled to compute the F factors using the full neutrino density matrix")
#endif
#ifdef NO_INTERPOLATION
		call addToLog("[precompiler] Compiled without interpolations for lepton densities and other quantities")
#endif
#ifdef DO_MUONS
		call addToLog("[precompiler] Compiled with contributions from muons")
#endif
#ifdef NO_NUE_ANNIHILATION
		call addToLog("[precompiler] Compiled without contributions from nunu<->ee annihilation to collision integrals")
#endif
#ifdef TESTSPEED
		call addToLog("[precompiler] Will execute speed test")
#endif
	end subroutine initialOperations

	subroutine endConfig
		call ini_file_close()
		call rename(trim(inilogname), trim(outputFolder)//'/ini.log')
		call addToLog("[config] Read configuration from ini file: complete.")
	end subroutine endConfig

	subroutine checkIniFile(ininame)
		character(len=*), intent(in) :: ininame

		if(num_args.eq.0) &
			call criticalError("You are not passing a configuration file, please provide one.")

		call addToLog("[config] reading additional configuration from "//trim(ininame))
        inilogname=trim(ininame)//".log"
		call ini_file_open(trim(ininame), trim(inilogname))
	end subroutine checkIniFile

	subroutine checkOutputFolder
		logical :: file_exist

		inquire(file=trim(outputFolder)//"/resume.dat", exist=file_exist)
		if (file_exist .and. .not. force_replace) then
			call criticalError("[config] resume.dat already exists! No need to proceed.")
		end if
		call addToLog("[config] Writing to: "//trim(outputFolder)//"/...")
		call system('mkdir -p '//trim(outputFolder))
	end subroutine checkOutputFolder

	subroutine checkNumThreads
		use omp_lib

		if (num_threads .gt.0) &
			CALL OMP_SET_NUM_THREADS(num_threads)
	end subroutine checkNumThreads

	subroutine allocateStuff()
		integer :: nf, ix

		if (has_offdiagonal()) then
			flavNumSqu = flavorNumber**2
		else
			flavNumSqu = flavorNumber
		end if

		nf = flavorNumber
		allocate(nuMasses(nf), nuFactor(nf), sterile(nf), Gs(nf))
		allocate(mixMat(nf,nf), mixMatInv(nf,nf))
		allocate(nuMassesMat(nf,nf), leptonDensities(nf,nf))
		call allocateCmplxMat(nuDensities)
		allocate(dampTermMatrixCoeffNue(nf,nf))
		allocate(dampTermMatrixCoeffNunu(nf,nf))
		allocate(dampingSetZero(nf,nf))
		allocate(GL_mat(nf,nf), GR_mat(nf,nf), GLR_vec(2, nf,nf))
		allocate(massSplittings(nf))
		allocate(mixingAngles(nf,nf))
		allocate(nsi_epsilon(nf,nf))

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

	subroutine setDampingFactors
		integer :: ix, iy
		character(len=300) :: tmpstr

		write(*,*) "[collint] Diagonal terms for collision integrals are set to zero: ", collint_diagonal_zero
		write(*,*) "[collint] Use dampings instead of full integrals for nue off-diagonal contributions: ", collint_offdiag_damping
		if (collint_damping_type.eq.0) then
			write(*,*) "[collint] off-diagonal contributions are zero"
		else if (collint_damping_type.eq.1) then
			write(*,*) "[collint] off-diagonal contributions are computed using Bennett:2020zkv expressions"
		else if (collint_damping_type.eq.2) then
			write(*,*) "[collint] off-diagonal contributions are computed using McKellar:1992ja expressions"
		else
			call criticalError("Invalid value for 'collint_damping_type', it can be [0, 1, 2].")
		end if

		call setDampingFactorCoeffs

		if (damping_read_zero) then
			do ix=1, flavorNumber
				do iy=ix+1, flavorNumber
					if (dampingSetZero(ix, iy)) then
						dampTermMatrixCoeffNue(ix, iy) = 0.d0
						dampTermMatrixCoeffNunu(ix, iy) = 0.d0
					end if
				end do
			end do
		end if

		if (collint_od_no_nue) &
			dampTermMatrixCoeffNue = 0.d0
		if (collint_od_no_nunu) &
			dampTermMatrixCoeffNunu = 0.d0
		if (collint_damping_type.eq.0) then
			dampTermMatrixCoeffNue = 0.d0
			dampTermMatrixCoeffNunu = 0.d0
		endif
		write(*,*)"[config] Damping factors (Nue):"
		call printMat(dampTermMatrixCoeffNue)
		write(*,*)"[config] Damping factors (NuNu):"
		call printMat(dampTermMatrixCoeffNunu)

		write(*,*)"[config] Damping factors done."
	end subroutine setDampingFactors

	subroutine init_matrices
		real(dl), dimension(:), allocatable :: diag_el
		integer :: ix, iy
		real(dl) :: fdm, fdme
		real(dl) :: eps
		character(len=300) :: tmpstr

		allocate(diag_el(flavorNumber))
		nsi_epsilon = 0.d0

		!GL
		diag_el = gLmt
		diag_el(1) = gLe
		do ix=1, flavorNumber
			if (sterile(ix)) diag_el(ix) = 0.d0
		end do
		call createDiagMat(GL_mat, flavorNumber, diag_el)
#ifndef DO_TESTS
#ifdef FULL_F_AB
		do ix=1, flavorNumber
			do iy=ix, flavorNumber
				write(tmpstr,"('nsi_GL_',2I1)") ix, iy
				eps = read_ini_real(trim(tmpstr), 0.d0)
				GL_mat(ix, iy) = GL_mat(ix, iy) + eps
				GL_mat(iy, ix) = GL_mat(ix, iy)
				nsi_epsilon(ix, iy) = nsi_epsilon(ix, iy) + eps
				nsi_epsilon(iy, ix) = nsi_epsilon(ix, iy)
			end do
		end do
#endif
#endif
		write(*,*) "G_L:"
		call printMat(GL_mat)
		GLR_vec(1,:,:) = GL_mat

		!GR
		diag_el = gRemt
		do ix=1, flavorNumber
			if (sterile(ix)) diag_el(ix) = 0.d0
		end do
		call createDiagMat(GR_mat, flavorNumber, diag_el)
#ifndef DO_TESTS
#ifdef FULL_F_AB
		do ix=1, flavorNumber
			do iy=ix, flavorNumber
				write(tmpstr,"('nsi_GR_',2I1)") ix, iy
				eps = read_ini_real(trim(tmpstr), 0.d0)
				GR_mat(ix, iy) = GR_mat(ix, iy) + eps
				GR_mat(iy, ix) = GR_mat(ix, iy)
				nsi_epsilon(ix, iy) = nsi_epsilon(ix, iy) + eps
				nsi_epsilon(iy, ix) = nsi_epsilon(ix, iy)
			end do
		end do
#endif
#endif
		write(*,*) "G_R:"
		call printMat(GR_mat)
		GLR_vec(2,:,:) = GR_mat
#ifdef FULL_F_AB
		write(*,*) "NSI epsilon:"
		call printMat(nsi_epsilon)
#endif

		call setDampingFactors

		!lepton densities
		leptonDensities = 0.d0
		nuDensities%re = 0.d0
		nuDensities%im = 0.d0

		!identity (Matrix complex)
		call createIdentityMat(idMat, flavorNumber)

		!nu density matrix
		allocate(nuDensMatVec(Ny), nuDensMatVecFD(Ny))
		ntotrho = Ny*(flavNumSqu) !independent elements in nuDensity(y)
		ntot = ntotrho + 2 ! + w,z
		ixw = ntotrho+1
		ixz = ntotrho+2

		allocate(nuDensVec(ntot))
		if (save_BBN .and. trim(outputFolder).ne."")&
			call openFile(3978, trim(outputFolder)//"/y_grid.dat", .true.)
		if(save_fd .and. trim(outputFolder).ne."")&
			call openFile(3154, trim(outputFolder)//"/fd.dat", .true.)
		do ix=1, Ny
			call allocateCmplxMat(nuDensMatVec(ix))
			call allocateCmplxMat(nuDensMatVecFD(ix))
			nuDensMatVec(ix)%y = y_arr(ix)
			nuDensMatVecFD(ix)%y = y_arr(ix)
			nuDensMatVec(ix)%re(:,:) = 0.d0
			nuDensMatVec(ix)%im(:,:) = 0.d0
			fdme = fermiDirac(y_arr(ix)/z_in)
			fdm = fermiDirac(y_arr(ix))
			if (save_BBN .and. trim(outputFolder).ne."")&
				write(3978, multidblfmt) y_arr(ix)
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
		if(save_BBN .and. trim(outputFolder).ne."")&
			close(3978)

		deallocate(diag_el)

		!quantities used to save intermediate steps
		allocate(intermediateSteps%Heff(Ny))
		allocate(intermediateSteps%commutator(Ny))
		allocate(intermediateSteps%colltermsNue(Ny))
		allocate(intermediateSteps%colltermsNunu(Ny))
		allocate(intermediateSteps%yvec(ntot))
		allocate(intermediateSteps%ydot(ntot))
		do ix=1, Ny
			call allocateCmplxMat(intermediateSteps%Heff(ix))
			call allocateCmplxMat(intermediateSteps%commutator(ix))
			call allocateCmplxMat(intermediateSteps%colltermsNue(ix))
			call allocateCmplxMat(intermediateSteps%colltermsNunu(ix))
		end do
	end subroutine init_matrices

	subroutine init_fermions
#ifndef NO_INTERPOLATION
		call init_interp_dme2_e
#endif
		call fermions(1)%initialize("electrons", .true., 1.d0, 1d3)
		electrons => fermions(1)
#ifdef DO_MUONS
		call fermions(2)%initialize("muons", .false., m_mu_o_m_e, x_muon_cut)
		muons => fermions(2)
#endif
#ifndef NO_INTERPOLATION
		call init_interp_jkyg12
#endif
	end subroutine init_fermions

	subroutine finish_y_arrays
		integer :: ix

		if (allocated(dy_arr)) &
			deallocate(dy_arr)
		if (allocated(fy_arr)) &
			deallocate(fy_arr)
		if (allocated(feq_arr)) &
			deallocate(feq_arr)
		allocate(dy_arr(Ny), fy_arr(Ny), feq_arr(Ny))

		do ix=1, Ny-1
			dy_arr(ix) = y_arr(ix+1) - y_arr(ix)
		end do
		dy_arr(Ny) = 0.d0

		feq_arr = 0.d0
		do ix=1, Ny
			feq_arr(ix) = fermiDirac(y_arr(ix))
		end do
	end subroutine

	subroutine setInterp
		allocate(interp_xvec(interp_nx), interp_yvec(interp_ny), interp_zvec(interp_nz), interp_xozvec(interp_nxz))
		interp_xvec = logspace(interp_logx_in, logx_fin, interp_nx)
		interp_yvec = logspace(interp_logy_min, interp_logy_max, interp_ny)
		interp_zvec = linspace(interp_zmin, interp_zmax, interp_nz)
		low_xoz = very_early_x/interp_zmax
		interp_xozvec = logspace(log10(low_xoz), logx_fin, interp_nxz)
	end subroutine setInterp

	subroutine setXYVectors
		character(len=300) :: tmparg, tmpstr
		real(dl), dimension(:), allocatable :: fake

		z_in=1.d0
		if (x_in.lt.very_early_x) then
			write(tmparg, multidblfmt) very_early_x
			call criticalError("initial x is too early: x_in cannot be smaller than "//tmparg)
		end if
		allocate(x_arr(Nx))
		logx_in  = log10(x_in)
		logx_fin = log10(x_fin)
		x_arr = logspace(logx_in, logx_fin, Nx)

		if (use_gauss_laguerre) then
			write(tmpstr,"('[config] Configuring Gauss-Laguerre with Ny=',I3,' and y_max=',"//dblfmt//")") Ny, y_max
			call addToLog(trim(tmpstr))
			call get_GLq_vectors(Ny, y_arr, w_gl_arr, w_gl_arr2, .true., 3, y_max)
		else
			allocate(y_arr(Ny))
			y_arr = linspace(y_min, y_max, Ny)
			write(tmpstr,"('[config] Using Ny=',I3)") Ny
			call addToLog(trim(tmpstr))
		end if

		call finish_y_arrays
		call get_GLq_vectors(N_opt_y, opt_y, opt_y_w, fake, .true., 2, opt_y_cut)
		call get_GLq_vectors(N_opt_xoz, opt_xoz, opt_xoz_w, fake, .true., 2, opt_xoz_cut)
	end subroutine setXYVectors

	subroutine setNuProperties
		integer :: ix
		character(len=300) :: tmparg, tmpstr

		do ix=1, flavorNumber
			if(sterile(ix)) then
				Gs(ix) = 0.d0
			else
				Gs(ix) = 1.d0
			endif
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
	end subroutine setNuProperties

	subroutine initConfig()
		character(len=300) :: logname

		if (verbose>0) write(*,*) '[config] init configuration'
		logname = "ini.log"
		num_args = command_argument_count()

		allocate(args(num_args))
		do ixa = 1, num_args
			call get_command_argument(ixa,args(ixa))
		end do
		if(num_args.gt.1) then
			logname = trim(args(2))
		end if
		call initialOperations(logname)
		call checkIniFile(args(1))

		call readBasicParams
		call checkOutputFolder

		call checkNumThreads

		call readVerbosity

		call readXY
		call setXYVectors

		call readPrecision
		call setInterp

		call readOutputConfig

		call allocateStuff

		call readCollint
		call readFTQED

		call readNuDef
		call setNuProperties

		call init_fermions
		call zin_solver

		call readNuMixing

		call setMixingMatrix
		call setMassMatrix
		call init_matrices

		call endConfig
	end subroutine initconfig
end module fpConfig
