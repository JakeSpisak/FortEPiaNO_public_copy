module fortepianoWrapper
	use fpversion
	use precision
	use variables
	use fpCosmology
	use fpConfig
	use fpEquations
	use fpOutput
	implicit None

	contains

	!general
	function get_version()
		character (len=5) :: get_version
		get_version=version
	end function get_version

	!functions for input parameters
	subroutine pass_basic_params( &
		fn, &
		fr, &
		outfolder, &
		nthr &
	)
		integer, intent(in) :: fn, nthr
		logical, intent(in) :: fr
		character(len=*), intent(in) :: outfolder

		flavorNumber = fn
		force_replace = fr
		if (trim(outfolder)/="") outputFolder=outfolder
		num_threads = nthr

		call writelog_int("flavorNumber", flavorNumber)
		call writelog_logical("force_replace", force_replace)
		call writelog_char("outputFolder", trim(outputFolder))
		call writelog_int("num_threads", num_threads)
	end subroutine pass_basic_params

	subroutine pass_collint( &
		cdt, &
		cdz, &
		cod, &
		cend, &
		cnnd, &
		ceno, &
		cnno, &
		dsz &
	)
		integer, intent(in) :: cdt
		logical, intent(in) :: cdz, cod, cend, cnnd, ceno, cnno
		logical, dimension(:,:), intent(in) :: dsz
		integer :: ix, iy
		character(len=16) :: tmpstr

		collint_damping_type = cdt
		collint_diagonal_zero = cdz
		collint_offdiag_damping = cod
		collint_d_no_nue = cend
		collint_d_no_nunu = cnnd
		collint_od_no_nue = ceno
		collint_od_no_nunu = cnno

		call writelog_int("collint_damping_type", collint_damping_type)
		call writelog_logical("collint_diagonal_zero", collint_diagonal_zero)
		call writelog_logical("collint_offdiag_damping", collint_offdiag_damping)
		call writelog_logical("collint_d_no_nue", collint_d_no_nue)
		call writelog_logical("collint_d_no_nunu", collint_d_no_nunu)
		call writelog_logical("collint_od_no_nue", collint_od_no_nue)
		call writelog_logical("collint_od_no_nunu", collint_od_no_nunu)

		damping_read_zero = .false.
		do ix=1, flavorNumber
			do iy=ix+1, flavorNumber
				dampingSetZero(ix,iy) = dsz(ix, iy)
				write(tmpstr, "('dampingSetZero',2I1)") ix, iy
				call writelog_logical(tmpstr, dampingSetZero(ix,iy))
			end do
		end do
	end subroutine pass_collint

	subroutine pass_ftqed( &
		ftc, &
		ftl, &
		ft3, &
		fte &
	)
		logical, intent(in) :: ftc, ftl, ft3, fte
		ftqed_temperature_corr = ftc
		ftqed_log_term = ftl
		ftqed_ord3 = ft3
		ftqed_e_mth_leptondens = fte

		call writelog_logical("ftqed_temperature_corr", ftqed_temperature_corr)
		call writelog_logical("ftqed_log_term", ftqed_log_term)
		call writelog_logical("ftqed_ord3", ftqed_ord3)
		call writelog_logical("ftqed_e_mth_leptondens", ftqed_e_mth_leptondens)
	end subroutine pass_ftqed

	subroutine pass_low_reheating( &
		T &
	)
		real(kind=kind(0.d0)), intent(in) :: T

		Trh = T

		call writelog_real("Trh", Trh)
	end subroutine pass_low_reheating

	subroutine pass_nu_def( &
		fac, &
		ster &
	)
		real(kind=kind(0.d0)), intent(in), dimension(:) :: fac
		logical, intent(in), dimension(:) :: ster
		integer :: ix
		character(len=10) :: tmpstr

		do ix=1, flavorNumber
			nuFactor(ix) = fac(ix)
			write(tmpstr, "('nuFactor',I1)") ix
			call writelog_real(trim(tmpstr), nuFactor(ix))
			sterile(ix) = ster(ix)
			write(tmpstr, "('sterile',I1)") ix
			call writelog_logical(trim(tmpstr), sterile(ix))
		end do
	end subroutine pass_nu_def

	subroutine pass_nu_mixing( &
		gssq, &
		dm, &
		ith &
	)
		logical, intent(in) :: gssq
		real(kind=kind(0.d0)), dimension(:), intent(in) :: dm
		real(kind=kind(0.d0)), dimension(:,:), intent(in) :: ith
		real(kind=kind(0.d0)), dimension(:,:), allocatable :: th
		integer :: ix, iy
		character(len=16) :: tmpstr

		allocate(th(flavorNumber, flavorNumber))
		th=0.d0
		giveSinSq = gssq
		massSplittings = 0.d0
		massSplittings(2) = dm(2)
		if (flavorNumber .gt. 2) then
			massSplittings(3) = dm(3)
		end if
		do ix = 4, flavorNumber
			massSplittings(ix) = dm(ix)
		end do
		do ix = 2, flavorNumber
			massSplittings(ix) = dm(ix)
			write(tmpstr, "('dm',2I1)") ix, 1
			call writelog_real(trim(tmpstr), massSplittings(ix))
		end do
		if (giveSinSq) then
			th = asin(sqrt(ith))
		else
			th = ith
		end if
		mixingAngles(1,2) = th(1, 2)
		if (flavorNumber .gt. 2) then
			mixingAngles(1,3) = th(1, 3)
			mixingAngles(2,3) = th(2, 3)
		end if
		do ix=4, flavorNumber
			do iy=1, ix-1
				mixingAngles(iy, ix) = th(iy, ix)
			end do
		end do
		do ix=2, flavorNumber
			do iy=1, ix-1
				mixingAngles(iy, ix) = th(iy, ix)
				write(tmpstr, "('theta',2I1)") iy, ix
				call writelog_real(trim(tmpstr), mixingAngles(iy,ix))
			end do
		end do
	end subroutine pass_nu_mixing

	subroutine pass_output_config( &
		chk, &
		sbbn, &
		sfd, &
		seee, &
		snf, &
		snd, &
		snn, &
		sw, &
		sz, &
		si &
	)
		logical, intent(in) :: chk, sbbn, sfd, seee, snf, snd, snn, sw, sz, si

		checkpoint = chk
		save_BBN = sbbn
		save_fd = sfd
		save_energy_entropy_evolution = seee
		save_Neff = snf
		save_nuDens_evolution = snd
		save_number_evolution = snn
		save_w_evolution = sw
		save_z_evolution = sz
		intermediateSteps%output = si

		call writelog_logical("checkpoint", checkpoint)
		call writelog_logical("save_BBN", save_BBN)
		call writelog_logical("save_fd", save_fd)
		call writelog_logical("save_energy_entropy_evolution", save_energy_entropy_evolution)
		call writelog_logical("save_Neff", save_Neff)
		call writelog_logical("save_nuDens_evolution", save_nuDens_evolution)
		call writelog_logical("save_number_evolution", save_number_evolution)
		call writelog_logical("save_w_evolution", save_w_evolution)
		call writelog_logical("save_z_evolution", save_z_evolution)
		call writelog_logical("intermediateSteps%output", intermediateSteps%output)
	end subroutine pass_output_config

	subroutine pass_precision( &
		mi, &
		tjk, &
		daz, &
		dad, &
		dao, &
		dr &
	)
		integer, intent(in) :: mi
		real(kind=kind(0.d0)), intent(in) :: tjk, daz, dad, dao, dr
    
		maxiter = mi
		toler_jkyg = tjk
		dlsoda_atol_z = daz
		dlsoda_atol_d = dad
		dlsoda_atol_o = dao
		dlsoda_rtol = dr

		interp_nx = interp_nx0
		interp_nz = interp_nz0
		interp_nxz = interp_nxz0
		interp_zmin = interp_zmin0
		interp_zmax = interp_zmax0

		call writelog_int("maxiter", maxiter)
		call writelog_real("toler_jkyg", toler_jkyg)
		call writelog_real("dlsoda_atol_z", dlsoda_atol_z)
		call writelog_real("dlsoda_atol_d", dlsoda_atol_d)
		call writelog_real("dlsoda_atol_o", dlsoda_atol_o)
		call writelog_real("dlsoda_rtol", dlsoda_rtol)
		call writelog_int("interp_nx", interp_nx)
		call writelog_int("interp_nz", interp_nz)
		call writelog_int("interp_nxz", interp_nxz)
		call writelog_real("interp_zmin", interp_zmin)
		call writelog_real("interp_zmax", interp_zmax)
	end subroutine pass_precision

	subroutine pass_verbosity( &
		ver, &
		Npd &
	)
		integer, intent(in) :: ver
		real(kind=kind(0.d0)), intent(in) :: Npd

		verbose = ver
		Nprintderivs = Npd

		call writelog_int("verbose", verbose)
		call writelog_real("Nprintderivs", Nprintderivs)
	end subroutine pass_verbosity

	subroutine pass_xy( &
		xi, &
		xf, &
		iNx, &
		ugl, &
		iNy, &
		ym, &
		yx &
	)
		real(kind=kind(0.d0)), intent(in) :: xi, xf, ym, yx
		integer, intent(in) :: iNx, iNy
		logical, intent(in) :: ugl
        
		x_in  = xi
		x_fin = xf
		Nx = iNx
		use_gauss_laguerre = ugl
		Ny = iNy
		y_min = ym
		y_max = yx

		call writelog_real("x_in", x_in)
		call writelog_real("x_fin", x_fin)
		call writelog_int("Nx", Nx)
		call writelog_logical("use_gauss_laguerre", use_gauss_laguerre)
		call writelog_int("Ny", Ny)
		call writelog_real("y_min", y_min)
		call writelog_real("y_max", y_max)
	end subroutine pass_xy

	!functions for config
	subroutine w_allocate_stuff
		call allocateStuff
	end subroutine w_allocate_stuff

	subroutine w_check_num_threads
		call checkNumThreads
	end subroutine w_check_num_threads

	subroutine w_check_output_folder
		call checkOutputFolder
	end subroutine w_check_output_folder

	subroutine w_end_config
		call endConfig
	end subroutine w_end_config

	subroutine w_init_fermions
		call init_fermions
	end subroutine w_init_fermions

	subroutine w_init_matrices
		call init_matrices
	end subroutine w_init_matrices

	subroutine w_initial_operations( &
		logname &
	)
		character(len=*), intent(in) :: logname
		call initialOperations(logname)
	end subroutine w_initial_operations

	subroutine w_set_interp
		call setInterp
	end subroutine w_set_interp

	subroutine w_set_low_reheating
		call setLowReheating
	end subroutine w_set_low_reheating

	subroutine w_set_mass_matrix
		call setMassMatrix
	end subroutine w_set_mass_matrix

	subroutine w_set_mixing_matrix
		call setMixingMatrix
	end subroutine w_set_mixing_matrix

	subroutine w_set_nu_properties
		call setNuProperties
	end subroutine w_set_nu_properties

	subroutine w_set_xy_vectors
		call setXYVectors
	end subroutine w_set_xy_vectors

	!functions for const
	subroutine store_nu_dens_mat(vec)
		real(kind=kind(0.d0)), dimension(:), intent(in) :: vec
		call vec_2_densMat(vec)
	end subroutine store_nu_dens_mat

	!functions for cosmology
	function w_photon_number_density(z) result(r)
		real(kind=kind(0.d0)), intent(in) :: z
		real(kind=kind(0.d0)) :: r
		r=photonNumberDensity(z)
	end function w_photon_number_density

	function w_photon_density(z) result(r)
		real(kind=kind(0.d0)), intent(in) :: z
		real(kind=kind(0.d0)) :: r
		r=photonDensity(z)
	end function w_photon_density

	!functions for equations
	subroutine w_solver
		call solver
	end subroutine w_solver

	subroutine w_zin_solver
		call zin_solver
	end subroutine w_zin_solver

	!functions for output
	subroutine w_final_results
		call finalresults
	end subroutine w_final_results
end module fortepianoWrapper
