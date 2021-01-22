program prepare_interpolations
	use fpConfig
	implicit none
	real(dl), dimension(:), allocatable :: fake

#ifdef NO_INTERPOLATION
	call criticalError("cannot prepare values for interpolation if NO_INTERPOLATION is active!")
#endif
	logFile = "prepare_interpolations.log"
	call openLogFile
	interp_nx = interp_nx0
	interp_nz = interp_nz0
	interp_nxz = interp_nxz0
	allocate(interp_xvec(interp_nx), interp_yvec(interp_ny), interp_zvec(interp_nz), interp_xozvec(interp_nxz))
	interp_xvec = logspace(logx_in, logx_fin, interp_nx)
	interp_yvec = logspace(log10(y_min), log10(y_max), interp_ny)
	interp_zvec = linspace(interp_zmin, interp_zmax, interp_nz)
	interp_xozvec = logspace(log10(x_in/interp_zmax), logx_fin, interp_nxz)
	call get_GLq_vectors(N_opt_y, opt_y, opt_y_w, fake, .true., 2, opt_y_cut)
	call get_GLq_vectors(N_opt_xoz, opt_xoz, opt_xoz_w, fake, .true., 2, opt_xoz_cut)

	tests_interpolations = .false.
	ftqed_temperature_corr = .true.
	ftqed_ord3 = .true.
	ftqed_log_term = .false.
	ftqed_e_mth_leptondens = .true.
	call init_fermions

	ftqed_temperature_corr = .true.
	ftqed_ord3 = .false.
	ftqed_log_term = .false.
	ftqed_e_mth_leptondens = .false.
	call init_fermions

	ftqed_temperature_corr = .false.
	ftqed_ord3 = .false.
	ftqed_log_term = .false.
	ftqed_e_mth_leptondens = .false.
	call init_fermions

	call addToLog("preparation of values for interpolations completed!")
end program
