program prepare_interpolations
	use fpConfig
	implicit none
	real(dl), dimension(:), allocatable :: fake
	real(dl), dimension(2) :: startx
	integer :: ix

#ifdef NO_INTERPOLATION
	call criticalError("cannot prepare values for interpolation if NO_INTERPOLATION is active!")
#endif
	logFile = "prepare_interpolations.log"
	call openLogFile
	interp_nx = interp_nx0
	interp_nz = interp_nz0
	interp_nxz = interp_nxz0
	interp_zmin = interp_zmin0
	interp_zmax = interp_zmax0
	toler_jkyg = 1.d-7
	maxiter = 100
	tot_factor_active_nu = 3.0
#ifndef XIN
#define XIN 0.001d0
#endif
#ifndef XFIN
#define XFIN 35.0d0
#endif
#ifndef YMIN
#define YMIN 0.01d0
#endif
#ifndef YMAX
#define YMAX 20.d0
#endif
#ifndef STARTX
#define STARTX 0.001d0
#endif
	x_in    = XIN
	x_fin   = XFIN
	logx_in  = log10(x_in)
	logx_fin = log10(x_fin)
	y_min = YMIN
	y_max = YMAX
	startx = (/STARTX, very_early_x/)
	allocate(interp_xvec(interp_nx), interp_yvec(interp_ny), interp_zvec(interp_nz), interp_xozvec(interp_nxz))
	interp_xvec = logspace(logx_in, logx_fin, interp_nx)
	interp_yvec = logspace(log10(y_min), log10(y_max), interp_ny)
	interp_zvec = linspace(interp_zmin, interp_zmax, interp_nz)
	call get_GLq_vectors(N_opt_y, opt_y, opt_y_w, fake, .true., 2, opt_y_cut)
	call get_GLq_vectors(N_opt_xoz, opt_xoz, opt_xoz_w, fake, .true., 2, opt_xoz_cut)
	do ix = 1, 2
		low_xoz = startx(ix)/interp_zmax
		interp_xozvec = logspace(log10(low_xoz), logx_fin, interp_nxz)

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
	end do

	call addToLog("preparation of values for interpolations completed!")
end program
