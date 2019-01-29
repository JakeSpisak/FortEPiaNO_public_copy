subroutine assert_double(testname, num1, num2, tol)
	use precision
	implicit none
	character(len=*) :: testname
	real(dl) :: num1, num2, tol

	if (abs(num1 - num2) .lt. tol) then
		write(*, fmt="(a)", advance="no") "."
	else
		print *, testname, " failed"
		write (*,"(*(E17.9))") num1, num2, num1 - num2, tol
		call exit()
	end if
end subroutine assert_double

subroutine assert_double_rel(testname, num1, num2, tol)
	use precision
	implicit none
	character(len=*) :: testname
	real(dl) :: num1, num2, tol

	if (abs((num1 - num2)/num1) .lt. tol) then
		write(*, fmt="(a)", advance="no") "."
	else
		print *, testname, " failed"
		write (*,"(*(E17.9))") num1, num2, (num1 - num2)/num1, tol
		call exit()
	end if
end subroutine assert_double_rel

subroutine assert_double_rel_verb(testname, num1, num2, tol)
	use precision
	implicit none
	character(len=*) :: testname
	real(dl) :: num1, num2, tol

	if (abs((num1 - num2)/num1) .lt. tol) then
		write(*, fmt="(a,E14.7)") testname, abs((num1 - num2)/num1)
	else
		print *, testname, " failed"
		write (*,"(*(E17.9))") num1, num2, (num1 - num2)/num1, tol
		call exit()
	end if
end subroutine assert_double_rel_verb

program tests
	use precision
	use ndConfig
	use ndErrors
	use ndEquations
	implicit none

	call openLogFile
	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "Initializations"
	call do_tests_initialization
	call init_interp_jkyg12
	call init_interp_ElDensity
	call init_interp_dme2_e
	call allocate_interpNuDens

	write(*,*) ""
	write(*,"(a)") "Starting tests"
	call do_tests_cosmology
	call do_tests_JKYG
	call do_tests_dzodx
	call do_tests_dme2
	call do_tests_Di
	call do_tests_Pi_ij
	call do_f_ann_sc_re_tests_eq
	call do_f_ann_sc_re_tests_full
	call do_tests_coll_int

	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "all tests were successfull!"

	write(*,*) ""
	write(*,"(a)") "now doing some timing tests..."
	call do_timing_tests()
	call closeLogFile

	contains

	subroutine do_tests_initialization
		integer :: ix
		checkpoint = .true.
		maxiter = 100
		toler   = 1.d-5
		toler_dme2 = 1.d-5
		toler_jkyg = 1.d-7
		toler_ed = 1.d-4
		dlsoda_atol = 1.d-6
		dlsoda_rtol = 1.d-6
		Nx = 100
		Ny = 100
		allocate(x_arr(Nx), y_arr(Ny), logy_arr(Ny))
		allocate(dy_arr(Ny), fy_arr(Ny))
		x_in    = 0.01d0
		x_fin   = 40.d0
		logx_in  = log10(x_in)
		logx_fin = log10(x_fin)
		x_arr = logspace(logx_in, logx_fin, Nx)
		y_min   = 0.01d0
		y_max   = 20.0d0
		logy_min = log10(y_min)
		logy_max = log10(y_max)
#ifdef LOGY
		y_arr = logspace(logy_min, logy_max, Ny)
		logy_arr = log10(y_arr)
#else
		y_arr = linspace(y_min, y_max, Ny)
		logy_arr = log10(y_arr)
#endif
		do ix=1, Ny-1
			dy_arr(ix) = y_arr(ix+1) - y_arr(ix)
		end do
		dy_arr(Ny) = 0.d0
		z_in = 1.00003d0
		m_lightest   = 0.0d0
		massOrdering = .true.
		!read cosmo parameters
		hubbleParam = i_HubbleParam
		photonTemperatureToday = i_photonTempToday
		!settings for collisional
		collision_offdiag = 1
		dme2_temperature_corr = .true.
		flavorNumber = 3
		flavNumSqu = flavorNumber**2
		call allocateStuff
		do ix=1, flavorNumber
			nuFactor(ix) = 1.d0
			sterile(ix) = 0.d0
		end do
		theta12      = i_theta12
		dm12         = i_dm12
		if (flavorNumber .gt. 2) then
			theta13      = i_theta13
			theta23      = i_theta23
			dm23         = i_dm23
			deltaCP13    = i_deltaCP13
		end if
		call setMixingMatrix()
		call setMassMatrix()
		call init_matrices
		allocate(interp_xvec(interp_nx), interp_zvec(interp_nz), interp_xozvec(interp_nx))
		interp_xvec = logspace(logx_in, logx_fin, interp_nx)
		interp_zvec = linspace(interp_zmin, interp_zmax, interp_nz)
		interp_xozvec = logspace(log10(x_in/interp_zmax), logx_fin, interp_nx)
	end subroutine do_tests_initialization

	subroutine do_tests_cosmology
		real(dl), dimension(:), allocatable :: ndmv_re
		integer :: i, iy

		allocate(ndmv_re(Ny))
		do i=1, flavorNumber
			do iy=1, Ny
				ndmv_re(iy) = 1.d0*i
				nuDensMatVec(iy)%re(i, i) = 1.d0*i
			end do
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
		end do

		write(*,*) ""
		write(*,"(a)") "Cosmology (35 tests)"
		call assert_double_rel("elDensF test 1", electronDensityFull(1.d0, 1.d0), 1.06283d0, 1d-4)
		call assert_double_rel("elDensF test 2", electronDensityFull(0.076d0, 1.32d0), 3.49493d0, 1d-4)
		call assert_double_rel("elDensF test 3", electronDensityFull(1.d1, 1.2d0), 0.0377723d0, 1d-4)
		call assert_double_rel("elDensF test 4", electronDensityFull(2.d1, 1.2d0), 0.0000421964d0, 1d-4)
		call assert_double_rel("elDensF test 5", electronDensityFull(3.d1, 1.2d0), 2.61468d-8, 5d-3)
		call assert_double_rel("elDens test 1", electronDensity(1.d0, 1.d0), 1.06283d0, 1d-4)
		call assert_double_rel("elDens test 2", electronDensity(0.076d0, 1.32d0), 3.49493d0, 1d-4)
		call assert_double_rel("elDens test 3", electronDensity(1.d1, 1.2d0), 0.0377723d0, 1d-2)
		call assert_double_rel("elDens test 4", electronDensity(2.d1, 1.2d0), 0.0000421964d0, 2d-2)
		call assert_double_rel("elDens test 5", electronDensity(3.d1, 1.2d0), 2.61468d-8, 5d-2)
		call assert_double_rel("photDens test 1", photonDensity(1.002d0), 0.66325322d0, 1d-7)
		call assert_double_rel("photDens test 2", photonDensity(1.34d0), 2.12142498d0, 1d-7)

		call assert_double_rel("nuDens test 1", nuDensity(1.d0, 1), 1.15145d0, 1d-4)
		call assert_double_rel("nuDens test 2", nuDensity(1.076d0, 1), 1.54346d0, 1d-4)
		call assert_double_rel("nuDens test 3", nuDensity(1.32d0, 1), 3.49577d0, 3d-4)
		call assert_double_rel("nuDens test 4", nuDensity(1.37d0, 2), 2.d0*4.05629d0, 5d-4)
		call assert_double_rel("nuDens test 5", nuDensity(1.003d0, 3), 3.d0*1.16533d0, 1d-4)

		call assert_double_rel("nuDensLin test 1", nuDensityLin(1.d0, 1), 1.15145d0, 1d-4)
		call assert_double_rel("nuDensLin test 2", nuDensityLin(1.076d0, 1), 1.54346d0, 1d-4)
		call assert_double_rel("nuDensLin test 3", nuDensityLin(1.32d0, 1), 3.49577d0, 2d-4)
		call assert_double_rel("nuDensLin test 4", nuDensityLin(1.37d0, 2), 2.d0*4.05629d0, 5d-4)
		call assert_double_rel("nuDensLin test 5", nuDensityLin(1.003d0, 3), 3.d0*1.16533d0, 1d-4)

		call assert_double_rel("nuDensLinEq test 1", nuDensityLinEq(1.d0), 1.15145d0, 1d-4)
		call assert_double_rel("nuDensLinEq test 2", nuDensityLinEq(1.37d0), 4.055034d0, 1d-4)

		call assert_double_rel("allNuDensLin test 1", allNuDensity(1.d0), 6*1.15145d0, 1d-3)
		call assert_double_rel("allNuDensLin test 2", allNuDensity(1.076d0), 6*1.54346d0, 1d-3)
		call assert_double_rel("allNuDensLin test 3", allNuDensity(1.32d0), 6*3.49577d0, 1d-3)

		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVec(iy)%re(i, i) = 1.d0
			end do
		end do
		call assert_double("radDens test 1", radDensity(0.7d0, 1.04d0), 6.11141d0, 1d-3)
		call assert_double("radDens test 2", radDensity(1.d0, 1.04d0), 6.06205d0, 1d-3)
		call assert_double("radDens test 3", radDensity(1.d0, 1.24d0), 12.309d0, 3d-3)

		call assert_double_rel("Neff test 1", Neff_from_rho_z(zid), 3.0d0, 1d-5)
		call assert_double_rel("Neff test 2", Neff_from_rho_z(1.39975d0), 3.011d0, 1d-4)
		do iy=1, Ny
			nuDensMatVec(iy)%re(1, 1) = 1.2d0
		end do
		call assert_double_rel("Neff test 3", Neff_from_rho_z(zid), 3.2d0, 1d-5)
		do iy=1, Ny
			nuDensMatVec(iy)%re(1, 1) = 1.00920d0
			nuDensMatVec(iy)%re(2, 2) = 1.00392d0
			nuDensMatVec(iy)%re(3, 3) = 1.00392d0
		end do
		call assert_double("Neff test 4", Neff_from_rho_z(1.397843d0), 3.045d0, 1d-3)
		do iy=1, Ny
			nuDensMatVec(iy)%re(1, 1) = 1.00699d0
			nuDensMatVec(iy)%re(2, 2) = 1.00511d0
			nuDensMatVec(iy)%re(3, 3) = 1.00519d0
		end do
		call assert_double("Neff test 5", Neff_from_rho_z(1.39779d0), 3.045d0, 1d-3)

		deallocate(ndmv_re)
	end subroutine do_tests_cosmology

	subroutine do_tests_JKYG
		real(dl), dimension(2) :: res
		write(*,*) ""
		write(*,"(a)") "JKYG (27 tests)"
		call assert_double_rel("J test 1", J_funcFull(0.01d0), 0.1666641d0, 1d-6)
		call assert_double_rel("J test 2", J_funcFull(1.d0), 0.1437972d0, 1d-6)
		call assert_double_rel("J test 3", J_funcFull(5.d0), 0.01339351d0, 1d-6)
		call assert_double_rel("J' test 1", JprimeFull(0.01d0), -0.0005065951d0, 1d-6)
		call assert_double_rel("J' test 2", JprimeFull(1.d0), -0.04125339d0, 1d-6)
		call assert_double_rel("J' test 3", JprimeFull(5.d0), -0.01015141d0, 1d-6)
		call assert_double_rel("K test 1", K_funcFull(0.01d0), 0.083319d0, 1d-6)
		call assert_double_rel("K test 2", K_funcFull(1.d0), 0.0550046d0, 1d-6)
		call assert_double_rel("K test 3", K_funcFull(5.d0), 0.00204432d0, 1d-6)
		call assert_double_rel("K' test 1", KprimeFull(0.01d0), -0.00262052d0, 1d-6)
		call assert_double_rel("K' test 2", KprimeFull(1.d0), -0.03378793d0, 1d-6)
		call assert_double_rel("K' test 3", KprimeFull(5.d0), -0.001860974d0, 1d-6)
		call assert_double_rel("Y test 1", Y_funcFull(0.01d0), 2.302883d0, 1d-6)
		call assert_double_rel("Y test 2", Y_funcFull(1.d0), 2.070647d0, 1d-6)
		call assert_double_rel("Y test 3", Y_funcFull(5.d0), 0.3145333d0, 1d-6)

		res = G12_funcFull(0.01d0)
		call assert_double_rel("G1 test 1", res(1), -0.0000658825d0, 3d-5)
		call assert_double_rel("G2 test 1", res(2), -0.0095518d0, 3d-5)
		res = G12_funcFull(1.d0)
		call assert_double_rel("G1 test 2", res(1), -0.00115846d0, 1d-5)
		call assert_double_rel("G2 test 2", res(2), -0.00806502d0, 1d-5)
		res = G12_funcFull(5.d0)
		call assert_double_rel("G1 test 3", res(1), -0.000108111d0, 1d-5)
		call assert_double_rel("G2 test 3", res(2), -0.000945107d0, 1d-5)

		res = dzodxcoef_interp_func(0.01d0)
		call assert_double_rel("A test 1", res(1), 0.00044351d0, 1d-5)
		call assert_double_rel("B test 1", res(2), 0.0140361d0, 1d-5)
		res = dzodxcoef_interp_func(1.d0)
		call assert_double_rel("A test 2", res(1), 0.0404956d0, 1d-5)
		call assert_double_rel("B test 2", res(2), 0.0143827d0, 1d-5)
		res = dzodxcoef_interp_func(5.d0)
		call assert_double_rel("A test 3", res(1), 0.034036d0, 1d-5)
		call assert_double_rel("B test 3", res(2), 0.0257897d0, 1d-5)
	end subroutine do_tests_JKYG

	subroutine do_tests_dzodx
		integer, parameter :: n=901
		real(dl), dimension(n) :: ydot
		integer :: m

		ydot = 0.d0
		do m=1, Ny
			ydot((m-1)*flavNumSqu + 1) = cos(0.02d0*y_arr(m))
			ydot((m-1)*flavNumSqu + 2) = y_arr(m)/20.d0
			ydot((m-1)*flavNumSqu + 3) = 1.d0
		end do
		write(*,*) ""
		write(*,"(a)") "dz/dx functions (X tests)"
		call dz_o_dx(0.01d0, 1.d0, ydot, n)
		call assert_double("dz_o_dx test 1", ydot(n), -0.17511d0, 1d-5)
		call dz_o_dx(1.d0/1.1d0, 1.1d0, ydot, n)
		call assert_double("dz_o_dx test 2", ydot(n), -0.159145d0, 5d-3)
		call dz_o_dx(6.d0, 1.2d0, ydot, n)
		call assert_double("dz_o_dx test 3", ydot(n), -0.359906d0, 1d-5)

		call dz_o_dx_lin(0.01d0, 1.d0, ydot, n)
		call assert_double("dz_o_dx_lin test 1", ydot(n), -0.17511d0, 1d-5)
		call dz_o_dx_lin(1.d0/1.1d0, 1.1d0, ydot, n)
		call assert_double("dz_o_dx_lin test 2", ydot(n), -0.159145d0, 5d-3)
		call dz_o_dx_lin(6.d0, 1.2d0, ydot, n)
		call assert_double("dz_o_dx_lin test 3", ydot(n), -0.359906d0, 1d-5)
	end subroutine do_tests_dzodx

	subroutine do_tests_dme2
		write(*,*) ""
		write(*,"(a)") "dme2 (14 tests)"
		call assert_double("dme2F test 1", dme2_electronFull(0.05d0, 0.d0, 1.0003d0), 0.02292d0, 1d-5)
		call assert_double("dme2F test 2", dme2_electronFull(0.05d0, 100.d0, 1.0003d0), 0.02292d0, 1d-5)
		call assert_double("dme2F test 3", dme2_electronFull(0.5d0, 0.d0, 1.1d0), 0.0283696d0, 1d-5)
		call assert_double("dme2F test 4", dme2_electronFull(1.23d0, 0.d0, 1.198d0), 0.0321547d0, 1d-5)
		call assert_double("dme2F test 5", dme2_electronFull(7.6d0, 0.d0, 1.3d0), 0.0260756d0, 1d-5)
		call assert_double("dme2F test 6", dme2_electronFull(35.d0, 0.d0, 1.39d0), 0.0295293d0, 1d-5)
		call assert_double("dme2 test 1", dme2_electron(0.05d0, 0.d0, 1.0003d0), 0.02292d0, 1d-5)
		call assert_double("dme2 test 2", dme2_electron(0.05d0, 100.d0, 1.0003d0), 0.02292d0, 1d-5)
		call assert_double("dme2 test 3", dme2_electron(0.5d0, 0.d0, 1.1d0), 0.0283696d0, 1d-5)
		call assert_double("dme2 test 4", dme2_electron(1.23d0, 0.d0, 1.198d0), 0.0321547d0, 1d-5)
		call assert_double("dme2 test 5", dme2_electron(7.6d0, 0.d0, 1.3d0), 0.0260756d0, 1d-5)
		call assert_double("dme2 test 6", dme2_electron(35.d0, 0.d0, 1.39d0), 0.0295293d0, 1d-5)
		call assert_double("Ebare_i_dme test 1", Ebare_i_dme(0.3d0, 0.4d0, 1.44d0), 1.3d0, 1d-7)
		call assert_double("Ebare_i_dme test 2", Ebare_i_dme(3.d0, 7.d0, 22.d0), 8.944272d0, 1d-7)
	end subroutine do_tests_dme2

	subroutine do_tests_Di
		write(*,*) ""
		write(*,"(a)") "D_i functions (18 tests)"
		call assert_double("D1 test 1", D1_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.4d0, 1d-7)
		call assert_double("D1 test 2", D1_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.4d0, 1d-7)
		call assert_double("D1 test 3", D1_full(0.01d0,5.d0,2.6d0,2.41d0), 0.04d0, 1d-7)
		call assert_double("D1 test 4", D1_full(10.03d0,5.4d0,8.8d0,6.63d0), 21.6d0, 1d-4)
		call assert_double("D1 test 5", D1_full(10.d0,2.d0,5.3d0,6.7d0), 8.d0, 1d-7)
		call assert_double("D1 test 6", D1_full(2.d0,10.d0,6.7d0,5.3d0), 8.d0, 1d-7)
		call assert_double("D2 test 1", D2_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), -0.00133333d0, 1d-7)
		call assert_double("D2 test 2", D2_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), -0.0213333d0, 1d-7)
		call assert_double("D2 test 3", D2_full(0.01d0,5.d0,2.6d0,2.41d0), -1.333333333d-6, 1d-12)
		call assert_double("D2 test 4", D2_full(10.03d0,5.4d0,8.8d0,6.63d0), -209.952d0, 1d-4)
		call assert_double_rel("D2 test 5", D2_full(10.d0,2.d0,5.3d0,6.7d0), -10.666667d0, 1d-7)
		call assert_double_rel("D2 test 6", D2_full(2.d0,10.d0,6.7d0,5.3d0), -10.666667d0, 1d-7)
		call assert_double("D3 test 1", D3_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.000192d0, 1d-7)
		call assert_double("D3 test 2", D3_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.000192d0, 1d-7)
		call assert_double("D3 test 3", D3_full(0.01d0,5.d0,2.6d0,2.41d0), 2.50454d-5, 1d-10)
		call assert_double_rel("D3 test 4", D3_full(10.03d0,5.4d0,8.8d0,6.63d0), 22692.22d0, 1d-7)
		call assert_double_rel("D3 test 5", D3_full(10.d0,2.d0,5.3d0,6.7d0), 918.2933d0, 1d-7)
		call assert_double_rel("D3 test 6", D3_full(2.d0,10.d0,6.7d0,5.3d0), 918.2933d0, 1d-7)
	end subroutine do_tests_Di

	subroutine do_tests_Pi_ij
		real(dl), dimension(2) :: temp_v2
		write(*,*) ""
		write(*,"(a)") "Pi(yi,yj) functions (24 tests)"
		call assert_double("Pi_1_12 test 1", PI1_12_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.00933333d0, 1d-7)
		call assert_double("Pi_1_12 test 2", PI1_12_full(0.4d0, 0.3d0, 0.2d0, 0.1d0), 0.0893333d0, 1d-7)
		call assert_double_rel("Pi_1_12 test 3", PI1_12_full(10.d0, 2.d0, 5.3d0, 6.7d0), 170.66667d0, 1d-7)
		call assert_double_rel("Pi_1_12 test 4", PI1_12_full(2.d0, 10.d0, 6.7d0, 5.3d0), 170.66667d0, 1d-7)

		call assert_double_rel("Pi_1_13 test 1", PI1_13_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.010666667d0, 1d-7)
		call assert_double_rel("Pi_1_13 test 2", PI1_13_full(0.4d0, 0.3d0, 0.2d0, 0.1d0), 0.010666667d0, 1d-7)
		call assert_double_rel("Pi_1_13 test 3", PI1_13_full(10.d0, 2.d0, 5.3d0, 6.7d0), 96.53333d0, 1d-7)
		call assert_double_rel("Pi_1_13 test 4", PI1_13_full(2.d0, 10.d0, 6.7d0, 5.3d0), 96.53333d0, 1d-7)

		temp_v2 = PI2_nn_f(0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.3d0, 0.4d0)
		call assert_double("Pi_2_14 test 1", temp_v2(1), 0.00267733d0, 1d-7)
		call assert_double("Pi_2_13 test 1", temp_v2(2), 0.000810667d0, 1d-7)
		temp_v2 = PI2_nn_f(0.4d0, 0.3d0, 0.2d0, 0.1d0, 0.2d0, 0.1d0)
		call assert_double("Pi_2_14 test 2", temp_v2(1), 0.00267733d0, 1d-7)
		call assert_double("Pi_2_13 test 2", temp_v2(2), 0.000810667d0, 1d-7)
		temp_v2 = PI2_nn_f(10.d0, 2.d0, 5.3d0, 6.7d0, 5.3d0, 6.7d0)
		call assert_double_rel("Pi_2_14 test 3", temp_v2(1), 1978.88d0, 1d-7)
		call assert_double_rel("Pi_2_13 test 3", temp_v2(2), 3293.01333d0, 1d-7)
		temp_v2 = PI2_nn_f(2.d0, 10.d0, 6.7d0, 5.3d0, 6.7d0, 5.3d0)
		call assert_double_rel("Pi_2_14 test 4", temp_v2(1), 1978.88d0, 1d-7)
		call assert_double_rel("Pi_2_13 test 4", temp_v2(2), 3293.01333d0, 1d-7)

		temp_v2 = PI2_ne_f(0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.2d0, 0.4d0)
		call assert_double("Pi_2_14 test 1", temp_v2(1), 0.00267733d0, 1d-7)
		call assert_double("Pi_2_12 test 1", temp_v2(2), 0.00427733d0, 1d-7)
		temp_v2 = PI2_ne_f(0.4d0, 0.3d0, 0.2d0, 0.1d0, 0.3d0, 0.1d0)
		call assert_double("Pi_2_14 test 2", temp_v2(1), 0.00267733d0, 1d-7)
		call assert_double("Pi_2_12 test 2", temp_v2(2), 0.00427733d0, 1d-7)
		temp_v2 = PI2_ne_f(10.d0, 2.d0, 5.3d0, 6.7d0, 2.d0, 6.7d0)
		call assert_double("Pi_2_14 test 3", temp_v2(1), 1978.88d0, 1d-7)
		call assert_double("Pi_2_12 test 3", temp_v2(2), 9420.8d0, 1d-7)
		temp_v2 = PI2_ne_f(2.d0, 10.d0, 6.7d0, 5.3d0, 10.d0, 5.3d0)
		call assert_double("Pi_2_14 test 4", temp_v2(1), 1978.88d0, 1d-7)
		call assert_double("Pi_2_12 test 4", temp_v2(2), 9420.8d0, 1d-7)
	end subroutine do_tests_Pi_ij

	subroutine do_f_ann_sc_re_tests_eq
		integer :: ix, iy
		real(dl) :: fdA, fdB, f1, f2, f3
		type(cmplxMatNN) :: nA, nB
		real(dl), dimension(3) :: tmparr
		character(len=300) :: tmparg

		f1 = fermiDirac(0.3d0)
		f2 = fermiDirac(0.4d0)
		f3 = fermiDirac(0.1d0)
		write(*,*) ""
		write(*,"(a)") "F_ann_re functions at equilibrium (144 tests)"
		call allocateCmplxMat(nA)
		call allocateCmplxMat(nB)
		fdA = fermiDirac(0.1d0)
		fdB = fermiDirac(0.2d0)
		nA%re = 0.d0
		nA%im = 0.d0
		nB%re = 0.d0
		nB%im = 0.d0
		do ix=1, flavorNumber
			nA%re(ix,ix) = fdA
			nB%re(ix,ix) = fdB
		end do
		!F_ab_ann_re(0.d0, 1.d0, n1,n2, y3, y4, 2,2, ix, iy) !L=1, R=2 !ix, iy in 1...3
		!FRReq[0.1, 0.2, 0.3, 0.4]={{-0.00259399, 0., 0.}, {0., -0.00259399, 0.}, {0., 0., -0.00259399}}
		tmparr = (/-0.00259399, -0.00259399, -0.00259399/)
		do ix=1, flavorNumber
			write(tmparg,"('FRR_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.1, 0.2, 0.3, 0.4]={{-0.0259319, 0., 0.}, {0., -0.00350124, 0.}, {0., 0., -0.00350124}}
		tmparr = (/-0.0259319, -0.00350124, -0.00350124/)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
		tmparr = (/-0.00820165, 0.00301367, 0.00301367/)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
		tmparr = (/-0.00820165, 0.00301367, 0.00301367/)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!second series
		fdA = fermiDirac(0.4d0)
		do ix=1, flavorNumber
			nA%re(ix,ix) = fdA
		end do
		!FRReq[0.4, 0.2, 0.3, 0.1]={{0.00129054, 0., 0.}, {0., 0.00129054, 0.}, {0., 0., 0.00129054}}
		tmparr = (/0.00129054, 0.00129054, 0.00129054/)
		do ix=1, flavorNumber
			write(tmparg,"('FRR_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,2, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,2, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.4, 0.2, 0.3, 0.1]={{0.0129014, 0., 0.}, {0., 0.00174191, 0.}, {0., 0., 0.00174191}}
		tmparr = (/0.0129014, 0.00174191, 0.00174191/)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,1, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,1, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
		tmparr = (/0.00408041, -0.00149933, -0.00149933/)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,1, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,1, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
		tmparr = (/0.00408041, -0.00149933, -0.00149933/)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,2, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,2, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,2, iy, ix), 0.d0, 1d-7)
			end do
		end do


		f1 = fermiDirac(0.2d0)
		f2 = fermiDirac(0.4d0)
		f3 = fermiDirac(0.1d0)
		write(*,*)
		write(*,"(a)") "F_sc_re functions at equilibrium (144 tests)"
		fdA = fermiDirac(0.1d0)
		fdB = fermiDirac(0.3d0)
		do ix=1, flavorNumber
			nA%re(ix,ix) = fdA
			nB%re(ix,ix) = fdB
		end do
		!FRReq[0.1, 0.2, 0.3, 0.4]={{-0.00259399, 0., 0.}, {0., -0.00259399, 0.}, {0., 0., -0.00259399}}
		tmparr = (/-0.00259399, -0.00259399, -0.00259399/)
		do ix=1, flavorNumber
			write(tmparg,"('FRR_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 2,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 2,2, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 2,2, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.1, 0.2, 0.3, 0.4]={{-0.0259319, 0., 0.}, {0., -0.00350124, 0.}, {0., 0., -0.00350124}}
		tmparr = (/-0.0259319, -0.00350124, -0.00350124/)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 1,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 1,1, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 1,1, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
		tmparr = (/-0.00820165, 0.00301367, 0.00301367/)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 2,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 2,1, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 2,1, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
		tmparr = (/-0.00820165, 0.00301367, 0.00301367/)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 1,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 1,2, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 1,2, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 1,2, iy, ix), 0.d0, 1d-7)
			end do
		end do

		fdA = fermiDirac(0.4d0)
		do ix=1, flavorNumber
			nA%re(ix,ix) = fdA
		end do
		!FRReq[0.4, 0.2, 0.3, 0.1]={{0.00129054, 0., 0.}, {0., 0.00129054, 0.}, {0., 0., 0.00129054}}
		tmparr = (/0.00129054, 0.00129054, 0.00129054/)
		do ix=1, flavorNumber
			write(tmparg,"('FRR_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 2,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 2,2, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 2,2, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.4, 0.2, 0.3, 0.1]={{0.0129014, 0., 0.}, {0., 0.00174191, 0.}, {0., 0., 0.00174191}}
		tmparr = (/0.0129014, 0.00174191, 0.00174191/)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 1,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 1,1, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 1,1, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
		tmparr = (/0.00408041, -0.00149933, -0.00149933/)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 2,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 2,1, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 2,1, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
		tmparr = (/0.00408041, -0.00149933, -0.00149933/)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 1,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 1,2, ix, iy), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f3, 1,2, iy, ix), 0.d0, 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f3, 1,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		call deallocateCmplxMat(nA)
		call deallocateCmplxMat(nB)
	end subroutine do_f_ann_sc_re_tests_eq

	subroutine do_f_ann_sc_re_tests_full
		integer :: a, b, ix, iy
		real(dl) :: fdA, fdB, f1, f2, f3
		type(cmplxMatNN) :: nA, nB
		character(len=300) :: tmparg
		real(dl), dimension(3) :: tmparrA, tmparrB
		real(dl), dimension(3, 3) :: tmpmatA, tmpmatB

		f1 = fermiDirac(0.1d0)
		f2 = fermiDirac(0.7d0)
		f3 = fermiDirac(1.9d0)
		write(*,*)
		write(*,"(a)") "F_ann functions, empty rho (144 tests)"
		call allocateCmplxMat(nA)
		call allocateCmplxMat(nB)
		nA%re = 0.d0
		nA%im = 0.d0
		nB%re = 0.d0
		nB%im = 0.d0
		do a=1, 2
			do b=1, 2
				do ix=1, flavorNumber
					iy=ix
					write(tmparg,"('F_ann',2I1,' test 1 empty rho ',2I1)") a, b, ix, iy
					call assert_double(trim(tmparg)//" im", F_ab_ann_im(nA, nB, f1, f2, a, b, ix, iy), 0.0d0, 1d-7)
					do iy=ix+1, flavorNumber
						write(tmparg,"('F_ann',2I1,' test 1 empty rho ',2I1)") a, b, ix, iy
						call assert_double(trim(tmparg)//" re", F_ab_ann_re(nA, nB, f1, f2, a, b, ix, iy), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_ann_im(nA, nB, f1, f2, a, b, ix, iy), 0.0d0, 1d-7)
						write(tmparg,"('F_ann',2I1,' test 1 empty rho ',2I1)") a, b, iy, ix
						call assert_double(trim(tmparg)//" re", F_ab_ann_re(nA, nB, f1, f2, a, b, iy, ix), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_ann_im(nA, nB, f1, f2, a, b, iy, ix), 0.0d0, 1d-7)
					end do
				end do
			end do
		end do
!		LL
		tmparrA = (/0.1685832, 0.0227616, 0.0227616/)
		do ix=1, flavorNumber
			write(tmparg,"('F_ann11 test 1 empty rho ',2I1)") ix,ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, ix,ix), tmparrA(ix), 1d-7)
		end do
!		RR
		tmparrA = (/0.0168635, 0.0168635, 0.0168635/)
		do ix=1, flavorNumber
			write(tmparg,"('F_ann22 test 1 empty rho ',2I1)") ix,ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, ix,ix), tmparrA(ix), 1d-7)
		end do
!		LR and RL
		tmparrA = (/0.0533189, -0.0195919, -0.0195919/)
		do ix=1, flavorNumber
			write(tmparg,"('F_ann21 test 1 empty rho ',2I1)") ix,ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, ix,ix), tmparrA(ix), 1d-7)
			write(tmparg,"('F_ann12 test 1 empty rho ',2I1)") ix,ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, ix,ix), tmparrA(ix), 1d-7)
		end do

		do a=1, 2
			do b=1, 2
				do ix=1, flavorNumber
					iy=ix
					write(tmparg,"('F_ann',2I1,' test 2 empty rho ',2I1)") a, b, ix, iy
					call assert_double(trim(tmparg)//" im", F_ab_ann_im(nB, nA, f2, f3, a, b, ix, iy), 0.0d0, 1d-7)
					do iy=ix+1, flavorNumber
						write(tmparg,"('F_ann',2I1,' test 2 empty rho ',2I1)") a, b, ix, iy
						call assert_double(trim(tmparg)//" re", F_ab_ann_re(nB, nA, f2, f3, a, b, ix, iy), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_ann_im(nB, nA, f2, f3, a, b, ix, iy), 0.0d0, 1d-7)
						write(tmparg,"('F_ann',2I1,' test 2 empty rho ',2I1)") a, b, iy, ix
						call assert_double(trim(tmparg)//" re", F_ab_ann_re(nB, nA, f2, f3, a, b, iy, ix), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_ann_im(nB, nA, f2, f3, a, b, iy, ix), 0.0d0, 1d-7)
					end do
				end do
			end do
		end do
!		LL
		tmparrA = (/0.046175, 0.00623441, 0.00623441/)
		do ix=1, flavorNumber
			write(tmparg,"('F_ann11 test 2 empty rho ',2I1)") ix,ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 1,1, ix,ix), tmparrA(ix), 1d-7)
		end do
!		RR
		tmparrA = (/0.00461893, 0.00461893, 0.00461893/)
		do ix=1, flavorNumber
			write(tmparg,"('F_ann22 test 2 empty rho ',2I1)") ix,ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 2,2, ix,ix), tmparrA(ix), 1d-7)
		end do
!		LR and RL
		tmparrA = (/0.0146041, -0.00536622, -0.00536622/)
		do ix=1, flavorNumber
			write(tmparg,"('F_ann21 test 2 empty rho ',2I1)") ix,ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 2,1, ix,ix), tmparrA(ix), 1d-7)
			write(tmparg,"('F_ann12 test 2 empty rho ',2I1)") ix,ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 1,2, ix,ix), tmparrA(ix), 1d-7)
		end do

		write(*,*)
		write(*,"(a)") "F_sc functions, empty rho (144 tests)"
		do a=1, 2
			do b=1, 2
				do ix=1, flavorNumber
					iy=ix
					write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, ix, iy
					call assert_double(trim(tmparg)//" re", F_ab_sc_re(nA, f1, nB, f2, a, b, ix, iy), 0.0d0, 1d-7)
					call assert_double(trim(tmparg)//" im", F_ab_sc_im(nA, f1, nB, f2, a, b, ix, iy), 0.0d0, 1d-7)
					do iy=ix+1, flavorNumber
						write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, ix, iy
						call assert_double(trim(tmparg)//" re", F_ab_sc_re(nA, f1, nB, f2, a, b, ix, iy), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_sc_im(nA, f1, nB, f2, a, b, ix, iy), 0.0d0, 1d-7)
						write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, iy, ix
						call assert_double(trim(tmparg)//" re", F_ab_sc_re(nA, f1, nB, f2, a, b, iy, ix), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_sc_im(nA, f1, nB, f2, a, b, iy, ix), 0.0d0, 1d-7)
					end do
				end do
			end do
		end do
		do a=1, 2
			do b=1, 2
				do ix=1, flavorNumber
					iy=ix
					write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, ix, iy
					call assert_double(trim(tmparg)//" re", F_ab_sc_re(nB, f2, nA, f3, a, b, ix, iy), 0.0d0, 1d-7)
					call assert_double(trim(tmparg)//" im", F_ab_sc_im(nB, f2, nA, f3, a, b, ix, iy), 0.0d0, 1d-7)
					do iy=ix+1, flavorNumber
						write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, ix, iy
						call assert_double(trim(tmparg)//" re", F_ab_sc_re(nB, f2, nA, f3, a, b, ix, iy), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_sc_im(nB, f2, nA, f3, a, b, ix, iy), 0.0d0, 1d-7)
						write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, iy, ix
						call assert_double(trim(tmparg)//" re", F_ab_sc_re(nB, f2, nA, f3, a, b, iy, ix), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_sc_im(nB, f2, nA, f3, a, b, iy, ix), 0.0d0, 1d-7)
					end do
				end do
			end do
		end do

		write(*,*)
		write(*,"(a)") "F_ann functions, full rho (144 tests)"
		!rhoA = {{0.9, 0.011 + 0.0001 I, -0.03 - 0.004 I}, {0.011 + 0.0001 I, 0.86, 0.001 I}, {-0.03 - 0.004 I, 0.001 I, 0.96}};
		nA%re(1,:) = (/0.9d0, 0.011d0, -0.03d0/)
		nA%re(2,:) = (/0.011d0, 0.86d0, 0.d0/)
		nA%re(3,:) = (/-0.03d0, 0.d0, 0.96d0/)
		nA%im(1,:) = (/0.d0, 0.0001d0, -0.004d0/)
		nA%im(2,:) = (/-0.0001d0, 0.d0, 0.001d0/)
		nA%im(3,:) = (/0.004d0, -0.001d0, 0.d0/)
		!rhoB = {{1.23, 0.1 - 0.08 I, 0.008 + 0.007 I}, {0.1 - 0.08 I, 1.7, -0.02}, {0.008 + 0.007 I, -0.02, 0.77}};
		nB%re(1,:) = (/1.23d0, 0.1d0, 0.008d0/)
		nB%re(2,:) = (/0.1d0, 1.7d0, -0.02d0/)
		nB%re(3,:) = (/0.008d0, -0.02d0, 0.77d0/)
		nB%im(1,:) = (/0.d0, -0.08d0, 0.007d0/)
		nB%im(2,:) = (/0.08d0, 0.d0, 0.d0/)
		nB%im(3,:) = (/-0.007d0, 0.d0, 0.d0/)
		!RR
		tmpmatA(1,:) = (/-0.0419512, -0.00402987, 0.000838691/)
		tmpmatA(2,:) = (/-0.00402987, -0.0565448, 0.000740187/)
		tmpmatA(3,:) = (/0.000838691, 0.000740187, -0.0275819/)
		tmpmatB(1,:) = (/0., 0.00279858, -0.000103477/)
		tmpmatB(2,:) = (/-0.00279858, 0., -0.0000142409/)
		tmpmatB(3,:) = (/0.000103477, 0.0000142409, 0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann22 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LL
		tmpmatA(1,:) = (/-0.4191489, 0.0101562, 0.00807072/)
		tmpmatA(2,:) = (/0.0101562, -0.0762081, 0.000864564/)
		tmpmatA(3,:) = (/0.00807072, 0.000864564, -0.0372565/)
		tmpmatB(1,:) = (/0., -0.0103262, 0.00186556/)
		tmpmatB(2,:) = (/0.0103262, 0., -0.000160603/)
		tmpmatB(3,:) = (/-0.00186556, 0.000160603, 0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann11 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LR
		tmpmatA(1,:) = (/-0.1326407, -0.00381177, 0.00112954/)
		tmpmatA(2,:) = (/-0.00381177, 0.0656931, -0.00085994/)
		tmpmatA(3,:) = (/0.00112954, -0.00085994, 0.0320443/)
		tmpmatB(1,:) = (/0., +0.00273011, -0.0000740348/)
		tmpmatB(2,:) = (/-0.00273011, 0., +0.0000165449/)
		tmpmatB(3,:) = (/+0.0000740348, -0.0000165449, 0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann12 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!RL
		tmpmatA(1,:) = (/-0.132567, -0.00399017, 0.00115456/)
		tmpmatA(2,:) = (/-0.00399017, 0.0655955, -0.000744167/)
		tmpmatA(3,:) = (/0.00115456, -0.000744167, 0.0320682/)
		tmpmatB(1,:) = (/0., +0.00287272, -0.0000508025/)
		tmpmatB(2,:) = (/-0.00287272, 0., +0.000138238/)
		tmpmatB(3,:) = (/+0.0000508025, -0.000138238, 0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann21 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!second series
		!RR
		tmpmatA(1,:) = (/-0.0689959, -0.00652399, 0.00140441/)
		tmpmatA(2,:) = (/-0.00652399, -0.0914345, 0.00121473/)
		tmpmatA(3,:) = (/0.00140441, 0.00121473, -0.0459115/)
		tmpmatB(1,:) = (/0., +0.00441142, -0.000161176/)
		tmpmatB(2,:) = (/-0.00441142, +0., +2.68659e-6/)
		tmpmatB(3,:) = (/+0.000161176, -2.68659e-6, +0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann22 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 2,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nB, nA, f2, f3, 2,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LL
		tmpmatA(1,:) = (/-0.689097, -0.0282022, -0.00945936/)
		tmpmatA(2,:) = (/-0.0282022, -0.1230982, 0.00126491/)
		tmpmatA(3,:) = (/-0.00945936, 0.00126491, -0.0620467/)
		tmpmatB(1,:) = (/+0., +0.0255293, -0.00317523/)
		tmpmatB(2,:) = (/-0.0255293, +0., -0.000390201/)
		tmpmatB(3,:) = (/+0.00317523, +0.000390201, 0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann11 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 1,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nB, nA, f2, f3, 1,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LR
		tmpmatA(1,:) = (/-0.2181504, -0.0071315, 0.000594222/)
		tmpmatA(2,:) = (/-0.0071315, 0.1062276, -0.00141126/)
		tmpmatA(3,:) = (/0.000594222, -0.00141126, 0.0533395/)
		tmpmatB(1,:) = (/0., +0.00460214, -0.00024319/)
		tmpmatB(2,:) = (/-0.00460214, 0., -3.12125e-6/)
		tmpmatB(3,:) = (/+0.00024319, +3.12125e-6, +0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann12 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 1,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nB, nA, f2, f3, 1,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!RL
		tmpmatA(1,:) = (/-0.2179453, -0.00641341, 0.00228429/)
		tmpmatA(2,:) = (/-0.00641341, 0.1059558, -0.00108876/)
		tmpmatA(3,:) = (/0.00228429, -0.00108876, 0.0534062/)
		tmpmatB(1,:) = (/0., +0.00461794, -0.000014447/)
		tmpmatB(2,:) = (/-0.00461794, 0., +0.000335862/)
		tmpmatB(3,:) = (/+0.000014447, -0.000335862, 0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann21 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 2,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nB, nA, f2, f3, 2,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do

		write(*,*)
		write(*,"(a)") "F_sc functions, full rho (144 tests)"
		!RR
		tmpmatA(1,:) = (/0.0093345,0.00309005,0.000821132/)
		tmpmatA(2,:) = (/0.00309005,0.0248957,-0.000671452/)
		tmpmatA(3,:) = (/0.000821132,-0.000671452,-0.00692823/)
		tmpmatB(1,:) = (/0.,-0.00257023,+0.000305504/)
		tmpmatB(2,:) = (/+0.00257023,0.,-0.0000359033/)
		tmpmatB(3,:) = (/-0.000305504,+0.0000359033,0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc22 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 2,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 2,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LL
		tmpmatA(1,:) = (/0.0931435,-0.0126826,0.00174399/)
		tmpmatA(2,:) = (/-0.0126826,0.0335189,-0.000806577/)
		tmpmatA(3,:) = (/0.00174399,-0.000806577,-0.00933077/)
		tmpmatB(1,:) = (/+0.,+0.00943206,-0.000486463/)
		tmpmatB(2,:) = (/-0.00943206,+0.,+0.0000563555/)
		tmpmatB(3,:) = (/+0.000486463,-0.0000563555,0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc11 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 1,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 1,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LR
		tmpmatA(1,:) = (/0.0295137,0.00292836,0.000605503/)
		tmpmatA(2,:) = (/0.00292836,-0.0289235,0.000780085/)
		tmpmatA(3,:) = (/0.000605503,0.000780085,0.00804914/)
		tmpmatB(1,:) = (/0.,-0.00251947,+0.000283677/)
		tmpmatB(2,:) = (/+0.00251947,0.,+0.0000417121/)
		tmpmatB(3,:) = (/-0.000283677,-0.0000417121,+0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc12 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 1,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 1,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!RL
		tmpmatA(1,:) = (/0.0294591,0.00306062,0.000586954/)
		tmpmatA(2,:) = (/0.00306062,-0.0288512,0.000694255/)
		tmpmatA(3,:) = (/0.000586954,0.000694255,0.00803138/)
		tmpmatB(1,:) = (/+0.,-0.0026252,+0.000266453/)
		tmpmatB(2,:) = (/+0.0026252,0.,-0.0000485076/)
		tmpmatB(3,:) = (/-0.000266453,+0.0000485076,+0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc21 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 2,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 2,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!second series
		!RR
		tmpmatA(1,:) = (/-0.00570605,-0.000732471,-0.00101446/)
		tmpmatA(2,:) = (/-0.000732471,-0.0129257,0.000196912/)
		tmpmatA(3,:) = (/-0.00101446,0.000196912,0.00109681/)
		tmpmatB(1,:) = (/0.,+0.000956151,-0.000198153/)
		tmpmatB(2,:) = (/-0.000956151,+0.,+6.56285e-6/)
		tmpmatB(3,:) = (/+0.000198153,-6.56285e-6,+0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc22 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nB, f2, nA, f3, 2,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nB, f2, nA, f3, 2,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LL
		tmpmatA(1,:) = (/-0.057286,-0.00823654,0.00289993/)
		tmpmatA(2,:) = (/-0.00823654,-0.0175649,0.00040623/)
		tmpmatA(3,:) = (/0.00289993,0.00040623,0.00150947/)
		tmpmatB(1,:) = (/+0.,+0.00522933,+6.31568e-6/)
		tmpmatB(2,:) = (/-0.00522933,+0.,+0.000156488/)
		tmpmatB(3,:) = (/-6.31568e-6,-0.000156488,0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc11 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nB, f2, nA, f3, 1,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nB, f2, nA, f3, 1,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LR
		tmpmatA(1,:) = (/-0.0180413,-0.000504739,-0.000710753/)
		tmpmatA(2,:) = (/-0.000504739,0.0150169,-0.00022877/)
		tmpmatA(3,:) = (/-0.000710753,-0.00022877,-0.00127426/)
		tmpmatB(1,:) = (/0.,+0.000884656,-0.00016741/)
		tmpmatB(2,:) = (/-0.000884656,0.,-7.62464e-6/)
		tmpmatB(3,:) = (/+0.00016741,+7.62464e-6,+0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc12 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nB, f2, nA, f3, 1,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nB, f2, nA, f3, 1,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!RL
		tmpmatA(1,:) = (/-0.0181182,-0.000773921,-0.00134429/)
		tmpmatA(2,:) = (/-0.000773921,0.0151188,-0.000349659/)
		tmpmatA(3,:) = (/-0.00134429,-0.000349659,-0.00129926/)
		tmpmatB(1,:) = (/+2.7336e-21,+0.000878732,-0.000253156/)
		tmpmatB(2,:) = (/-0.000878732,-1.09344e-20,-0.000134696/)
		tmpmatB(3,:) = (/+0.000253156,+0.000134696,+0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc21 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nB, f2, nA, f3, 2,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nB, f2, nA, f3, 2,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		call deallocateCmplxMat(nA)
		call deallocateCmplxMat(nB)
	end subroutine do_f_ann_sc_re_tests_full

	subroutine do_tests_coll_int
		real(dl) :: x,z, y1
		type(coll_args) :: collArgs
		integer :: i, iy!,j, k
		real(dl) :: errr1,errr2, res1,res2,res3,res4, cf, ref
		INTEGER :: IFAIL, ITRANS, N, NPTS, NRAND
		real(dl) ::  VK(2)
		real(dl), dimension(:), allocatable :: ndmv_re
		real(dl), dimension(3) :: tmparrS, tmparrA
		character(len=300) :: tmparg

		x=0.05d0
		y1=1.2d0
		z=1.06d0
		npts=1
		nrand=1
		n=2

		call allocateCmplxMat(collArgs%n)
		allocate(ndmv_re(Ny))

		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVec(iy)%re(i, i) = 1.d0*i
				ndmv_re(iy) = nuDensMatVec(iy)%re(i, i)
			end do
#ifdef LOGY
			call interpNuDens%re(i, i)%replace(Ny, logy_arr, ndmv_re)
#else
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
#endif
		end do

		collArgs%x = x
		collArgs%z = z
		collArgs%y1 = y1
		collArgs%dme2 = 0.d0!dme2_electron(x, 0.d0, z)
		collArgs%n%re = 0.d0
		collArgs%n%im = 0.d0
		collArgs%n%y = y1
		collArgs%n%logy = log10(y1)
		do i=1, flavorNumber
			collArgs%n%re(i,i) = 1.d0 / fermiDirac(y1/z)
		end do

		write(*,*) ""
		write(*,"(a)") "Collision integrals (X tests)"

		!first series of sc and ann tests
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		vk=(/5.2d0,2.35132498560248d0/)
		ref = -4.8611d0
		res2=coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
		call assert_double_rel("test coll sc 3 sing 1", res2, ref, 2d-4)
		res1=coll_nue_4_sc_int_re(n, vk, collArgs)
		call assert_double_rel("test coll sc 4 sing 1", res1, ref, 2d-4)
		collArgs%ix1 = 2
		collArgs%ix2 = 2
		ref = -1.0161576d0
		res2 = coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
		call assert_double_rel("test coll sc 3 sing 2", res2, ref, 2d-4)
		res1 = coll_nue_4_sc_int_re(n, vk, collArgs)
		call assert_double_rel("test coll sc 4 sing 2", res1, ref, 2d-4)
		collArgs%ix1 = 3
		collArgs%ix2 = 3
		ref = -0.993358d0
		res2 = coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
		call assert_double_rel("test coll sc 3 sing 3", res2, ref, 2d-4)
		res1=coll_nue_4_sc_int_re(n, vk, collArgs)
		call assert_double_rel("test coll sc 4 sing 3", res1, ref, 2d-4)

		tmparrS = (/-160.602d0, -29.083d0, -23.845d0/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			ifail=0
			itrans=0
			call D01GCF(n,coll_nue_4_sc_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
			write(tmparg,"('test coll sc 4 - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res1, tmparrS(i), 3d-3)
			res2 = integrate_coll_int_3(coll_nue_3_sc_int_re, collArgs)
			write(tmparg,"('test coll sc 3 - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrS(i), 3d-3)
!			write(*,"(I2,*(E17.9))") i, res1, res2, tmparrS(i)
		end do

		collArgs%ix1 = 1
		collArgs%ix2 = 1
		vk=(/2.24741107392278d0, 3.0d0/)
		ref = -2.75887d0
		res2=coll_nue_3_ann_int_re(21, 3.d0, collArgs)
		call assert_double_rel("test coll ann 3 sing 1", res2, ref, 1d-4)
		res1=coll_nue_4_ann_int_re(n, vk, collArgs)
		call assert_double_rel("test coll ann 4 sing 1", res1, ref, 1d-4)
		collArgs%ix1 = 2
		collArgs%ix2 = 2
		ref = -0.981776d0
		res2 = coll_nue_3_ann_int_re(21, 3.d0, collArgs)
		call assert_double_rel("test coll ann 3 sing 2", res2, ref, 1d-4)
		res1=coll_nue_4_ann_int_re(n, vk, collArgs)
		call assert_double_rel("test coll ann 4 sing 2", res1, ref, 1d-4)
		collArgs%ix1 = 3
		collArgs%ix2 = 3
		ref = -1.47266d0
		res2 = coll_nue_3_ann_int_re(21, 3.d0, collArgs)
		call assert_double_rel("test coll ann 3 sing 3", res2, ref, 1d-4)
		res1=coll_nue_4_ann_int_re(n, vk, collArgs)
		call assert_double_rel("test coll ann 4 sing 3", res1, ref, 1d-4)

		tmparrA = (/-40.745d0, -17.299d0, -25.948d0/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			ifail=0
			itrans=0
			call D01GCF(n,coll_nue_4_ann_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
			write(tmparg,"('test coll ann 4 - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res1, tmparrA(i), 2d-1)
			res2 = integrate_coll_int_3(coll_nue_3_ann_int_re, collArgs)
			write(tmparg,"('test coll ann 3 - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrA(i), 2d-2)
!			write(*,"(*(E17.9))") res1, res2, tmparrA(i)
		end do

		!second series of sc and ann tests
		y1 = 13.2d0
		z = 1.3d0
		do i=1, flavorNumber
			collArgs%n%re(i,i) = 1.d0 / fermiDirac(y1/z)
			do iy=1, Ny
				nuDensMatVec(iy)%re(i, i) = 1.d0
				ndmv_re(iy) = nuDensMatVec(iy)%re(i, i)
			end do
#ifdef LOGY
			call interpNuDens%re(i, i)%replace(Ny, logy_arr, ndmv_re)
#else
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
#endif
		end do
		collArgs%x = 0.5d0
		collArgs%z = z
		collArgs%y1 = y1
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		vk=(/5.2d0,14.3669d0/)
		ref = -3642.11d0
		res2=coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
		call assert_double_rel("test coll sc 3 sing 4", res2, ref, 2d-2)
		res1=coll_nue_4_sc_int_re(n, vk, collArgs)
		call assert_double_rel("test coll sc 4 sing 4", res1, ref, 2d-2)
		collArgs%ix1 = 2
		collArgs%ix2 = 2
		ref = -780.84d0
		res2 = coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
		call assert_double_rel("test coll sc 3 sing 5", res2, ref, 2d-2)
		res1 = coll_nue_4_sc_int_re(n, vk, collArgs)
		call assert_double_rel("test coll sc 4 sing 5", res1, ref, 2d-2)
		collArgs%ix1 = 3
		collArgs%ix2 = 3
		ref = -780.84d0
		res2 = coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
		call assert_double_rel("test coll sc 3 sing 6", res2, ref, 2d-2)
		res1=coll_nue_4_sc_int_re(n, vk, collArgs)
		call assert_double_rel("test coll sc 4 sing 6", res1, ref, 2d-2)

		tmparrS = (/-607434d0, -130051d0, -130051d0/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			ifail=0
			itrans=0
			call D01GCF(n,coll_nue_4_sc_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
			write(tmparg,"('test coll sc 4 b - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res1, tmparrS(i), 0.025d0)
			res2 = integrate_coll_int_3(coll_nue_3_sc_int_re, collArgs)
			write(tmparg,"('test coll sc 3 b - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrS(i), 0.025d0)
!			write(*,"(I2,*(E17.9))") i, res1, res2, tmparrS(i)
		end do

		collArgs%ix1 = 1
		collArgs%ix2 = 1
		vk=(/14.1982d0, 3.0d0/)
		ref = -1880.98d0
		res2=coll_nue_3_ann_int_re(21, 3.d0, collArgs)
		call assert_double_rel("test coll ann 3 sing 4", res2, ref, 3d-3)
		res1=coll_nue_4_ann_int_re(n, vk, collArgs)
		call assert_double_rel("test coll ann 4 sing 4", res1, ref, 3d-3)
		collArgs%ix1 = 2
		collArgs%ix2 = 2
		ref = -263.819d0
		res2 = coll_nue_3_ann_int_re(21, 3.d0, collArgs)
		call assert_double_rel("test coll ann 3 sing 5", res2, ref, 3d-3)
		res1=coll_nue_4_ann_int_re(n, vk, collArgs)
		call assert_double_rel("test coll ann 4 sing 5", res1, ref, 3d-3)
		collArgs%ix1 = 3
		collArgs%ix2 = 3
		ref = -263.819d0
		res2 = coll_nue_3_ann_int_re(21, 3.d0, collArgs)
		call assert_double_rel("test coll ann 3 sing 6", res2, ref, 3d-3)
		res1=coll_nue_4_ann_int_re(n, vk, collArgs)
		call assert_double_rel("test coll ann 4 sing 6", res1, ref, 3d-3)

		tmparrA = (/-151024d0, -31971.3d0, -31971.3d0/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			ifail=0
			itrans=0
			call D01GCF(n,coll_nue_4_ann_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
			write(tmparg,"('test coll ann 4 - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res1, tmparrA(i), 3d-3)
			res2 = integrate_coll_int_3(coll_nue_3_ann_int_re, collArgs)
			write(tmparg,"('test coll ann 3 - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrA(i), 1d-3)
!			write(*,"(*(E17.9))") res1, res2, tmparrA(i)
		end do
		call criticalError("stop")
		deallocate(ndmv_re)
	end subroutine do_tests_coll_int

	subroutine test_speed_coll_int
		real(dl) :: x,z, y1
		type(coll_args) :: collArgs
		integer :: i, ix, Npt!,j, k
		real(dl) :: errr1,errr2, res1,res2,res3,res4, cf
		INTEGER :: IFAIL, ITRANS, N, NPTS, NRAND
		real(dl) ::  VK(2)
		real(dl), dimension(:), allocatable :: ndmv_re
		real(8) :: timer1

		x=0.05d0
		y1=1.2d0
		z=1.06d0
		npts=1
		nrand=1
		n=2
		Npt=1000!number of calculations for each comparison

		call allocateCmplxMat(collArgs%n)
		allocate(ndmv_re(Ny))

		do i=1, flavorNumber
			do ix=1, Ny
				nuDensMatVec(ix)%re(i, i) = 1.d0*i * fermiDirac(y_arr(ix)/z)
				ndmv_re(ix) = nuDensMatVec(ix)%re(i, i)
			end do
#ifdef LOGY
			call interpNuDens%re(i, i)%replace(Ny, logy_arr, ndmv_re)
#else
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
#endif
		end do

		collArgs%x = x
		collArgs%z = z
		collArgs%y1 = y1
		collArgs%dme2 = 0.d0!dme2_electron(x, 0.d0, z)
		collArgs%n%re = 0.d0
		collArgs%n%im = 0.d0
		do i=1, flavorNumber
			collArgs%n%re(i,i) = 1.d0
		end do
		collArgs%ix1 = 1
		collArgs%ix2 = 1

		if (timing_tests) then
			write (*,*) "[interactions] timing 2D integrals..."
			call tic(timer1)
			do ix=1, Npt
				call random_number(x)
				call random_number(z)
				collArgs%x = (x_fin-x_in)*x + x_in
				collArgs%z = 0.4d0*z + z_in
				ifail=0
				itrans=0
				res2 = integrate_coll_int_3(coll_nue_3_sc_int_re, collArgs)
			end do
			call toc(timer1, "<sc re semilin>")
			call tic(timer1)
			do ix=1, Npt
				call random_number(x)
				call random_number(z)
				collArgs%x = (x_fin-x_in)*x + x_in
				collArgs%z = 0.4d0*z + z_in
				ifail=0
				itrans=0
				call D01GCF(n,coll_nue_4_sc_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
			end do
			call toc(timer1, "<sc re D01GCF>")

			call tic(timer1)
			do ix=1, Npt
				call random_number(x)
				call random_number(z)
				collArgs%x = (x_fin-x_in)*x + x_in
				collArgs%z = 0.4d0*z + z_in
				ifail=0
				itrans=0
				res2 = integrate_coll_int_3(coll_nue_3_ann_int_re, collArgs)
			end do
			call toc(timer1, "<ann re semilin>")
			call tic(timer1)
			do ix=1, Npt
				call random_number(x)
				call random_number(z)
				collArgs%x = (x_fin-x_in)*x + x_in
				collArgs%z = 0.4d0*z + z_in
				ifail=0
				itrans=0
				call D01GCF(n,coll_nue_4_ann_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
			end do
			call toc(timer1, "<ann re D01GCF>")
		end if
		deallocate(ndmv_re)
		call addToLog("[interactions] ...done!")
	end subroutine test_speed_coll_int

	subroutine do_timing_tests
		timing_tests = .true.
		call test_dzodx_speed
		call test_nuDens_speed
		call init_interp_ElDensity
		call init_interp_dme2_e
		call init_interp_FD
		call init_interp_d123
		call test_speed_coll_int
	end subroutine do_timing_tests

end program tests
