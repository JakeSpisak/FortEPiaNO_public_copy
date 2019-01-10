subroutine assert_double(testname, num1, num2, tol)
	use precision
	implicit none
	character(len=*) :: testname
	real(dl) :: num1, num2, tol

	if (abs(num1 - num2) .lt. tol) then
		write(*, fmt="(a)", advance="no") "."
	else
		print *, testname, " failed"
		write (*,"(*(E15.7))") num1, num2, num1 - num2, tol
		call exit()
	end if
end subroutine assert_double

program tests
	use precision
	use ndConfig
	use ndErrors
	use ndEquations
	implicit none

	integer :: ix
	real(dl) :: fdA, fdB
	real(dl), dimension(2) :: temp_v2
	type(cmplxMatNN) :: nA, nB


	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "Initializations"
	checkpoint = .true.
	maxiter = 100
	toler   = 1.d-5
	dlsoda_atol = 1.d-6
	dlsoda_rtol = 1.d-6
	Nx = 100
	Ny = 100
	allocate(x_arr(Nx), y_arr(Ny), logy_arr(Ny))
	x_in    = 0.01d0
	x_fin   = 40.d0
	logx_in  = log10(x_in)
	logx_fin = log10(x_fin)
	x_arr = logspace(logx_in, logx_fin, Nx)
	y_min   = 0.01d0
	y_max   = 20.0d0
#ifdef LOGY
	logy_min = log10(y_min)
	logy_max = log10(y_max)
	y_arr = logspace(logy_min, logy_max, Ny)
	logy_arr = log10(y_arr)
#else
	y_arr = linspace(y_min, y_max, Ny)
	logy_arr = log10(y_arr)
#endif
	z_in = 1.00003d0
	m_lightest   = 0.0d0
	massOrdering = .true.
	!read cosmo parameters
	hubbleParam = i_HubbleParam
	photonTemperatureToday = i_photonTempToday
	!settings for collisional
	collision_offdiag = 1
	dme2_temperature_corr = .true.
	only_1a_1s = .false.
	flavorNumber = 3
	if (collision_offdiag.ne.3) then
		flavNumSqu = flavorNumber**2
	else
		flavNumSqu = flavorNumber
	end if
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


	write(*,*) ""
	write(*,"(a)") "Starting tests"


	write(*,*) ""
	write(*,"(a)") "D_i functions (12 tests)"
	call assert_double("D1 test 1", D1_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.4d0, 1d-7)
	call assert_double("D1 test 2", D1_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.4d0, 1d-7)
	call assert_double("D1 test 3", D1_full(0.01d0,5.d0,2.6d0,2.41d0), 0.04d0, 1d-7)
	call assert_double("D1 test 4", D1_full(10.03d0,5.4d0,8.8d0,6.63d0), 21.6d0, 1d-4)
	call assert_double("D2 test 1", D2_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), -0.00133333d0, 1d-7)
	call assert_double("D2 test 2", D2_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), -0.0213333d0, 1d-7)
	call assert_double("D2 test 3", D2_full(0.01d0,5.d0,2.6d0,2.41d0), -1.333333333d-6, 1d-12)
	call assert_double("D2 test 4", D2_full(10.03d0,5.4d0,8.8d0,6.63d0), -209.952d0, 1d-4)
	call assert_double("D3 test 1", D3_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.000192d0, 1d-7)
	call assert_double("D3 test 2", D3_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.000192d0, 1d-7)
	call assert_double("D3 test 3", D3_full(0.01d0,5.d0,2.6d0,2.41d0), 2.50454d-5, 1d-10)
	call assert_double("D3 test 4", D3_full(10.03d0,5.4d0,8.8d0,6.63d0), 22692.2d0, 3d-1)


	write(*,*) ""
	write(*,"(a)") "Pi(yi,yj) functions (12 tests)"
	call assert_double("Pi_1_12 test 1", PI1_12_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.00933333d0, 1d-7)
	call assert_double("Pi_1_12 test 2", PI1_12_full(0.4d0, 0.3d0, 0.2d0, 0.1d0), 0.0893333d0, 1d-7)
	call assert_double("Pi_1_13 test 1", PI1_13_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.0106667d0, 1d-7)
	call assert_double("Pi_1_13 test 2", PI1_13_full(0.4d0, 0.3d0, 0.2d0, 0.1d0), 0.0106667d0, 1d-7)
	temp_v2 = PI2_nn_f(0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.3d0, 0.4d0)
	call assert_double("Pi_2_14 test 1", temp_v2(1), 0.00267733d0, 1d-7)
	call assert_double("Pi_2_13 test 1", temp_v2(2), 0.000810667d0, 1d-7)
	temp_v2 = PI2_nn_f(0.4d0, 0.3d0, 0.2d0, 0.1d0, 0.2d0, 0.1d0)
	call assert_double("Pi_2_14 test 2", temp_v2(1), 0.00267733d0, 1d-7)
	call assert_double("Pi_2_13 test 2", temp_v2(2), 0.000810667d0, 1d-7)
	temp_v2 = PI2_ne_f(0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.2d0, 0.4d0)
	call assert_double("Pi_2_14 test 1", temp_v2(1), 0.00267733d0, 1d-7)
	call assert_double("Pi_2_12 test 1", temp_v2(2), 0.00427733d0, 1d-7)
	temp_v2 = PI2_ne_f(0.4d0, 0.3d0, 0.2d0, 0.1d0, 0.3d0, 0.1d0)
	call assert_double("Pi_2_14 test 2", temp_v2(1), 0.00267733d0, 1d-7)
	call assert_double("Pi_2_12 test 2", temp_v2(2), 0.00427733d0, 1d-7)


	write(*,*) ""
	write(*,"(a)") "F_ann_re functions (72 tests)"
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
	call assert_double("FRR_ann test 1 - 1,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,2, 1, 1), -0.00259399d0, 1d-7)
	call assert_double("FRR_ann test 1 - 1,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,2, 1, 2), 0.d0, 1d-7)
	call assert_double("FRR_ann test 1 - 1,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,2, 1, 3), 0.d0, 1d-7)
	call assert_double("FRR_ann test 1 - 2,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,2, 2, 1), 0.d0, 1d-7)
	call assert_double("FRR_ann test 1 - 2,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,2, 2, 2), -0.00259399d0, 1d-7)
	call assert_double("FRR_ann test 1 - 2,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,2, 2, 3), 0.d0, 1d-7)
	call assert_double("FRR_ann test 1 - 3,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,2, 3, 1), 0.d0, 1d-7)
	call assert_double("FRR_ann test 1 - 3,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,2, 3, 2), 0.d0, 1d-7)
	call assert_double("FRR_ann test 1 - 3,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,2, 3, 3), -0.00259399d0, 1d-7)
	!FLLeq[0.1, 0.2, 0.3, 0.4]={{-0.0259319, 0., 0.}, {0., -0.00350124, 0.}, {0., 0., -0.00350124}}
	call assert_double("FLL_ann test 1 - 1,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,1, 1, 1), -0.0259319d0, 1d-7)
	call assert_double("FLL_ann test 1 - 1,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,1, 1, 2), 0.d0, 1d-7)
	call assert_double("FLL_ann test 1 - 1,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,1, 1, 3), 0.d0, 1d-7)
	call assert_double("FLL_ann test 1 - 2,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,1, 2, 1), 0.d0, 1d-7)
	call assert_double("FLL_ann test 1 - 2,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,1, 2, 2), -0.00350124d0, 1d-7)
	call assert_double("FLL_ann test 1 - 2,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,1, 2, 3), 0.d0, 1d-7)
	call assert_double("FLL_ann test 1 - 3,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,1, 3, 1), 0.d0, 1d-7)
	call assert_double("FLL_ann test 1 - 3,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,1, 3, 2), 0.d0, 1d-7)
	call assert_double("FLL_ann test 1 - 3,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,1, 3, 3), -0.00350124d0, 1d-7)
	!FRLeq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
	call assert_double("FRL_ann test 1 - 1,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,1, 1, 1), -0.00820165d0, 1d-7)
	call assert_double("FRL_ann test 1 - 1,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,1, 1, 2), 0.d0, 1d-7)
	call assert_double("FRL_ann test 1 - 1,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,1, 1, 3), 0.d0, 1d-7)
	call assert_double("FRL_ann test 1 - 2,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,1, 2, 1), 0.d0, 1d-7)
	call assert_double("FRL_ann test 1 - 2,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,1, 2, 2), 0.00301367d0, 1d-7)
	call assert_double("FRL_ann test 1 - 2,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,1, 2, 3), 0.d0, 1d-7)
	call assert_double("FRL_ann test 1 - 3,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,1, 3, 1), 0.d0, 1d-7)
	call assert_double("FRL_ann test 1 - 3,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,1, 3, 2), 0.d0, 1d-7)
	call assert_double("FRL_ann test 1 - 3,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 2,1, 3, 3), 0.00301367d0, 1d-7)
	!FLReq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
	call assert_double("FLR_ann test 1 - 1,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,2, 1, 1), -0.00820165d0, 1d-7)
	call assert_double("FLR_ann test 1 - 1,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,2, 1, 2), 0.d0, 1d-7)
	call assert_double("FLR_ann test 1 - 1,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,2, 1, 3), 0.d0, 1d-7)
	call assert_double("FLR_ann test 1 - 2,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,2, 2, 1), 0.d0, 1d-7)
	call assert_double("FLR_ann test 1 - 2,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,2, 2, 2), 0.00301367d0, 1d-7)
	call assert_double("FLR_ann test 1 - 2,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,2, 2, 3), 0.d0, 1d-7)
	call assert_double("FLR_ann test 1 - 3,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,2, 3, 1), 0.d0, 1d-7)
	call assert_double("FLR_ann test 1 - 3,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,2, 3, 2), 0.d0, 1d-7)
	call assert_double("FLR_ann test 1 - 3,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.4d0, 1,2, 3, 3), 0.00301367d0, 1d-7)
	fdA = fermiDirac(0.4d0)
	do ix=1, flavorNumber
		nA%re(ix,ix) = fdA
	end do
	!FRReq[0.4, 0.2, 0.3, 0.1]={{0.00129054, 0., 0.}, {0., 0.00129054, 0.}, {0., 0., 0.00129054}}
	call assert_double("FRR_ann test 2 - 1,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,2, 1, 1), 0.00129054d0, 1d-7)
	call assert_double("FRR_ann test 2 - 1,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,2, 1, 2), 0.d0, 1d-7)
	call assert_double("FRR_ann test 2 - 1,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,2, 1, 3), 0.d0, 1d-7)
	call assert_double("FRR_ann test 2 - 2,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,2, 2, 1), 0.d0, 1d-7)
	call assert_double("FRR_ann test 2 - 2,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,2, 2, 2), 0.00129054d0, 1d-7)
	call assert_double("FRR_ann test 2 - 2,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,2, 2, 3), 0.d0, 1d-7)
	call assert_double("FRR_ann test 2 - 3,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,2, 3, 1), 0.d0, 1d-7)
	call assert_double("FRR_ann test 2 - 3,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,2, 3, 2), 0.d0, 1d-7)
	call assert_double("FRR_ann test 2 - 3,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,2, 3, 3), 0.00129054d0, 1d-7)
	!FLLeq[0.4, 0.2, 0.3, 0.1]={{0.0129014, 0., 0.}, {0., 0.00174191, 0.}, {0., 0., 0.00174191}}
	call assert_double("FLL_ann test 2 - 1,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,1, 1, 1), 0.0129014d0, 1d-7)
	call assert_double("FLL_ann test 2 - 1,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,1, 1, 2), 0.d0, 1d-7)
	call assert_double("FLL_ann test 2 - 1,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,1, 1, 3), 0.d0, 1d-7)
	call assert_double("FLL_ann test 2 - 2,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,1, 2, 1), 0.d0, 1d-7)
	call assert_double("FLL_ann test 2 - 2,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,1, 2, 2), 0.00174191d0, 1d-7)
	call assert_double("FLL_ann test 2 - 2,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,1, 2, 3), 0.d0, 1d-7)
	call assert_double("FLL_ann test 2 - 3,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,1, 3, 1), 0.d0, 1d-7)
	call assert_double("FLL_ann test 2 - 3,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,1, 3, 2), 0.d0, 1d-7)
	call assert_double("FLL_ann test 2 - 3,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,1, 3, 3), 0.00174191d0, 1d-7)
	!FRLeq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
	call assert_double("FRL_ann test 2 - 1,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,1, 1, 1), 0.00408041d0, 1d-7)
	call assert_double("FRL_ann test 2 - 1,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,1, 1, 2), 0.d0, 1d-7)
	call assert_double("FRL_ann test 2 - 1,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,1, 1, 3), 0.d0, 1d-7)
	call assert_double("FRL_ann test 2 - 2,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,1, 2, 1), 0.d0, 1d-7)
	call assert_double("FRL_ann test 2 - 2,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,1, 2, 2), -0.00149933d0, 1d-7)
	call assert_double("FRL_ann test 2 - 2,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,1, 2, 3), 0.d0, 1d-7)
	call assert_double("FRL_ann test 2 - 3,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,1, 3, 1), 0.d0, 1d-7)
	call assert_double("FRL_ann test 2 - 3,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,1, 3, 2), 0.d0, 1d-7)
	call assert_double("FRL_ann test 2 - 3,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 2,1, 3, 3), -0.00149933d0, 1d-7)
	!FLReq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
	call assert_double("FLR_ann test 2 - 1,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,2, 1, 1), 0.00408041d0, 1d-7)
	call assert_double("FLR_ann test 2 - 1,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,2, 1, 2), 0.d0, 1d-7)
	call assert_double("FLR_ann test 2 - 1,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,2, 1, 3), 0.d0, 1d-7)
	call assert_double("FLR_ann test 2 - 2,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,2, 2, 1), 0.d0, 1d-7)
	call assert_double("FLR_ann test 2 - 2,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,2, 2, 2), -0.00149933d0, 1d-7)
	call assert_double("FLR_ann test 2 - 2,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,2, 2, 3), 0.d0, 1d-7)
	call assert_double("FLR_ann test 2 - 3,1", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,2, 3, 1), 0.d0, 1d-7)
	call assert_double("FLR_ann test 2 - 3,2", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,2, 3, 2), 0.d0, 1d-7)
	call assert_double("FLR_ann test 2 - 3,3", F_ab_ann_re(0.d0, 1.d0, nA,nB, 0.3d0, 0.1d0, 1,2, 3, 3), -0.00149933d0, 1d-7)


	write(*,*) ""
	write(*,"(a)") "F_sc_re functions (72 tests)"
	fdA = fermiDirac(0.1d0)
	fdB = fermiDirac(0.3d0)
	do ix=1, flavorNumber
		nA%re(ix,ix) = fdA
		nB%re(ix,ix) = fdB
	end do
	!FRReq[0.1, 0.2, 0.3, 0.4]={{-0.00259399, 0., 0.}, {0., -0.00259399, 0.}, {0., 0., -0.00259399}}
	call assert_double("FRR_sc test 1 - 1,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,2, 1, 1), -0.00259399d0, 1d-7)
	call assert_double("FRR_sc test 1 - 1,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,2, 1, 2), 0.d0, 1d-7)
	call assert_double("FRR_sc test 1 - 1,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,2, 1, 3), 0.d0, 1d-7)
	call assert_double("FRR_sc test 1 - 2,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,2, 2, 1), 0.d0, 1d-7)
	call assert_double("FRR_sc test 1 - 2,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,2, 2, 2), -0.00259399d0, 1d-7)
	call assert_double("FRR_sc test 1 - 2,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,2, 2, 3), 0.d0, 1d-7)
	call assert_double("FRR_sc test 1 - 3,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,2, 3, 1), 0.d0, 1d-7)
	call assert_double("FRR_sc test 1 - 3,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,2, 3, 2), 0.d0, 1d-7)
	call assert_double("FRR_sc test 1 - 3,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,2, 3, 3), -0.00259399d0, 1d-7)
	!FLLeq[0.1, 0.2, 0.3, 0.4]={{-0.0259319, 0., 0.}, {0., -0.00350124, 0.}, {0., 0., -0.00350124}}
	call assert_double("FLL_sc test 1 - 1,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,1, 1, 1), -0.0259319d0, 1d-7)
	call assert_double("FLL_sc test 1 - 1,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,1, 1, 2), 0.d0, 1d-7)
	call assert_double("FLL_sc test 1 - 1,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,1, 1, 3), 0.d0, 1d-7)
	call assert_double("FLL_sc test 1 - 2,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,1, 2, 1), 0.d0, 1d-7)
	call assert_double("FLL_sc test 1 - 2,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,1, 2, 2), -0.00350124d0, 1d-7)
	call assert_double("FLL_sc test 1 - 2,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,1, 2, 3), 0.d0, 1d-7)
	call assert_double("FLL_sc test 1 - 3,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,1, 3, 1), 0.d0, 1d-7)
	call assert_double("FLL_sc test 1 - 3,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,1, 3, 2), 0.d0, 1d-7)
	call assert_double("FLL_sc test 1 - 3,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,1, 3, 3), -0.00350124d0, 1d-7)
	!FRLeq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
	call assert_double("FRL_sc test 1 - 1,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,1, 1, 1), -0.00820165d0, 1d-7)
	call assert_double("FRL_sc test 1 - 1,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,1, 1, 2), 0.d0, 1d-7)
	call assert_double("FRL_sc test 1 - 1,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,1, 1, 3), 0.d0, 1d-7)
	call assert_double("FRL_sc test 1 - 2,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,1, 2, 1), 0.d0, 1d-7)
	call assert_double("FRL_sc test 1 - 2,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,1, 2, 2), 0.00301367d0, 1d-7)
	call assert_double("FRL_sc test 1 - 2,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,1, 2, 3), 0.d0, 1d-7)
	call assert_double("FRL_sc test 1 - 3,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,1, 3, 1), 0.d0, 1d-7)
	call assert_double("FRL_sc test 1 - 3,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,1, 3, 2), 0.d0, 1d-7)
	call assert_double("FRL_sc test 1 - 3,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 2,1, 3, 3), 0.00301367d0, 1d-7)
	!FLReq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
	call assert_double("FLR_sc test 1 - 1,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,2, 1, 1), -0.00820165d0, 1d-7)
	call assert_double("FLR_sc test 1 - 1,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,2, 1, 2), 0.d0, 1d-7)
	call assert_double("FLR_sc test 1 - 1,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,2, 1, 3), 0.d0, 1d-7)
	call assert_double("FLR_sc test 1 - 2,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,2, 2, 1), 0.d0, 1d-7)
	call assert_double("FLR_sc test 1 - 2,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,2, 2, 2), 0.00301367d0, 1d-7)
	call assert_double("FLR_sc test 1 - 2,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,2, 2, 3), 0.d0, 1d-7)
	call assert_double("FLR_sc test 1 - 3,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,2, 3, 1), 0.d0, 1d-7)
	call assert_double("FLR_sc test 1 - 3,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,2, 3, 2), 0.d0, 1d-7)
	call assert_double("FLR_sc test 1 - 3,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.4d0, 1,2, 3, 3), 0.00301367d0, 1d-7)

	fdA = fermiDirac(0.4d0)
	do ix=1, flavorNumber
		nA%re(ix,ix) = fdA
	end do
	!FRReq[0.4, 0.2, 0.3, 0.1]={{0.00129054, 0., 0.}, {0., 0.00129054, 0.}, {0., 0., 0.00129054}}
	call assert_double("FRR_sc test 2 - 1,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,2, 1, 1), 0.00129054d0, 1d-7)
	call assert_double("FRR_sc test 2 - 1,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,2, 1, 2), 0.d0, 1d-7)
	call assert_double("FRR_sc test 2 - 1,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,2, 1, 3), 0.d0, 1d-7)
	call assert_double("FRR_sc test 2 - 2,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,2, 2, 1), 0.d0, 1d-7)
	call assert_double("FRR_sc test 2 - 2,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,2, 2, 2), 0.00129054d0, 1d-7)
	call assert_double("FRR_sc test 2 - 2,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,2, 2, 3), 0.d0, 1d-7)
	call assert_double("FRR_sc test 2 - 3,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,2, 3, 1), 0.d0, 1d-7)
	call assert_double("FRR_sc test 2 - 3,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,2, 3, 2), 0.d0, 1d-7)
	call assert_double("FRR_sc test 2 - 3,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,2, 3, 3), 0.00129054d0, 1d-7)
	!FLLeq[0.4, 0.2, 0.3, 0.1]={{0.0129014, 0., 0.}, {0., 0.00174191, 0.}, {0., 0., 0.00174191}}
	call assert_double("FLL_sc test 2 - 1,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,1, 1, 1), 0.0129014d0, 1d-7)
	call assert_double("FLL_sc test 2 - 1,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,1, 1, 2), 0.d0, 1d-7)
	call assert_double("FLL_sc test 2 - 1,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,1, 1, 3), 0.d0, 1d-7)
	call assert_double("FLL_sc test 2 - 2,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,1, 2, 1), 0.d0, 1d-7)
	call assert_double("FLL_sc test 2 - 2,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,1, 2, 2), 0.00174191d0, 1d-7)
	call assert_double("FLL_sc test 2 - 2,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,1, 2, 3), 0.d0, 1d-7)
	call assert_double("FLL_sc test 2 - 3,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,1, 3, 1), 0.d0, 1d-7)
	call assert_double("FLL_sc test 2 - 3,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,1, 3, 2), 0.d0, 1d-7)
	call assert_double("FLL_sc test 2 - 3,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,1, 3, 3), 0.00174191d0, 1d-7)
	!FRLeq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
	call assert_double("FRL_sc test 2 - 1,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,1, 1, 1), 0.00408041d0, 1d-7)
	call assert_double("FRL_sc test 2 - 1,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,1, 1, 2), 0.d0, 1d-7)
	call assert_double("FRL_sc test 2 - 1,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,1, 1, 3), 0.d0, 1d-7)
	call assert_double("FRL_sc test 2 - 2,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,1, 2, 1), 0.d0, 1d-7)
	call assert_double("FRL_sc test 2 - 2,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,1, 2, 2), -0.00149933d0, 1d-7)
	call assert_double("FRL_sc test 2 - 2,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,1, 2, 3), 0.d0, 1d-7)
	call assert_double("FRL_sc test 2 - 3,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,1, 3, 1), 0.d0, 1d-7)
	call assert_double("FRL_sc test 2 - 3,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,1, 3, 2), 0.d0, 1d-7)
	call assert_double("FRL_sc test 2 - 3,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 2,1, 3, 3), -0.00149933d0, 1d-7)
	!FLReq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
	call assert_double("FLR_sc test 2 - 1,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,2, 1, 1), 0.00408041d0, 1d-7)
	call assert_double("FLR_sc test 2 - 1,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,2, 1, 2), 0.d0, 1d-7)
	call assert_double("FLR_sc test 2 - 1,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,2, 1, 3), 0.d0, 1d-7)
	call assert_double("FLR_sc test 2 - 2,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,2, 2, 1), 0.d0, 1d-7)
	call assert_double("FLR_sc test 2 - 2,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,2, 2, 2), -0.00149933d0, 1d-7)
	call assert_double("FLR_sc test 2 - 2,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,2, 2, 3), 0.d0, 1d-7)
	call assert_double("FLR_sc test 2 - 3,1", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,2, 3, 1), 0.d0, 1d-7)
	call assert_double("FLR_sc test 2 - 3,2", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,2, 3, 2), 0.d0, 1d-7)
	call assert_double("FLR_sc test 2 - 3,3", F_ab_sc_re(0.d0, 1.d0, nA, 0.2d0, nB, 0.1d0, 1,2, 3, 3), -0.00149933d0, 1d-7)

	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "all tests were successfull!"
end program tests
