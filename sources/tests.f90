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

	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "Initializations"
	call do_tests_initialization

	write(*,*) ""
	write(*,"(a)") "Starting tests"
	call do_tests_Di
	call do_tests_Pi_ij
	call do_f_ann_sc_re_tests_eq
	call do_f_ann_sc_re_tests_full

	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "all tests were successfull!"

	write(*,*) ""
	write(*,"(a)") "now doing some timing tests..."
	call do_timing_tests()

	contains

	subroutine do_tests_initialization
		integer :: ix
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
		logy_min = log10(y_min)
		logy_max = log10(y_max)
#ifdef LOGY
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
	end subroutine do_tests_initialization

	subroutine do_tests_Di
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
	end subroutine do_tests_Di

	subroutine do_tests_Pi_ij
		real(dl), dimension(2) :: temp_v2
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
			call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.1, 0.2, 0.3, 0.4]={{-0.0259319, 0., 0.}, {0., -0.00350124, 0.}, {0., 0., -0.00350124}}
		tmparr = (/-0.0259319, -0.00350124, -0.00350124/)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
		tmparr = (/-0.00820165, 0.00301367, 0.00301367/)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
		tmparr = (/-0.00820165, 0.00301367, 0.00301367/)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, iy, ix), 0.d0, 1d-7)
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
			call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 2,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 2,2, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 2,2, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.4, 0.2, 0.3, 0.1]={{0.0129014, 0., 0.}, {0., 0.00174191, 0.}, {0., 0., 0.00174191}}
		tmparr = (/0.0129014, 0.00174191, 0.00174191/)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 1,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 1,1, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 1,1, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
		tmparr = (/0.00408041, -0.00149933, -0.00149933/)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 2,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 2,1, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 2,1, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
		tmparr = (/0.00408041, -0.00149933, -0.00149933/)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 1,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 1,2, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_ann_re(nA, nB, f1, f3, 1,2, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_ann_im(nA, nB, f1, f3, 1,2, iy, ix), 0.d0, 1d-7)
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
			call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 2,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 2,2, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 2,2, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.1, 0.2, 0.3, 0.4]={{-0.0259319, 0., 0.}, {0., -0.00350124, 0.}, {0., 0., -0.00350124}}
		tmparr = (/-0.0259319, -0.00350124, -0.00350124/)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 1,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 1,1, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 1,1, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
		tmparr = (/-0.00820165, 0.00301367, 0.00301367/)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 2,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 2,1, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 2,1, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.1, 0.2, 0.3, 0.4]={{-0.00820165, 0., 0.}, {0., 0.00301367, 0.}, {0., 0., 0.00301367}}
		tmparr = (/-0.00820165, 0.00301367, 0.00301367/)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 1,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 1,2, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f2, 1,2, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f2, 1,2, iy, ix), 0.d0, 1d-7)
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
			call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 2,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 2,2, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 2,2, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.4, 0.2, 0.3, 0.1]={{0.0129014, 0., 0.}, {0., 0.00174191, 0.}, {0., 0., 0.00174191}}
		tmparr = (/0.0129014, 0.00174191, 0.00174191/)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 1,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 1,1, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 1,1, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
		tmparr = (/0.00408041, -0.00149933, -0.00149933/)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 2,1, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 2,1, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 2,1, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.4, 0.2, 0.3, 0.1]={{0.00408041, 0., 0.}, {0., -0.00149933, 0.}, {0., 0., -0.00149933}}
		tmparr = (/0.00408041, -0.00149933, -0.00149933/)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 1,2, ix, ix), tmparr(ix), 1d-7)
			call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 1,2, ix, iy), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(tmparg//"re", F_ab_sc_re(nA, f1, nB, f3, 1,2, iy, ix), 0.d0, 1d-7)
				call assert_double(tmparg//"im", F_ab_sc_im(nA, f1, nB, f3, 1,2, iy, ix), 0.d0, 1d-7)
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
		nA%im(2,:) = (/0.0001d0, 0.d0, 0.001d0/)
		nA%im(3,:) = (/-0.004d0, 0.001d0, 0.d0/)
		!rhoB = {{1.23, 0.1 - 0.08 I, 0.008 + 0.007 I}, {0.1 - 0.08 I, 1.7, -0.02}, {0.008 + 0.007 I, -0.02, 0.77}};
		nB%re(1,:) = (/1.23d0, 0.1d0, 0.008d0/)
		nB%re(2,:) = (/0.1d0, 1.7d0, -0.02d0/)
		nB%re(3,:) = (/0.008d0, -0.02d0, 0.77d0/)
		nB%im(1,:) = (/0.d0, -0.08d0, 0.007d0/)
		nB%im(2,:) = (/-0.08d0, 0.d0, 0.d0/)
		nB%im(3,:) = (/0.007d0, 0.d0, 0.d0/)
		!RR
		tmpmatA(1,:) = (/-0.0419526, -0.00402972, 0.000838691/)
		tmpmatA(2,:) = (/-0.00402972, -0.0565451, 0.000746814/)
		tmpmatA(3,:) = (/0.000838691, 0.000746814, -0.027583/)
		tmpmatB(1,:) = (/0.0000229817, 0.00279841, - 0.000103477/)
		tmpmatB(2,:) = (/0.00279841, 0.0000183936, - 0.0000638581/)
		tmpmatB(3,:) = (/-0.000103477, -0.0000638581, 5.41474e-6/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann22 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LL
		tmpmatA(1,:) = (/-0.4191434, 0.0101556, 0.00807072/)
		tmpmatA(2,:) = (/0.0101556, -0.0762068, 0.000840218/)
		tmpmatA(3,:) = (/0.00807072, 0.000840218, -0.0372522/)
		tmpmatB(1,:) = (/- 0.0000844192, - 0.0103256, 0.00186556/)
		tmpmatB(2,:) = (/- 0.0103256, -0.0000654895, 0.0000216574/)
		tmpmatB(3,:) = (/0.00186556, 0.0000216574, - 0.0000178139/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann11 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LR
		tmpmatA(1,:) = (/-0.1326454, -0.00381131, 0.00112954/)
		tmpmatA(2,:) = (/-0.00381131,  0.0656934, -0.00086764/)
		tmpmatA(3,:) = (/0.00112954, -0.00086764, 0.0320456/)
		tmpmatB(1,:) = (/0.0000726631, 0.00272959, - 0.0000740348/)
		tmpmatB(2,:) = (/0.00272959, - 0.0000213695, 0.0000741896/)
		tmpmatB(3,:) = (/- 0.0000740348, 0.0000741896, - 6.29078e-6/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann12 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!RL
		tmpmatA(1,:) = (/-0.1325653, -0.00399034, 0.00115456/)
		tmpmatA(2,:) = (/-0.00399034, 0.0655944, -0.000723211/)
		tmpmatA(3,:) = (/0.00115456, -0.000723211, 0.0320646/)
		tmpmatB(1,:) = (/- 0.0000266998, 0.00287291, - 0.0000508025/)
		tmpmatB(2,:) = (/0.00287291, 0.0000563695, - 0.0000186415/)
		tmpmatB(3,:) = (/- 0.0000508025, - 0.0000186415, 0.0000153332/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann21 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!second series
		!RR
		tmpmatA(1,:) = (/-0.0690001, -0.00652359, 0.00140441/)
		tmpmatA(2,:) = (/-0.00652359, -0.0914355, 0.00123319/)
		tmpmatA(3,:) = (/0.00140441, 0.00123319, -0.0459147/)
		tmpmatB(1,:) = (/0.0000640169, 0.00441096, -0.000161176/)
		tmpmatB(2,:) = (/0.00441096, 0.0000512366, -0.000135525/)
		tmpmatB(3,:) = (/-0.000161176, -0.000135525, 0.0000150831/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann22 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 2,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nB, nA, f2, f3, 2,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LL
		tmpmatA(1,:) = (/-0.6890818, -0.0282017, -0.00945936/)
		tmpmatA(2,:) = (/-0.0282017, -0.1230948, 0.00119709/)
		tmpmatA(3,:) = (/-0.00945936, 0.00119709, -0.0620348/)
		tmpmatB(1,:) = (/-0.000235155, 0.0255287, -0.00317523/)
		tmpmatB(2,:) = (/0.0255287, -0.000182425, 0.000117498/)
		tmpmatB(3,:) = (/-0.00317523, 0.000117498, -0.0000496218/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann11 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 1,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nB, nA, f2, f3, 1,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LR
		tmpmatA(1,:) = (/-0.2181635, -0.00713197, 0.000594222/)
		tmpmatA(2,:) = (/-0.00713197, 0.1062286, -0.00143271/)
		tmpmatA(3,:) = (/0.000594222, -0.00143271, 0.0533432/)
		tmpmatB(1,:) = (/0.000202408, 0.00460267, - 0.00024319/)
		tmpmatB(2,:) = (/0.00460267, -0.000059526, 0.000157452/)
		tmpmatB(3,:) = (/- 0.00024319, 0.000157452, -0.0000175234/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann12 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, 1,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nB, nA, f2, f3, 1,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!RL
		tmpmatA(1,:) = (/-0.2179405, -0.00641388, 0.00228429/)
		tmpmatA(2,:) = (/-0.00641388, 0.1059528, -0.00103039/)
		tmpmatA(3,:) = (/0.00228429, -0.00103039, 0.053396/)
		tmpmatB(1,:) = (/- 0.0000743741, 0.00461848, - 0.000014447/)
		tmpmatB(2,:) = (/0.00461848, 0.000157021, - 0.000101135/)
		tmpmatB(3,:) = (/- 0.000014447, - 0.000101135, 0.0000427115/)
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
		tmpmatA(1,:) = (/0.0093356, 0.00308994, 0.000821132/)
		tmpmatA(2,:) = (/0.00308994, 0.0248959, -0.000676366/)
		tmpmatA(3,:) = (/0.000821132, -0.000676366, -0.00692737/)
		tmpmatB(1,:) = (/-0.000017038, -0.00257011, 0.000305504/)
		tmpmatB(2,:) = (/-0.00257011, -0.0000136365, 8.81481e-7/)
		tmpmatB(3,:) = (/0.000305504, 8.81481e-7, -4.01434e-6/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc22 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 2,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 2,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LL
		tmpmatA(1,:) = (/0.0931395, -0.0126822, 0.00174399/)
		tmpmatA(2,:) = (/-0.0126822, 0.033518, -0.000788527/)
		tmpmatA(3,:) = (/0.00174399, -0.000788527, -0.00933392/)
		tmpmatB(1,:) = (/0.000062586, 0.00943161, -0.000486463/)
		tmpmatB(2,:) = (/0.00943161, 0.0000485521, -0.0000787672/)
		tmpmatB(3,:) = (/-0.000486463, -0.0000787672, 0.0000132067/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc11 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 1,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 1,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LR
		tmpmatA(1,:) = (/0.0295172, 0.00292802, 0.000605503/)
		tmpmatA(2,:) = (/0.00292802, -0.0289238, 0.000785794/)
		tmpmatA(3,:) = (/0.000605503, 0.000785794, 0.00804814/)
		tmpmatB(1,:) = (/-0.0000538704, -0.00251908, 0.000283677/)
		tmpmatB(2,:) = (/-0.00251908, 0.0000158427, -1.02409e-6/)
		tmpmatB(3,:) = (/0.000283677, -1.02409e-6, 4.66381e-6/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc12 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 1,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 1,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!RL
		tmpmatA(1,:) = (/0.0294578, 0.00306075, 0.000586954/)
		tmpmatA(2,:) = (/0.00306075, -0.0288504, 0.000678718/)
		tmpmatA(3,:) = (/0.000586954, 0.000678718, 0.0080341/)
		tmpmatB(1,:) = (/0.0000197945, -0.00262534, 0.000266453/)
		tmpmatB(2,:) = (/-0.00262534, -0.0000417908, 0.0000677983/)
		tmpmatB(3,:) = (/0.000266453, +0.0000677983, -0.0000113676/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc21 test 1 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, f1, nB, f2, 2,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, f1, nB, f2, 2,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!second series
		!RR
		tmpmatA(1,:) = (/-0.0057045, -0.000732622, -0.00101446/)
		tmpmatA(2,:) = (/-0.000732622, -0.0129254, 0.000189991/)
		tmpmatA(3,:) = (/-0.00101446, 0.000189991, 0.00109802/)
		tmpmatB(1,:) = (/-0.0000239973, 0.000956323, -0.000198153/)
		tmpmatB(2,:) = (/0.000956323, -0.0000192065, 0.0000583729/)
		tmpmatB(3,:) = (/-0.000198153, 0.0000583729, -5.65404e-6/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc22 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nB, f2, nA, f3, 2,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nB, f2, nA, f3, 2,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LL
		tmpmatA(1,:) = (/-0.0572917, -0.00823674, 0.00289993/)
		tmpmatA(2,:) = (/-0.00823674, -0.0175661, 0.000431652/)
		tmpmatA(3,:) = (/0.00289993, 0.000431652, 0.00150503/)
		tmpmatB(1,:) = (/0.00008815, 0.00522956, 6.31568e-6/)
		tmpmatB(2,:) = (/0.00522956, 0.0000683837,-0.0000338275/)
		tmpmatB(3,:) = (/6.31568e-6,-0.0000338275, 0.0000186012/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc11 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nB, f2, nA, f3, 1,1, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nB, f2, nA, f3, 1,1, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!LR
		tmpmatA(1,:) = (/-0.0180364, -0.000504563, -0.000710753/)
		tmpmatA(2,:) = (/-0.000504563, 0.0150165, -0.000220729/)
		tmpmatA(3,:) = (/-0.000710753, -0.000220729, -0.00127566/)
		tmpmatB(1,:) = (/-0.0000758744,0.000884456,-0.00016741/)
		tmpmatB(2,:) = (/0.000884456,0.0000223139,-0.0000678169/)
		tmpmatB(3,:) = (/-0.00016741,-0.0000678169, 6.5688e-6/)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_sc12 test 2 full rho ',2I1)") ix,iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nB, f2, nA, f3, 1,2, ix,iy), tmpmatA(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nB, f2, nA, f3, 1,2, ix,iy), tmpmatB(ix,iy), 1d-7)
			end do
		end do
		!RL
		tmpmatA(1,:) = (/-0.01812, -0.000773746, -0.00134429/)
		tmpmatA(2,:) = (/-0.000773746, 0.0151199, -0.000371541/)
		tmpmatA(3,:) = (/-0.00134429, -0.000371541, -0.00129544/)
		tmpmatB(1,:) = (/0.0000278798,0.000878531,-0.000253156/)
		tmpmatB(2,:) = (/0.000878531,-0.0000588607,0.0000291168/)
		tmpmatB(3,:) = (/-0.000253156,0.0000291168,-0.0000160108/)
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

	subroutine do_timing_tests
		timing_tests = .true.
		call init_interp_ElDensity
		call init_interp_dme2_e
		call init_interp_FD
		call init_interp_d123
	end subroutine do_timing_tests
end program tests
