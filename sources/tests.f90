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

subroutine assert_double_verb(testname, num1, num2, tol)
	use precision
	implicit none
	character(len=*) :: testname
	real(dl) :: num1, num2, tol

	if (abs(num1 - num2) .lt. tol) then
		write(*, fmt="(a,'  ',E14.7)") testname, num1 - num2
	else
		print *, testname, " failed"
		write (*,"(*(E17.9))") num1, num2, num1 - num2, tol
		call exit()
	end if
end subroutine assert_double_verb

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
		write(*, fmt="(a,'  ',E14.7)") testname, (num1 - num2)/num1
	else
		print *, testname, " failed"
		write (*,"(*(E17.9))") num1, num2, (num1 - num2)/num1, tol
		call exit()
	end if
end subroutine assert_double_rel_verb

program tests
	use precision
	use variables
	use ndConfig
	use ndErrors
	use ndEquations
	use ndstuff
	implicit none

	call openLogFile
	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "Initializations"
	call do_tests_initialization
	call init_interp_jkyg12
	call init_interp_dme2_e
	call init_interp_ElDensity
	call allocate_interpNuDens

	write(*,*) ""
	write(*,"(a)") "Starting tests"
	call do_basic_tests
	call do_test_linearized_integrals
	call do_test_commutator
	call do_test_nu_matrices
	call do_tests_cosmology
	call do_tests_JKYG
	call do_tests_dzodx
	call do_tests_dme2
	call do_tests_Di
	call do_tests_Pi_ij
	call do_f_ann_sc_re_tests_eq
	call do_f_ann_sc_re_tests_full
	call do_tests_coll_int
	call do_test_drho_dx
	call do_test_collision_terms
	call do_test_damping_factors
	call do_test_zin

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
		outputFolder = ""
		maxiter = 100
		toler   = 1.d-5
		toler_dme2 = 1.d-5
		toler_jkyg = 1.d-7
		toler_ed = 1.d-4
		dlsoda_atol = 1.d-6
		dlsoda_rtol = 1.d-6
		Nx = 100
		Ny = 100
		interp_nx = interp_nx0
		interp_nz = interp_nz0
		interp_nxz = interp_nxz0
		interp_zmin = interp_zmin0
		interp_zmax = interp_zmax0
		allocate(x_arr(Nx), y_arr(Ny))
		allocate(dy_arr(Ny), fy_arr(Ny))
		x_in    = 0.001d0
		x_fin   = 40.d0
		logx_in  = log10(x_in)
		logx_fin = log10(x_fin)
		x_arr = logspace(logx_in, logx_fin, Nx)
		y_min = 0.01d0
		y_max = 20.0d0
		y_cen = 0.01d0
		Nylog = 1
		y_arr = linspace(y_min, y_max, Ny)
		do ix=1, Ny-1
			dy_arr(ix) = y_arr(ix+1) - y_arr(ix)
		end do
		dy_arr(Ny) = 0.d0
		m_lightest   = 0.0d0
		collision_offdiag = 1
		dme2_temperature_corr = .true.
		flavorNumber = 3
		flavNumSqu = flavorNumber**2
		call allocateStuff
		do ix=1, flavorNumber
			nuFactor(ix) = 1.d0
			sterile(ix) = .false.
		end do
		tot_factor_active_nu = 3.0
		tot_factor_nu = 0.d0
		theta12      = 0.5840d0
		dm21         = 7.53d-05
		if (flavorNumber .gt. 2) then
			theta13      = 0.1485d0
			theta23      = 0.7954d0
			dm31         = 0.0025153d0
			deltaCP13    = 0.d0
		end if
		z_in=1.0000575
		save_w_evolution = .true.
		save_nuDens_evolution = .true.
		call setMixingMatrix()
		call setMassMatrix()
		call init_matrices
		allocate(interp_xvec(interp_nx), interp_zvec(interp_nz), interp_xozvec(interp_nxz))
		interp_xvec = logspace(logx_in, logx_fin, interp_nx)
		interp_zvec = linspace(interp_zmin, interp_zmax, interp_nz)
		interp_xozvec = logspace(log10(x_in/interp_zmax), logx_fin, interp_nxz)
	end subroutine do_tests_initialization

	subroutine do_basic_tests
		real(dl), dimension(:), allocatable :: tarr
		integer :: N, ix
		character(len=30) :: tmparg

		N = 12
		allocate(tarr(N))
		tarr = linspace(0.0d0, 11.d0, N)
		do ix=1, N
			write(tmparg, '("tarr lin ",I2)') ix
			call assert_double(tmparg, tarr(ix), (ix-1)*1.d0, 1d-7)
		end do
		tarr = logspace(-5d0, 6.d0, N)
		do ix=1, N
			write(tmparg, '("tarr log ",I2)') ix
			call assert_double_rel(tmparg, tarr(ix), 10.d0**(ix-6), 1d-7)
		end do
		tarr = loglinspace(0.01d0, 1.d0, 10.d0, N, 3)
		call assert_double_rel("tarr linlog 1", tarr(1), 0.01d0, 1d-7)
		call assert_double_rel("tarr linlog 2", tarr(2), 0.1d0, 1d-7)
		do ix=3, N
			write(tmparg, '("tarr linlog ",I2)') ix
			call assert_double_rel(tmparg, tarr(ix), (ix-2)*1.d0, 1d-7)
		end do

		call assert_double_rel("y_arr linlog 1", y_arr(1), 0.01d0, 1d-7)
		call assert_double_rel("y_arr linlog 2", y_arr(2), 0.21191919191919d0, 1d-7)
	end subroutine do_basic_tests

	subroutine do_test_linearized_integrals
		integer :: ia, ib
		real(dl), dimension(:), allocatable :: fy1_arr
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(fy1_arr(Ny))
		allocate(fy2_arr(Ny, Ny))

		write(*,*) ""
		write(*,"(a)") "integral_linearized (4 tests)"

		do ia=1, Ny
			fy1_arr(ia) = 0.35d0*y_arr(ia) + 11.41d0
		end do
		call assert_double_rel("intLin1D 1", integral_linearized_1d(Ny, dy_arr, fy1_arr), 298.08588d0, 1d-7)
		do ia=1, Ny
			fy1_arr(ia) = 0.35d0/(exp(y_arr(ia))+1.d0) + 0.2d0
		end do
		call assert_double_rel("intLin1D 2", integral_linearized_1d(Ny, dy_arr, fy1_arr), 4.23886d0, 1d-4)

		do ia=1, Ny
			do ib=1, Ny
				fy2_arr(ia,ib) = 0.35d0*y_arr(ia) - 1.96d0*y_arr(ib) + 11.41d0
			end do
		end do
		call assert_double_rel("intLin2D 1", integral_linearized_2d(Ny, Ny, dy_arr, dy_arr, fy2_arr), -1877.341249d0, 1d-7)
		do ia=1, Ny
			do ib=1, Ny
				fy2_arr(ia,ib) = 0.35d0/(exp(y_arr(ia))+1.d0) + 1.96d0/(exp(y_arr(ib))+1.d0) + 0.2d0
			end do
		end do
		call assert_double_rel("intLin2D 2", integral_linearized_2d(Ny, Ny, dy_arr, dy_arr, fy2_arr), 111.967d0, 0.003d0)

		deallocate(fy1_arr)
		deallocate(fy2_arr)
	end subroutine do_test_linearized_integrals

	subroutine do_test_commutator
		real(dl), dimension(:,:), allocatable :: m1, m2, m3, res
		integer :: i,j
		character(len=300) :: tmparg

		allocate(m1(4,4), m2(4,4), res(4,4))
		m1(1,:) = (/0.2, 0.6, 3.01, 0.3/)
		m1(2,:) = (/-0.12, 0.045, -0.05, 0./)
		m1(3,:) = (/0., 0.11, 0.22, -0.11/)
		m1(4,:) = (/1., 4., -0.03, 0.04/)
		m2(1,:) = (/1., 2., -1., -2./)
		m2(2,:) = (/0.5, -0.6, 0.7, -0.8/)
		m2(3,:) = (/0.01, 0.66, 0.55, 0.31/)
		m2(4,:) = (/0.98, 0.4, -0.78, 0.4/)

		write(*,*)
		write(*,"(a)") "Commutator & anticommutator (32 tests)"
		call Commutator(m1, m2, m3)
		res(1,:) = (/2.8641, 9.5666, -1.1085, -0.1569/)
		res(2,:) = (/0.53, 2.55, -1.589, 0.1475/)
		res(3,:) = (/-0.2834, -1.301, 0.175, -0.0187/)
		res(4,:) = (/2.4909, -2.524, -0.9939, -5.5891/)
		do i=1,4
			do j=1,4
				write(tmparg,"('commutator ',2I1)") i,j
				call assert_double_rel(trim(tmparg), m3(i,j), res(i,j), 2d-7)
			end do
		end do
		call Anticommutator(m1, m2, m3)
		res(1,:) = (/-1.2159, -5.2734, 4.3915, 0.5031/)
		res(2,:) = (/-0.726, -3.15, 1.837, 0.2295/)
		res(3,:) = (/0.1822, 1.3714, 0.3926, -0.1089/)
		res(4,:) = (/3.5869, 1.7164, 4.4985, -4.7975/)
		do i=1,4
			do j=1,4
				write(tmparg,"('anticommutator ',2I1)") i,j
				call assert_double_rel(trim(tmparg), m3(i,j), res(i,j), 2d-7)
			end do
		end do
		deallocate(m1, m2, m3, res)
	end subroutine do_test_commutator

	subroutine do_test_nu_matrices
		real(dl), dimension(3,3) :: m
		real(dl), dimension(:,:), allocatable :: ide
		integer :: i,j
		character(len=300) :: tmparg

		write(*,*)
		write(*,"(a)") "Neutrino and lepton matrices (63 tests)"

		m(1,:) = (/0.825082, 0.54529713, 0.1479548/)
		m(2,:) = (/-0.47410446, 0.525726198, 0.7062839/)
		m(3,:) = (/0.30735086, -0.6528882, 0.69229506/)
		do i=1,3
			do j=1,3
				write(tmparg,"('mixing matrix ',2I1)") i,j
				call assert_double_rel(trim(tmparg), mixMat(i,j), m(i,j), 1d-7)
			end do
		end do
		do i=1,3
			do j=1,3
				write(tmparg,"('inverse mixing matrix ',2I1)") i,j
				call assert_double_rel(trim(tmparg), mixMatInv(i,j), m(j,i), 1d-7)
			end do
		end do
		call createIdentityMat(ide, 3)
		m = matmul(mixMat, mixMatInv)
		do i=1,3
			do j=1,3
				write(tmparg,"('mixing matrix product A ',2I1)") i,j
				call assert_double(trim(tmparg), m(i,j), ide(i,j), 1d-7)
			end do
		end do
		m = matmul(mixMatInv, mixMat)
		do i=1,3
			do j=1,3
				write(tmparg,"('mixing matrix product B ',2I1)") i,j
				call assert_double(trim(tmparg), m(i,j), ide(i,j), 1d-7)
			end do
		end do
		deallocate(ide)

		m(1,:) = (/0.00007745186685918483, 0.00028443083990052604, 0.00023082995114815183/)
		m(2,:) = (/0.00028443083990052604, 0.001275536536721503, 0.0012040271446586716/)
		m(3,:) = (/0.00023082995114815185, 0.0012040271446586716, 0.001237611596419312/)
		do i=1,3
			do j=1,3
				write(tmparg,"('mass matrix ',2I1)") i,j
				call assert_double(trim(tmparg), nuMassesMat(i,j), m(i,j), 1d-7)
			end do
		end do

		leptonDensities=0.d0
		call updateLeptonDensities(0.076d0, 1.32d0)
		m(1,:) = (/-0.0026760938/1.22, 0., 0./)
		m(2,:) = (/0.,-2.58816e-6/1.22,0./)
		m(3,:) = (/0.,0.,0./)
		do i=1,3
			do j=1,3
				write(tmparg,"('lepton matrix A ',2I1)") i,j
				if (m(i,j).gt.1d-16) then
					call assert_double_rel(trim(tmparg), leptonDensities(i,j), m(i,j), 1d-6)
				else
					call assert_double(trim(tmparg), leptonDensities(i,j), m(i,j), 1d-7)
				end if
			end do
		end do
		call updateLeptonDensities(2.d1, 1.2d0)
		m(1,:) = (/-1.7572758d-23/0.22d0, 0.d0, 0.d0/)
		m(2,:) = (/0.,0.,0./)
		m(3,:) = (/0.,0.,0./)
		do i=1,3
			do j=1,3
				write(tmparg,"('lepton matrix B ',2I1)") i,j
				call assert_double(trim(tmparg), leptonDensities(i,j), m(i,j), 5d-26)
			end do
		end do
	end subroutine do_test_nu_matrices

	subroutine do_tests_cosmology
		real(dl), dimension(:), allocatable :: ndmv_re
		integer :: i, iy

		allocate(ndmv_re(Ny))
		do i=1, flavorNumber
			do iy=1, Ny
				ndmv_re(iy) = 1.d0*i
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i
			end do
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
		end do

		write(*,*) ""
		write(*,"(a)") "Cosmology (35 tests)"
		call assert_double_rel("elDensF test 1", electronDensityFull(1.d0, 1.d0), 1.06102d0, 1d-4)
		call assert_double_rel("elDensF test 2", electronDensityFull(0.076d0, 1.32d0), 3.48762d0, 1d-4)
		call assert_double_rel("elDensF test 3", electronDensityFull(1.d1, 1.2d0), 0.0377464d0, 1d-4)
		call assert_double_rel("elDensF test 4", electronDensityFull(2.d1, 1.2d0), 0.0000421797d0, 1d-4)
		call assert_double_rel("elDensF test 5", electronDensityFull(3.d1, 1.2d0), 2.61396d-8, 5d-3)
		call assert_double_rel("elDens test 1", electronDensity(1.d0, 1.d0), 1.06102d0, 1d-4)
		call assert_double_rel("elDens test 2", electronDensity(0.076d0, 1.32d0), 3.48762d0, 1d-4)
		call assert_double_rel("elDens test 3", electronDensity(1.d1, 1.2d0), 0.0377464d0, 1d-2)
		call assert_double_rel("elDens test 4", electronDensity(2.d1, 1.2d0), 0.0000421797d0, 2d-2)
		call assert_double_rel("elDens test 5", electronDensity(3.d1, 1.2d0), 2.61396d-8, 5d-2)
		call assert_double_rel("photDens test 1", photonDensity(1.002d0), 0.66325322d0, 1d-7)
		call assert_double_rel("photDens test 2", photonDensity(1.34d0), 2.12142498d0, 1d-7)

		call assert_double_rel("nuDens test 1", nuDensity(1.d0, 1), 0.575727d0, 1d-4)
		call assert_double_rel("nuDens test 2", nuDensity(1.076d0, 1), 0.575727d0, 1d-4)
		call assert_double_rel("nuDens test 3", nuDensity(1.32d0, 1), 0.575727d0, 3d-4)
		call assert_double_rel("nuDens test 4", nuDensity(1.37d0, 2), 2.d0*0.575727d0, 5d-4)
		call assert_double_rel("nuDens test 5", nuDensity(1.003d0, 3), 3.d0*0.575727d0, 1d-4)

		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy) / 1.d0)
			end do
		end do
		call assert_double_rel("nuDensLin test 1", nuDensityLin(1), 0.575727d0, 1d-4)
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy) / 1.076d0)
			end do
		end do
		call assert_double_rel("nuDensLin test 2", nuDensityLin(1), 0.5d0*1.54346d0, 1d-4)
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy) / 1.32d0)
			end do
		end do
		call assert_double_rel("nuDensLin test 3", nuDensityLin(1), 0.5d0*3.49577d0, 2d-4)
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy) / 1.37d0)
			end do
		end do
		call assert_double_rel("nuDensLin test 4", nuDensityLin(2), 0.5d0*2.d0*4.05629d0, 5d-4)
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy) / 1.003d0)
			end do
		end do
		call assert_double_rel("nuDensLin test 5", nuDensityLin(3), 0.5d0*3.d0*1.16533d0, 1d-4)

		call assert_double_rel("nuDensLinEq test 1", nuDensityLinEq(1.d0), 0.575727d0, 1d-4)
		call assert_double_rel("nuDensLinEq test 2", nuDensityLinEq(1.37d0), 2.02814d0, 5d-4)

		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy))
			end do
		end do
		call assert_double_rel("allNuDensLin test 1", allNuDensity(), 6*0.575727d0, 1d-3)

		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
			end do
		end do
		call assert_double("radDens test 1", radDensity(0.7d0, 1.04d0), 3.79519d0, 1d-3)
		call assert_double("radDens test 2", radDensity(1.d0, 1.04d0), 3.74595d0, 1d-3)
		call assert_double("radDens test 3", radDensity(1.d0, 1.24d0), 5.86425d0, 3d-3)

		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
			end do
		end do
		call assert_double_rel("Neff test 1", Neff_from_rho_z(1.d0, zid), 3.0d0, 1d-5)
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
			end do
		end do
		call assert_double_rel("Neff test 2", Neff_from_rho_z(1.d0, 1.39975d0), 3.011d0, 1d-4)
		do iy=1, Ny
			nuDensMatVecFD(iy)%re(1, 1) = 1.2d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(2, 2) = 1.d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(3, 3) = 1.d0 * fermiDirac(y_arr(iy))
		end do
		call assert_double_rel("Neff test 3", Neff_from_rho_z(1.d0, zid), 3.2d0, 1d-5)
		do iy=1, Ny
			nuDensMatVecFD(iy)%re(1, 1) = 1.00920d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(2, 2) = 1.00392d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(3, 3) = 1.00392d0 * fermiDirac(y_arr(iy))
		end do
		call assert_double("Neff test 4", Neff_from_rho_z(1.d0, 1.397843d0), 3.045d0, 1d-3)
		do iy=1, Ny
			nuDensMatVecFD(iy)%re(1, 1) = 1.00699d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(2, 2) = 1.00511d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(3, 3) = 1.00519d0 * fermiDirac(y_arr(iy))
		end do
		call assert_double("Neff test 5", Neff_from_rho_z(1.d0, 1.39779d0), 3.045d0, 1d-3)

		call assert_double("muonDensF test 1", muonDensityFull(1.d0, 1.d0), 0.d0, 1d-15)
		call assert_double_rel("muonDensF test 2", muonDensityFull(0.1d0, 1.32d0), 0.000146007d0, 1d-4)
		call assert_double_rel("muonDensF test 3", muonDensityFull(0.01d0, 1.2d0), 1.86472d0, 1d-4)
		call assert_double_rel("muonDensF test 4", muonDensityFull(3.d-3, 1.2d0), 2.3396d0, 1d-4)
		call assert_double_rel("muonDensF test 5", muonDensityFull(1.d-3, 1.d0), 1.14785d0, 1d-4)
		call assert_double("muonDensF test 6", muonDensityFull(4.d0, 1.3d0), 0.d0, 1d-15)

		call assert_double("muonDens test 1", muonDensity(1.d0, 1.d0), 0.d0, 1d-15)
		call assert_double_rel("muonDens test 2", muonDensity(1.d-1, 1.32d0), 0.000146007d0, 1.1d-2)
		call assert_double_rel("muonDens test 3", muonDensity(1.d-2, 1.2d0), 1.86472d0, 1d-4)
		call assert_double_rel("muonDens test 4", muonDensity(3.d-3, 1.2d0), 2.3396d0, 1d-4)
		call assert_double_rel("muonDens test 5", muonDensity(1.d-3, 1.d0), 1.14785d0, 1d-4)
		call assert_double("muonDens test 6", muonDensity(4.d0, 1.3d0), 0.d0, 1d-15)

		do iy=1, Ny
			nuDensMatVecFD(iy)%re(1, 1) = 1.d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(2, 2) = 1.d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(3, 3) = 1.d0 * fermiDirac(y_arr(iy))
		end do
		call assert_double("radDens test 4", radDensity(0.01d0, 1.24d0), 8.16007d0, 1d-3)
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
		integer :: n
		real(dl), dimension(:), allocatable :: ydot
		integer :: m

		n=ntot
		allocate(ydot(n))
		ydot = 0.d0
		do m=1, Ny
			ydot((m-1)*flavNumSqu + 1) = cos(0.02d0*y_arr(m))
			ydot((m-1)*flavNumSqu + 2) = y_arr(m)/20.d0
			ydot((m-1)*flavNumSqu + 3) = 1.d0
		end do
		write(*,*) ""
		write(*,"(a)") "dz/dx functions (6 tests)"
		call dz_o_dx(0.01d0, 1.d0, ydot, n)
		call assert_double_rel("dz_o_dx test 1", ydot(n), -0.1751102d0, 4d-6)
		call dz_o_dx(1.1d0, 1.1d0, ydot, n)
		call assert_double_rel("dz_o_dx test 2", ydot(n), -0.0946571d0, 5d-6)
		call dz_o_dx(6.d0, 1.2d0, ydot, n)
		call assert_double_rel("dz_o_dx test 3", ydot(n), -0.15262978d0, 4d-6)

		call dz_o_dx_lin(0.01d0, 1.2d0, 1.d0, ydot, n)
		call assert_double_rel("dz_o_dx_lin test 1a", ydot(n), -0.1751102d0, 4d-6)
		call assert_double_rel("dz_o_dx_lin test 1b", ydot(n-1), -0.10615d0, 4d-6)
		call dz_o_dx_lin(1.1d0, 1.1d0, 1.1d0, ydot, n)
		call assert_double_rel("dz_o_dx_lin test 2a", ydot(n), -0.0946571d0, 3d-6)
		call assert_double_rel("dz_o_dx_lin test 2b", ydot(n-1), -0.137812d0, 3d-6)
		call dz_o_dx_lin(6.d0, 1.d0, 1.2d0, ydot, n)
		call assert_double_rel("dz_o_dx_lin test 3a", ydot(n), -0.15262978d0, 2d-6)
		call assert_double_rel("dz_o_dx_lin test 3b", ydot(n-1), -0.183428d0, 2d-6)
		deallocate(ydot)
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
		write(*,"(a)") "D_i functions (48 tests)"
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

		call assert_double("D2c test 1", D2_cases(0.1d0, 0.2d0, 0.3d0, 0.4d0), -0.00133333d0, 1d-7)
		call assert_double("D2c test 2", D2_cases(0.4d0, 0.2d0, 0.3d0, 0.1d0), -0.0213333d0, 1d-7)
		call assert_double("D2c test 3", D2_cases(0.01d0,5.d0,2.6d0,2.41d0), -1.333333333d-6, 1d-12)
		call assert_double("D2c test 4", D2_cases(10.03d0,5.4d0,8.8d0,6.63d0), -209.952d0, 1d-4)
		call assert_double_rel("D2c test 5", D2_cases(10.d0,2.d0,5.3d0,6.7d0), -10.666667d0, 1d-7)
		call assert_double_rel("D2c test 6", D2_cases(2.d0,10.d0,6.7d0,5.3d0), -10.666667d0, 1d-7)
		call assert_double("D3c test 1", D3_cases(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.000192d0, 1d-7)
		call assert_double("D3c test 2", D3_cases(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.000192d0, 1d-7)
		call assert_double("D3c test 3", D3_cases(0.01d0,5.d0,2.6d0,2.41d0), 2.50454d-5, 1d-10)
		call assert_double_rel("D3c test 4", D3_cases(10.03d0,5.4d0,8.8d0,6.63d0), 22692.22d0, 1d-7)
		call assert_double_rel("D3c test 5", D3_cases(10.d0,2.d0,5.3d0,6.7d0), 918.2933d0, 1d-7)
		call assert_double_rel("D3c test 6", D3_cases(2.d0,10.d0,6.7d0,5.3d0), 918.2933d0, 1d-7)

		call assert_double("D1p test 1", D1_bis(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.4d0, 1d-7)
		call assert_double("D1p test 2", D1_bis(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.4d0, 1d-7)
		call assert_double("D1p test 3", D1_bis(0.01d0,5.d0,2.6d0,2.41d0), 0.04d0, 1d-7)
		call assert_double("D1p test 4", D1_bis(10.03d0,5.4d0,8.8d0,6.63d0), 21.6d0, 1d-4)
		call assert_double("D1p test 5", D1_bis(10.d0,2.d0,5.3d0,6.7d0), 8.d0, 1d-7)
		call assert_double("D1p test 6", D1_bis(2.d0,10.d0,6.7d0,5.3d0), 8.d0, 1d-7)
		call assert_double("D2p test 1", D2_bis(0.1d0, 0.2d0, 0.3d0, 0.4d0), -0.00133333d0, 1d-7)
		call assert_double("D2p test 2", D2_bis(0.4d0, 0.2d0, 0.3d0, 0.1d0), -0.0213333d0, 1d-7)
		call assert_double("D2p test 3", D2_bis(0.01d0,5.d0,2.6d0,2.41d0), -1.333333333d-6, 1d-12)
		call assert_double("D2p test 4", D2_bis(10.03d0,5.4d0,8.8d0,6.63d0), -209.952d0, 1d-4)
		call assert_double_rel("D2p test 5", D2_bis(10.d0,2.d0,5.3d0,6.7d0), -10.666667d0, 1d-7)
		call assert_double_rel("D2p test 6", D2_bis(2.d0,10.d0,6.7d0,5.3d0), -10.666667d0, 1d-7)
		call assert_double("D3p test 1", D3_bis(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.000192d0, 1d-7)
		call assert_double("D3p test 2", D3_bis(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.000192d0, 1d-7)
		call assert_double("D3p test 3", D3_bis(0.01d0,5.d0,2.6d0,2.41d0), 2.50454d-5, 1d-10)
		call assert_double_rel("D3p test 4", D3_bis(10.03d0,5.4d0,8.8d0,6.63d0), 22692.22d0, 1d-7)
		call assert_double_rel("D3p test 5", D3_bis(10.d0,2.d0,5.3d0,6.7d0), 918.2933d0, 1d-7)
		call assert_double_rel("D3p test 6", D3_bis(2.d0,10.d0,6.7d0,5.3d0), 918.2933d0, 1d-7)
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
		real(dl) :: fdA, fdB, f1, f2, f3, r1, r2
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

		write(*,*)
		write(*,"(a)") "F_ann functions, final tests (36 tests)"
		!rhoA = {{0.9, 0.011 + 0.0001 I, -0.03 - 0.004 I}, {0.011 + 0.0001 I, 0.86, 0.001 I}, {-0.03 - 0.004 I, 0.001 I, 0.96}};
		nA%re(1,:) = (/0.9d0, 0.011d0, -0.03d0/)
		nA%re(2,:) = (/0.011d0, 0.86d0, 0.d0/)
		nA%re(3,:) = (/-0.03d0, 0.d0, 0.96d0/)
		nA%im(1,:) = (/0.d0, 0.0001d0, -0.004d0/)
		nA%im(2,:) = (/-0.0001d0, 0.d0, 0.001d0/)
		nA%im(3,:) = (/0.004d0, -0.001d0, 0.d0/)
		!rhoB = rhoA;
		nB = nA
		!RR
		tmpmatA(1,:) = (/-0.023346, -0.000578945, 0.00164337/)
		tmpmatA(2,:) = (/-0.000578945, -0.0212209, 7.06991e-6/)
		tmpmatA(3,:) = (/0.00164337, 7.06991e-6, -0.0266302/)
		tmpmatB(1,:) = (/0., -5.90586e-6, 0.00021888/)
		tmpmatB(2,:) = (/5.90586e-6, 0., -0.0000530457/)
		tmpmatB(3,:) = (/-0.00021888, 0.0000530457, 0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann22 test 1 full rho ',2I1)") ix,iy
				call assert_double_rel(trim(tmparg)//"re", F_ab_ann_re(nA, nB, 0.1d0, 0.7d0, 2,2, ix,iy), tmpmatA(ix,iy), 1d-5)
				if (abs(tmpmatB(ix,iy)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, 0.1d0, 0.7d0, 2,2, ix,iy), tmpmatB(ix,iy), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", F_ab_ann_im(nA, nB, 0.1d0, 0.7d0, 2,2, ix,iy), tmpmatB(ix,iy), 1d-5)
				end if
			end do
		end do
		!LL
		tmpmatA(1,:) = (/0.0580013, -0.00422802, 0.0113064/)
		tmpmatA(2,:) = (/-0.00422802, 0.0104606, -0.0000779103/)
		tmpmatA(3,:) = (/0.0113064, -0.0000779103, 0.00354312/)
		tmpmatB(1,:) = (/0., -0.0000361964, 0.00150834/)
		tmpmatB(2,:) = (/0.0000361964, 0., -0.0000807178/)
		tmpmatB(3,:) = (/-0.00150834, 0.0000807178, 0./)
		do ix=1, flavorNumber
			do iy=1, flavorNumber
				write(tmparg,"('F_ann22 test 1 full rho ',2I1)") ix,iy
				call assert_double_rel(trim(tmparg)//"re", F_ab_sc_re(nA, 0.1d0, nB, 0.7d0, 1,1, ix,iy), tmpmatA(ix,iy), 1d-5)
				r2 = F_ab_sc_im(nA, 0.1d0, nB, 0.7d0, 1,1, ix,iy)
				if (abs(tmpmatB(ix,iy)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", r2, tmpmatB(ix,iy), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", r2, tmpmatB(ix,iy), 1d-5)
				end if
			end do
		end do
		call deallocateCmplxMat(nA)
		call deallocateCmplxMat(nB)
	end subroutine do_f_ann_sc_re_tests_full

	subroutine do_tests_coll_int
		real(dl) :: x,z,dme2
		type(coll_args) :: collArgs
		integer :: i, j, iy1, iy
		real(dl) :: errr1,errr2, res1,res2,res3,res4, cf, ref
		INTEGER :: IFAIL, ITRANS, N, NPTS, NRAND
		real(dl) ::  VK(2)
		real(dl), dimension(:), allocatable :: ndmv_re, tmpfy2_arr
		real(dl), dimension(3) :: tmparrS, tmparrA, tmperr, tmperr3,tmperr4
		real(dl), dimension(3, 3) :: tmpmatA, tmpmatB
		character(len=300) :: tmparg
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(fy2_arr(Ny, Ny))

		write(*,*) ""
		write(*,"(a)") "Collision integrals (57 tests)"

		x=0.05d0
		iy1=7 !1.22151515151515
		z=1.06d0
		dme2=0.1d0
		npts=6
		nrand=1
		n=2

		allocate(ndmv_re(Ny))

		do iy=1, Ny
			nuDensMatVecFD(iy)%re = 0.d0
			nuDensMatVecFD(iy)%im = 0.d0
		end do
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy))
				ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, i)
			end do
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
			nuDensMatVecFD(iy1)%re(i, i) = 1.d0
		end do

		collArgs%x = x
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%y2 = 0.d0
		collArgs%y3 = 0.d0
		collArgs%y4 = 0.d0
		collArgs%dme2 = dme2

		!first series of fake sc and ann tests
		vk=(/5.2d0,2.3613731184d0/)
		tmparrS = (/-5.0770022d0, -1.0769122d0, -1.0577857d0/)
		tmparrA = (/1d-6, 1d-6, 1d-6/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i
			res2=coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
			call assert_double_rel("test coll sc 3 fake "//trim(tmparg), res2, tmparrS(i), tmparrA(i))
			res1=coll_nue_4_sc_int_re(n, vk, collArgs)
			call assert_double_rel("test coll sc 4 fake "//trim(tmparg), res1, tmparrS(i), tmparrA(i))
!			write(*,"(I2,*(E17.9))") i, res1, res2, tmparrS(i)
		end do

		vk=(/2.229999170d0, 3.0d0/)
		tmparrS = (/-2.3278997d0, -0.798282d0, -1.1974237d0/)
		tmparrA = (/1d-6, 1d-6, 1d-6/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i
			res2=coll_nue_3_ann_int_re(21, 3.d0, collArgs)
			call assert_double_rel("test coll ann 3 fake "//trim(tmparg), res2, tmparrS(i), tmparrA(i))
			res1=coll_nue_4_ann_int_re(n, vk, collArgs)
			call assert_double_rel("test coll ann 4 fake "//trim(tmparg), res1, tmparrS(i), tmparrA(i))
!			write(*,"(I2,*(E17.9))") i, res1, res2, tmparrS(i)
		end do

		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy))
				ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, i)
			end do
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
		end do

		!first series of sc and ann tests
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		do i=1, Ny
			do j=1, Ny
				fy2_arr(i,j) = coll_nue_3_sc_int_re(i, y_arr(j), collArgs)
			end do
		end do
		write(*,*) ""
		allocate(tmpfy2_arr(Ny))
		call openFile(987, "test_sc3_re.dat", .true.)
		do i=1, Ny
			tmpfy2_arr = fy2_arr(i,:)
			write(987,multidblfmt) tmpfy2_arr
		end do
		close(987)
		do i=1, Ny
			do j=1, Ny
				fy2_arr(i,j) = coll_nue_3_ann_int_re(i, y_arr(j), collArgs)
			end do
		end do
		call openFile(987, "test_ann3_re.dat", .true.)
		do i=1, Ny
			tmpfy2_arr = fy2_arr(i,:)
			write(987,multidblfmt) tmpfy2_arr
		end do
		close(987)

		vk=(/5.2d0,2.3613731184d0/)
		tmparrS = (/-0.170915d0, -0.1904502d0, -0.46065747d0/)
		tmparrA = (/1d-6, 1d-6, 1d-6/)
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		res2=coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
		call assert_double("test coll sc 3 sing 1", res2, tmparrS(1), tmparrA(1))
		res1=coll_nue_4_sc_int_re(n, vk, collArgs)
		call assert_double("test coll sc 4 sing 1", res1, tmparrS(1), tmparrA(1))
		do i=2, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i
			res2=coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
			call assert_double_rel("test coll sc 3 sing "//trim(tmparg), res2, tmparrS(i), tmparrA(i))
			res1=coll_nue_4_sc_int_re(n, vk, collArgs)
			call assert_double_rel("test coll sc 4 sing "//trim(tmparg), res1, tmparrS(i), tmparrA(i))
		end do

		write(*,*) ""
		tmparrS = (/-2.8769d0, -3.4185d0, -8.39213d0/)
		tmperr3 = (/2d-3, 2d-3, 2d-3/)
		tmperr4 = (/2d-2, 2d-2, 2d-2/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			ifail=0
			itrans=0
!			call D01GCF(n,coll_nue_4_sc_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
!			write(tmparg,"('test coll sc 4 - ',2I1)") i, i
!			call assert_double_rel_verb(trim(tmparg), res1, tmparrS(i), tmperr4(i))
			res2 = integrate_coll_int_3(coll_nue_3_sc_int_re, collArgs)
			write(tmparg,"('test coll sc 3 - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrS(i), tmperr3(i))
!			write(*,"(I2,*(E17.9))") i, res1, res2, tmparrS(i)
		end do

		vk=(/2.229999170d0, 3.0d0/)
		tmparrS = (/0.18420472d0, -0.2786386d0, -0.76934608d0/)
		tmperr = (/1d-6, 1d-6, 1d-6/)
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i
			res2=coll_nue_3_ann_int_re(21, 3.d0, collArgs)
			call assert_double_rel("test coll ann 3 sing "//trim(tmparg), res2, tmparrS(i), tmperr(i))
			res1=coll_nue_4_ann_int_re(n, vk, collArgs)
			call assert_double_rel("test coll ann 4 sing "//trim(tmparg), res1, tmparrS(i), tmperr(i))
		end do

		write(*,*) ""
		tmparrA = (/2.81889d0, -4.67305d0, -12.9279d0/)
		tmperr3 = (/0.02d0, 0.02d0, 0.02d0/)
		tmperr4 = (/0.2d0, 0.1d0, 0.07d0/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			ifail=0
			itrans=0
!			call D01GCF(n,coll_nue_4_ann_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
!			write(tmparg,"('test coll ann 4 - ',2I1)") i, i
!			call assert_double_rel_verb(trim(tmparg), res1, tmparrA(i), tmperr4(i))
			res2 = integrate_coll_int_3(coll_nue_3_ann_int_re, collArgs)
			write(tmparg,"('test coll ann 3 - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrA(i), tmperr3(i))
!			write(*,"(I2,*(E17.9))") i, res1, res2, tmparrA(i)
		end do

		!second series of sc and ann tests
		iy1 = 67
		z = 1.3d0
		collArgs%y1 = y_arr(iy1)
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
				ndmv_re(iy) = nuDensMatVecFD(iy)%re(i, i)
			end do
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
		end do
		collArgs%x = 0.5d0
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%ix1 = 3
		collArgs%ix2 = 3
		do i=1, Ny
			do j=1, Ny
				fy2_arr(i,j) = coll_nue_3_sc_int_re(i, y_arr(j), collArgs)
			end do
		end do
		call openFile(987, "test_sc3b_re.dat", .true.)
		do i=1, Ny
			tmpfy2_arr = fy2_arr(i,:)
			write(987,multidblfmt) tmpfy2_arr
		end do
		close(987)
		do i=1, Ny
			do j=1, Ny
				fy2_arr(i,j) = coll_nue_3_ann_int_re(i, y_arr(j), collArgs)
			end do
		end do
		call openFile(987, "test_ann3b_re.dat", .true.)
		do i=1, Ny
			tmpfy2_arr = fy2_arr(i,:)
			write(987,multidblfmt) tmpfy2_arr
		end do
		close(987)

		vk=(/5.2d0,14.50977264d0/)
		tmparrS = (/0.044973361d0, 0.0096546081d0, 0.0096546081d0/)
		tmperr = (/1d-6, 1d-6, 1d-6/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i+3
			res2=coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
			call assert_double_rel("test coll sc 3 sing "//trim(tmparg), res2, tmparrS(i), tmperr(i))
			res1=coll_nue_4_sc_int_re(n, vk, collArgs)
			call assert_double_rel("test coll sc 4 sing "//trim(tmparg), res1, tmparrS(i), tmperr(i))
!			write(*,"(I2,*(E17.9))") i, res1, res2, tmparrS(i)
		end do

		write(*,*) ""
		tmparrS = (/2.13185d0, 0.458106d0, 0.458106d0/)
		tmperr3 = (/3d-4, 4d-4, 4d-4/)
		tmperr4 = (/1d-6, 1d-6, 1d-6/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			ifail=0
			itrans=0
!			call D01GCF(n,coll_nue_4_sc_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
!			write(tmparg,"('test coll sc 4 b - ',2I1)") i, i
!			call assert_double_rel_verb(trim(tmparg), res1, tmparrS(i), tmperr4(i))
			res2 = integrate_coll_int_3(coll_nue_3_sc_int_re, collArgs)
			write(tmparg,"('test coll sc 3 b - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrS(i), tmperr3(i))
			write(*,"(I2,*(E17.9))") i, res1, res2, tmparrS(i)
		end do

		collArgs%ix1 = 1
		collArgs%ix2 = 1
		vk=(/14.31505386d0, 3.0d0/)
		tmparrS = (/0.068397661d0, 0.00955919314d0, 0.00955919314d0/)
		tmperr = (/1d-6, 1d-6, 1d-6/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i+3
			res2=coll_nue_3_ann_int_re(21, 3.d0, collArgs)
			call assert_double_rel("test coll ann 3 sing "//trim(tmparg), res2, tmparrS(i), tmperr(i))
			res1=coll_nue_4_ann_int_re(n, vk, collArgs)
			call assert_double_rel("test coll ann 4 sing "//trim(tmparg), res1, tmparrS(i), tmperr(i))
!			write(*,"(I2,*(E17.9))") i, res1, res2, tmparrS(i)
		end do

		write(*,*) ""
		tmparrA = (/5.49883d0, 1.16187d0, 1.16187d0/)
		tmperr3 = (/2d-3, 3d-3, 3d-3/)
		tmperr4 = (/1d-5, 1d-5, 1d-5/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			ifail=0
			itrans=0
!			call D01GCF(n,coll_nue_4_ann_int_re, region, npts, vk, nrand,itrans,res1,ERRr1,ifail, collArgs)
!			write(tmparg,"('test coll ann 4 b - ',2I1)") i, i
!			call assert_double_rel_verb(trim(tmparg), res1, tmparrA(i), tmperr4(i))
			res2 = integrate_coll_int_3(coll_nue_3_ann_int_re, collArgs)
			write(tmparg,"('test coll ann 3 b - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrA(i), tmperr3(i))
!			write(*,"(I2,*(E17.9))") i, res1, res2, tmparrA(i)
		end do

		x=0.05d0
		iy1=7 !1.22151515151515
		z=1.06d0
		write(*,*)
		write(*,"(a)") "Collision integrands im part and off diagonal (36 tests)"
		collArgs%x = x
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%dme2 = dme2

		do iy=1, Ny
			nuDensMatVecFD(iy)%re = 0.d0
			nuDensMatVecFD(iy)%im = 0.d0
		end do
		nuDensMatVecFD(7)%re(1, :) = (/1.d0*fermiDirac(y_arr(7)), 0.1d0, 0.01d0/)
		nuDensMatVecFD(7)%re(2, :) = (/0.1d0, 1.05d0*fermiDirac(y_arr(7)), 0.03d0/)
		nuDensMatVecFD(7)%re(3, :) = (/0.01d0, 0.03d0, 1.1d0*fermiDirac(y_arr(7))/)
		nuDensMatVecFD(7)%im(1, :) = (/0.d0, -0.1d0, 0.03d0/)
		nuDensMatVecFD(7)%im(2, :) = (/0.1d0, 0.d0, -0.02d0/)
		nuDensMatVecFD(7)%im(3, :) = (/-0.03d0, 0.02d0, 0.d0/)
		nuDensMatVecFD(21)%re(1, :) = (/1.1d0*fermiDirac(y_arr(21)), 0.d0, 0.12d0/)
		nuDensMatVecFD(21)%re(2, :) = (/0.d0, 1.4d0*fermiDirac(y_arr(21)), 0.04d0/)
		nuDensMatVecFD(21)%re(3, :) = (/0.12d0, 0.04d0, 1.7d0*fermiDirac(y_arr(21))/)
		nuDensMatVecFD(21)%im(1, :) = (/0.d0, -0.05d0, 0.d0/)
		nuDensMatVecFD(21)%im(2, :) = (/0.05d0, 0.0d0, 0.01d0/)
		nuDensMatVecFD(21)%im(3, :) = (/0.d0, -0.01d0, 0.d0/)

		tmpmatA(1,:) = (/0.0347964, -0.371928, -1.7873/)
		tmpmatA(2,:) = (/-0.371928, 0.103022, 0.539575/)
		tmpmatA(3,:) = (/-1.7873, 0.539575, 0.0851156/)
		tmpmatB(1,:) = (/0., 1.12312, -0.110645/)
		tmpmatB(2,:) = (/-1.12312, 0., 0.260258/)
		tmpmatB(3,:) = (/0.110645, -0.260258, 0./)
		tmparrA(:) = (/1d-5, 1d-5, 1d-5/)
		tmparrS(:) = (/1d-5, 1d-5, 1d-5/)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('coll sc integrand full ',2I1)") i,j
				collArgs%ix1 = i
				collArgs%ix2 = j
				res1 = coll_nue_3_sc_int_re(21, 5.2d0, collArgs)
				res2 = coll_nue_3_sc_int_im(21, 5.2d0, collArgs)
				call assert_double_rel(trim(tmparg)//"re", res1, tmpmatA(i,j), tmparrA(i))
				if (i.eq.j) then
					call assert_double(trim(tmparg)//"im", res2, tmpmatB(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", res2, tmpmatB(i,j), tmparrS(i))
				end if
			end do
		end do

		tmpmatA(1,:) = (/0.374939, -0.146968, 1.06745/)
		tmpmatA(2,:) = (/-0.146968, 0.16068, -0.0485291/)
		tmpmatA(3,:) = (/1.06745, -0.0485291, -0.028232/)
		tmpmatB(1,:) = (/0., -0.230122, -0.0597533/)
		tmpmatB(2,:) = (/0.230122, 0., 0.204294/)
		tmpmatB(3,:) = (/0.0597533, -0.204294, 0./)
		tmparrA(:) = (/1d-5, 1d-5, 1d-5/)
		tmparrS(:) = (/1d-5, 1d-5, 1d-5/)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('coll ann integrand full ',2I1)") i,j
				collArgs%ix1 = i
				collArgs%ix2 = j
				res1 = coll_nue_3_ann_int_re(21, 3.d0, collArgs)
				res2 = coll_nue_3_ann_int_im(21, 3.d0, collArgs)
				call assert_double_rel(trim(tmparg)//"re", res1, tmpmatA(i,j), tmparrA(i))
				if (i.eq.j) then
					call assert_double(trim(tmparg)//"im", res2, tmpmatB(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", res2, tmpmatB(i,j), tmparrS(i))
				end if
			end do
		end do
!		call criticalError("stop")
		deallocate(ndmv_re)
	end subroutine do_tests_coll_int

	pure real(dl) function fakecollint1(a, b, o)
		use variables
		integer, intent(in) :: a
		real(dl), intent(in) :: b
		type(coll_args), intent(in) :: o
		fakecollint1=1.d0
	end function
	pure real(dl) function fakecollint0(a, b, o)
		use variables
		integer, intent(in) :: a
		real(dl), intent(in) :: b
		type(coll_args), intent(in) :: o
		fakecollint0=0.d0
	end function
	pure real(dl) function fakecollinty(a, b, o)
		use variables
		integer, intent(in) :: a
		real(dl), intent(in) :: b
		type(coll_args), intent(in) :: o
		fakecollinty=1.d4*b**2
	end function

	subroutine do_test_collision_terms
		real(dl) :: x,z,dme2
		integer :: i, j, iy1, iy
		real(dl) :: y,res1,res2
		type(coll_args) :: collArgs
		real(dl), dimension(3) :: tmparrS, tmparrA
		real(dl), dimension(3, 3) :: tmpmatA, tmpmatB
		character(len=300) :: tmparg
		type(cmplxMatNN) :: cts

		call allocateCmplxMat(cts)

		x = 0.75d0
		iy1 = 7 !1.22151515151515
		z = 1.186d0
		dme2 = 0.1d0
		write(*,*)
		write(*,"(a)") "Collision_terms (36 tests)"
		collArgs%ix1 = 1
		collArgs%ix1 = 1
		collArgs%x = x
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%y2 = 0.d0
		collArgs%y3 = 0.d0
		collArgs%y4 = 0.d0
		collArgs%dme2 = dme2

		!fake function
		tmpmatA(1,:) = (/348548, 348548, 348548/)
		tmpmatA(2,:) = (/348548, 348548, 348548/)
		tmpmatA(3,:) = (/348548, 348548, 348548/)
		tmpmatB(1,:) = (/0., 0.26128, 0.26128/)
		tmpmatB(2,:) = (/-0.26128, 0., 0.26128/)
		tmpmatB(3,:) = (/-0.26128, -0.26128, 0./)
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollinty, fakecollint1)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('collision_terms f ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				if (abs(tmpmatB(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), tmparrS(i))
				end if
			end do
		end do

		x=0.05d0
		iy1 = 7 !1.22151515151515
!		z=1.06d0
		dme2=0.1d0
!		iy1 = 67 !13.3366666666667
		z = 1.3d0
		collArgs%ix1 = 1
		collArgs%ix1 = 1
		collArgs%x = x
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%y2 = 0.d0
		collArgs%y3 = 0.d0
		collArgs%y4 = 0.d0
		collArgs%dme2 = dme2
		do iy=1, Ny
			nuDensMatVecFD(iy)%re(:,:) = 0.d0
			nuDensMatVecFD(iy)%im(:,:) = 0.d0
			do i=1, flavorNumber
				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
			end do
		end do

		!real params
!		do iy=1, Ny
!			y = y_arr(iy)
!			nuDensMatVecFD(iy)%re(1,:) = (/1.d0*fermiDirac(y), 0.d0, 0.d0/)
!			nuDensMatVecFD(iy)%re(2,:) = (/0.d0, 1.d0*fermiDirac(y), 0.d0/)
!			nuDensMatVecFD(iy)%re(3,:) = (/0.d0, 0.d0, 1.d0*fermiDirac(y)/)
!!			nuDensMatVecFD(iy)%im(1,:) = (/0.d0, -0.001d0, 0.003d0/)
!!			nuDensMatVecFD(iy)%im(2,:) = (/0.001d0, 0.d0, -0.002d0/)
!!			nuDensMatVecFD(iy)%im(3,:) = (/-0.003d0, 0.002d0, 0.d0/)
!			nuDensMatVecFD(iy)%im(:,:) = 0.d0
!		end do

		collArgs%ix1=1
		collArgs%ix2=1
!		print *,collArgs
!		write(*,multidblfmt) integrate_coll_int_3(coll_nue_3_sc_int_re, collArgs),integrate_coll_int_3(coll_nue_3_ann_int_re, collArgs)
		tmpmatA(1,:) = (/-149.474, 0., 0./)
		tmpmatA(2,:) = (/0., -38.1733, 0./)
		tmpmatA(3,:) = (/0., 0., -38.1733/)
		tmpmatB(1,:) = (/0., 0., 0./)
		tmpmatB(2,:) = (/0., 0., 0./)
		tmpmatB(3,:) = (/0., 0., 0./)
		tmparrA(:) = (/0.08d0, 0.05d0, 0.05d0/)
		tmparrS(:) = (/0.00001d0, 0.00001d0, 0.00001d0/)
		cts = get_collision_terms(collArgs, coll_nue_3_int_re, coll_nue_3_int_im)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
!		write(*,multidblfmt)cts
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('collision_terms r ',2I1)") i,j
				if (abs(tmpmatA(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				end if
				if (abs(tmpmatB(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), tmparrS(i))
				end if
			end do
		end do
	end subroutine do_test_collision_terms

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

		allocate(ndmv_re(Ny))

		do i=1, flavorNumber
			do ix=1, Ny
				nuDensMatVecFD(ix)%re(i, i) = 1.d0*i * fermiDirac(y_arr(ix)/z)
				ndmv_re(ix) = nuDensMatVecFD(ix)%re(i, i)
			end do
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
		end do

		collArgs%x = x
		collArgs%z = z
		collArgs%y1 = y1
		collArgs%dme2 = 0.d0!dme2_electron(x, 0.d0, z)
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		collArgs%iy = 12

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
		
		call tic(timer1)
		do ix=1, 10*Npt
			call random_number(x)
			call random_number(z)
			collArgs%x = (x_fin-x_in)*x + x_in
			collArgs%z = 0.4d0*z + z_in
			ifail=0
			itrans=0
			res1 = integrate_coll_int_3(coll_nue_3_ann_int_re, collArgs)
			res2 = integrate_coll_int_3(coll_nue_3_sc_int_re, collArgs)
			res2 = integrate_coll_int_3(coll_nue_3_int_re, collArgs)
		end do
		call toc(timer1, "<reset>")
		call tic(timer1)
		do ix=1, 10*Npt
			call random_number(x)
			call random_number(z)
			collArgs%x = (x_fin-x_in)*x + x_in
			collArgs%z = 0.4d0*z + z_in
			ifail=0
			itrans=0
			res1 = integrate_coll_int_3(coll_nue_3_ann_int_re, collArgs)
			res2 = integrate_coll_int_3(coll_nue_3_sc_int_re, collArgs)
		end do
		call toc(timer1, "<sep ann, sc>")
		call tic(timer1)
		do ix=1, 10*Npt
			call random_number(x)
			call random_number(z)
			collArgs%x = (x_fin-x_in)*x + x_in
			collArgs%z = 0.4d0*z + z_in
			ifail=0
			itrans=0
			res2 = integrate_coll_int_3(coll_nue_3_int_re, collArgs)
		end do
		call toc(timer1, "<sum ann+sc>")
		
		deallocate(ndmv_re)
		call addToLog("[interactions] ...done!")
	end subroutine test_speed_coll_int

	subroutine do_test_drho_dx
		real(dl) :: x,z,fd, dme2, sqrtraddens
		type(cmplxMatNN) :: res, outp
		character(len=300) :: tmparg
		integer :: i, j, iy
		real(dl), dimension(:,:,:), allocatable :: GLR_vectmp

		allocate(GLR_vectmp(2, flavorNumber, flavorNumber))
		GLR_vectmp = GLR_vec
		GLR_vec = 0.d0
		write(*,*) ""
		write(*,"(a)") "d rho/d x [without collision_terms or with fake ones] (72 tests)"
		call allocateCmplxMat(res)
		x = 0.06d0
		z = 1.23d0
		dme2 = 0.1d0
		call updateLeptonDensities(x,z)
		iy = 12 !2.231111111111111
		do i=1, Ny
			fd = fermiDirac(y_arr(i))
			nuDensMatVecFD(i)%re(1, :) = (/fd, 0.01d0, 0.04d0/)
			nuDensMatVecFD(i)%re(2, :) = (/0.01d0, fd, -0.3d0/)
			nuDensMatVecFD(i)%re(3, :) = (/0.04d0, -0.3d0, fd/)
			nuDensMatVecFD(i)%im(1, :) = (/0.0, -0.02, -0.1/)
			nuDensMatVecFD(i)%im(2, :) = (/0.02, 0.0, 0.1/)
			nuDensMatVecFD(i)%im(3, :) = (/0.1, -0.1, 0.0/)
		end do
		sqrtraddens = sqrt(radDensity(x,z))

		fd = fermiDirac(y_arr(iy))
		res%re(1,:) = (/605.541d0/fd, 15531.6d0, 73346.3d0/)
		res%re(2,:) = (/15531.6d0, -2653.78d0/fd, -604.533d0/)
		res%re(3,:) = (/73346.3d0, -604.533d0, 2048.24d0/fd/)
		res%im(1,:) = (/0., 8489.27, 30142.1/)
		res%im(2,:) = (/-8489.27, 0., -865.415/)
		res%im(3,:) = (/-30142.1, 865.415, 0./)
		call drhoy_dx_fullMat(outp, x, z, iy, dme2, sqrtraddens, fakecollint0, fakecollint0)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('drho/dx a ',2I1)") i,j
				if (abs(res%re(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-4)
				end if
				if (abs(res%im(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 2d-4)
				end if
			end do
		end do

		fd = fermiDirac(y_arr(iy))
		res%re(1,:) = (/1394.47d0/fd, 16320.5d0, 74135.2d0/)
		res%re(2,:) = (/16320.5d0, -1864.85d0/fd, 184.352d0/)
		res%re(3,:) = (/74135.2d0, 184.352d0, 2837.16d0/fd/)
		res%im(1,:) = (/0., 9278.19, 30931./)
		res%im(2,:) = (/-9278.19, 0., -76.52251/)
		res%im(3,:) = (/-30931., 76.52251, 0./)
		call drhoy_dx_fullMat(outp,x,z,iy, dme2, sqrtraddens, fakecollint1, fakecollint1)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('drho/dx b ',2I1)") i,j
				if (abs(res%re(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-4)
				end if
				if (abs(res%im(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-4)
				end if
			end do
		end do

		x = 1.76d0
		z = 1.31d0
		call updateLeptonDensities(x,z)
		iy = 34 !6.67333333333333
		do i=1, Ny
			fd = fermiDirac(y_arr(i))
			nuDensMatVecFD(i)%re(1, :) = (/fd, 0.01d0, 0.04d0/)
			nuDensMatVecFD(i)%re(2, :) = (/0.01d0, fd, -0.3d0/)
			nuDensMatVecFD(i)%re(3, :) = (/0.04d0, -0.3d0, fd/)
			nuDensMatVecFD(i)%im(1, :) = (/0.0, -0.02, -0.1/)
			nuDensMatVecFD(i)%im(2, :) = (/0.02, 0.0, 0.1/)
			nuDensMatVecFD(i)%im(3, :) = (/0.1, -0.1, 0.0/)
		end do
		sqrtraddens = sqrt(radDensity(x,z))

		fd = fermiDirac(y_arr(iy))
		res%re(1,:) = (/164543.d0/fd, 346799.d0, 481933.d0/)
		res%re(2,:) = (/346799.d0, -721106.d0/fd, -83688.4d0/)
		res%re(3,:) = (/481933.d0, -83688.4d0, 556564.d0/fd/)
		res%im(1,:) = (/0., 369988., 411122./)
		res%im(2,:) = (/-369988., 0., 6601.2/)
		res%im(3,:) = (/-411122., -6601.2, 0./)
		call drhoy_dx_fullMat(outp,x,z,iy, dme2, sqrtraddens, fakecollint0, fakecollint0)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('drho/dx c ',2I1)") i,j
				if (abs(res%re(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-4)
				end if
				if (abs(res%im(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-4)
				end if
			end do
		end do

		res%re(1,:) = (/164693.d0/fd, 346950.d0, 482083.d0/)
		res%re(2,:) = (/346950.d0, -720956.d0/fd, -83538.3d0/)
		res%re(3,:) = (/482083.d0, -83538.3d0, 556714.d0/fd/)
		res%im(1,:) = (/0., 370138., 411272./)
		res%im(2,:) = (/-370138., 0., 6751.29/)
		res%im(3,:) = (/-411272., -6751.29, 0./)
		call drhoy_dx_fullMat(outp,x,z,iy, dme2, sqrtraddens, fakecollinty, fakecollinty)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('drho/dx d ',2I1)") i,j
				if (abs(res%re(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-4)
				end if
				if (abs(res%im(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-4)
				end if
			end do
		end do
		GLR_vec = GLR_vectmp
		deallocate(GLR_vectmp)
	end subroutine do_test_drho_dx

	subroutine do_test_damping_factors
		real(dl) :: x,z,dme2
		integer :: i, j, ix, iy1, iy
		real(dl) :: y,res1,res2
		type(coll_args) :: collArgs
		real(dl), dimension(3) :: tmparrS, tmparrA
		real(dl), dimension(3, 3) :: tmpmatA, tmpmatB
		character(len=300) :: tmparg
		type(cmplxMatNN) :: cts

		call allocateCmplxMat(cts)

		collision_offdiag = 2

		x = 0.75d0
		iy1 = 7 !1.22151515151515
		z = 1.186d0
		dme2 = 0.1d0
		write(*,*)
		write(*,"(a)") "Damping terms (28 tests)"
		collArgs%x = x
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%y2 = 0.d0
		collArgs%y3 = 0.d0
		collArgs%y4 = 0.d0
		collArgs%dme2 = dme2
		do iy=1, Ny
			y = y_arr(iy)
			nuDensMatVecFD(iy)%re(1,:) = (/1.d0*fermiDirac(y), 10.d0, 33.d0/)
			nuDensMatVecFD(iy)%re(2,:) = (/10.d0, 1.d0*fermiDirac(y), 45.d0/)
			nuDensMatVecFD(iy)%re(3,:) = (/33.d0, 45.d0, 1.d0*fermiDirac(y)/)
			nuDensMatVecFD(iy)%im(1,:) = (/0.d0, -0.001d0, 0.003d0/)
			nuDensMatVecFD(iy)%im(2,:) = (/0.001d0, 0.d0, -0.002d0/)
			nuDensMatVecFD(iy)%im(3,:) = (/-0.003d0, 0.002d0, 0.d0/)
		end do

		!2+1
		sterile(3) = .true.
		call setDampingFactorCoeffs
		tmpmatA(1,:) = (/0., -0.928646, -6.56674/)
		tmpmatA(2,:) = (/0., 0., -7.45108/)
		tmpmatA(3,:) = (/0., 0., 0./)
		tmpmatB(1,:) = (/0., 0.0000928646, -0.000596976/)
		tmpmatB(2,:) = (/0., 0., 0.000331159/)
		tmpmatB(3,:) = (/0., 0., 0./)
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollinty, fakecollint1)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=i+1, flavorNumber
				write(tmparg,"('damping term 2+1 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				if (abs(tmpmatB(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), tmparrS(i))
				end if
			end do
		end do

!		!3+0
		sterile(3) = .false.
		call setDampingFactorCoeffs
		tmpmatA(1,:) = (/0., -0.928646, -3.06453/)
		tmpmatA(2,:) = (/0., 0., -1.76139/)
		tmpmatA(3,:) = (/0., 0., 0./)
		tmpmatB(1,:) = (/0., 0.0000928646, -0.000278594/)
		tmpmatB(2,:) = (/0., 0., 0.0000782838/)
		tmpmatB(3,:) = (/0., 0., 0./)
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollinty, fakecollint1)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=i+1, flavorNumber
				write(tmparg,"('damping term 3+0 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				if (abs(tmpmatB(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), tmparrS(i))
				end if
			end do
		end do

		flavorNumber = 2
		flavNumSqu = flavorNumber**2
		call deallocateStuff
		call allocateStuff
		!1+1
		sterile(2) = .true.
		call setDampingFactorCoeffs
		tmpmatA(1,:) = (/0., -1.98992, 0./)
		tmpmatA(2,:) = (/0., 0., 0./)
		tmpmatA(3,:) = (/0., 0., 0./)
		tmpmatB(1,:) = (/0., 0.000198992, 0./)
		tmpmatB(2,:) = (/0., 0., 0./)
		tmpmatB(3,:) = (/0., 0., 0./)
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollinty, fakecollint1)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=i+1, flavorNumber
				write(tmparg,"('damping term 1+1 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				if (abs(tmpmatB(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), tmparrS(i))
				end if
			end do
		end do

		!2+0
		sterile(2) = .false.
		call setDampingFactorCoeffs
		tmpmatA(1,:) = (/0., -0.928646, 0./)
		tmpmatA(2,:) = (/0., 0., 0./)
		tmpmatA(3,:) = (/0., 0., 0./)
		tmpmatB(1,:) = (/0., 0.0000928646, 0./)
		tmpmatB(2,:) = (/0., 0., 0./)
		tmpmatB(3,:) = (/0., 0., 0./)
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollinty, fakecollint1)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=i+1, flavorNumber
				write(tmparg,"('damping term 2+0 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				if (abs(tmpmatB(i,j)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), tmparrS(i))
				end if
			end do
		end do

		flavorNumber = 4
		flavNumSqu = flavorNumber**2
		call deallocateStuff
		call allocateStuff
		call setDampingFactorCoeffs
		do ix=1, Ny
			deallocate(nuDensMatVecFD(ix)%re, nuDensMatVecFD(ix)%im)
			allocate(nuDensMatVecFD(ix)%re(flavorNumber,flavorNumber), nuDensMatVecFD(ix)%im(flavorNumber,flavorNumber))
		end do
		do iy=1, Ny
			y = y_arr(iy)
			nuDensMatVecFD(iy)%re(1,:) = (/1.d0*fermiDirac(y), 10.d0, 33.d0, 12.d0/)
			nuDensMatVecFD(iy)%re(2,:) = (/10.d0, 1.d0*fermiDirac(y), 45.d0, 8.d0/)
			nuDensMatVecFD(iy)%re(3,:) = (/33.d0, 45.d0, 1.d0*fermiDirac(y), 23.d0/)
			nuDensMatVecFD(iy)%re(4,:) = (/12.d0, 8.d0, 23.d0, 1.d0*fermiDirac(y)/)
			nuDensMatVecFD(iy)%im(1,:) = (/0.d0, -0.001d0, 0.003d0, 0.1d0/)
			nuDensMatVecFD(iy)%im(2,:) = (/0.001d0, 0.d0, -0.002d0, -0.2d0/)
			nuDensMatVecFD(iy)%im(3,:) = (/-0.003d0, 0.002d0, 0.d0, 0.5d0/)
			nuDensMatVecFD(iy)%im(4,:) = (/-0.1d0, 0.2d0, -0.5d0, 0.d0/)
		end do
		!3+1
		sterile = .false.
		sterile(4) = .true.
		call setDampingFactorCoeffs
		tmpmatA(1,:) = (/-0.928646, -3.06453, -2.38791/)
		tmpmatA(2,:) = (/0., -1.76139, -1.32464/)
		tmpmatA(3,:) = (/0., 0., -3.80833/)
		tmpmatB(1,:) = (/0.0000928646, -0.000278594, -0.0198992/)
		tmpmatB(2,:) = (/0., 0.0000782838, 0.0331159/)
		tmpmatB(3,:) = (/0., 0., -0.0827898/)
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollinty, fakecollint1)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=i+1, flavorNumber
				write(tmparg,"('damping term 3+1 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i, j-1), tmparrA(i))
				if (abs(tmpmatB(i,j-1)).lt.1d-7) then
					call assert_double(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i, j-1), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i, j-1), tmparrS(i))
				end if
			end do
		end do

		!restore initial settings
		flavorNumber = 3
		flavNumSqu = flavorNumber**2
		call deallocateStuff
		call allocateStuff
		call setDampingFactorCoeffs
		do ix=1, Ny
			deallocate(nuDensMatVecFD(ix)%re, nuDensMatVecFD(ix)%im)
			allocate(nuDensMatVecFD(ix)%re(flavorNumber,flavorNumber), nuDensMatVecFD(ix)%im(flavorNumber,flavorNumber))
		end do
		collision_offdiag = 1
	end subroutine do_test_damping_factors

	subroutine do_test_zin
		write(*,*)
		write(*,"(a)") "z_in solver (4 tests)"
		dme2_temperature_corr = .false.
		x_in=0.05d0
		z_in=0.d0
		call zin_solver
		call assert_double("z_in test 1", z_in-1.d0, 0.575d-04, 1d-6)
		x_in=1d-3
		call zin_solver
		call assert_double("z_in test 2", z_in-1.d0, 2.3d-8, 1d-9)
		dme2_temperature_corr = .true.
		x_in=0.05d0
		z_in=0.d0
		call zin_solver
		call assert_double("z_in test 1", z_in-1.d0, 0.56d-04, 1d-6)
		x_in=1d-3
		call zin_solver
		call assert_double("z_in test 2", z_in-1.d0, 7.7d-8, 1d-9)
	end subroutine do_test_zin

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
