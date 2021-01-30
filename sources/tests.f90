program tests
	use precision
	use variables
	use fpConfig
	use fpErrors
	use fpEquations
	use fpCosmology
	use fpStuff
	use sgTestUtils
	use ftqed
	implicit none

	integer, parameter :: fu = 89237, fv = 89238, fw = 89239
	character(len=1), dimension(2), parameter :: chLR=(/'L','R'/)

	call openLogFile
	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "Initializations"
	call do_test_initialization
	call init_fermions
	call allocate_interpNuDens

	call do_basic_tests
	call do_test_NC_integrals
	call do_test_commutator
	call do_test_JKG
	call do_test_dme2
	call do_test_cosmology
	call do_test_nu_matrices
	call do_test_dzodx
	call do_test_dme2
	call do_test_Di
	call do_test_Pi_ij
	call do_f_ann_sc_re_tests_eq
	call do_f_ann_sc_re_tests_full
	call do_test_F_nu
	call do_test_interp_nudens
	call do_test_collint_nunu
	call do_test_coll_int
	call do_test_drho_dx
	call do_test_collision_terms
	call do_test_damping_factors
	call do_test_zin
	call do_test_damping_bennett
	call do_test_GL
	call do_test_matterPotential
	call do_test_diagonalization

	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "all tests were successfull!"

!	write(*,*) ""
!	write(*,"(a)") "now doing some timing tests..."
!	call do_timing_tests()

	call closeLogFile

	contains

	subroutine do_test_initialization
		real(dl), dimension(:), allocatable :: fake
		integer :: ix

		checkpoint = .true.
		outputFolder = ""
		maxiter = 100
		toler_jkyg = 1.d-7
		dlsoda_atol_z = 1.d-6
		dlsoda_atol_d = 1.d-6
		dlsoda_atol_o = 1.d-6
		dlsoda_rtol = 1.d-6
		Nx = 100
		Ny = 100
		interp_nx = interp_nx0
		interp_nz = interp_nz0
		interp_nxz = interp_nxz0
		interp_zmin = interp_zmin0
		interp_zmax = interp_zmax0
		allocate(x_arr(Nx), y_arr(Ny))
		x_in    = 0.001d0
		x_fin   = 35.d0
		logx_in  = log10(x_in)
		logx_fin = log10(x_fin)
		x_arr = logspace(logx_in, logx_fin, Nx)
		y_min = 0.01d0
		y_max = 20.0d0
		use_gauss_laguerre = .false.
		y_arr = linspace(y_min, y_max, Ny)
		call finish_y_arrays
		call get_GLq_vectors(N_opt_y, opt_y, opt_y_w, fake, .true., 2, opt_y_cut)
		call get_GLq_vectors(N_opt_xoz, opt_xoz, opt_xoz_w, fake, .true., 2, opt_xoz_cut)
		collint_damping_type = 2
		collint_diagonal_zero = .false.
		collint_offdiag_damping = .false.
		damping_read_zero = .false.
		ftqed_temperature_corr = .true.
		ftqed_log_term = .false.
		ftqed_ord3 = .false.
		ftqed_e_mth_leptondens = .false.
		flavorNumber = 3
		flavNumSqu = flavorNumber**2
		call allocateStuff
		do ix=1, flavorNumber
			nuFactor(ix) = 1.d0
			sterile(ix) = .false.
		end do
		tot_factor_active_nu = 3.0
		tot_factor_nu = 0.d0
		mixingAngles(1,2) = 0.5840d0
		massSplittings(1) = 0.d0
		massSplittings(2) = 7.53d-05
		if (flavorNumber .gt. 2) then
			mixingAngles(1,3) = 0.1485d0
			mixingAngles(2,3) = 0.7954d0
			massSplittings(3) = 0.0025283d0
		end if
		z_in=1.0000575
		collint_d_no_nue = .false.
		collint_d_no_nunu = .false.
		collint_od_no_nue = .false.
		collint_od_no_nunu = .false.
		save_fd = .true.
		save_energy_entropy_evolution = .true.
		save_Neff = .true.
		save_nuDens_evolution = .true.
		save_number_evolution = .true.
		save_w_evolution = .true.
		save_z_evolution = .true.
		call setMixingMatrix()
		call setMassMatrix()
		call init_matrices
		allocate(interp_xvec(interp_nx), interp_yvec(interp_ny), interp_zvec(interp_nz), interp_xozvec(interp_nxz))
		interp_xvec = logspace(logx_in, logx_fin, interp_nx)
		interp_yvec = logspace(log10(y_min), log10(y_max), interp_ny)
		interp_zvec = linspace(interp_zmin, interp_zmax, interp_nz)
		low_xoz = x_in/interp_zmax
		interp_xozvec = logspace(log10(low_xoz), logx_fin, interp_nxz)
		call get_GLq_vectors(N_opt_y, opt_y, opt_y_w, fake, .true., 2, opt_y_cut)
		call get_GLq_vectors(N_opt_xoz, opt_xoz, opt_xoz_w, fake, .true., 2, opt_xoz_cut)
	end subroutine do_test_initialization

	subroutine do_basic_tests
		real(dl), dimension(:), allocatable :: tarr
		integer :: N, ix
		character(len=30) :: tmparg

		call printTestBlockName("basic tests")

		collint_damping_type = 0
		collint_offdiag_damping = .false.
		call assert_logical("has_offdiagonal 1", has_offdiagonal(), .true.)
		collint_damping_type = 2
		collint_offdiag_damping = .true.
		call assert_logical("has_offdiagonal 2", has_offdiagonal(), .true.)
		collint_damping_type = 0
		call assert_logical("has_offdiagonal 3", has_offdiagonal(), .false.)
		collint_damping_type = 2
		collint_offdiag_damping = .false.

		N = 12
		allocate(tarr(N))
		tarr = linspace(0.0d0, 11.d0, N)
		do ix=1, N
			write(tmparg, "('tarr lin ',I2)") ix
			call assert_double(tmparg, tarr(ix), (ix-1)*1.d0, 1d-7)
		end do
		tarr = logspace(-5d0, 6.d0, N)
		do ix=1, N
			write(tmparg, "('tarr log ',I2)") ix
			call assert_double_rel(tmparg, tarr(ix), 10.d0**(ix-6), 1d-7)
		end do
		tarr = geomspace(1.d-5, 1.d6, N)
		do ix=1, N
			write(tmparg, "('tarr geom ',I2)") ix
			call assert_double_rel(tmparg, tarr(ix), 10.d0**(ix-6), 1d-7)
		end do
		tarr = loglinspace(0.01d0, 1.d0, 10.d0, N, 3)
		call assert_double_rel("tarr linlog 1", tarr(1), 0.01d0, 1d-7)
		call assert_double_rel("tarr linlog 2", tarr(2), 0.1d0, 1d-7)
		do ix=3, N
			write(tmparg, "('tarr linlog ',I2)") ix
			call assert_double_rel(tmparg, tarr(ix), (ix-2)*1.d0, 1d-7)
		end do

		call assert_double_rel("y_arr linlog 1", y_arr(1), 0.01d0, 1d-7)
		call assert_double_rel("y_arr linlog 2", y_arr(2), 0.21191919191919d0, 1d-7)

		call assert_double_rel("FD 1", fermiDirac(0.d0), 0.5d0, 1d-7)
		call assert_double_rel("FD 2", fermiDirac(5.d0), 0.0066928509242848554d0, 1d-7)
		call assert_double_rel("FD 3", fermiDirac(20.d0), 2.0611536181902037d-09, 1d-7)
		call assert_double_rel("f_eq saved A", feq_arr(1), fermiDirac(y_arr(1)), 1d-3)
		call assert_double_rel("f_eq saved B", feq_arr(Ny), fermiDirac(y_arr(Ny)), 1d-3)
		call printTotalTests
		call resetTestCounter
	end subroutine do_basic_tests

	subroutine do_test_NC_integrals
		integer :: ia, ib
		real(dl), dimension(:), allocatable :: fy1_arr
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(fy1_arr(Ny))
		allocate(fy2_arr(Ny, Ny))

		call printTestBlockName("integral_linearized")

		do ia=1, Ny
			fy1_arr(ia) = 0.35d0*y_arr(ia) + 11.41d0
		end do
		call assert_double_rel("intLin1D 1", integral_NC_1d(Ny, dy_arr, fy1_arr), 298.08588d0, 1d-7)
		do ia=1, Ny
			fy1_arr(ia) = 0.35d0/(exp(y_arr(ia))+1.d0) + 0.2d0
		end do
		call assert_double_rel("intLin1D 2", integral_NC_1d(Ny, dy_arr, fy1_arr), 4.23886d0, 1d-4)

		do ia=1, Ny
			do ib=1, Ny
				fy2_arr(ia,ib) = 0.35d0*y_arr(ia) - 1.96d0*y_arr(ib) + 11.41d0
			end do
		end do
		call assert_double_rel("intLin2D 1", integral_NC_2d(Ny, Ny, dy_arr, dy_arr, fy2_arr), -1877.341249d0, 1d-7)
		do ia=1, Ny
			do ib=1, Ny
				fy2_arr(ia,ib) = 0.35d0/(exp(y_arr(ia))+1.d0) + 1.96d0/(exp(y_arr(ib))+1.d0) + 0.2d0
			end do
		end do
		call assert_double_rel("intLin2D 2", integral_NC_2d(Ny, Ny, dy_arr, dy_arr, fy2_arr), 111.967d0, 0.003d0)

		call assert_double_rel("lnP test 1", &
			integrate_ftqed_ln(deltaP_ln_integrand, 1.d0, 1.d0), 0.542467d0, 1d-3)
		call assert_double_rel("lnP test 2", &
			integrate_ftqed_ln(deltaP_ln_integrand, 0.01d0, 1.d0), 2.01577d0, 5d-3)
		call assert_double_rel("lnP test 3", &
			integrate_ftqed_ln(deltaP_ln_integrand, 0.2d0, 1.3d0), 3.05627d0, 1d-3)

		call assert_double_rel("lnG1 test 1", &
			integrate_ftqed_ln(G1_ln_integrand, 0.1d0, 10.d0), -403.3d0, 1d-3)
		call assert_double_rel("lnG1 test 2", &
			integrate_ftqed_ln(G1_ln_integrand, 1.5d0, 1.5d0), -1.998d0, 3d-3)
		call assert_double_rel("lnG1 test 3", &
			integrate_ftqed_ln(G1_ln_integrand, 0.5d0, 0.9d0), -1.901d0, 5d-4)

		call assert_double_rel("lnG2 test 1", &
			integrate_ftqed_ln(G2_ln_integrand, 0.1d0, 10.d0), 2.01577d0, 5d-2)
		call assert_double_rel("lnG2 test 2", &
			integrate_ftqed_ln(G2_ln_integrand, 1.5d0, 1.5d0), 2.073d0, 5d-2)
		call assert_double_rel("lnG2 test 3", &
			integrate_ftqed_ln(G2_ln_integrand, 0.5d0, 0.9d0), 2.204d0, 1d-2)

		deallocate(fy1_arr)
		deallocate(fy2_arr)
		call printTotalTests
		call resetTestCounter
	end subroutine do_test_NC_integrals

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

		call printTestBlockName("Commutator & anticommutator")
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
		call printTotalTests
		call resetTestCounter
	end subroutine do_test_commutator

	subroutine do_test_nu_matrices
		real(dl), dimension(3) :: er
		real(dl), dimension(3,3) :: m, nr, ni
		real(dl), dimension(:,:), allocatable :: ide
		integer :: i,j, iy
		character(len=300) :: tmparg

		call printTestBlockName("Neutrino and lepton matrices")

		open(unit=fu, file="test_outputs/mixmat.dat", status="old")
		do i=1, 3
			read (fu, *) m(i,:)
		end do
		close(fu)
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

		open(unit=fu, file="test_outputs/massmat.dat", status="old")
		do i=1, 3
			read (fu, *) m(i,:)
		end do
		close(fu)
		do i=1,3
			do j=1,3
				write(tmparg,"('mass matrix ',2I1)") i,j
				call assert_double(trim(tmparg), nuMassesMat(i,j), m(i,j), 1d-7)
			end do
		end do

		leptonDensities=0.d0
		nuDensities%re=0.d0
		nuDensities%im=0.d0
		do iy=1, Ny
			nuDensMatVecFD(iy)%re = 0.d0
			nuDensMatVecFD(iy)%im = 0.d0
		end do
		!A
		call updateMatterDensities(0.076d0, 1.32d0)
		open(unit=fu, file="test_outputs/leptmatA.dat", status="old")
		do i=1, 3
			read (fu, *) m(i,:)
		end do
		close(fu)
		er = (/5d-5,1d-3,0.d0/)
		do i=1,3
			do j=1,3
				write(tmparg,"('lepton matrix A ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg), leptonDensities(i,j), m(i,j), 1d-16, er(i))
				write(tmparg,"('nu density matrix re A ',2I1)") i,j
				call assert_double(trim(tmparg), nuDensities%re(i,j), 0.d0, 1d-7)
				write(tmparg,"('nu density matrix im A ',2I1)") i,j
				call assert_double(trim(tmparg), nuDensities%im(i,j), 0.d0, 1d-7)
			end do
		end do
		!B
		do iy=1, Ny
			nuDensMatVecFD(iy)%re = fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%im = 0.01d0 * fermiDirac(y_arr(iy))
		end do
		call updateMatterDensities(0.0176d0, 1.d0)
		write(*,*)""
		open(unit=fu, file="test_outputs/leptmatB.dat", status="old")
		do i=1, 3
			read (fu, *) m(i,:)
		end do
		close(fu)
		open(unit=fu, file="test_outputs/nudmB_re.dat", status="old")
		do i=1, 3
			read (fu, *) nr(i,:)
		end do
		close(fu)
		open(unit=fu, file="test_outputs/nudmB_im.dat", status="old")
		do i=1, 3
			read (fu, *) ni(i,:)
		end do
		close(fu)
		er = (/5d-5,1d-4,0.d0/)
		do i=1,3
			do j=1,3
				write(tmparg,"('lepton matrix B ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg), leptonDensities(i,j), m(i,j), 1d-16, er(i))
				write(tmparg,"('nu density matrix re B ',2I1)") i,j
				call assert_double(trim(tmparg), nuDensities%re(i,j), nr(i,j), 2d-6)
				write(tmparg,"('nu density matrix im B ',2I1)") i,j
				call assert_double(trim(tmparg), nuDensities%im(i,j), ni(i,j), 1d-7)
			end do
		end do
		!C
		do iy=1, Ny
			nuDensMatVecFD(iy)%re = 0.01d0*exp(-y_arr(iy))
			nuDensMatVecFD(iy)%im = 0.4d0*exp(-y_arr(iy))
		end do
		call updateMatterDensities(2.d1, 1.2d0)
		open(unit=fu, file="test_outputs/leptmatC.dat", status="old")
		do i=1, 3
			read (fu, *) m(i,:)
		end do
		close(fu)
		open(unit=fu, file="test_outputs/nudmC_re.dat", status="old")
		do i=1, 3
			read (fu, *) nr(i,:)
		end do
		close(fu)
		open(unit=fu, file="test_outputs/nudmC_im.dat", status="old")
		do i=1, 3
			read (fu, *) ni(i,:)
		end do
		close(fu)
		do i=1,3
			do j=1,3
				write(tmparg,"('lepton matrix C ',2I1)") i,j
				call assert_double(trim(tmparg), leptonDensities(i,j), m(i,j), 1d-25)
				write(tmparg,"('nu density matrix re C ',2I1)") i,j
				call assert_double(trim(tmparg), nuDensities%re(i,j), nr(i,j), 1d-21)
				write(tmparg,"('nu density matrix im C ',2I1)") i,j
				call assert_double(trim(tmparg), nuDensities%im(i,j), ni(i,j), 1d-19)
			end do
		end do
		call printTotalTests
		call resetTestCounter
	end subroutine do_test_nu_matrices

	subroutine do_test_cosmology
		real(dl), dimension(:), allocatable :: ndmv_re
		integer :: i, ix, iy, n
		real(dl) :: x, z, r, rn
		real(dl), dimension(8) :: ve1, ve2
		character(len=100) :: tmpstr

		allocate(ndmv_re(Ny))
		do iy=1, Ny
			nuDensMatVecFD(iy)%re = 0.d0
			nuDensMatVecFD(iy)%im = 0.d0
		end do
		do i=1, flavorNumber
			do iy=1, Ny
				ndmv_re(iy) = 1.d0*i
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i
			end do
			call interpNuDens%re(i, i)%replace(Ny, y_arr, ndmv_re)
		end do

		call printTestBlockName("Cosmology")

		n=8
		ve1=1d-6
		ve1(7)=2e-6
		ve1(8)=1e-5
		ve2=(/1d-4, 1d-4, 1d-4, 1d-4, 1d-3, 1.5d-3, 3.5d-3, 3d-3/)
		open(unit=fu, file="test_outputs/elDens.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("elDensF test "//trim(tmpstr), electrons%energyDensityFull(x, z, .false.), r, ve1(i))
			call assert_double_rel("elDens test "//trim(tmpstr), electrons%energyDensity(x, z, .false.), r, ve2(i))
		end do
		close(fu)
		ve2=(/1d-4, 1d-4, 1d-4, 1d-4, 1.1d-3, 2.d-3, 3.5d-3, 3d-3/)
		open(unit=fu, file="test_outputs/elPress.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("elPressF test "//trim(tmpstr), electrons%pressureFull(x, z, .false.), r, ve1(i))
			call assert_double_rel("elPress test "//trim(tmpstr), electrons%pressure(x, z, .false.), r, ve2(i))
		end do
		close(fu)
		ve2=(/1d-4, 1d-4, 1d-4, 1d-4, 1d-3, 1.5d-3, 3.5d-3, 3d-3/)
		open(unit=fu, file="test_outputs/elEntropy.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("elEntropy test "//trim(tmpstr), electrons%entropy(x, z), r, ve2(i))
		end do
		close(fu)
		open(unit=fu, file="test_outputs/elNumDens.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("elNumDens test "//trim(tmpstr), electrons%numberDensity(x, z, .false.), r, ve1(i))
		end do
		close(fu)

		n=7
		ve1=1d-6
		ve1(7)=2e-6
		ve1(8)=1e-5
		ve2=(/1d-4, 5d-3, 1d-4, 1d-4, 1d-4, 1d-4, 1d-4, 1d-4/)
		open(unit=fu, file="test_outputs/muDens.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel_safe("muDensF test "//trim(tmpstr), muons%energyDensityFull(x, z, .false.), r, 1d-20, ve1(i))
			call assert_double_rel_safe("muDens test "//trim(tmpstr), muons%energyDensity(x, z, .false.), r, 1d-20, ve2(i))
		end do
		close(fu)
		ve2=(/1d-4, 5d-3, 1d-4, 1d-4, 1d-4, 1d-4, 1d-4, 1d-4/)
		open(unit=fu, file="test_outputs/muPress.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel_safe("muPressF test "//trim(tmpstr), muons%pressureFull(x, z, .false.), r, 1d-20, ve1(i))
			call assert_double_rel_safe("muPress test "//trim(tmpstr), muons%pressure(x, z, .false.), r, 1d-20, ve2(i))
		end do
		close(fu)
		ve2=(/1d-4, 5d-3, 1d-4, 1d-4, 1d-4, 1d-4, 1d-4, 1d-4/)
		open(unit=fu, file="test_outputs/muEntropy.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel_safe("muEntropy test "//trim(tmpstr), muons%entropy(x, z), r, 1d-20, ve2(i))
		end do
		close(fu)
		open(unit=fu, file="test_outputs/muNumDens.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel_safe("muNumDens test "//trim(tmpstr), muons%numberDensity(x, z, .false.), r, 1d-20, ve1(i))
		end do
		close(fu)

		n=5
		ve1=1d-7
		open(unit=fu, file="test_outputs/photDens.dat", status="old")
		do i=1, n
			read (fu, *) z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("photDens test "//trim(tmpstr), photonDensity(z), r, ve1(i))
		end do
		close(fu)
		open(unit=fu, file="test_outputs/photEntr.dat", status="old")
		do i=1, n
			read (fu, *) z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("photEntropy test "//trim(tmpstr), photonEntropy(z), r, ve1(i))
		end do
		close(fu)
		open(unit=fu, file="test_outputs/photNumDens.dat", status="old")
		do i=1, n
			read (fu, *) z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("photNumDens test "//trim(tmpstr), photonNumberDensity(z), r, ve1(i))
		end do
		close(fu)

		n=5
		ve1=1d-4
		ve2=(/1.,1.,1.,2.,3.,0.,0.,0./)
		open(unit=fu, file="test_outputs/nuDens.dat", status="old")
		do i=1, n
			read (fu, *) z, r, rn
			write(tmpstr, "(I1)") i
			call assert_double_rel("nuDens test "//trim(tmpstr), nuDensity(z,int(ve2(i))), ve2(i)*r, ve1(i))
		end do
		close(fu)

		ve1=(/1d-4,1d-4,2d-4,5d-4,1d-4,0.d0,0.d0,0.d0/)
		open(unit=fu, file="test_outputs/nuDensEq.dat", status="old")
		do i=1, n
			read (fu, *) z, r, rn
			do ix=1, flavorNumber
				do iy=1, Ny
				nuDensMatVecFD(iy)%re(ix, ix) = 1.d0*ix * fermiDirac(y_arr(iy) / z)
				end do
			end do
			write(tmpstr, "(I1)") i
			call assert_double_rel("nuDensNC test "//trim(tmpstr), nuDensityNC(int(ve2(i)), int(ve2(i))), ve2(i)*r, ve1(i))
			call assert_double_rel("nuNumDensNC test "//trim(tmpstr), nuNumberDensityNC(int(ve2(i)), int(ve2(i))), ve2(i)*rn, ve1(i))
		end do
		close(fu)

		call assert_double("nuDensNC test 6", nuDensityNC(1, 2), 0.0d0, 1d-7)
		call assert_double("nuDensNC test 7", nuDensityNC(1, 2, .false.), 0.0d0, 1d-7)
		open(unit=fu, file="test_outputs/nuDens.dat", status="old")
		read (fu, *) z, r, rn
		close(fu)
		do iy=1, Ny
			nuDensMatVecFD(iy)%re(1, 2) = 2.d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%im(1, 2) = 0.1d0 * fermiDirac(y_arr(iy))
		end do
		call assert_double_rel("nuDensNC test 8", nuDensityNC(1, 2), 2.d0*r, 1d-4)
		call assert_double_rel("nuDensNC test 9", nuDensityNC(1, 2, .false.), 0.1d0*r, 1d-4)
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy))
			end do
		end do
		call assert_double_rel("allnuDensNC test 1", allNuDensity(), 6*r, 1d-3)

		ftqed_log_term = .false.
		ftqed_ord3 = .false.
		n=5
		ve1=(/1d-6,1d-6,1d-6,1d-6,1.d6,0.d0,0.d0,0.d0/)
		open(unit=fu, file="test_outputs/ftqed_dR2.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("dRhoft o2 test "//trim(tmpstr), deltaRhoTot_em(x, z), r, ve1(i))
		end do
		close(fu)
		ftqed_ord3 = .true.
		open(unit=fu, file="test_outputs/ftqed_dRo23.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("dRhoft o23 test "//trim(tmpstr), deltaRhoTot_em(x, z), r, ve1(i))
		end do
		close(fu)
		ftqed_log_term = .true.
		ftqed_ord3 = .false.
		ve1=(/3d-5,1d-6,1d-5,4d-5,3d-5,0.d0,0.d0,0.d0/)
		open(unit=fu, file="test_outputs/ftqed_dRo2l.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("dRhoft o2l test "//trim(tmpstr), deltaRhoTot_em(x, z), r, ve1(i))
		end do
		close(fu)
		ftqed_ord3 = .true.
		open(unit=fu, file="test_outputs/ftqed_dRo23l.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("dRhoft o23l test "//trim(tmpstr), deltaRhoTot_em(x, z), r, ve1(i))
		end do
		close(fu)
		ftqed_log_term = .false.
		ftqed_ord3 = .false.
		ve1=(/1d-6,1d-6,1d-6,1d-6,1d-6,0.d0,0.d0,0.d0/)
		open(unit=fu, file="test_outputs/ftqed_dP2.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("dPft o2 test "//trim(tmpstr), deltaPTot_em(x, z), r, ve1(i))
		end do
		close(fu)
		ftqed_ord3 = .true.
		open(unit=fu, file="test_outputs/ftqed_dPo23.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("dPft o23 test "//trim(tmpstr), deltaPTot_em(x, z), r, ve1(i))
		end do
		close(fu)
		ftqed_log_term = .true.
		ftqed_ord3 = .false.
		ve1=(/2d-5,1d-6,4d-6,4d-5,2d-5,0.d0,0.d0,0.d0/)
		open(unit=fu, file="test_outputs/ftqed_dPo2l.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("dPft o2l test "//trim(tmpstr), deltaPTot_em(x, z), r, ve1(i))
		end do
		close(fu)
		ftqed_ord3 = .true.
		open(unit=fu, file="test_outputs/ftqed_dPo23l.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("dPft o23l test "//trim(tmpstr), deltaPTot_em(x, z), r, ve1(i))
		end do
		close(fu)
		ftqed_log_term = .false.
		ftqed_ord3 = .false.

		ftqed_temperature_corr = .false.
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
			end do
		end do
		n=4
		ve1=1d-5
		open(unit=fu, file="test_outputs/radDens.dat", status="old")
		do i=1, n
			read (fu, *) x, z, r
			write(tmpstr, "(I1)") i
			call assert_double_rel("radDens test "//trim(tmpstr), totalRadiationDensity(x, z), r, ve1(i))
		end do
		close(fu)

		!check Neff values from some previous papers
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
			end do
		end do
		call assert_double_rel("Neff test 1", Neff_from_rho_z(zid), 3.0d0, 1d-5)
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
			end do
		end do
		call assert_double_rel("Neff test 2", Neff_from_rho_z(1.39975d0), 3.011d0, 1d-4)
		do iy=1, Ny
			nuDensMatVecFD(iy)%re(1, 1) = 1.2d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(2, 2) = 1.d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(3, 3) = 1.d0 * fermiDirac(y_arr(iy))
		end do
		call assert_double_rel("Neff test 3", Neff_from_rho_z(zid), 3.2d0, 1d-5)
		do iy=1, Ny
			nuDensMatVecFD(iy)%re(1, 1) = 1.00920d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(2, 2) = 1.00392d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(3, 3) = 1.00392d0 * fermiDirac(y_arr(iy))
		end do
		call assert_double("Neff test 4", Neff_from_rho_z(1.397843d0), 3.045d0, 1d-3)
		do iy=1, Ny
			nuDensMatVecFD(iy)%re(1, 1) = 1.00699d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(2, 2) = 1.00511d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(3, 3) = 1.00519d0 * fermiDirac(y_arr(iy))
		end do
		call assert_double("Neff test 5", Neff_from_rho_z(1.39779d0), 3.045d0, 1d-3)

		do iy=1, Ny
			nuDensMatVecFD(iy)%re(1, 1) = 1.d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(2, 2) = 1.d0 * fermiDirac(y_arr(iy))
			nuDensMatVecFD(iy)%re(3, 3) = 1.d0 * fermiDirac(y_arr(iy))
		end do
		deallocate(ndmv_re)

		ftqed_temperature_corr = .true.

		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
			end do
		end do

		open(unit=fu, file="test_outputs/ftqed_dR2.dat", status="old")
		read (fu, *) x, z, rn
		close(fu)
		open(unit=fu, file="test_outputs/radDens.dat", status="old")
		read (fu, *) x, z, r
		close(fu)
		call assert_double("radDens test ftqed", totalRadiationDensity(x, z), r+rn, 1d-5)

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_cosmology

	subroutine do_test_JKG
		real(dl), dimension(2) :: res
		real(dl), dimension(4) :: v1, v2, zs
		real(dl) :: x, z, r1, r2
		integer :: ix, iy
		character(len=300) :: tmparg

		call printTestBlockName("JKG")

		do iy=0,4,2
			write(tmparg,"('test_outputs/jkg_j',I1,'.dat')") iy
			open(unit=fu, file=trim(tmparg), status="old")
			do ix=1,4
				read (fu, *) x,r1
				write(tmparg,"('J',I1,' test ',I1)") iy,ix
				call assert_double(trim(tmparg), J_funcFull(x, iy), r1, 1d-7)
			end do
			close(fu)
		end do
		do iy=0,4,2
			write(tmparg,"('test_outputs/jkg_jp',I1,'.dat')") iy
			open(unit=fu, file=trim(tmparg), status="old")
			do ix=1,4
				read (fu, *) x,r1
				write(tmparg,"('Jp',I1,' test ',I1)") iy,ix
				call assert_double(trim(tmparg), JprimeFull(x, iy), r1, 1d-7)
			end do
			close(fu)
		end do
		do iy=0,2,2
			write(tmparg,"('test_outputs/jkg_k',I1,'.dat')") iy
			open(unit=fu, file=trim(tmparg), status="old")
			do ix=1,4
				read (fu, *) x,r1
				write(tmparg,"('K',I1,' test ',I1)") iy,ix
				call assert_double(trim(tmparg), K_funcFull(x, iy), r1, 1d-7)
			end do
			close(fu)
		end do
		do iy=0,2,2
			write(tmparg,"('test_outputs/jkg_kp',I1,'.dat')") iy
			open(unit=fu, file=trim(tmparg), status="old")
			do ix=1,4
				read (fu, *) x,r1
				write(tmparg,"('Kp',I1,' test ',I1)") iy,ix
				call assert_double(trim(tmparg), KprimeFull(x, iy), r1, 1d-7)
			end do
			close(fu)
		end do

		open(unit=fu, file="test_outputs/jkg_elcontr.dat", status="old")
		do ix=1,3
			read (fu, *) x,r1,r2
			res = electrons%dzodx_terms(x)
			write(tmparg,"(I1)") ix
			call assert_double("elContr test "//trim(tmparg)//"a", res(1), r1, 1d-7)
			call assert_double("elContr test "//trim(tmparg)//"b", res(2), r2, 1d-7)
		end do
		close(fu)
		open(unit=fu, file="test_outputs/jkg_mucontr.dat", status="old")
		do ix=1,2
			read (fu, *) x,r1,r2
			res = muons%dzodx_terms(x)
			write(tmparg,"(I1)") ix
			call assert_double_rel_safe("elContr test "//trim(tmparg)//"a", res(1), r1, 1d-15,1d-7)
			call assert_double_rel_safe("elContr test "//trim(tmparg)//"b", res(2), r2, 1d-15,1d-7)
		end do
		close(fu)

		zs = (/1., 10., 1.5, 0.1/)
		open(unit=fu, file="test_outputs/jkg_g12_o2.dat", status="old")
		do ix=1,4
			read (fu, *) x,r1,r2
			res = G12_funcFull(x*zs(ix), zs(ix))
			write(tmparg,"(I1)") ix
			call assert_double("G1 o2 test "//trim(tmparg), res(1), r1, 1d-7)
			call assert_double("G2 o2 test "//trim(tmparg), res(2), r2, 1d-7)
		end do
		close(fu)

		ftqed_log_term = .true.
		v1 = (/3d-4, 2d-5, 2d-4, 4d-5/)
		v2 = (/1d-6, 1d-5, 5d-5, 2d-5/)
		open(unit=fu, file="test_outputs/jkg_g12_o2l.dat", status="old")
		do ix=1,4
			read (fu, *) x,z,r1,r2
			res = G12_funcFull(x, z)
			write(tmparg,"(I1)") ix
			call assert_double_rel("G1 o2l test "//trim(tmparg), res(1), r1, v1(ix))
			call assert_double_rel("G2 o2l test "//trim(tmparg), res(2), r2, v2(ix))
		end do
		close(fu)
		ftqed_log_term = .false.

		ftqed_ord3 = .true.
		v1 = (/1d-4, 1d-5, 1d-5, 1d-5/)
		v2 = (/1d-6, 1d-5, 1d-5, 1d-5/)
		open(unit=fu, file="test_outputs/jkg_g12_o23.dat", status="old")
		do ix=1,4
			read (fu, *) x,z,r1,r2
			res = G12_funcFull(x,z)
			write(tmparg,"(I1)") ix
			call assert_double("G1 o23 test "//trim(tmparg), res(1), r1, v1(ix))
			call assert_double("G2 o23 test "//trim(tmparg), res(2), r2, v2(ix))
		end do
		close(fu)

		ftqed_log_term = .true.
		v1 = (/1d-4, 1d-5, 1d-5, 1d-5/)
		v2 = (/1d-6, 1d-5, 1d-5, 1d-5/)
		open(unit=fu, file="test_outputs/jkg_g12_o23l.dat", status="old")
		do ix=1,4
			read (fu, *) x,z,r1,r2
			res = G12_funcFull(x, z)
			write(tmparg,"(I1)") ix
			call assert_double("G1 o23l test "//trim(tmparg), res(1), r1, v1(ix))
			call assert_double("G2 o23l test "//trim(tmparg), res(2), r2, v2(ix))
		end do
		close(fu)
		ftqed_log_term = .false.
		ftqed_ord3 = .false.

#ifndef NO_INTERPOLATION
		open(unit=fu, file="test_outputs/jkg_ab.dat", status="old")
		do ix=1,3
			read (fu, *) x,r1,r2
			res = dzodxcoef_interp_func(x)
			write(tmparg,"(I1)") ix
			call assert_double("A test "//trim(tmparg), res(1), r1, 1d-5)
			call assert_double("B test "//trim(tmparg), res(2), r2, 1d-5)
		end do
		close(fu)
#endif

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_JKG

	subroutine do_test_dzodx
		integer :: n
		real(dl), dimension(:), allocatable :: ydot
		integer :: m
		real(dl) :: x, z, r, rn
		integer :: ix
		character(len=300) :: tmparg
		real(dl), dimension(3) :: v1

		n=ntot
		allocate(ydot(n))
		ydot = 0.d0
		do m=1, Ny
			ydot((m-1)*flavNumSqu + 1) = cos(0.02d0*y_arr(m))
			ydot((m-1)*flavNumSqu + 2) = y_arr(m)/20.d0
			ydot((m-1)*flavNumSqu + 3) = 1.d0
		end do
		call printTestBlockName("dz/dx functions")

#ifndef NO_INTERPOLATION
		open(unit=fu, file="test_outputs/dzodx_g.dat", status="old")
		v1=(/8d-3,7d-6,6d-6/)
		do ix=1,3
			read (fu, *) x,z,r
			write(tmparg,"(I1)") ix
			call dz_o_dx_old(x, z, ydot, n)
			call assert_double_rel("dz_o_dx old test "//trim(tmparg), ydot(n), r, v1(ix))
		end do
		close(fu)
#endif

		open(unit=fu, file="test_outputs/dzodx_g.dat", status="old")
		open(unit=fv, file="test_outputs/dzodx_n.dat", status="old")
		v1=(/5d-5,6d-6,6d-6/)
		do ix=1,3
			read (fu, *) x,z,r
			read (fv, *) x,z,rn
			write(tmparg,"(I1)") ix
			call dz_o_dx(x, 1.2d0, z, ydot, n)
			call assert_double_rel("dz_o_dx test "//trim(tmparg), ydot(n), r, v1(ix))
			call assert_double_rel("dw_o_dx test "//trim(tmparg), ydot(n-1), rn, v1(ix))
		end do
		close(fu)
		close(fv)

		deallocate(ydot)
		call printTotalTests
		call resetTestCounter
	end subroutine do_test_dzodx

	subroutine do_test_dme2
		real(dl) :: x, y, z, r
		integer :: ix
		character(len=300) :: tmparg
		real(dl), dimension(6) :: v1

		call printTestBlockName("dme2")

		ftqed_log_term=.false.
		open(unit=fu, file="test_outputs/ftqed_dme2.dat", status="old")
		do ix=1,6
			read (fu, *) x,y,z,r
			write(tmparg,"(I1)") ix
			call assert_double_rel("dme2F test "//trim(tmparg), dme2_electronFull(x, y, z), r, 1d-6)
#ifndef NO_INTERPOLATION
			call assert_double_rel("dme2 test "//trim(tmparg), dme2_electron(x, y, z), r, 1d-5)
			call assert_double_rel("dme2nl test "//trim(tmparg), dme2_nolog(x, z), r, 1d-5)
#endif
		end do
		close(fu)

		call assert_double("Ebare_i_dme test 1", Ebare_i_dme(0.3d0, 0.4d0, 1.44d0), sqrt(0.3d0**2+0.4d0**2+1.44d0), 1d-7)
		call assert_double("Ebare_i_dme test 2", Ebare_i_dme(3.d0, 7.d0, 22.d0), sqrt(3.d0**2+7.d0**2+22.d0), 1d-7)

		ftqed_log_term=.true.
		write(*,*)
		write(*,*) "now with log term in dme2"
		open(unit=fu, file="test_outputs/ftqed_dme2l.dat", status="old")
		v1=(/6d-5,1d-5,6d-4,7d-4,5d-4,1d-5/)
		do ix=1,6
			read (fu, *) x,y,z,r
			write(tmparg,"(I1)") ix
			call assert_double_rel("dme2F w log test "//trim(tmparg), dme2_electronFull(x, y, z), r, v1(ix))
		end do
		close(fu)
		open(unit=fu, file="test_outputs/ftqed_dme2.dat", status="old")
		read (fu, *) x,y,z,r
		close(fu)
		call assert_double_rel("dme2F logt test 1", dme2_electronFull(x, 0.d0, z, .false.), r, 1d-6)
		call assert_double_rel("dme2F logt test 2", dme2_electronFull(x, 1.d0, z, .false.), r, 1d-6)
		call assert_double_rel("dme2F logt test 3", dme2_electronFull(x, 10.d0, z, .false.), r, 1d-6)
		ftqed_log_term=.false.

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_dme2

	subroutine do_test_Di
		real(dl) :: y1, y2, y3, y4, d1, d2, d3
		integer :: ix
		character(len=100) :: tmpstr

		call printTestBlockName("D_i functions")

		do ix=1,6
			write(tmpstr, "('test_outputs/Di_',I1,'.dat')") ix
			open(unit=fu, file=trim(tmpstr), status="old")
			read (fu, *) y1, y2, y3, y4, d1, d2, d3
			close(fu)
			write(tmpstr, "(I1)") ix
			call assert_double_rel("D1  test "//trim(tmpstr), D1_full(y1, y2, y3, y4),  d1, 1d-7)
			call assert_double_rel("D1p test "//trim(tmpstr), D1_bis(y1, y2, y3, y4),   d1, 1d-7)
			call assert_double_rel("D2  test "//trim(tmpstr), D2_full(y1, y2, y3, y4),  d2, 1d-7)
			call assert_double_rel("D2c test "//trim(tmpstr), D2_cases(y1, y2, y3, y4), d2, 1d-7)
			call assert_double_rel("D2p test "//trim(tmpstr), D2_bis(y1, y2, y3, y4),   d2, 1d-7)
			call assert_double_rel("D3  test "//trim(tmpstr), D3_full(y1, y2, y3, y4),  d3, 1d-7)
			call assert_double_rel("D3c test "//trim(tmpstr), D3_cases(y1, y2, y3, y4), d3, 1d-7)
			call assert_double_rel("D3p test "//trim(tmpstr), D3_bis(y1, y2, y3, y4),   d3, 1d-7)
		end do

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_Di

	subroutine do_test_Pi_ij
		real(dl), dimension(2) :: temp_v2
		real(dl) :: y1,y2,y3,y4,x,dm,r1,r2
		integer :: ix
		character(len=100) :: tmpstr

		call printTestBlockName("Pi(yi,yj) functions")

		do ix=1,4
			write(tmpstr, "('test_outputs/Pi1_',I1,'.dat')") ix
			open(unit=fu, file=trim(tmpstr), status="old")
			read (fu, *) y1, y2, y3, y4, r1, r2
			close(fu)
			write(tmpstr, "(I1)") ix
			call assert_double_rel("Pi_1_12 test "//trim(tmpstr), PI1_12_full(y1,y2,y3,y4), r1, 1d-7)
			call assert_double_rel("Pi_1_13 test "//trim(tmpstr), PI1_13_full(y1,y2,y3,y4), r2, 1d-7)
		end do

		do ix=1,4
			write(tmpstr, "('test_outputs/Pi2nn_',I1,'.dat')") ix
			open(unit=fu, file=trim(tmpstr), status="old")
			read (fu, *) y1, y2, y3, y4, x, dm, r1, r2
			close(fu)
			write(tmpstr, "(I1)") ix
			temp_v2 = PI2_nn_f(y1,y2,y3,y4,Ebare_i_dme(x, y3, dm),Ebare_i_dme(x, y4, dm))
			call assert_double_rel("Pi_2_14a test "//trim(tmpstr), temp_v2(1), r1, 1d-7)
			call assert_double_rel("Pi_2_13  test "//trim(tmpstr), temp_v2(2), r2, 1d-7)
		end do

		do ix=1,4
			write(tmpstr, "('test_outputs/Pi2ne_',I1,'.dat')") ix
			open(unit=fu, file=trim(tmpstr), status="old")
			read (fu, *) y1, y2, y3, y4, x, dm, r1, r2
			close(fu)
			write(tmpstr, "(I1)") ix
			temp_v2 = PI2_ne_f(y1,y2,y3,y4,Ebare_i_dme(x, y2, dm),Ebare_i_dme(x, y4, dm))
			call assert_double_rel("Pi_2_14s test "//trim(tmpstr), temp_v2(1), r1, 1d-7)
			call assert_double_rel("Pi_2_12  test "//trim(tmpstr), temp_v2(2), r2, 1d-7)
		end do

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_Pi_ij

	subroutine do_f_ann_sc_re_tests_eq
		integer :: ix, iy
		real(dl) :: fdA, fdB, f1, f2, f3
		real(dl), dimension(3,3) :: m
		type(cmplxMatNN) :: nA, nB
		character(len=300) :: tmparg

		f1 = fermiDirac(0.3d0)
		f2 = fermiDirac(0.4d0)
		f3 = fermiDirac(0.1d0)
		call printTestBlockName("F_ann_re functions at equilibrium")
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
		!FRReq[0.1, 0.2, 0.3, 0.4]
		open(unit=fu, file="test_outputs/FRR_eq_1.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FRR_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,2, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.1, 0.2, 0.3, 0.4]
		open(unit=fu, file="test_outputs/FLL_eq_1.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,1, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.1, 0.2, 0.3, 0.4]
		open(unit=fu, file="test_outputs/FRL_eq_1.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 2,1, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.1, 0.2, 0.3, 0.4]
		open(unit=fu, file="test_outputs/FLR_eq_1.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_ann eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_ann eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_ann eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, 1,2, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, 1,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!second series
		fdA = fermiDirac(0.4d0)
		do ix=1, flavorNumber
			nA%re(ix,ix) = fdA
		end do
		!FRReq[0.4, 0.2, 0.3, 0.1]
		open(unit=fu, file="test_outputs/FRR_eq_2.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FRR_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,2, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,2, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,2, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.4, 0.2, 0.3, 0.1]
		open(unit=fu, file="test_outputs/FLL_eq_2.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,1, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,1, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,1, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.4, 0.2, 0.3, 0.1]
		open(unit=fu, file="test_outputs/FRL_eq_2.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,1, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,1, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 2,1, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.4, 0.2, 0.3, 0.1]
		open(unit=fu, file="test_outputs/FLR_eq_2.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_ann eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,2, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_ann eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,2, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_ann eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f3, 1,2, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f3, 1,2, iy, ix), 0.d0, 1d-7)
			end do
		end do


		f1 = fermiDirac(0.2d0)
		f2 = fermiDirac(0.4d0)
		f3 = fermiDirac(0.1d0)
		call printTestBlockName("F_sc_re functions at equilibrium")
		fdA = fermiDirac(0.1d0)
		fdB = fermiDirac(0.3d0)
		do ix=1, flavorNumber
			nA%re(ix,ix) = fdA
			nB%re(ix,ix) = fdB
		end do
		!FRReq[0.1, 0.2, 0.3, 0.4]
		open(unit=fu, file="test_outputs/FRR_eq_1.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FRR_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 2,2, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 2,2, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 2,2, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.1, 0.2, 0.3, 0.4]
		open(unit=fu, file="test_outputs/FLL_eq_1.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 1,1, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 1,1, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 1,1, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.1, 0.2, 0.3, 0.4]
		open(unit=fu, file="test_outputs/FRL_eq_1.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 2,1, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 2,1, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 2,1, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.1, 0.2, 0.3, 0.4]
		open(unit=fu, file="test_outputs/FLR_eq_1.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_sc eq test 1 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 1,2, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_sc eq test 1 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 1,2, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_sc eq test 1 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, 1,2, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, 1,2, iy, ix), 0.d0, 1d-7)
			end do
		end do

		fdA = fermiDirac(0.4d0)
		do ix=1, flavorNumber
			nA%re(ix,ix) = fdA
		end do
		!FRReq[0.4, 0.2, 0.3, 0.1]
		open(unit=fu, file="test_outputs/FRR_eq_2.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FRR_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 2,2, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 2,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRR_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 2,2, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 2,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRR_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 2,2, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 2,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLLeq[0.4, 0.2, 0.3, 0.1]
		open(unit=fu, file="test_outputs/FLL_eq_2.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FLL_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 1,1, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 1,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLL_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 1,1, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 1,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLL_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 1,1, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 1,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FRLeq[0.4, 0.2, 0.3, 0.1]
		open(unit=fu, file="test_outputs/FRL_eq_2.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FRL_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 2,1, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 2,1, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FRL_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 2,1, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 2,1, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FRL_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 2,1, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 2,1, iy, ix), 0.d0, 1d-7)
			end do
		end do
		!FLReq[0.4, 0.2, 0.3, 0.1]
		open(unit=fu, file="test_outputs/FLR_eq_2.dat", status="old")
		do ix=1, 3
			read (fu, *) m(ix,:)
		end do
		close(fu)
		do ix=1, flavorNumber
			write(tmparg,"('FLR_sc eq test 2 - ',2I1)") ix, ix
			call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 1,2, ix, ix), m(ix,ix), 1d-7)
			call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 1,2, ix, ix), 0.d0, 1d-7)
			do iy=ix+1, flavorNumber
				write(tmparg,"('FLR_sc eq test 2 - ',2I1)") ix, iy
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 1,2, ix, iy), m(ix,iy), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 1,2, ix, iy), 0.d0, 1d-7)
				write(tmparg,"('FLR_sc eq test 2 - ',2I1)") iy, ix
				call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f3, 1,2, iy, ix), m(iy,ix), 1d-7)
				call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f3, 1,2, iy, ix), 0.d0, 1d-7)
			end do
		end do
		call deallocateCmplxMat(nA)
		call deallocateCmplxMat(nB)
		call printTotalTests
		call resetTestCounter
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
		call printTestBlockName("F_ann functions, empty rho")
		call allocateCmplxMat(nA)
		call allocateCmplxMat(nB)
		nA%re = 0.d0
		nA%im = 0.d0
		nB%re = 0.d0
		nB%im = 0.d0
		!imaginary parts and off-diagonal are zero
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
		!diagonal
		do a=1,2
			do b=1,2
				open(unit=fu, file="test_outputs/Fann_er1_"//chLR(a)//chLR(b)//"_re.dat", status="old")
				do ix=1, 3
					read (fu, *) tmpmatA(ix,:)
				end do
				close(fu)
				do ix=1, flavorNumber
					write(tmparg,"('F_ann test 1 "//chLR(a)//chLR(b)//" empty rho ',2I1)") ix,ix
					call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, a,b, ix,ix), tmpmatA(ix,ix), 1d-7)
				end do
			end do
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
		do a=1,2
			do b=1,2
				open(unit=fu, file="test_outputs/Fann_er2_"//chLR(a)//chLR(b)//"_re.dat", status="old")
				do ix=1, 3
					read (fu, *) tmpmatA(ix,:)
				end do
				close(fu)
				do ix=1, flavorNumber
						write(tmparg,"('F_ann test 2 "//chLR(a)//chLR(b)//" empty rho ',2I1)") ix,ix
						call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, a,b, ix,ix), tmpmatA(ix,ix), 1d-7)
				end do
			end do
		end do

		call printTestBlockName("F_sc functions, empty rho")
		do a=1, 2
			do b=1, 2
				do ix=1, flavorNumber
					iy=ix
					write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, ix, iy
					call assert_double(trim(tmparg)//" re", F_ab_sc_re(nA, nB, f1, f2, a, b, ix, iy), 0.0d0, 1d-7)
					call assert_double(trim(tmparg)//" im", F_ab_sc_im(nA, nB, f1, f2, a, b, ix, iy), 0.0d0, 1d-7)
					do iy=ix+1, flavorNumber
						write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, ix, iy
						call assert_double(trim(tmparg)//" re", F_ab_sc_re(nA, nB, f1, f2, a, b, ix, iy), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_sc_im(nA, nB, f1, f2, a, b, ix, iy), 0.0d0, 1d-7)
						write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, iy, ix
						call assert_double(trim(tmparg)//" re", F_ab_sc_re(nA, nB, f1, f2, a, b, iy, ix), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_sc_im(nA, nB, f1, f2, a, b, iy, ix), 0.0d0, 1d-7)
					end do
				end do
			end do
		end do
		do a=1, 2
			do b=1, 2
				do ix=1, flavorNumber
					iy=ix
					write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, ix, iy
					call assert_double(trim(tmparg)//" re", F_ab_sc_re(nB, nA, f2, f3, a, b, ix, iy), 0.0d0, 1d-7)
					call assert_double(trim(tmparg)//" im", F_ab_sc_im(nB, nA, f2, f3, a, b, ix, iy), 0.0d0, 1d-7)
					do iy=ix+1, flavorNumber
						write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, ix, iy
						call assert_double(trim(tmparg)//" re", F_ab_sc_re(nB, nA, f2, f3, a, b, ix, iy), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_sc_im(nB, nA, f2, f3, a, b, ix, iy), 0.0d0, 1d-7)
						write(tmparg,"('F_sc',2I1,' test 1 empty rho ',2I1)") a, b, iy, ix
						call assert_double(trim(tmparg)//" re", F_ab_sc_re(nB, nA, f2, f3, a, b, iy, ix), 0.0d0, 1d-7)
						call assert_double(trim(tmparg)//" im", F_ab_sc_im(nB, nA, f2, f3, a, b, iy, ix), 0.0d0, 1d-7)
					end do
				end do
			end do
		end do

		call printTestBlockName("F_ann functions, full rho")
		!rhoA
		open(unit=fu, file="test_outputs/Fann_fr12_rhoA_re.dat", status="old")
		open(unit=fv, file="test_outputs/Fann_fr12_rhoA_im.dat", status="old")
		do ix=1, 3
			read (fu, *) nA%re(ix,:)
			read (fv, *) nA%im(ix,:)
		end do
		close(fu)
		close(fv)
		!rhoB
		open(unit=fu, file="test_outputs/Fann_fr12_rhoB_re.dat", status="old")
		open(unit=fv, file="test_outputs/Fann_fr12_rhoB_im.dat", status="old")
		do ix=1, 3
			read (fu, *) nB%re(ix,:)
			read (fv, *) nB%im(ix,:)
		end do
		close(fu)
		close(fv)
		!first series
		do a=1,2
			do b=1,2
				open(unit=fu, file="test_outputs/Fann_fr1_"//chLR(a)//chLR(b)//"_re.dat", status="old")
				open(unit=fv, file="test_outputs/Fann_fr1_"//chLR(a)//chLR(b)//"_im.dat", status="old")
				do ix=1, 3
					read (fu, *) tmpmatA(ix,:)
					read (fv, *) tmpmatB(ix,:)
				end do
				close(fu)
				close(fv)
				do ix=1, flavorNumber
					do iy=1, flavorNumber
						write(tmparg,"('F_ann test 1 "//chLR(a)//chLR(b)//" full rho ',2I1)") ix,iy
						call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, a,b, ix,iy), tmpmatA(ix,iy), 1d-7)
						call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, a,b, ix,iy), tmpmatB(ix,iy), 1d-7)
					end do
				end do
			end do
		end do
		!second series
		do a=1,2
			do b=1,2
				open(unit=fu, file="test_outputs/Fann_fr2_"//chLR(a)//chLR(b)//"_re.dat", status="old")
				open(unit=fv, file="test_outputs/Fann_fr2_"//chLR(a)//chLR(b)//"_im.dat", status="old")
				do ix=1, 3
					read (fu, *) tmpmatA(ix,:)
					read (fv, *) tmpmatB(ix,:)
				end do
				close(fu)
				close(fv)
				do ix=1, flavorNumber
					do iy=1, flavorNumber
					write(tmparg,"('F_ann test 2 "//chLR(a)//chLR(b)//" full rho ',2I1)") ix,iy
					call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, a,b, ix,iy), tmpmatA(ix,iy), 1d-7)
					call assert_double(trim(tmparg)//"im", F_ab_ann_im(nB, nA, f2, f3, a,b, ix,iy), tmpmatB(ix,iy), 1d-7)
					end do
				end do
			end do
		end do

		call printTestBlockName("F_sc functions, full rho")
		!first series
		do a=1,2
			do b=1,2
				open(unit=fu, file="test_outputs/Fsc_fr1_"//chLR(a)//chLR(b)//"_re.dat", status="old")
				open(unit=fv, file="test_outputs/Fsc_fr1_"//chLR(a)//chLR(b)//"_im.dat", status="old")
				do ix=1, 3
					read (fu, *) tmpmatA(ix,:)
					read (fv, *) tmpmatB(ix,:)
				end do
				close(fu)
				close(fv)
				do ix=1, flavorNumber
					do iy=1, flavorNumber
						write(tmparg,"('F_sc test 1 "//chLR(a)//chLR(b)//" full rho ',2I1)") ix,iy
						call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, a,b, ix,iy), tmpmatA(ix,iy), 1d-7)
						call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, a,b, ix,iy), tmpmatB(ix,iy), 1d-7)
					end do
				end do
			end do
		end do
		!second series
		do a=1,2
			do b=1,2
				open(unit=fu, file="test_outputs/Fsc_fr2_"//chLR(a)//chLR(b)//"_re.dat", status="old")
				open(unit=fv, file="test_outputs/Fsc_fr2_"//chLR(a)//chLR(b)//"_im.dat", status="old")
				do ix=1, 3
					read (fu, *) tmpmatA(ix,:)
					read (fv, *) tmpmatB(ix,:)
				end do
				close(fu)
				close(fv)
				do ix=1, flavorNumber
					do iy=1, flavorNumber
						write(tmparg,"('F_sc test 2 "//chLR(a)//chLR(b)//" full rho ',2I1)") ix,iy
						call assert_double(trim(tmparg)//"re", F_ab_sc_re(nB, nA, f2, f3, a,b, ix,iy), tmpmatA(ix,iy), 1d-7)
						call assert_double(trim(tmparg)//"im", F_ab_sc_im(nB, nA, f2, f3, a,b, ix,iy), tmpmatB(ix,iy), 1d-7)
					end do
				end do
			end do
		end do

		call printTestBlockName("F_ann functions, final tests")
		!rhoA
		open(unit=fu, file="test_outputs/Fann_fr34_rhoA_re.dat", status="old")
		open(unit=fv, file="test_outputs/Fann_fr34_rhoA_im.dat", status="old")
		do ix=1, 3
			read (fu, *) nA%re(ix,:)
			read (fv, *) nA%im(ix,:)
		end do
		close(fu)
		close(fv)
		!rhoB
		open(unit=fu, file="test_outputs/Fann_fr34_rhoB_re.dat", status="old")
		open(unit=fv, file="test_outputs/Fann_fr34_rhoB_im.dat", status="old")
		do ix=1, 3
			read (fu, *) nB%re(ix,:)
			read (fv, *) nB%im(ix,:)
		end do
		close(fu)
		close(fv)
		!annihilation
		!first series
		do a=1,2
			do b=1,2
				open(unit=fu, file="test_outputs/Fann_fr3_"//chLR(a)//chLR(b)//"_re.dat", status="old")
				open(unit=fv, file="test_outputs/Fann_fr3_"//chLR(a)//chLR(b)//"_im.dat", status="old")
				do ix=1, 3
					read (fu, *) tmpmatA(ix,:)
					read (fv, *) tmpmatB(ix,:)
				end do
				close(fu)
				close(fv)
				do ix=1, flavorNumber
					do iy=1, flavorNumber
						write(tmparg,"('F_ann test 3 "//chLR(a)//chLR(b)//" full rho ',2I1)") ix,iy
						call assert_double(trim(tmparg)//"re", F_ab_ann_re(nA, nB, f1, f2, a,b, ix,iy), tmpmatA(ix,iy), 1d-7)
						call assert_double(trim(tmparg)//"im", F_ab_ann_im(nA, nB, f1, f2, a,b, ix,iy), tmpmatB(ix,iy), 1d-7)
					end do
				end do
			end do
		end do
		!second series
		do a=1,2
			do b=1,2
				open(unit=fu, file="test_outputs/Fann_fr4_"//chLR(a)//chLR(b)//"_re.dat", status="old")
				open(unit=fv, file="test_outputs/Fann_fr4_"//chLR(a)//chLR(b)//"_im.dat", status="old")
				do ix=1, 3
					read (fu, *) tmpmatA(ix,:)
					read (fv, *) tmpmatB(ix,:)
				end do
				close(fu)
				close(fv)
				do ix=1, flavorNumber
					do iy=1, flavorNumber
					write(tmparg,"('F_ann test 4 "//chLR(a)//chLR(b)//" full rho ',2I1)") ix,iy
					call assert_double(trim(tmparg)//"re", F_ab_ann_re(nB, nA, f2, f3, a,b, ix,iy), tmpmatA(ix,iy), 1d-7)
					call assert_double(trim(tmparg)//"im", F_ab_ann_im(nB, nA, f2, f3, a,b, ix,iy), tmpmatB(ix,iy), 1d-7)
					end do
				end do
			end do
		end do
		!scattering
		!first series
		do a=1,2
			do b=1,2
				open(unit=fu, file="test_outputs/Fsc_fr3_"//chLR(a)//chLR(b)//"_re.dat", status="old")
				open(unit=fv, file="test_outputs/Fsc_fr3_"//chLR(a)//chLR(b)//"_im.dat", status="old")
				do ix=1, 3
					read (fu, *) tmpmatA(ix,:)
					read (fv, *) tmpmatB(ix,:)
				end do
				close(fu)
				close(fv)
				do ix=1, flavorNumber
					do iy=1, flavorNumber
						write(tmparg,"('F_sc test 3 "//chLR(a)//chLR(b)//" full rho ',2I1)") ix,iy
						call assert_double(trim(tmparg)//"re", F_ab_sc_re(nA, nB, f1, f2, a,b, ix,iy), tmpmatA(ix,iy), 1d-7)
						call assert_double(trim(tmparg)//"im", F_ab_sc_im(nA, nB, f1, f2, a,b, ix,iy), tmpmatB(ix,iy), 1d-7)
					end do
				end do
			end do
		end do
		!second series
		do a=1,2
			do b=1,2
				open(unit=fu, file="test_outputs/Fsc_fr4_"//chLR(a)//chLR(b)//"_re.dat", status="old")
				open(unit=fv, file="test_outputs/Fsc_fr4_"//chLR(a)//chLR(b)//"_im.dat", status="old")
				do ix=1, 3
					read (fu, *) tmpmatA(ix,:)
					read (fv, *) tmpmatB(ix,:)
				end do
				close(fu)
				close(fv)
				do ix=1, flavorNumber
					do iy=1, flavorNumber
						write(tmparg,"('F_sc test 4 "//chLR(a)//chLR(b)//" full rho ',2I1)") ix,iy
						call assert_double(trim(tmparg)//"re", F_ab_sc_re(nB, nA, f2, f3, a,b, ix,iy), tmpmatA(ix,iy), 1d-7)
						call assert_double(trim(tmparg)//"im", F_ab_sc_im(nB, nA, f2, f3, a,b, ix,iy), tmpmatB(ix,iy), 1d-7)
					end do
				end do
			end do
		end do
		call deallocateCmplxMat(nA)
		call deallocateCmplxMat(nB)
		call printTotalTests
		call resetTestCounter
	end subroutine do_f_ann_sc_re_tests_full

	subroutine do_test_coll_int
		real(dl) :: x,z,dme2
		type(coll_args) :: collArgs
		integer :: i, j, iy1, iy, oi
		real(dl) :: errr1,errr2, res1,res2,res3,res4, cf, ref
		real(dl) :: y1,y2a,y3a,y4a,y2s,y3s,y4s
		real(dl), dimension(3) :: tmparrS, tmparrA, tmperr, tmperr3,tmperr4
		real(dl), dimension(3, 3) :: tmpmatA, tmpmatB
		character(len=300) :: tmparg
		real(dl), dimension(:,:), allocatable :: fy2_arr

		allocate(fy2_arr(Ny, Ny))

		call printTestBlockName("Collision integrals")

		x=0.05d0
		iy1=7 !1.22151515151515
		oi=21
		z=1.06d0
		dme2=0.1d0
		y1 = y_arr(iy1)
		y2s = 5.2d0
		y3s = y_arr(oi)
		y4s = sqrt((y1+Ebare_i_dme(x, y2s, dme2)-y3s)**2-x*x-dme2)
		y2a = y_arr(oi)
		y3a = 3.d0
		y4a = sqrt((y1+y2a-Ebare_i_dme(x, y3a, dme2))**2-x*x-dme2)

		do iy=1, Ny
			nuDensMatVecFD(iy)%re = 0.d0
			nuDensMatVecFD(iy)%im = 0.d0
		end do
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy))
			end do
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
		open(unit=fu, file="test_outputs/collint_sc_fake.dat", status="old")
		read (fu, *) tmparrS(1:3)
		close(fu)
		tmparrA = (/1d-6, 1d-6, 1d-6/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i
			res2=coll_nue_sc_int(oi, y2s, collArgs, F_ab_sc_re)
			call assert_double_rel("test coll sc fake "//trim(tmparg), res2, tmparrS(i), tmparrA(i))
		end do

		open(unit=fu, file="test_outputs/collint_ann_fake.dat", status="old")
		read (fu, *) tmparrS(1:3)
		close(fu)
		tmparrA = (/1d-6, 1d-6, 1d-6/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i
			res2=coll_nue_ann_int(oi, y3a, collArgs, F_ab_ann_re)
			call assert_double_rel("test coll ann fake "//trim(tmparg), res2, tmparrS(i), tmparrA(i))
		end do

		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0*i * fermiDirac(y_arr(iy))
			end do
		end do

		!first series of sc and ann tests
		write(*,*) ""
		open(unit=fu, file="test_outputs/collint_sc_sing1.dat", status="old")
		read (fu, *) tmparrS(1:3)
		close(fu)
		tmparrA = (/1d-6, 1d-6, 1d-6/)
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i
			res2=coll_nue_sc_int(oi, y2s, collArgs, F_ab_sc_re)
			call assert_double_rel("test coll sc sing1 "//trim(tmparg), res2, tmparrS(i), tmparrA(i))
		end do

		open(unit=fu, file="test_outputs/collint_ann_sing1.dat", status="old")
		read (fu, *) tmparrS(1:3)
		close(fu)
		tmperr = (/1d-6, 1d-6, 1d-6/)
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i
			res2=coll_nue_ann_int(oi, y3a, collArgs, F_ab_ann_re)
			call assert_double_rel("test coll ann sing1 "//trim(tmparg), res2, tmparrS(i), tmperr(i))
		end do

		blocking=.false.
		write(*,*) ""
		open(unit=fu, file="test_outputs/coll_sc_A.dat", status="old")
		read (fu, *) tmparrS(1:3)
		close(fu)
		tmperr3 = (/2d-3, 2d-3, 2d-3/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			res2 = integrate_collint_nue_NC(coll_nue_sc_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
			write(tmparg,"('test coll sc A ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrS(i), tmperr3(i))
		end do

		write(*,*) ""
		open(unit=fu, file="test_outputs/coll_ann_A.dat", status="old")
		read (fu, *) tmparrA(1:3)
		close(fu)
		tmperr3 = (/2d-3, 2d-3, 2d-3/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			res2 = integrate_collint_nue_NC(coll_nue_ann_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
			write(tmparg,"('test coll ann A ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrA(i), tmperr3(i))
		end do
		blocking=.true.

		!second series of sc and ann tests
		iy1 = 67
		z = 1.3d0
		collArgs%y1 = y_arr(iy1)
		do i=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
			end do
		end do
		collArgs%x = 0.5d0
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%ix1 = 3
		collArgs%ix2 = 3
		do i=1, Ny
			do j=1, Ny
				fy2_arr(i,j) = coll_nue_sc_int(i, y_arr(j), collArgs, F_ab_sc_re)
			end do
		end do

		open(unit=fu, file="test_outputs/collint_sc_sing2.dat", status="old")
		read (fu, *) tmparrS(1:3)
		close(fu)
		tmperr = (/1d-6, 1d-6, 1d-6/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i
			res2=coll_nue_sc_int(oi, y2s, collArgs, F_ab_sc_re)
			call assert_double_rel("test coll sc sing2 "//trim(tmparg), res2, tmparrS(i), tmperr(i))
		end do

		write(*,*) ""
		open(unit=fu, file="test_outputs/coll_sc_B.dat", status="old")
		read (fu, *) tmparrS(1:3)
		close(fu)
		tmperr3 = (/2d-4, 2d-4, 2d-4/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			res2 = integrate_collint_nue_NC(coll_nue_sc_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
			write(tmparg,"('test coll sc B - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrS(i), tmperr3(i))
		end do

		collArgs%ix1 = 1
		collArgs%ix2 = 1
		open(unit=fu, file="test_outputs/collint_ann_sing2.dat", status="old")
		read (fu, *) tmparrS(1:3)
		close(fu)
		tmperr = (/1d-6, 1d-6, 1d-6/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			write(tmparg,"(I1)") i
			res2=coll_nue_ann_int(oi, y3a, collArgs, F_ab_ann_re)
			call assert_double_rel("test coll ann sing2 "//trim(tmparg), res2, tmparrS(i), tmperr(i))
		end do

		write(*,*) ""
		open(unit=fu, file="test_outputs/coll_ann_B.dat", status="old")
		read (fu, *) tmparrA(1:3)
		close(fu)
		tmperr3 = (/1.1d-4, 1.1d-4, 1.1d-4/)
		do i=1, flavorNumber
			collArgs%ix1 = i
			collArgs%ix2 = i
			res2 = integrate_collint_nue_NC(coll_nue_ann_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
			write(tmparg,"('test coll ann B - ',2I1)") i, i
			call assert_double_rel_verb(trim(tmparg), res2, tmparrA(i), tmperr3(i))
		end do

		x=0.05d0
		iy1=7 !1.22151515151515
		z=1.06d0
		call printTestBlockName("Collision integrands, im and off-diag")
		collArgs%x = x
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%dme2 = dme2

		do iy=1, Ny
			nuDensMatVecFD(iy)%re = 0.d0
			nuDensMatVecFD(iy)%im = 0.d0
		end do
		open(unit=fu, file="test_outputs/collint_imoff_rhoiy1_re.dat", status="old")
		open(unit=fv, file="test_outputs/collint_imoff_rhoiy1_im.dat", status="old")
		do i=1, 3
			read (fu, *) nuDensMatVecFD(iy1)%re(i,:)
			nuDensMatVecFD(iy1)%re(i, i) = nuDensMatVecFD(iy1)%re(i, i)*fermiDirac(y_arr(iy1))
			read (fv, *) nuDensMatVecFD(iy1)%im(i,:)
		end do
		close(fu)
		close(fv)
		open(unit=fu, file="test_outputs/collint_imoff_rhooi_re.dat", status="old")
		open(unit=fv, file="test_outputs/collint_imoff_rhooi_im.dat", status="old")
		do i=1, 3
			read (fu, *) nuDensMatVecFD(oi)%re(i,:)
			nuDensMatVecFD(oi)%re(i, i) = nuDensMatVecFD(oi)%re(i, i)*fermiDirac(y_arr(oi))
			read (fv, *) nuDensMatVecFD(oi)%im(i,:)
		end do
		close(fu)
		close(fv)

		open(unit=fu, file="test_outputs/collint_sc_imoff_integrand_re.dat", status="old")
		open(unit=fv, file="test_outputs/collint_sc_imoff_integrand_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		tmparrA(:) = (/1d-5, 1d-5, 1d-5/)
		tmparrS(:) = (/1d-5, 1d-5, 1d-5/)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('collint sc imoff integrand ',2I1)") i,j
				collArgs%ix1 = i
				collArgs%ix2 = j
				res1 = coll_nue_sc_int(21, 5.2d0, collArgs, F_ab_sc_re)
				res2 = coll_nue_sc_int(21, 5.2d0, collArgs, F_ab_sc_im)
				call assert_double_rel(trim(tmparg)//"re", res1, tmpmatA(i,j), tmparrA(i))
				if (i.eq.j) then
					call assert_double(trim(tmparg)//"im", res2, tmpmatB(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", res2, tmpmatB(i,j), tmparrS(i))
				end if
			end do
		end do

		open(unit=fu, file="test_outputs/collint_ann_imoff_integrand_re.dat", status="old")
		open(unit=fv, file="test_outputs/collint_ann_imoff_integrand_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		tmparrA(:) = (/1d-5, 1d-5, 1d-5/)
		tmparrS(:) = (/1d-5, 1d-5, 1d-5/)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('collint ann imoff integrand ',2I1)") i,j
				collArgs%ix1 = i
				collArgs%ix2 = j
				res1 = coll_nue_ann_int(21, 3.d0, collArgs, F_ab_ann_re)
				res2 = coll_nue_ann_int(21, 3.d0, collArgs, F_ab_ann_im)
				call assert_double_rel(trim(tmparg)//"re", res1, tmpmatA(i,j), tmparrA(i))
				if (i.eq.j) then
					call assert_double(trim(tmparg)//"im", res2, tmpmatB(i,j), 1d-7)
				else
					call assert_double_rel(trim(tmparg)//"im", res2, tmpmatB(i,j), tmparrS(i))
				end if
			end do
		end do
		call printTotalTests
		call resetTestCounter
	end subroutine do_test_coll_int

	pure real(dl) function fakecollintnue1(a, b, o, F_ab_ann, F_ab_sc)
		use variables
		integer, intent(in) :: a
		real(dl), intent(in) :: b
		type(coll_args), intent(in) :: o
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		fakecollintnue1=1.d0
	end function
	pure real(dl) function fakecollintnue0(a, b, o, F_ab_ann, F_ab_sc)
		use variables
		integer, intent(in) :: a
		real(dl), intent(in) :: b
		type(coll_args), intent(in) :: o
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		fakecollintnue0=0.d0
	end function
	pure real(dl) function fakecollintnuey(a, b, o, F_ab_ann, F_ab_sc)
		use variables
		integer, intent(in) :: a
		real(dl), intent(in) :: b
		type(coll_args), intent(in) :: o
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		fakecollintnuey=1.d4*b**2
	end function
	pure real(dl) function fakecollintnuefuncA(a, b, o, F_ab_ann, F_ab_sc)
		use variables
		integer, intent(in) :: a
		real(dl), intent(in) :: b
		type(coll_args), intent(in) :: o
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		fakecollintnuefuncA=o%z*o%z/o%x/o%x*nuDensMatVecFD(a)%re(o%ix1, o%ix1)*fermiDirac(b)/o%y1/o%y1*0.3*(y_arr(a)+b)
	end function
	pure real(dl) function fakecollintnuefuncB(a, b, o, F_ab_ann, F_ab_sc)
		use variables
		integer, intent(in) :: a
		real(dl), intent(in) :: b
		type(coll_args), intent(in) :: o
		procedure (F_annihilation) :: F_ab_ann
		procedure (F_scattering) :: F_ab_sc
		fakecollintnuefuncB=o%z*o%z/o%x/o%x*nuDensMatVecFD(a)%re(o%ix1, o%ix1)*fermiDirac(b)/o%y1/o%y1*0.3*(sqrt(y_arr(a))+sqrt(b))
	end function
	pure real(dl) function fakecollintnunu1(a, b, o, F_nu_sc, F_nu_pa)
		use variables
		integer, intent(in) :: a, b
		type(coll_args), intent(in) :: o
		procedure (Fnunu) :: F_nu_sc, F_nu_pa
		fakecollintnunu1=1.d0
	end function
	pure real(dl) function fakecollintnunu0(a, b, o, F_nu_sc, F_nu_pa)
		use variables
		integer, intent(in) :: a, b
		type(coll_args), intent(in) :: o
		procedure (Fnunu) :: F_nu_sc, F_nu_pa
		fakecollintnunu0=0.d0
	end function
	pure real(dl) function fakecollintnunuy(a, b, o, F_nu_sc, F_nu_pa)
		use variables
		integer, intent(in) :: a, b
		type(coll_args), intent(in) :: o
		procedure (Fnunu) :: F_nu_sc, F_nu_pa
		fakecollintnunuy=1.d4*b**2
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
		call printTestBlockName("Collision_terms")
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		collArgs%x = x
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%y2 = 0.d0
		collArgs%y3 = 0.d0
		collArgs%y4 = 0.d0
		collArgs%dme2 = dme2

		!fake function
		open(unit=fu, file="test_outputs/collterms_fake_x22_re.dat", status="old")
		open(unit=fv, file="test_outputs/collterms_fake_x22_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		tmparrA(:) = (/1d-4,1d-4,1d-4/)
		tmparrS(:) = (/1d-4,1d-4,1d-4/)
		cts = get_collision_terms(collArgs, fakecollintnuey, fakecollintnunu0)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('collision_terms fake x22 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7, tmparrS(i))
			end do
		end do
		open(unit=fu, file="test_outputs/collterms_fake_x22_zerooffdiag_re.dat", status="old")
		open(unit=fv, file="test_outputs/collterms_fake_x22_zerooffdiag_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		collint_damping_type = 0
		collint_offdiag_damping = .true.
		cts = get_collision_terms(collArgs, fakecollintnuey, fakecollintnunu0)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('collision_terms fake x22 zero off diagonal ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), 1d-7, tmparrA(i))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7, tmparrS(i))
			end do
		end do
		open(unit=fu, file="test_outputs/collterms_fake_x22_zerodiag_re.dat", status="old")
		open(unit=fv, file="test_outputs/collterms_fake_x22_zerodiag_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		collint_damping_type = 2
		collint_diagonal_zero = .true.
		collint_offdiag_damping = .false.
		cts = get_collision_terms(collArgs, fakecollintnuey, fakecollintnunu0)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('collision_terms zd ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), 1d-7, tmparrA(i))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7, tmparrS(i))
			end do
		end do
		collint_diagonal_zero = .false.

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

		collArgs%ix1=1
		collArgs%ix2=1
		open(unit=fu, file="test_outputs/collterms_real_diag_re.dat", status="old")
		open(unit=fv, file="test_outputs/collterms_real_diag_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		tmparrA(:) = (/2d-4,5d-4,5d-4/)
		tmparrS(:) = (/1d-5,1d-5,1d-5/)
		cts = get_collision_terms(collArgs, coll_nue_int, fakecollintnunu0)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('collision_terms real diag ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), 1d-7, tmparrA(i))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7, tmparrS(i))
			end do
		end do

!		x=0.75d0
!		iy1 = 7 !1.22151515151515
!!		z=1.06d0
!		dme2=0.1d0
!!		iy1 = 67 !13.3366666666667
!		z = 1.186d0
!		collArgs%ix1 = 1
!		collArgs%ix1 = 1
!		collArgs%x = x
!		collArgs%z = z
!		collArgs%iy = iy1
!		collArgs%y1 = y_arr(iy1)
!		collArgs%y2 = 0.d0
!		collArgs%y3 = 0.d0
!		collArgs%y4 = 0.d0
!		collArgs%dme2 = dme2
!		do iy=1, Ny
!			nuDensMatVecFD(iy)%re(:,:) = 0.d0
!			nuDensMatVecFD(iy)%im(:,:) = 0.d0
!			do i=1, flavorNumber
!				nuDensMatVecFD(iy)%re(i, i) = 1.d0 * fermiDirac(y_arr(iy))
!			end do
!		end do
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
!		collArgs%ix1=1
!		collArgs%ix2=1
!		open(unit=fu, file="test_outputs/collterms_real_full_re.dat", status="old")
!		open(unit=fv, file="test_outputs/collterms_real_full_im.dat", status="old")
!		do i=1, 3
!			read (fu, *) tmpmatA(i,:)
!			read (fv, *) tmpmatB(i,:)
!		end do
!		close(fu)
!		close(fv)
!		tmparrA(:) = (/0.08d0, 0.05d0, 0.05d0/)
!		tmparrS(:) = (/0.00001d0, 0.00001d0, 0.00001d0/)
!		cts = get_collision_terms(collArgs, coll_nue_int, fakecollintnunu0)
!		cts%re(:,:) = cts%re(:,:) * overallFactor
!		cts%im(:,:) = cts%im(:,:) * overallFactor
!		do i=1, flavorNumber
!			do j=1, flavorNumber
!				write(tmparg,"('collision_terms real diag ',2I1)") i,j
!				call assert_double_rel_safe(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), 1d-7, tmparrA(i))
!				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7, tmparrS(i))
!			end do
!		end do
		call printTotalTests
		call resetTestCounter
	end subroutine do_test_collision_terms

	subroutine do_test_drho_dx
		real(dl) :: x,z,fd, dme2, sqrtraddens
		type(cmplxMatNN) :: res, outp
		character(len=300) :: tmparg
		integer :: i, j, iy
		real(dl), dimension(:,:,:), allocatable :: GLR_vectmp

		ftqed_temperature_corr = .false.
		allocate(GLR_vectmp(2, flavorNumber, flavorNumber))
		GLR_vectmp = GLR_vec
		GLR_vec = 0.d0
		call printTestBlockName("d rho/d x [w/o real collision_terms]")
		call allocateCmplxMat(res)
		x = 0.06d0
		z = 1.23d0
		dme2 = 0.1d0
		call updateMatterDensities(x, z)
		nuDensities%re=0.d0
		nuDensities%im=0.d0
		iy = 12 !2.231111111111111
		do i=1, Ny
			fd = fermiDirac(y_arr(i))
			open(unit=fu, file="test_outputs/drhodx_rho_re.dat", status="old")
			open(unit=fv, file="test_outputs/drhodx_rho_im.dat", status="old")
			do j=1, 3
				read (fu, *) nuDensMatVecFD(i)%re(j,:)
				read (fv, *) nuDensMatVecFD(i)%im(j,:)
			end do
			close(fu)
			close(fv)
			do j=1,3
				nuDensMatVecFD(i)%re(j,j) = nuDensMatVecFD(i)%re(j,j)*fd
			end do
		end do
		sqrtraddens = sqrt(totalRadiationDensity(x,z))

		fd = fermiDirac(y_arr(iy))
		open(unit=fu, file="test_outputs/drhodx_A_re.dat", status="old")
		open(unit=fv, file="test_outputs/drhodx_A_im.dat", status="old")
		do j=1, 3
			read (fu, *) res%re(j,:)
			read (fv, *) res%im(j,:)
		end do
		close(fu)
		close(fv)
		do j=1,3
			res%re(j,j) = res%re(j,j)/fd
		end do
		call drhoy_dx_fullMat(outp, x, 1.d0, z, iy, dme2, sqrtraddens, fakecollintnue0, fakecollintnunu0)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('drho/dx A ',2I1)") i,j
				if ((i.eq.2 .and. j.eq.3) .or. (i.eq.3 .and. j.eq.2)) then
					call assert_double_rel(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-3)
				else
					call assert_double_rel_safe(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-7, 1d-4)
				end if
				if ((i.eq.2 .and. j.eq.3) .or. (i.eq.3 .and. j.eq.2)) then
					call assert_double_rel(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 2d-3)
				else
					call assert_double_rel_safe(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-7, 2d-4)
				end if
			end do
		end do

		fd = fermiDirac(y_arr(iy))
		open(unit=fu, file="test_outputs/drhodx_B_re.dat", status="old")
		open(unit=fv, file="test_outputs/drhodx_B_im.dat", status="old")
		do j=1, 3
			read (fu, *) res%re(j,:)
			read (fv, *) res%im(j,:)
		end do
		close(fu)
		close(fv)
		do j=1,3
			res%re(j,j) = res%re(j,j)/fd
		end do
		call drhoy_dx_fullMat(outp,x,1.d0, z,iy, dme2, sqrtraddens, fakecollintnue1, fakecollintnunu0)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('drho/dx B ',2I1)") i,j
				if ((i.eq.2 .and. j.eq.3) .or. (i.eq.3 .and. j.eq.2)) then
					call assert_double_rel(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 2d-3)
				else
					call assert_double_rel_safe(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-7, 1d-4)
				end if
				if ((i.eq.2 .and. j.eq.3) .or. (i.eq.3 .and. j.eq.2)) then
					call assert_double_rel(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1.5d-2)
				else
					call assert_double_rel_safe(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-7, 1d-4)
				end if
			end do
		end do

		x = 1.76d0
		z = 1.31d0
		call updateMatterDensities(x, z)
		nuDensities%re=0.d0
		nuDensities%im=0.d0
		iy = 34 !6.67333333333333
		do i=1, Ny
			fd = fermiDirac(y_arr(i))
			open(unit=fu, file="test_outputs/drhodx_rho_re.dat", status="old")
			open(unit=fv, file="test_outputs/drhodx_rho_im.dat", status="old")
			do j=1, 3
				read (fu, *) nuDensMatVecFD(i)%re(j,:)
				read (fv, *) nuDensMatVecFD(i)%im(j,:)
			end do
			close(fu)
			close(fv)
			do j=1,3
				nuDensMatVecFD(i)%re(j,j) = nuDensMatVecFD(i)%re(j,j)*fd
			end do
		end do
		sqrtraddens = sqrt(totalRadiationDensity(x,z))

		fd = fermiDirac(y_arr(iy))
		open(unit=fu, file="test_outputs/drhodx_C_re.dat", status="old")
		open(unit=fv, file="test_outputs/drhodx_C_im.dat", status="old")
		do j=1, 3
			read (fu, *) res%re(j,:)
			read (fv, *) res%im(j,:)
		end do
		close(fu)
		close(fv)
		do j=1,3
			res%re(j,j) = res%re(j,j)/fd
		end do
		call drhoy_dx_fullMat(outp,x,1.d0, z,iy, dme2, sqrtraddens, fakecollintnue0, fakecollintnunu0)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('drho/dx C ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-7, 1d-4)
				call assert_double_rel_safe(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-7, 1d-4)
			end do
		end do

		open(unit=fu, file="test_outputs/drhodx_D_re.dat", status="old")
		open(unit=fv, file="test_outputs/drhodx_D_im.dat", status="old")
		do j=1, 3
			read (fu, *) res%re(j,:)
			read (fv, *) res%im(j,:)
		end do
		close(fu)
		close(fv)
		do j=1,3
			res%re(j,j) = res%re(j,j)/fd
		end do
		call drhoy_dx_fullMat(outp,x,1.d0,z,iy, dme2, sqrtraddens, fakecollintnuey, fakecollintnunu0)
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('drho/dx D ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", outp%re(i,j), res%re(i,j), 1d-7, 1d-4)
				call assert_double_rel_safe(trim(tmparg)//"im", outp%im(i,j), res%im(i,j), 1d-7, 1d-4)
			end do
		end do
		GLR_vec = GLR_vectmp
		deallocate(GLR_vectmp)
		ftqed_temperature_corr = .true.
		call printTotalTests
		call resetTestCounter
	end subroutine do_test_drho_dx

	subroutine do_test_damping_factors
		real(dl) :: x,z,dme2
		integer :: i, j, ix, iy1, iy
		real(dl) :: y,res1,res2, fd
		type(coll_args) :: collArgs
		real(dl), dimension(3) :: tmparrS, tmparrA
		real(dl), dimension(3, 3) :: tmpmatA, tmpmatB
		character(len=300) :: tmparg
		type(cmplxMatNN) :: cts

		call allocateCmplxMat(cts)

		collint_damping_type = 2
		collint_diagonal_zero = .false.
		collint_offdiag_damping = .true.

		x = 0.75d0
		iy1 = 7 !1.22151515151515
		z = 1.186d0
		dme2 = 0.1d0
		call printTestBlockName("Damping terms - McKellar")
		collArgs%x = x
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%y2 = 0.d0
		collArgs%y3 = 0.d0
		collArgs%y4 = 0.d0
		collArgs%dme2 = dme2
		do iy=1, Ny
			fd = fermiDirac(y_arr(iy))
			open(unit=fu, file="test_outputs/damping_testrho3_re.dat", status="old")
			open(unit=fv, file="test_outputs/damping_testrho3_im.dat", status="old")
			do j=1, 3
				read (fu, *) nuDensMatVecFD(iy)%re(j,:)
				read (fv, *) nuDensMatVecFD(iy)%im(j,:)
			end do
			close(fu)
			close(fv)
			do j=1,3
				nuDensMatVecFD(iy)%re(j,j) = nuDensMatVecFD(iy)%re(j,j)*fd
			end do
		end do

		!2+1
		sterile(3) = .true.
		call setDampingFactors
		open(unit=fu, file="test_outputs/damping_2+1_re.dat", status="old")
		open(unit=fv, file="test_outputs/damping_2+1_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		tmpmatA = tmpmatA * z**4
		tmpmatB = tmpmatB * z**4
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollintnuey, fakecollintnunu0)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=i+1, flavorNumber
				write(tmparg,"('damping term 2+1 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7, tmparrS(i))
			end do
		end do

!		!3+0
		sterile(3) = .false.
		call setDampingFactors
		open(unit=fu, file="test_outputs/damping_3+0_re.dat", status="old")
		open(unit=fv, file="test_outputs/damping_3+0_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		tmpmatA = tmpmatA * z**4
		tmpmatB = tmpmatB * z**4
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollintnuey, fakecollintnunu0)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=i+1, flavorNumber
				write(tmparg,"('damping term 3+0 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7, tmparrS(i))
			end do
		end do

		flavorNumber = 2
		flavNumSqu = flavorNumber**2
		call deallocateStuff
		call allocateStuff
		!1+1
		sterile(2) = .true.
		call setDampingFactors
		open(unit=fu, file="test_outputs/damping_1+1_re.dat", status="old")
		open(unit=fv, file="test_outputs/damping_1+1_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		tmpmatA = tmpmatA * z**4
		tmpmatB = tmpmatB * z**4
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollintnuey, fakecollintnunu0)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=i+1, flavorNumber
				write(tmparg,"('damping term 1+1 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7, tmparrS(i))
			end do
		end do

		!2+0
		sterile(2) = .false.
		call setDampingFactors
		open(unit=fu, file="test_outputs/damping_2+0_re.dat", status="old")
		open(unit=fv, file="test_outputs/damping_2+0_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		tmpmatA = tmpmatA * z**4
		tmpmatB = tmpmatB * z**4
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollintnuey, fakecollintnunu0)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=i+1, flavorNumber
				write(tmparg,"('damping term 2+0 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i,j), tmparrA(i))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i,j), 1d-7, tmparrS(i))
			end do
		end do

		flavorNumber = 4
		flavNumSqu = flavorNumber**2
		call deallocateStuff
		call allocateStuff
		call setDampingFactors
		do ix=1, Ny
			deallocate(nuDensMatVecFD(ix)%re, nuDensMatVecFD(ix)%im)
			allocate(nuDensMatVecFD(ix)%re(flavorNumber,flavorNumber), nuDensMatVecFD(ix)%im(flavorNumber,flavorNumber))
		end do
		do iy=1, Ny
			fd = fermiDirac(y_arr(iy))
			open(unit=fu, file="test_outputs/damping_testrho4_re.dat", status="old")
			open(unit=fv, file="test_outputs/damping_testrho4_im.dat", status="old")
			do j=1, 4
				read (fu, *) nuDensMatVecFD(iy)%re(j,:)
				read (fv, *) nuDensMatVecFD(iy)%im(j,:)
			end do
			close(fu)
			close(fv)
			do j=1, 4
				nuDensMatVecFD(iy)%re(j,j) = nuDensMatVecFD(iy)%re(j,j)*fd
			end do
		end do
		!3+1
		sterile = .false.
		sterile(4) = .true.
		call setDampingFactors
		open(unit=fu, file="test_outputs/damping_3+1_re.dat", status="old")
		open(unit=fv, file="test_outputs/damping_3+1_im.dat", status="old")
		do i=1, 3
			read (fu, *) tmpmatA(i,:)
			read (fv, *) tmpmatB(i,:)
		end do
		close(fu)
		close(fv)
		tmpmatA = tmpmatA * z**4
		tmpmatB = tmpmatB * z**4
		tmparrA(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		tmparrS(:) = (/0.0001d0, 0.0001d0, 0.0001d0/)
		cts = get_collision_terms(collArgs, fakecollintnuey, fakecollintnunu0)
		cts%re(:,:) = cts%re(:,:) * overallFactor
		cts%im(:,:) = cts%im(:,:) * overallFactor
		do i=1, flavorNumber
			do j=i+1, flavorNumber
				write(tmparg,"('damping term 3+1 ',2I1)") i,j
				call assert_double_rel(trim(tmparg)//"re", cts%re(i,j), tmpmatA(i, j-1), tmparrA(i))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), tmpmatB(i, j-1), 1d-7, tmparrS(i))
			end do
		end do

		!restore initial settings
		flavorNumber = 3
		flavNumSqu = flavorNumber**2
		call deallocateStuff
		call allocateStuff
		call setDampingFactors
		do ix=1, Ny
			deallocate(nuDensMatVecFD(ix)%re, nuDensMatVecFD(ix)%im)
			allocate(nuDensMatVecFD(ix)%re(flavorNumber,flavorNumber), nuDensMatVecFD(ix)%im(flavorNumber,flavorNumber))
		end do
		collint_damping_type = 2
		collint_diagonal_zero = .false.
		collint_offdiag_damping = .false.
		sterile = .false.
		call printTotalTests
		call resetTestCounter
	end subroutine do_test_damping_factors

	subroutine do_test_GL
		real(dl), dimension(:), allocatable :: xa, wa, wa2, fx1, dx, ya
		real(dl), dimension(:,:), allocatable :: fx2a, fx2b
		real(dl) :: inta, intb, tol, r
		integer :: nix, nx, ix, iy, i, n, m
		character(len=300) :: tmparg
		integer, parameter :: nn=17
		real(dl), dimension(nn,3) :: tmparrA, tmparrB, tmperrA, tmperrB
		type(coll_args) :: collArgs
		real(dl), dimension(:), allocatable :: ydot
		real(dl) :: x, z, rn
		real(dl), dimension(9) :: ve1, ve2
		character(len=100) :: tmpstr

		call printTestBlockName("Gauss-Laguerre quadrature")

		use_gauss_laguerre = .true.
		open(unit=fu, file="test_outputs/gl_test_func.dat", status="old")
		read (fu, *) r
		close(fu)
		do nx=90, 10, -1
			call get_GLq_vectors(nx, xa, wa, wa2, .false., 3, 20.d0)
			call finish_y_arrays

			allocate(fx1(nx))
			do ix=1,nx
				fx1(ix) = fermiDirac(xa(ix))
			end do
			inta = integral_GL_1d(wa, fx1)
			write(tmparg, "('test GL quadrature on Fermi-Dirac, nx=',I2)") nx
			call assert_double_rel(trim(tmparg), inta, PISQ*PISQD15*7.d0/8.0, 3d-6)
			deallocate(fx1)

			allocate(fx2a(nx, nx), fx2b(nx, nx), ya(nx), dx(nx))
			ya = linspace(0.01d0, 20.d0, nx)
			do ix=1, nx-1
				dx(ix) = ya(ix+1) - ya(ix)
			end do
			do ix=1,nx
				do iy=1,nx
					fx2a(ix, iy) = fermiDirac(xa(ix))*fermiDirac(xa(iy))/(xa(ix)**2 * xa(iy))
!					fx2b(ix, iy) = fermiDirac(ya(ix))*fermiDirac(ya(iy))*(ya(ix) * ya(iy)**2)
				end do
			end do
			inta = integral_GL_2d(nx, wa, wa, fx2a)
!			intb = integral_NC_2d(nx, nx, dx, dx, fx2b)
			if (nx.gt.77) then
				tol=1d-5
			elseif (nx.gt.42) then
				tol=1d-4
			elseif (nx.gt.22) then
				tol=1d-3
			elseif (nx.gt.16) then
				tol=3d-3
			elseif (nx.gt.11) then
				tol=1d-2
			else
				tol=3d-2
			end if
			write(tmparg, "('test GL quadrature 2D on Fermi-Dirac, nx=',I2)") nx
			call assert_double_rel(trim(tmparg), inta, r, tol)
!			print *,nx, inta, intb
			deallocate(fx2a, fx2b, xa, wa, wa2, ya, dx)
		end do

		call printTestBlockName("GL quadrature of collision integrals")
		collArgs%x = 0.05d0
		collArgs%z = 1.06d0
		collArgs%y2 = 0.d0
		collArgs%y3 = 0.d0
		collArgs%y4 = 0.d0
		collArgs%dme2 = 0.1d0
		tmparrA(:,:) = 0.d0
		tmparrB(:,:) = 0.d0
		print*,"writing list of nix, Ny, y_arr(iy) with iy=int(Ny/3)"
		do nix=1, nn
			nx = nix*5+5
			call get_GLq_vectors(nx, xa, wa, wa2, .false., 3, 20.d0)
			call finish_y_arrays
			collArgs%iy = int(nx/3)
			print*,nix,nx, xa(collArgs%iy)
		end do

		tmperrA=1d-5
		tmperrB=1d-5
		!Ny=20
		tmperrA(3,:) = (/0.07,0.07,0.07/)
		tmperrB(3,:) = (/0.07,0.07,0.07/)
		!Ny=25
		tmperrA(4,:) = (/0.05,0.05,0.05/)
		tmperrB(4,:) = (/0.05,0.05,0.05/)
		!Ny=30
		tmperrA(5,:) = (/0.03,0.03,0.03/)
		tmperrB(5,:) = (/0.03,0.03,0.03/)
		!Ny=35
		tmperrA(6,:) = (/0.02,0.02,0.02/)
		tmperrB(6,:) = (/0.02,0.02,0.02/)
		!Ny=40
		tmperrA(7,:) = (/0.012,0.012,0.012/)
		tmperrB(7,:) = (/0.012,0.012,0.012/)
		!Ny=45
		tmperrA(8,:) = (/0.007,0.007,0.007/)
		tmperrB(8,:) = (/0.007,0.007,0.007/)
		!Ny=50
		tmperrA(9,:) = (/0.005,0.005,0.005/)
		tmperrB(9,:) = (/0.005,0.005,0.005/)
		!Ny=55
		tmperrA(10,:) = (/0.0025,0.0025,0.0025/)
		tmperrB(10,:) = (/0.0025,0.0025,0.0025/)
		!Ny=60
		tmperrA(11,:) = (/5e-4,5e-4,5e-4/)
		tmperrB(11,:) = (/5e-4,5e-4,5e-4/)
		!Ny=65
		tmperrA(12,:) = (/8d-4,8d-4,8d-4/)
		tmperrB(12,:) = (/8d-4,8d-4,8d-4/)
		!Ny=70
		tmperrA(13,:) = (/2d-3,2d-3,2d-3/)
		tmperrB(13,:) = (/2d-3,2d-3,2d-3/)
		!Ny=75
		tmperrA(14,:) = (/3d-3,3d-3,3d-3/)
		tmperrB(14,:) = (/3d-3,3d-3,3d-3/)
		!Ny=80
		tmperrA(15,:) = (/4d-3,4d-3,4d-3/)
		tmperrB(15,:) = (/4d-3,4d-3,4d-3/)
		!Ny=85
		tmperrA(16,:) = (/4d-3,4d-3,4d-3/)
		tmperrB(16,:) = (/4d-3,4d-3,4d-3/)
		!Ny=90
		tmperrA(17,:) = (/5d-3,5d-3,5d-3/)
		tmperrB(17,:) = (/5d-3,5d-3,5d-3/)
		open(unit=fu, file="test_outputs/GL_sc_int.dat", status="old")
		open(unit=fv, file="test_outputs/GL_ann_int.dat", status="old")
		do nix=1, nn
			read (fu, *) Ny, r, tmparrA(nix,:)
			if (Ny .ne. nix*5+5) &
				call criticalError("Ny does not match. The file 'GL_sc_int.dat' cannot be read")
			read (fv, *) Ny, r, tmparrB(nix,:)
			if (Ny .ne. nix*5+5) &
				call criticalError("Ny does not match. The file 'GL_ann_int.dat' cannot be read")
		end do
		close(fu)
		close(fv)
		do nix=nn,3,-1
			Ny=nix*5+5
			call get_GLq_vectors(Ny, y_arr, w_gl_arr, w_gl_arr2, .false., 3, 20.d0)
			call finish_y_arrays
!			do ix=1, Ny-1
!				dy_arr(ix) = y_arr(ix+1) - y_arr(ix)
!			end do
			collArgs%iy = int(Ny/3)
			collArgs%y1 = y_arr(collArgs%iy)
			write(*,*)Ny,y_arr(collArgs%iy)

			do iy=1, Ny
				nuDensMatVecFD(iy)%re = 0.d0
				nuDensMatVecFD(iy)%im = 0.d0
!				nuDensMatVecFD(iy)%y = y_arr(iy)
			end do
			do ix=1, flavorNumber
				do iy=1, Ny
					nuDensMatVecFD(iy)%re(ix, ix) = 1.d0*ix * fermiDirac(y_arr(iy))
				end do
!				nuDensMatVecFD(collArgs%iy)%re(ix, ix) = 1.d0
			end do

!			write(*,*) ""
			do ix=1, flavorNumber
				collArgs%ix1 = ix
				collArgs%ix2 = ix
				inta = integrate_collint_nue_GL(fakecollintnuefuncB, collArgs, F_ab_ann_re, F_ab_sc_re)
!				print*, integrate_collint_nue_GL(coll_nue_sc_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
				write(tmparg,"('test coll sc GL, N=',I2,' - ',2I1)") Ny, ix, ix
				call assert_double_rel_safe_verb(trim(tmparg), inta, tmparrA(nix,ix), 1d-30, tmperrA(nix,ix))
			end do

!			write(*,*) ""
			do ix=1, flavorNumber
				collArgs%ix1 = ix
				collArgs%ix2 = ix
				inta = integrate_collint_nue_GL(fakecollintnuefuncA, collArgs, F_ab_ann_re, F_ab_sc_re)
!				print*, integrate_collint_nue_GL(coll_nue_ann_int_w, collArgs, F_ab_ann_re, F_ab_sc_re)
				write(tmparg,"('test coll ann GL, N=',I2,' - ',2I1)") Ny, ix, ix
				call assert_double_rel_safe_verb(trim(tmparg), inta, tmparrB(nix,ix), 1d-30, tmperrB(nix,ix))
			end do
		end do

		call printTestBlockName("other applications of GL quadrature")
		Ny=50
		call get_GLq_vectors(Ny, y_arr, w_gl_arr, w_gl_arr2, .false., 3, 20.d0)
		call finish_y_arrays

		n=5
		ve2=(/1.,1.,1.,2.,3.,0.,0.,0.,0./)
		ve1=(/1d-4,1d-4,2d-4,5d-4,1d-4,0.d0,0.d0,0.d0,0.d0/)
		do ix=1, flavorNumber
			do iy=1, Ny
				nuDensMatVecFD(iy)%re(ix, ix) = 1.d0*ix * fermiDirac(y_arr(iy))
			end do
		end do
		open(unit=fu, file="test_outputs/nuDens.dat", status="old")
		do i=1, n
			read (fu, *) z, r, rn
			write(tmpstr, "(I1)") i
			call assert_double_rel("nuDensGL test "//trim(tmpstr), nuDensity(z,int(ve2(i))), ve2(i)*r, ve1(i))
		end do
		close(fu)
		open(unit=fu, file="test_outputs/nuDensEq.dat", status="old")
		do i=1, n
			read (fu, *) z, r, rn
			do ix=1, flavorNumber
				do iy=1, Ny
					nuDensMatVecFD(iy)%re(ix, ix) = 1.d0*ix * fermiDirac(y_arr(iy) / z)
				end do
			end do
			write(tmpstr, "(I1)") i
			call assert_double_rel("nuDensGL test "//trim(tmpstr), nuDensityGL(int(ve2(i)), int(ve2(i))), ve2(i)*r, ve1(i))
			call assert_double_rel("nuNumDensGL test "//trim(tmpstr), nuNumberDensityGL(int(ve2(i)), int(ve2(i))), ve2(i)*rn, ve1(i))
		end do
		close(fu)

		n=ntot
		allocate(ydot(n))
		nuFactor=0.d0
		nuFactor(1:3)=1.d0
		Ny=25
		call get_GLq_vectors(Ny, y_arr, w_gl_arr, w_gl_arr2, .false., 3, 20.d0)
		call finish_y_arrays
		ydot = 0.d0
		do m=1, Ny
			ydot((m-1)*flavNumSqu + 1) = 1.d0/y_arr(m)
			ydot((m-1)*flavNumSqu + 2) = y_arr(m)/20.d0
			ydot((m-1)*flavNumSqu + 3) = 1.d0
		end do

		feq_arr = 0.d0
		do ix=1, Ny
			feq_arr(ix) = fermiDirac(y_arr(ix))
		end do
		open(unit=fu, file="test_outputs/dzodx_g_A.dat", status="old")
		open(unit=fv, file="test_outputs/dzodx_n_A.dat", status="old")
		ve1=(/5d-5,6d-6,6d-6,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
		do ix=1,3
			read (fu, *) x,z,r
			read (fv, *) x,z,rn
			write(tmparg,"(I1)") ix
			call dz_o_dx(x, 1.2d0, z, ydot, n)
			call assert_double_rel("dz_o_dx GL test "//trim(tmparg), ydot(n), r, ve1(ix))
			call assert_double_rel("dw_o_dx GL test "//trim(tmparg), ydot(n-1), rn, ve1(ix))
		end do
		close(fu)
		close(fv)
		deallocate(ydot)
		call printTotalTests
		call resetTestCounter
	end subroutine do_test_GL

	subroutine do_test_zin
		call printTestBlockName("z_in solver")
		x_in=0.05d0
		z_in=0.d0
		call zin_solver
		call assert_double("z_in test 1", z_in-1.d0, 0.09788d0, 1d-4)
		x_in=1d-3
		call zin_solver
		call assert_double("z_in test 2", z_in-1.d0, 0.29017d-03, 1d-4)
		call printTotalTests
		call resetTestCounter
	end subroutine do_test_zin

	subroutine do_test_matterPotential
		type(cmplxMatNN) :: Heff
		complex(dl), dimension(maxFlavorNumber, maxFlavorNumber) :: Heffc
		character(len=300) :: tmparg
		real(dl), dimension(3,3) :: r1, r2
		integer :: i,j

		call printTestBlockName("matter potential, including neutrinos")

		!A
		open(unit=fu, file="test_outputs/matterpot_nmm.dat", status="old")
		open(unit=fv, file="test_outputs/matterpot_ldm.dat", status="old")
		open(unit=fw, file="test_outputs/matterpot_ndm1.dat", status="old")
		do j=1, 3
			read (fu, *) nuMassesMat(j,:)
			read (fv, *) leptonDensities(j,:)
			read (fw, *) nuDensities%re(j,:)
		end do
		close(fu)
		close(fv)
		close(fw)
		nuDensities%im = 0.d0
		Heff = H_eff(0.04d0)
		open(unit=fu, file="test_outputs/matterpot_t1.dat", status="old")
		do j=1, 3
			read (fu, *) r1(j,:)
		end do
		close(fu)
		r2 = 0.d0
		do i=1, flavorNumber
			do j=i, flavorNumber
				write(tmparg,"('H_eff A ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", Heff%re(i,j), r1(i,j), 1d-7, 1d-7)
				call assert_double(trim(tmparg)//"re", Heff%im(i,j), r2(i,j), 1d-7)
			end do
		end do
		Heffc = H_eff_cmplx(0.04d0)
		do i=1, flavorNumber
			do j=i, flavorNumber
				write(tmparg,"('H_eff_cmplx A ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", dble(Heffc(i,j)), r1(i,j), 1d-7, 1d-7)
				call assert_double_rel_safe(trim(tmparg)//"im", dimag(Heffc(i,j)), r2(i,j), 1d-7, 1d-7)
			end do
		end do
		!B
		open(unit=fu, file="test_outputs/matterpot_ndm2_re.dat", status="old")
		open(unit=fv, file="test_outputs/matterpot_ndm2_im.dat", status="old")
		do j=1, 3
			read (fu, *) nuDensities%re(j,:)
			read (fv, *) nuDensities%im(j,:)
		end do
		close(fu)
		close(fv)
		Heff = H_eff(0.7d0)
		open(unit=fu, file="test_outputs/matterpot_t2_re.dat", status="old")
		open(unit=fv, file="test_outputs/matterpot_t2_im.dat", status="old")
		do j=1, 3
			read (fu, *) r1(j,:)
			read (fv, *) r2(j,:)
		end do
		close(fu)
		close(fv)
		do i=1, flavorNumber
			do j=i, flavorNumber
				write(tmparg,"('H_eff B ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", Heff%re(i,j), r1(i,j), 1d-7, 1d-7)
				call assert_double_rel_safe(trim(tmparg)//"im", Heff%im(i,j), r2(i,j), 1d-7, 1d-7)
			end do
		end do
		Heffc = H_eff_cmplx(0.7d0)
		do i=1, flavorNumber
			do j=i, flavorNumber
				write(tmparg,"('H_eff_cmplx B ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", dble(Heffc(i,j)), r1(i,j), 1d-7, 1d-7)
				call assert_double_rel_safe(trim(tmparg)//"im", dimag(Heffc(i,j)), r2(i,j), 1d-7, 1d-7)
			end do
		end do

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_matterPotential

	subroutine do_test_diagonalization
		type(cmplxMatNN) :: m
		character(len=300) :: tmparg
		complex(dl), dimension(maxFlavorNumber, maxFlavorNumber) :: cmp1, cmp2
		real(dl), dimension(maxFlavorNumber) :: tmpvec, rv
		real(dl), dimension(3,3) :: r1, r2
		integer :: iy,i,j

		call printTestBlockName("diagonalization")

		open(unit=fu, file="test_outputs/diagonalization_nmm.dat", status="old")
		open(unit=fv, file="test_outputs/diagonalization_ldm.dat", status="old")
		do j=1, 3
			read (fu, *) nuMassesMat(j,:)
			read (fv, *) leptonDensities(j,:)
		end do
		close(fu)
		close(fv)
		open(unit=fu, file="test_outputs/diagonalization_ndm_re.dat", status="old")
		open(unit=fv, file="test_outputs/diagonalization_ndm_im.dat", status="old")
		do j=1, 3
			read (fu, *) nuDensities%re(j,:)
			read (fv, *) nuDensities%im(j,:)
		end do
		close(fu)
		close(fv)
		y_arr(1) = 1.

		tmpvec = 0.d0
		cmp2(:,:) = cmplx(0.d0, 0.d0)
		cmp1 = H_eff_cmplx(y_arr(1))
		call HEigensystem(flavorNumber, cmp1, flavorNumber, tmpvec, cmp2, flavorNumber, 0)
		open(unit=fu, file="test_outputs/diagonalization_eigenvalues.dat", status="old")
		read (fu, *) rv(:)
		close(fu)
		do i=1, maxFlavorNumber
			write(tmparg,"('HEigensystem ',I1)") i
			call assert_double_rel_safe(trim(tmparg)//"re", tmpvec(i), rv(i), 1d-7, 1d-6)
		end do

		do iy=1, Ny
		open(unit=fu, file="test_outputs/diagonalization_nm_re.dat", status="old")
		open(unit=fv, file="test_outputs/diagonalization_nm_im.dat", status="old")
		do j=1, 3
			read (fu, *) nuDensMatVecFD(iy)%re(j,:)
			read (fv, *) nuDensMatVecFD(iy)%im(j,:)
		end do
		close(fu)
		close(fv)
		end do

		m = rho_diag_mass(1)
		open(unit=fu, file="test_outputs/diagonalization_rho_mass_basis_re.dat", status="old")
		open(unit=fv, file="test_outputs/diagonalization_rho_mass_basis_im.dat", status="old")
		do j=1, 3
			read (fu, *) r1(j,:)
			read (fv, *) r2(j,:)
		end do
		close(fu)
		close(fv)
		do i=1, flavorNumber
			do j=i, flavorNumber
				write(tmparg,"('rho mass basis ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", m%re(i,j), r1(i,j), 1d-7, 1d-5)
				call assert_double(trim(tmparg)//"im", m%im(i,j), r2(i,j), 1d-7)
			end do
		end do

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_diagonalization

	subroutine do_test_damping_bennett
		real(dl) :: x,w,z,dme2,y1
		integer :: ix, iy1, iy, j
		real(dl) :: y,res1,res2
		type(coll_args) :: collArgs
		real(dl), dimension(3) :: tmparrS, tmparrA
		real(dl), dimension(3, 3) :: tmpmatA, tmpmatB
		real(dl), dimension(12) :: errv
		character(len=300) :: tmparg
		type(cmplxMatNN) :: cts

		call allocateCmplxMat(cts)

		call printTestBlockName("damping factors a la Bennett:2020zkv")

		Ny=50
		call get_GLq_vectors(Ny, y_arr, w_gl_arr, w_gl_arr2, .false., 3, 20.d0)
		call finish_y_arrays
		errv=(/2d-2,2d-2,2d-2,2d-2,2d-2,2d-2,1d-2,1d-2,1d-2,7d-2,0.d0,0.d0/)
		open(unit=fu, file="test_outputs/damping_bennett_sv.dat", status="old")
		open(unit=fv, file="test_outputs/damping_bennett_dy.dat", status="old")
		do j=1, 10
			read (fu, *) res1
			read (fv, *) res2
			write(tmparg, "('dy_damping_pi ',E10.3)") res1
			call assert_double_rel(trim(tmparg), dy_damping_pi(res1), res2, errv(j))
		end do
		close(fu)
		close(fv)

		Ny=100
		deallocate(y_arr)
		allocate(y_arr(Ny))
		y_arr = linspace(y_min, y_max, Ny)
		errv=(/1d-3,1d-3,2d-3,1d-3,1d-3,1d-3,1d-3,1d-3,3d-3,1d-3,1d-3,2d-3/)
		open(unit=fu, file="test_outputs/damping_bennett_sv.dat", status="old")
		open(unit=fv, file="test_outputs/damping_bennett_dy.dat", status="old")
		do j=1, 12
			read (fu, *) res1
			read (fv, *) res2
			write(tmparg, "('dy_damping_fit ',E10.3)") res1
			call assert_double_rel(trim(tmparg), dy_damping_fit(res1), res2, errv(j))
		end do
		close(fu)
		close(fv)

		call assert_double_rel("kappa A", -15.4485396d0, kappa_damp(12.d0, 0.44d0, 0.33d0), 1d-7)
		call assert_double_rel("kappa B", -0.12818209d0, kappa_damp(1.d0, 0.14d0, 0.01d0), 1d-7)

		call assert_double_rel("nunu_damp_integrand A", 3.948450853637d-6, nunu_damp_integrand(0.01d0, 3.d0, 5.d0), 1d-7)
		call assert_double_rel("nunu_damp_integrand B", -3.34556795185d0,  nunu_damp_integrand(0.01d0, 13.d0, 5.d0), 1d-7)
		call assert_double_rel("nunu_damp_integrand C", -37.225769332d0,   nunu_damp_integrand(0.01d0, 3.d0, 1.d0), 1d-7)
		call assert_double_rel("nunu_damp_integrand D", 366.6887058d0,     nunu_damp_integrand(10.d0, 3.d0, 5.d0), 1d-7)
		call assert_double_rel("nunu_damp_integrand E", -0.080526245965d0, nunu_damp_integrand(10.d0, 30.d0, 5.d0), 1d-7)
		call assert_double_rel("nunu_damp_integrand F", 328.64337311d0,    nunu_damp_integrand(10.d0, 3.d0, 1.d0), 1d-7)

		errv=(/2d-2,2d-2,2d-2,2d-2,2d-2,2d-2,1d-2,1d-2,1d-2,3d-2,2d-2,2d-2/)
		open(unit=fu, file="test_outputs/damping_bennett_sv.dat", status="old")
		open(unit=fv, file="test_outputs/damping_bennett_dy.dat", status="old")
		do j=1, 12
			read (fu, *) res1
			read (fv, *) res2
			write(tmparg, "('dy_damping ',E10.3)") res1
			call assert_double_rel(trim(tmparg), dy_damping(res1), res2, errv(j))
		end do
		close(fu)
		close(fv)

		collint_damping_type = 1
		collint_diagonal_zero = .false.
		collint_offdiag_damping = .true.
		call setDampingFactors

		call assert_double_rel("c Nue  1,2", dampTermMatrixCoeffNue (1,2), 0.357d0, 1d-2)
		call assert_double_rel("c Nue  1,3", dampTermMatrixCoeffNue (1,3), 0.357d0, 1d-2)
		call assert_double_rel("c Nue  2,3", dampTermMatrixCoeffNue (2,3), 0.1257d0, 1d-2)
		call assert_double_rel("c Nunu 1,2", dampTermMatrixCoeffNunu(1,2), 1.d0, 1d-2)
		call assert_double_rel("c Nunu 1,3", dampTermMatrixCoeffNunu(1,3), 1.d0, 1d-2)
		call assert_double_rel("c Nunu 2,3", dampTermMatrixCoeffNunu(2,3), 1.d0, 1d-2)

		! tests for comparing with complete terms
		x = 0.75d0
		iy1 = 7 !1.22151515151515
		y1 = y_arr(iy1)
		z = 1.186d0
		dme2 = 0.1d0
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
			nuDensMatVecFD(iy)%re(1,:) = (/1.1d0*fermiDirac(y), 10.d0, 33.d0/)
			nuDensMatVecFD(iy)%re(2,:) = (/10.d0, 1.1d0*fermiDirac(y), 46.d0/)
			nuDensMatVecFD(iy)%re(3,:) = (/33.d0, 46.d0, 1.1d0*fermiDirac(y)/)
			nuDensMatVecFD(iy)%im(1,:) = (/0.d0, -0.001d0, 0.003d0/)
			nuDensMatVecFD(iy)%im(2,:) = (/0.001d0, 0.d0, -0.002d0/)
			nuDensMatVecFD(iy)%im(3,:) = (/-0.003d0, 0.002d0, 0.d0/)
		end do

		cts = get_collision_terms(collArgs, fakecollintnuey, fakecollintnunu0)
		res1 = integrate_collint_nue_NC(fakecollintnuey, collArgs, F_ab_ann_re, F_ab_sc_re) &
			* collTermFactor/(y1**2 * x**4)
!		call printMat(cts%re)
		do ix=1, 3
			write(tmparg,"('damping Bennett:2020zkv d A',2I1)") ix,ix
			call assert_double_rel(trim(tmparg)//" re", cts%re(ix, ix), res1, 1d-4)
			do iy=ix+1, 3
				write(tmparg,"('damping Bennett:2020zkv od A',2I1)") ix,iy
				call assert_double_rel(trim(tmparg)//" re", cts%re(ix, iy), &
					-(dampTermMatrixCoeffNue(ix,iy)+dampTermMatrixCoeffNunu(ix,iy)) &
					* 2.d0 * dy_damping_fit(y1/z) * z**4 * y1**3 * nuDensMatVecFD(iy1)%re(ix,iy) &
					* collTermFactor/(y1**2 * x**4), &
					1d-4)
				call assert_double_rel(trim(tmparg)//" im", cts%im(ix, iy), &
					-(dampTermMatrixCoeffNue(ix,iy)+dampTermMatrixCoeffNunu(ix,iy)) &
					* 2.d0 * dy_damping_fit(y1/z) * z**4 * y1**3 * nuDensMatVecFD(iy1)%im(ix,iy) &
					* collTermFactor/(y1**2 * x**4), &
					1d-4)
			end do
		end do
		x = 4.75d0
		iy1 = 3
		y1 = y_arr(iy1)
		w = 1.234d0
		z = 1.386d0
		collArgs%x = x
		collArgs%w = w
		collArgs%z = z
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		cts = get_collision_terms(collArgs, fakecollintnuey, fakecollintnunu0)
		res1 = integrate_collint_nue_NC(fakecollintnuey, collArgs, F_ab_ann_re, F_ab_sc_re) &
			* collTermFactor/(y1**2 * x**4)
!		call printMat(cts%re)
		do ix=1, 3
			write(tmparg,"('damping Bennett:2020zkv d B',2I1)") ix,ix
			call assert_double_rel(trim(tmparg)//" re", cts%re(ix, ix), res1, 1d-4)
			do iy=ix+1, 3
				write(tmparg,"('damping Bennett:2020zkv od B',2I1)") ix,iy
				call assert_double_rel(trim(tmparg)//" re", cts%re(ix, iy), &
					-(dampTermMatrixCoeffNue(ix,iy)+dampTermMatrixCoeffNunu(ix,iy)) &
					* 2.d0 * dy_damping_fit(y1/z) * z**4 * y1**3 * nuDensMatVecFD(iy1)%re(ix,iy) &
					* collTermFactor/(y1**2 * x**4), &
					1d-4)
				call assert_double_rel(trim(tmparg)//" im", cts%im(ix, iy), &
					-(dampTermMatrixCoeffNue(ix,iy)+dampTermMatrixCoeffNunu(ix,iy)) &
					* 2.d0 * dy_damping_fit(y1/z) * z**4 * y1**3 * nuDensMatVecFD(iy1)%im(ix,iy) &
					* collTermFactor/(y1**2 * x**4), &
					1d-4)
			end do
		end do

		collint_damping_type = 2
		collint_diagonal_zero = .false.
		collint_offdiag_damping = .false.
		call setDampingFactors

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_damping_bennett

	subroutine do_test_F_nu
		integer :: i, j
		real(dl) :: res
		real(dl), dimension(3, 3) :: Fpr, Fpi, Fsr, Fsi
		character(len=300) :: tmparg
		type(cmplxMatNN) :: m1, m2, m3, m4

		call allocateCmplxMat(m1)
		call allocateCmplxMat(m2)
		call allocateCmplxMat(m3)
		call allocateCmplxMat(m4)

		call printTestBlockName("phase space functions for nunu")

#ifdef FULL_F_NU
		open(unit=fu, file="test_outputs/Fnunu_m1_re.dat", status="old")
		open(unit=fv, file="test_outputs/Fnunu_m1_im.dat", status="old")
		do j=1, 3
			read (fu, *) m1%re(j,:)
			read (fv, *) m1%im(j,:)
		end do
		close(fu)
		close(fv)
		open(unit=fu, file="test_outputs/Fnunu_m2_re.dat", status="old")
		open(unit=fv, file="test_outputs/Fnunu_m2_im.dat", status="old")
		do j=1, 3
			read (fu, *) m2%re(j,:)
			read (fv, *) m2%im(j,:)
		end do
		close(fu)
		close(fv)
		open(unit=fu, file="test_outputs/Fnunu_m3_re.dat", status="old")
		open(unit=fv, file="test_outputs/Fnunu_m3_im.dat", status="old")
		do j=1, 3
			read (fu, *) m3%re(j,:)
			read (fv, *) m3%im(j,:)
		end do
		close(fu)
		close(fv)
		open(unit=fu, file="test_outputs/Fnunu_m4_re.dat", status="old")
		open(unit=fv, file="test_outputs/Fnunu_m4_im.dat", status="old")
		do j=1, 3
			read (fu, *) m4%re(j,:)
			read (fv, *) m4%im(j,:)
		end do
		close(fu)
		close(fv)
		open(unit=fu, file="test_outputs/Fnunu_sc_full_re.dat", status="old")
		open(unit=fv, file="test_outputs/Fnunu_sc_full_im.dat", status="old")
		do j=1, 3
			read (fu, *) Fsr(j,:)
			read (fv, *) Fsi(j,:)
		end do
		close(fu)
		close(fv)
		open(unit=fu, file="test_outputs/Fnunu_pa_full_re.dat", status="old")
		open(unit=fv, file="test_outputs/Fnunu_pa_full_im.dat", status="old")
		do j=1, 3
			read (fu, *) Fpr(j,:)
			read (fv, *) Fpi(j,:)
		end do
		close(fu)
		close(fv)
#else
		open(unit=fu, file="test_outputs/Fnunu_d1.dat", status="old")
		open(unit=fv, file="test_outputs/Fnunu_d2.dat", status="old")
		do j=1, 3
			read (fu, *) m1%re(j,:)
			read (fv, *) m2%re(j,:)
		end do
		close(fu)
		close(fv)
		open(unit=fu, file="test_outputs/Fnunu_d3.dat", status="old")
		open(unit=fv, file="test_outputs/Fnunu_d4.dat", status="old")
		do j=1, 3
			read (fu, *) m3%re(j,:)
			read (fv, *) m4%re(j,:)
		end do
		close(fu)
		close(fv)
		open(unit=fu, file="test_outputs/Fnunu_sc_diag.dat", status="old")
		open(unit=fv, file="test_outputs/Fnunu_pa_diag.dat", status="old")
		do j=1, 3
			read (fu, *) Fsr(j,:)
			read (fv, *) Fpr(j,:)
		end do
		close(fu)
		close(fv)
		m1%im = 0.d0
		m2%im = 0.d0
		m3%im = 0.d0
		m4%im = 0.d0
		Fsi=0.d0
		Fpi=0.d0
#endif
		do i=1, 3
			do j=i,3
				write(tmparg,"('F_nu_sc ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", F_nu_sc_re(m1, m2, m3, m4, i, j), Fsr(i,j), 1d-7, 1d-4)
				call assert_double_rel_safe(trim(tmparg)//"im", F_nu_sc_im(m1, m2, m3, m4, i, j), Fsi(i,j), 1d-7, 1d-4)

				write(tmparg,"('F_nu_pa ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", F_nu_pa_re(m1, m2, m3, m4, i, j), Fpr(i,j), 1d-7, 1d-4)
				call assert_double_rel_safe(trim(tmparg)//"im", F_nu_pa_im(m1, m2, m3, m4, i, j), Fpi(i,j), 1d-7, 1d-4)
			end do
		end do

		call deallocateCmplxMat(m1)
		call deallocateCmplxMat(m2)
		call deallocateCmplxMat(m3)
		call deallocateCmplxMat(m4)

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_F_nu

	subroutine do_test_interp_nudens
		integer :: i, j
		real(dl), dimension(2, 2) :: ndr, ndi
		type(cmplxMatNN) :: nm
		type(cmplxMatNN), dimension(:), allocatable :: vdm
		character(len=300) :: tmparg
		character(len=1), dimension(3) :: fn

		call printTestBlockName("interpolation of nudens")
		allocate(vdm(3))
		open(unit=fu, file="test_outputs/interp_nunu_y.dat", status="old")
		read (fu, *) vdm(1)%y, vdm(2)%y, vdm(3)%y
		close(fu)
		fn=(/"A","B","C"/)
		do i=1, 3
			allocate(vdm(i)%re(2,2), vdm(i)%im(2,2))
			open(unit=fu, file="test_outputs/interp_nunu_"//fn(i)//"m_re.dat", status="old")
			open(unit=fv, file="test_outputs/interp_nunu_"//fn(i)//"m_im.dat", status="old")
			do j=1, 2
				read (fu, *) vdm(i)%re(j,:)
				read (fv, *) vdm(i)%im(j,:)
				vdm(i)%re(j,j) = vdm(i)%re(j,j)*fermiDirac(vdm(i)%y)
			end do
			close(fu)
			close(fv)
		end do
#ifdef RHO_OFFDIAG_INTERP_DIV_FD
		open(unit=fu, file="test_outputs/interp_nunu_t1_fd_re.dat", status="old")
		open(unit=fv, file="test_outputs/interp_nunu_t1_fd_im.dat", status="old")
#else
		open(unit=fu, file="test_outputs/interp_nunu_t1_re.dat", status="old")
		open(unit=fv, file="test_outputs/interp_nunu_t1_im.dat", status="old")
#endif
		do j=1, 2
			read (fu, *) ndr(j,:)
			read (fv, *) ndi(j,:)
			ndr(j,j) = ndr(j,j)*fermiDirac(1.5d0)
		end do
		close(fu)
		close(fv)
		nm = get_interpolated_nudens(vdm, 1.5d0, 2, 3)
		call assert_double("ndr A y", nm%y, 1.5d0, 1d-7)
		do i=1,2
			do j=1,2
				write(tmparg,"('ndr A ',2I1)") i,j
				call assert_double(trim(tmparg)//"re", nm%re(i,j), ndr(i,j), 1d-7)
				call assert_double(trim(tmparg)//"im", nm%im(i,j), ndi(i,j), 1d-7)
			end do
		end do
		call deallocateCmplxMat(nm)

#ifdef RHO_OFFDIAG_INTERP_DIV_FD
		open(unit=fu, file="test_outputs/interp_nunu_t2_fd_re.dat", status="old")
		open(unit=fv, file="test_outputs/interp_nunu_t2_fd_im.dat", status="old")
#else
		open(unit=fu, file="test_outputs/interp_nunu_t2_re.dat", status="old")
		open(unit=fv, file="test_outputs/interp_nunu_t2_im.dat", status="old")
#endif
		do j=1, 2
			read (fu, *) ndr(j,:)
			read (fv, *) ndi(j,:)
			ndr(j,j) = ndr(j,j)*fermiDirac(8.d0)
		end do
		close(fu)
		close(fv)
		nm = get_interpolated_nudens(vdm, 8.d0, 2, 3)
		call assert_double("ndr B y", nm%y, 8.d0, 1d-7)
		do i=1,2
			do j=1,2
				write(tmparg,"('ndr B ',2I1)") i,j
				call assert_double(trim(tmparg)//"re", nm%re(i,j), ndr(i,j), 1d-7)
				call assert_double(trim(tmparg)//"im", nm%im(i,j), ndi(i,j), 1d-7)
			end do
		end do
		call deallocateCmplxMat(nm)

		ndr=0.d0
		ndi=0.d0
		nm = get_interpolated_nudens(vdm, 0.1d0, 2, 3)
		call assert_double("ndr C y", nm%y, 0.1d0, 1d-7)
		do i=1,2
			do j=1,2
				write(tmparg,"('ndr C ',2I1)") i,j
				call assert_double(trim(tmparg)//"re", nm%re(i,j), ndr(i,j), 1d-7)
				call assert_double(trim(tmparg)//"im", nm%im(i,j), ndi(i,j), 1d-7)
			end do
		end do
		call deallocateCmplxMat(nm)

		ndr=0.d0
		ndi=0.d0
		nm = get_interpolated_nudens(vdm, 20.d0, 2, 3)
		call assert_double("ndr D y", nm%y, 20.d0, 1d-7)
		do i=1,2
			do j=1,2
				write(tmparg,"('ndr D ',2I1)") i,j
				call assert_double(trim(tmparg)//"re", nm%re(i,j), ndr(i,j), 1d-7)
				call assert_double(trim(tmparg)//"im", nm%im(i,j), ndi(i,j), 1d-7)
			end do
		end do
		call deallocateCmplxMat(nm)

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_interp_nudens

	subroutine do_test_collint_nunu
		integer :: i, j, iy1, iy2, iy3
		real(dl) :: y2, y3, y4, fsc, fpa, res1, res2, ey
		real(dl), dimension(3, 3) :: ndr, ndi, er, ei
		type(coll_args) :: collArgs
		type(cmplxMatNN) :: n4, cts
		real(dl), dimension(2) :: pi2_vec
		character(len=300) :: tmparg

		call printTestBlockName("collision integrals of nunu")
		collArgs%x = 0.05d0
		collArgs%z = 1.06d0
		collArgs%y2 = 0.d0
		collArgs%y3 = 0.d0
		collArgs%y4 = 0.d0
		collArgs%dme2 = 0.0d0

		Ny=50
		call get_GLq_vectors(Ny, y_arr, w_gl_arr, w_gl_arr2, .false., 3, 20.d0)
		call finish_y_arrays
		!print*,"y range (GL):",y_arr(1),y_arr(Ny)
		do j=1, Ny
			nuDensMatVecFD(j)%y=y_arr(j)
			nuDensMatVecFD(j)%re=0.d0
			nuDensMatVecFD(j)%im=0.d0
			nuDensMatVecFD(j)%re(1,1) = 1.0d0 * y_arr(j) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(2,2) = 1.1d0 * y_arr(j) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(3,3) = 0.9d0 * y_arr(j) * fermiDirac(y_arr(j))
		end do

		!A
		iy1=7
		iy2=2
		iy3=5
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		open(unit=fu, file="test_outputs/collint_nunu_int_A.dat", status="old")
		do j=1, 3
			read (fu, *) ndr(j,:)
		end do
		close(fu)
		ndi = 0.d0
		er = 2d-3
		ei = 1d-3
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				write(tmparg,"('nunu int A ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", coll_nunu_int(iy2, iy3, collArgs, F_nu_sc_re, F_nu_pa_re), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", coll_nunu_int(iy2, iy3, collArgs, F_nu_sc_im, F_nu_pa_im), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do

		!B
		iy1=10
		iy2=4
		iy3=7
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		open(unit=fu, file="test_outputs/collint_nunu_int_B.dat", status="old")
		do j=1, 3
			read (fu, *) ndr(j,:)
		end do
		close(fu)
		ndi = 0.d0
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				write(tmparg,"('nunu int B ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", coll_nunu_int(iy2, iy3, collArgs, F_nu_sc_re, F_nu_pa_re), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", coll_nunu_int(iy2, iy3, collArgs, F_nu_sc_im, F_nu_pa_im), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do

		!full integral
		!A
		iy1=10
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		open(unit=fu, file="test_outputs/collint_nunu_integral_GL_A.dat", status="old")
		do j=1, 3
			read (fu, *) ndr(j,:)
		end do
		close(fu)
		ndi = 0.d0
		er = 1d-3
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				write(tmparg,"('nunu integral A ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", integrate_collint_nunu_GL(coll_nunu_int, collArgs, F_nu_sc_re, F_nu_pa_re), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", integrate_collint_nunu_GL(coll_nunu_int, collArgs, F_nu_sc_im, F_nu_pa_im), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do

		!B
		iy1=5
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		open(unit=fu, file="test_outputs/collint_nunu_integral_GL_B.dat", status="old")
		do j=1, 3
			read (fu, *) ndr(j,:)
		end do
		close(fu)
		ndi = 0.d0
		er = 1.1d-3
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				write(tmparg,"('nunu integral B ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", integrate_collint_nunu_GL(coll_nunu_int, collArgs, F_nu_sc_re, F_nu_pa_re), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", integrate_collint_nunu_GL(coll_nunu_int, collArgs, F_nu_sc_im, F_nu_pa_im), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do

#ifdef FULL_F_NU
		!C
		iy1=7
		iy2=2
		iy3=5
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		do j=1, Ny
			nuDensMatVecFD(j)%y=y_arr(j)
			nuDensMatVecFD(j)%re=0.d0
			nuDensMatVecFD(j)%im=0.d0
			nuDensMatVecFD(j)%re(1, 1) = 1.d0 * y_arr(j) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(1, 2) = 0.1d0
			nuDensMatVecFD(j)%re(1, 3) = 0.d0
			nuDensMatVecFD(j)%re(2, 2) = 1.1d0 * y_arr(j) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(2, 3) = -0.1d0 * y_arr(j)
			nuDensMatVecFD(j)%re(3, 3) = 0.9d0 * y_arr(j) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(2, 1) = nuDensMatVecFD(j)%re(1, 2)
			nuDensMatVecFD(j)%re(3, 1) = nuDensMatVecFD(j)%re(1, 3)
			nuDensMatVecFD(j)%re(3, 2) = nuDensMatVecFD(j)%re(2, 3)
			nuDensMatVecFD(j)%im(1, 2) = 0.d0
			nuDensMatVecFD(j)%im(1, 3) = 0.1d0 * y_arr(j)
			nuDensMatVecFD(j)%im(2, 3) = 0.2d0
			nuDensMatVecFD(j)%im(2, 1) = - nuDensMatVecFD(j)%im(1, 2)
			nuDensMatVecFD(j)%im(3, 1) = - nuDensMatVecFD(j)%im(1, 3)
			nuDensMatVecFD(j)%im(3, 2) = - nuDensMatVecFD(j)%im(2, 3)
		end do

		open(unit=fu, file="test_outputs/collint_nunu_int_C_re.dat", status="old")
		open(unit=fv, file="test_outputs/collint_nunu_int_C_im.dat", status="old")
		do j=1, 3
			read (fu, *) ndr(j,:)
			read (fv, *) ndi(j,:)
		end do
		close(fu)
		close(fv)
		er(1,:) = (/1d-3,1d-3,1d-3/)
		er(2,:) = (/1d-3,1d-3,1d-3/)
		er(3,:) = (/1d-3,1d-3,1d-3/)
		ei(1,:) = (/1d-3,1d-3,1d-3/)
		ei(2,:) = (/1d-3,1d-3,1d-3/)
		ei(3,:) = (/1d-3,1d-3,1d-3/)
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				write(tmparg,"('nunu int C ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", coll_nunu_int(iy2, iy3, collArgs, F_nu_sc_re, F_nu_pa_re), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", coll_nunu_int(iy2, iy3, collArgs, F_nu_sc_im, F_nu_pa_im), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do

		!full integral
		iy1=10
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		open(unit=fu, file="test_outputs/collint_nunu_integral_GL_C_re.dat", status="old")
		open(unit=fv, file="test_outputs/collint_nunu_integral_GL_C_im.dat", status="old")
		do i=1, 3
			read (fu, *) ndr(i,:)
			read (fv, *) ndi(i,:)
		end do
		close(fu)
		close(fv)
		er(1,:) = (/1.1d-1, 1.0d-1, 1.0d-1/)
		er(2,:) = (/1.0d-1, 1.1d-1, 1.2d-1/)
		er(3,:) = (/1.0d-1, 1.2d-1, 1.1d-1/)
		ei(1,:) = (/1.0d-1, 1.1d-1, 1.2d-1/)
		ei(2,:) = (/1.1d-1, 1.0d-1, 1.0d-1/)
		ei(3,:) = (/1.2d-1, 1.0d-1, 1.0d-1/)
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				write(tmparg,"('nunu integral C ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", integrate_collint_nunu_GL(coll_nunu_int, collArgs, F_nu_sc_re, F_nu_pa_re), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", integrate_collint_nunu_GL(coll_nunu_int, collArgs, F_nu_sc_im, F_nu_pa_im), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do
#endif

		do j=1, Ny
			nuDensMatVecFD(j)%y=y_arr(j)
			nuDensMatVecFD(j)%re=0.d0
			nuDensMatVecFD(j)%im=0.d0
			nuDensMatVecFD(j)%re(1,1) = (1.0d0 ) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(2,2) = (1.0d0 ) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(3,3) = (1.0d0 ) * fermiDirac(y_arr(j))
		end do
		do j=0, 5
			if (j.eq.0) then
				collArgs%iy=1
			else
				collArgs%iy = j*10
			end if
			if (j.eq.5) then
				ey = 0.15d0
			elseif(j.eq.1) then
				ey=0.02d0
			else
				ey=1.d-2
			end if
			collArgs%y1 = y_arr(collArgs%iy)
			res1=integrate_collint_nunu_NC(coll_nunu_int, collArgs, F_nu_sc_da, F_nu_pa_da)/collArgs%y1**3/8.d0
			res2=dy_damping_fit(collArgs%y1)
			!print*,"s",collargs%y1,res1,res2
			write(tmparg,"('dy-sim GL',I1)") j
			call assert_double_rel_safe(trim(tmparg), res1, res2, 1d-7, ey)
		end do

		!back to linlog momenta -> NC method
		Ny=100
		y_arr = linspace(y_min, y_max, Ny)
		call finish_y_arrays
		do j=1, Ny
			nuDensMatVecFD(j)%y=y_arr(j)
			nuDensMatVecFD(j)%re=0.d0
			nuDensMatVecFD(j)%im=0.d0
			nuDensMatVecFD(j)%re(1,1) = 1.0d0 * y_arr(j) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(2,2) = 1.1d0 * y_arr(j) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(3,3) = 0.9d0 * y_arr(j) * fermiDirac(y_arr(j))
		end do

		!A
		iy1=20
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		open(unit=fu, file="test_outputs/collint_nunu_integral_NC_A.dat", status="old")
		do j=1, 3
			read (fu, *) ndr(j,:)
		end do
		close(fu)
		ndi = 0.d0
		er = 1d-3
		er(3,3) = 1.4d-2
		ei = 1d-3
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				write(tmparg,"('nunu integral NC A ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", integrate_collint_nunu_NC(coll_nunu_int, collArgs, F_nu_sc_re, F_nu_pa_re), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", integrate_collint_nunu_NC(coll_nunu_int, collArgs, F_nu_sc_im, F_nu_pa_im), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do

		!B
		iy1=4
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		open(unit=fu, file="test_outputs/collint_nunu_integral_NC_B.dat", status="old")
		do j=1, 3
			read (fu, *) ndr(j,:)
		end do
		close(fu)
		ndi = 0.d0
		er = 1d-3
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				write(tmparg,"('nunu integral NC B ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", integrate_collint_nunu_NC(coll_nunu_int, collArgs, F_nu_sc_re, F_nu_pa_re), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", integrate_collint_nunu_NC(coll_nunu_int, collArgs, F_nu_sc_im, F_nu_pa_im), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do

#ifdef FULL_F_NU
		!C
		do j=1, Ny
			nuDensMatVecFD(j)%y=y_arr(j)
			nuDensMatVecFD(j)%re=0.d0
			nuDensMatVecFD(j)%im=0.d0
			nuDensMatVecFD(j)%re(1, 1) = 1.d0 * y_arr(j) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(1, 2) = 0.1d0
			nuDensMatVecFD(j)%re(1, 3) = 0.d0
			nuDensMatVecFD(j)%re(2, 2) = 1.1d0 * y_arr(j) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(2, 3) = -0.1d0 * y_arr(j)
			nuDensMatVecFD(j)%re(3, 3) = 0.9d0 * y_arr(j) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(2, 1) = nuDensMatVecFD(j)%re(1, 2)
			nuDensMatVecFD(j)%re(3, 1) = nuDensMatVecFD(j)%re(1, 3)
			nuDensMatVecFD(j)%re(3, 2) = nuDensMatVecFD(j)%re(2, 3)
			nuDensMatVecFD(j)%im(1, 2) = 0.d0
			nuDensMatVecFD(j)%im(1, 3) = 0.1d0 * y_arr(j)
			nuDensMatVecFD(j)%im(2, 3) = 0.2d0
			nuDensMatVecFD(j)%im(2, 1) = - nuDensMatVecFD(j)%im(1, 2)
			nuDensMatVecFD(j)%im(3, 1) = - nuDensMatVecFD(j)%im(1, 3)
			nuDensMatVecFD(j)%im(3, 2) = - nuDensMatVecFD(j)%im(2, 3)
		end do

		iy1=7
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		open(unit=fu, file="test_outputs/collint_nunu_integral_NC_C_re.dat", status="old")
		open(unit=fv, file="test_outputs/collint_nunu_integral_NC_C_im.dat", status="old")
		do i=1, 3
			read (fu, *) ndr(i,:)
			read (fv, *) ndi(i,:)
		end do
		close(fu)
		close(fv)
		er = 1d-3
		ei = 1d-3
		do i=1, flavorNumber
			do j=1, flavorNumber
				collArgs%ix1 = i
				collArgs%ix2 = j
				write(tmparg,"('nunu integral NC C ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", integrate_collint_nunu_NC(coll_nunu_int, collArgs, F_nu_sc_re, F_nu_pa_re), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", integrate_collint_nunu_NC(coll_nunu_int, collArgs, F_nu_sc_im, F_nu_pa_im), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do
#endif

		!now test that get_collision_terms does what expected
		write(*,*)""
		call allocateCmplxMat(cts)
		collArgs%ix1 = 1
		collArgs%ix2 = 1
		collArgs%x = 0.05d0
		collArgs%z = 1.06d0
		collArgs%y2 = 0.d0
		collArgs%y3 = 0.d0
		collArgs%y4 = 0.d0
		collArgs%dme2 = 0.0d0
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)

		res1 = integrate_collint_nunu_NC(fakecollintnunu1, collArgs, F_nu_sc_re, F_nu_pa_re)
		cts = get_collision_terms(collArgs, fakecollintnue0, fakecollintnunu1)
		cts%re(:,:) = cts%re(:,:) * collArgs%y1**2 * collArgs%x**4 / collTermFactor
		cts%im(:,:) = cts%im(:,:) * collArgs%y1**2 * collArgs%x**4 / collTermFactor
#ifdef FULL_F_NU
		ndr = res1/4.d0
		ndi(1,:) = (/    0.d0,  res1/4., res1/4./)
		ndi(2,:) = (/-res1/4.,     0.d0, res1/4./)
		ndi(3,:) = (/-res1/4., -res1/4.,    0.d0/)
#else
		ndr(1,:) = (/res1/4.,     0.d0,    0.d0/)
		ndr(2,:) = (/   0.d0,  res1/4.,    0.d0/)
		ndr(3,:) = (/   0.d0,     0.d0, res1/4./)
		ndi = 0.d0
#endif
		er = 1d-7
		ei = 1d-7
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('collision_terms nunu A ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", cts%re(i,j), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do

		res1 = integrate_collint_nunu_NC(fakecollintnunu1, collArgs, F_nu_sc_re, F_nu_pa_re)
		res2 = integrate_collint_nue_NC(fakecollintnue1, collArgs, F_ab_sc_re, F_ab_ann_re)
		cts = get_collision_terms(collArgs, fakecollintnue1, fakecollintnunu1)
		cts%re(:,:) = cts%re(:,:) * collArgs%y1**2 * collArgs%x**4 / collTermFactor
		cts%im(:,:) = cts%im(:,:) * collArgs%y1**2 * collArgs%x**4 / collTermFactor
#ifdef FULL_F_NU
		ndr = res1/4.d0 + res2
		ndi(1,:) = (/         0.d0,  res1/4.+res2, res1/4.+res2/)
		ndi(2,:) = (/-res1/4.-res2,          0.d0, res1/4.+res2/)
		ndi(3,:) = (/-res1/4.-res2, -res1/4.-res2,         0.d0/)
#else
		ndr(1,:) = (/res1/4. + res2,            res2,           res2/)
		ndr(2,:) = (/          res2,  res1/4. + res2,           res2/)
		ndr(3,:) = (/          res2,            res2, res1/4. + res2/)
		ndi(1,:) = (/ 0.d0, +res2, res2/)
		ndi(2,:) = (/-res2,  0.d0, res2/)
		ndi(3,:) = (/-res2, -res2, 0.d0/)
#endif
		er = 1d-7
		ei = 1d-7
		do i=1, flavorNumber
			do j=1, flavorNumber
				write(tmparg,"('collision_terms nunu B ',2I1)") i,j
				call assert_double_rel_safe(trim(tmparg)//"re", cts%re(i,j), ndr(i,j), 1d-7, er(i,j))
				call assert_double_rel_safe(trim(tmparg)//"im", cts%im(i,j), ndi(i,j), 1d-7, ei(i,j))
			end do
		end do

		do j=1, Ny
			nuDensMatVecFD(j)%y=y_arr(j)
			nuDensMatVecFD(j)%re=0.d0
			nuDensMatVecFD(j)%im=0.d0
			nuDensMatVecFD(j)%re(1,1) = (1.0d0 ) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(2,2) = (1.0d0 ) * fermiDirac(y_arr(j))
			nuDensMatVecFD(j)%re(3,3) = (1.0d0 ) * fermiDirac(y_arr(j))
		end do
		do j=0, 5
			if (j.eq.0) then
				collArgs%iy=1
			else
				collArgs%iy = j*10
			end if
			collArgs%y1 = y_arr(collArgs%iy)
			res1=integrate_collint_nunu_NC(coll_nunu_int, collArgs, F_nu_sc_da, F_nu_pa_da)/collArgs%y1**3/8.d0
			res2=dy_damping_fit(collArgs%y1)
			write(tmparg,"('dy-sim GL',I1)") j
			call assert_double_rel_safe(trim(tmparg), res1, res2, 1d-7, 1d-2)
		end do

		collArgs%z = 1.1d0
		collArgs%y2 = 0.d0
		collArgs%y3 = 0.d0
		collArgs%y4 = 0.d0
		collArgs%dme2 = 0.0d0
		collArgs%iy = iy1
		collArgs%y1 = y_arr(iy1)
		collArgs%x = 0.06d0
		cts = get_collision_terms(collArgs, coll_nue_int, fakecollintnunu0)
		call printMat(cts%re)
		cts = get_collision_terms(collArgs, fakecollintnue0, coll_nunu_int)
		call printMat(cts%re)

		collArgs%x = 0.1d0
		cts = get_collision_terms(collArgs, coll_nue_int, fakecollintnunu0)
		call printMat(cts%re)
		cts = get_collision_terms(collArgs, fakecollintnue0, coll_nunu_int)
		call printMat(cts%re)

		call printTotalTests
		call resetTestCounter
	end subroutine do_test_collint_nunu

	subroutine do_timing_tests
		timing_tests = .true.
#ifndef NO_INTERPOLATION
		call test_dzodx_speed
#endif
		call test_nuDens_speed
		call time_electron_energyDensity

#ifndef NO_INTERPOLATION
		call init_interp_dme2_e
#endif
		call init_interp_FD
		call init_interp_d123
		call test_speed_coll_int
	end subroutine do_timing_tests

end program tests
