module fpversion
implicit none
character (len=5) :: version = "1.0.0"
end module

module Precision
implicit none

integer, parameter :: sp = KIND(1.0)
integer, parameter :: dl = KIND(1.d0)
integer, parameter :: qp = 16
end module Precision

module constants
	use precision
	implicit none

	real(dl), parameter :: one_third = 1.d0/3.d0
	real(dl), parameter :: four_thirds = 4.d0/3.d0
	real(dl), parameter :: one_sixth = 1.d0/6.d0
	real(dl), parameter :: SQRT2 = sqrt(2.d0)
	real(dl), parameter :: PI  =3.141592653589793238463d0
	real(dl), parameter :: PIx2=2.d0*PI
	real(dl), parameter :: PIx8D3=8.d0*PI/3.d0
	real(dl), parameter :: PID2=PI/2.d0
	real(dl), parameter :: PID3=PI/3.d0
	real(dl), parameter :: PISQ = PI*PI
	real(dl), parameter :: PISQx2=2.d0*PISQ
	real(dl), parameter :: PISQD2 = PISQ/2.d0
	real(dl), parameter :: PISQD15=PISQ/15.d0
	real(dl), parameter :: PICub = PI*PI*PI
	real(dl), parameter :: Gev2eV = 1.d9, Mev2eV = 1.d6
	real(dl), parameter :: zetaOfThree = 1.20205690315959428540d0
	real(dl), parameter :: zetaOfThreex2DPISQ = zetaOfThree * 2.d0 / PISQ

	real(dl), parameter :: zero = 0.0d0

	!most constants are from the PDG 2020
#ifndef SINSQTHW
#define SINSQTHW 0.23121
#endif
	!numbers used in the real calculation
	real(dl), parameter :: m_e = 0.51099895d0*Mev2eV!eV
	real(dl), parameter :: m_mu = 105.6583745d0*Mev2eV!eV
	real(dl), parameter :: sin2thW_Z = 0.23121d0
	real(dl), parameter :: m_W = 80.379d0*Gev2eV!eV
	real(dl), parameter :: planck_mass = 1.220890d19*Gev2eV
	real(dl), parameter :: m_e_sq = m_e**2
	real(dl), parameter :: m_e_cub = m_e**3
	real(dl), parameter :: m_mu_o_m_e = m_mu/m_e
	real(dl), parameter :: x_muon_cut = 0.5d0!do not compute mu densities above this value, to avoid float overflow
	real(dl), parameter :: G_F = 1.1663787d-5/(Gev2eV*Gev2eV)
	real(dl), parameter :: G_Fsq = G_F * G_F
	real(dl), parameter :: cos2thW_Z = 1.d0-sin2thW_Z
#ifdef GLR_ZERO_MOMENTUM
	!from 10.1016/j.ppnp.2013.03.004
	real(dl), parameter :: sin2thW = 0.23871d0
	real(dl), parameter :: gLe = 0.727d0
	real(dl), parameter :: gLmt = -0.273d0
	real(dl), parameter :: gRemt = 0.233d0
#else
	real(dl), parameter :: sin2thW = SINSQTHW
	real(dl), parameter :: gLe = sin2thW + 0.5d0
	real(dl), parameter :: gLmt = sin2thW - 0.5d0
	real(dl), parameter :: gRemt = sin2thW
#endif
	real(dl), parameter :: alpha_fine = 1.d0/137.035999084d0
	real(dl), parameter :: electron_charge = sqrt(4*PI*alpha_fine)
	real(dl), parameter :: electron_charge_sq = electron_charge ** 2
	real(dl), parameter :: electron_charge_cub = electron_charge ** 3

	integer,  parameter :: maxFlavorNumber = 6
	integer,  parameter :: i_flavorNumber = 3
	!from PDG 2020
	real(dl), parameter :: i_theta12 = 0.307d0
	real(dl), parameter :: i_theta13 = 0.0218d0
	real(dl), parameter :: i_theta23 = 0.545d0
	real(dl), parameter :: i_dm21 = 7.53d-05
	real(dl), parameter :: i_dm31 = 0.002453d0+i_dm21

	real(dl), parameter :: leptDensFactor = -2*SQRT2*G_F*m_e**6/(m_W**2)
	real(dl), parameter :: collTermFactor = G_Fsq/(8.d0*PICub) * m_e_cub
	real(dl), parameter :: overallFactor = planck_mass / sqrt(PIx8D3)
	real(dl), parameter :: dampTermFactor = 7.d0*PISQ*PISQ/135.d0

	real(dl), parameter :: zid = (11.d0/4.d0)**(1.d0/3.d0)

	character(len=5), parameter :: dblfmt = "E17.9"
	character(len=14), parameter :: multidblfmt = "(1P, *("//dblfmt//"))"
end module constants

module variables
	use precision
	use constants
	implicit none

	logical :: timing_tests = .false.
	character(len=300) :: outputFolder
	logical :: firstWrite = .true.
	logical :: firstPoint = .false.
	logical :: checkpoint = .false.
	logical :: force_replace = .false.

	integer :: collint_damping_type
	logical :: collint_offdiag_damping
	logical :: collint_diagonal_zero
	logical :: collint_d_no_nue, collint_d_no_nunu
	logical :: collint_od_no_nue, collint_od_no_nunu
	logical :: damping_read_zero
	logical :: ftqed_temperature_corr
	logical :: ftqed_log_term
	logical :: ftqed_ord3
	logical :: ftqed_e_mth_leptondens
	logical :: save_fd, save_Neff, save_nuDens_evolution, save_w_evolution, save_z_evolution
	logical :: save_energy_entropy_evolution
	logical :: save_number_evolution

	!variables that will be read from config file
	logical :: giveSinSq
	integer :: flavorNumber, flavNumSqu

	!complex matrix, it will host the neutrino density matrix.
	!Intended to have the relative shape correction with respect to the FD in the diagonal
	type cmplxMatNN
		real(dl), dimension(:,:), allocatable :: re, im
		real(dl) :: x, y, z
		logical :: a=.false. !allocated?
	end type cmplxMatNN

	type coll_args
		real(dl) :: y1, y2, y3, y4, x, w, z, dme2
		integer :: ix1, ix2, iy
	end type coll_args

	type derivatives_terms
		logical :: output
		real(dl) :: x, norm
		real(dl), dimension(:), allocatable :: yvec, ydot
		type(cmplxMatNN), dimension(:), allocatable :: Heff, commutator, collterms
	end type derivatives_terms

	type(derivatives_terms) :: intermediateSteps

	real(dl), dimension(:), allocatable :: nuMasses, nuFactor
	real(dl) :: tot_factor_active_nu, tot_factor_nu
	logical , dimension(:), allocatable :: sterile
	real(dl), dimension(:,:), allocatable :: mixMat, mixMatInv, nuMassesMat, leptonDensities
	type(cmplxMatNN) :: nuDensities
	real(dl), dimension(:,:), allocatable :: dampTermMatrixCoeffNue, dampTermMatrixCoeffNunu
	real(dl), dimension(:,:), allocatable :: GL_mat, GR_mat
	real(dl), dimension(:,:,:), allocatable :: GLR_vec
	real(dl), dimension(:,:), allocatable :: idMat
	real(dl), dimension(:), allocatable :: massSplittings
	real(dl), dimension(:,:), allocatable :: mixingAngles

	type(cmplxMatNN), dimension(:), allocatable :: nuDensMatVec, nuDensMatVecFD
	real(dl), dimension(:), allocatable :: nuDensVec
	integer :: ntot
	integer :: ntotrho

	!technical settings
	integer :: verbose = 1
	real(dl) :: Nprintderivs = 100.d0
	logical :: use_gauss_laguerre
	integer :: Nx, Ny
	real(dl) :: x_in, x_fin, y_min, y_max, z_in, logx_in, logx_fin
	real(dl), dimension(:), allocatable :: x_arr, y_arr
	real(dl), dimension(:), allocatable :: feq_arr
	real(dl), dimension(:), allocatable :: y_gl, w_gl, w_gl_arr, w_gl_arr2
	real(dl), dimension(:), allocatable :: dy_arr, fy_arr
	integer :: maxiter
	real(dl) :: toler_jkyg
	real(dl) :: dlsoda_rtol
	real(dl) :: dlsoda_atol_z, dlsoda_atol_d, dlsoda_atol_o

	integer, parameter :: N_opt_xoz = 63
	real(dl), parameter :: opt_xoz_cut = 30.d0
	real(dl), dimension(:), allocatable :: opt_xoz, opt_xoz_w

	integer, parameter :: N_opt_y = 110
	real(dl), parameter :: opt_y_cut = 100.d0
	real(dl), dimension(:), allocatable :: opt_y, opt_y_w

	!2D integrals for log terms of ftqed functions
	integer, parameter :: ln_2dint_Npts = 100 !cannot use too many nor too few
	real(dl), parameter :: ln_2dint_lower = 0.001d0
	real(dl), parameter :: ln_2dint_upper = 100.d0 !must be high enough to catch features of G1/G2
	real(dl), dimension(:), allocatable :: ln_2dint_y, ln_2dint_dy

	!used for 2d interpolation:
	integer, parameter :: interp_nx0 = 750, interp_nz0 = 250, interp_nxz0 = 1800
	integer :: interp_nx, interp_nz, interp_nxz
	integer, parameter :: interp_ny = 100
	logical :: tests_interpolations = .true.
	real(dl), parameter :: interp_logy_min = -2.d0
	real(dl), parameter :: interp_logy_max = 1.5d0
	real(dl), parameter :: interp_zmin0 = 0.9d0, interp_zmax0 = 1.5d0
	real(dl), parameter :: very_early_x=0.1d0*m_e/m_mu!initial temperature is 10 times the muon mass
	real(dl), parameter :: interp_logx_in=log10(very_early_x)
	real(dl) :: interp_zmin, interp_zmax
	real(dl), dimension(:), allocatable :: interp_xvec
	real(dl), dimension(:), allocatable :: interp_yvec
	real(dl), dimension(:), allocatable :: interp_zvec
	real(dl), dimension(:), allocatable :: interp_xozvec
	real(dl) :: low_xoz

	contains

	function get_interpolation_folder()
		character (len=300) :: get_interpolation_folder
		write(get_interpolation_folder, &
			"('interpolations/xi',SP 1P E9.2,'_xf',1P E8.2E1,'_yl',1P E8.2E1,'_yu',1P E8.2E1,'_zl',1P E8.2E1,'_zu',1P E8.2E1,'_xz',1P E8.2E1,'/')" &
		) &
			x_in, x_fin, y_min, y_max, interp_zmin, interp_zmax, low_xoz
		call system("mkdir -p "//trim(get_interpolation_folder))
	end function

	pure subroutine allocateCmplxMat(m)
		type(cmplxMatNN), intent(inout) :: m
		if (.not. m%a) then
			m%a=.true.
			if (.not.allocated(m%re)) &
				allocate(m%re(flavorNumber,flavorNumber))
			if (.not.allocated(m%im)) &
				allocate(m%im(flavorNumber,flavorNumber))
		end if
	end subroutine allocateCmplxMat

	elemental function has_offdiagonal()
		!decide if off-diagonal collision terms are present or not
		logical :: has_offdiagonal
		has_offdiagonal = .not.(collint_offdiag_damping .and. collint_damping_type.eq.0)
	end function has_offdiagonal

	pure subroutine matrix_to_vector(mat, k, vec)
		type(cmplxMatNN), intent(in) :: mat
		integer, intent(inout) :: k
		real(dl), dimension(:), intent(out) :: vec
		integer :: i,j

		do i=1, flavorNumber
			vec(k+i-1) = mat%re(i,i)
		end do
		k=k+flavorNumber
		if (has_offdiagonal()) then
			do i=1, flavorNumber-1
				do j=i+1, flavorNumber
					vec(k) = mat%re(i,j)
					vec(k+1) = mat%im(i,j)
					k=k+2
				end do
			end do
		end if
	end subroutine matrix_to_vector

	subroutine mat_2_vec(mat, N, vec)
		!save the real and imaginary elements of the input array of cmplxMatNN
		!into a 1d vector, with their appropriate order
		type(cmplxMatNN), dimension(:), allocatable, intent(in) :: mat
		integer, intent(in) :: N
		real(dl), dimension(:), intent(out) :: vec
		integer :: k,m

		k=1
		do m=1, N
			call matrix_to_vector(mat(m), k, vec)
		end do
	end subroutine mat_2_vec

	subroutine densMat_2_vec(vec)
		!save the real and imaginary elements of the neutrino density matrix
		!into a 1d vector, with their appropriate order
		real(dl), dimension(:), intent(out) :: vec

		call mat_2_vec(nuDensMatVec, Ny, vec)
	end subroutine densMat_2_vec

	pure subroutine vector_to_matrix(vec, k, mat)
		real(dl), dimension(:), intent(in) :: vec
		integer, intent(inout) :: k
		type(cmplxMatNN), intent(inout) :: mat !must be already allocated!
		integer :: i,j

		do i=1, flavorNumber
			mat%re(i,i) = vec(k+i-1)
			mat%im(i,i) = 0.d0
		end do
		k=k+flavorNumber
		if (has_offdiagonal()) then
			do i=1, flavorNumber-1
				do j=i+1, flavorNumber
					mat%re(i,j) = vec(k)
					mat%im(i,j) = vec(k+1)
					mat%re(j,i) = vec(k)
					mat%im(j,i) = -vec(k+1)
					k=k+2
				end do
			end do
		end if
	end subroutine vector_to_matrix

	subroutine vec_2_densMat(vec)
		!save the real and imaginary elements of the neutrino density matrix
		!from a 1d vector, with their appropriate order, into nuDensMatVec and nuDensMatVecFD
		real(dl), dimension(:), intent(in) :: vec
		integer :: i,j,k,m

		k=1
		do m=1, Ny
			call vector_to_matrix(vec, k, nuDensMatVec(m))
			nuDensMatVecFD(m)%re = nuDensMatVec(m)%re
			nuDensMatVecFD(m)%im = nuDensMatVec(m)%im
			do i=1, flavorNumber
				nuDensMatVecFD(m)%re(i,i) = (1.d0 + nuDensMatVec(m)%re(i,i)) * feq_arr(m)
			end do
		end do
	end subroutine vec_2_densMat
end module variables

module fpInterfaces1
	implicit None
	interface
		pure real(dl) function F_annihilation(n1, n2, f3, f4, a, b, i, j)
			use precision
			use variables
			type(cmplxMatNN), intent(in) :: n1, n2
			real(dl), intent(in) :: f3, f4
			integer, intent(in) :: a, b, i, j
		end function
	end interface
	interface
		pure real(dl) function F_scattering(n1, n3, f2, f4, a, b, i, j)
			use precision
			use variables
			type(cmplxMatNN), intent(in) :: n1, n3
			real(dl), intent(in) :: f2, f4
			integer, intent(in) :: a, b, i, j
		end function
	end interface
	interface
		pure real(dl) function Fnunu(n1, n2, n3, n4, i, j)
			use precision
			use variables
			type(cmplxMatNN), intent(in) :: n1, n2, n3, n4
			integer, intent(in) :: i, j
		end function
	end interface
	interface
		real(dl) function nuDensity_integrator(i1, i2, reim)
			use precision
			integer, intent(in) :: i1, i2
			logical, intent(in), optional :: reim
		end function
	end interface
	interface
		pure real(dl) function ftqed_ln_integrand(x, z, y, k)
			use precision
			real(dl), intent(in) :: x, y, z, k
		end function
	end interface
end module fpInterfaces1

module fpInterfaces2
	interface
		pure real(dl) function collision_integrand_nue(a, b, o, F_ab_ann, F_ab_sc)
			use precision
			use variables
			use fpInterfaces1
			procedure (F_annihilation) :: F_ab_ann
			procedure (F_scattering) :: F_ab_sc
			integer, intent(in) :: a
			real(dl), intent(in) :: b
			type(coll_args), intent(in) :: o
		end function
	end interface
	interface
		pure real(dl) function collision_integrand_nunu(a, b, o, F_nu_sc, F_nu_pa)
			use precision
			use variables
			use fpInterfaces1
			procedure (Fnunu) :: F_nu_sc, F_nu_pa
			integer, intent(in) :: a, b
			type(coll_args), intent(in) :: o
		end function
	end interface
end module fpInterfaces2

module fpInterfaces3
	interface
		pure real(dl) function collision_integrator_nue(f, obj, F_ab_ann, F_ab_sc)
			use precision
			use variables
			use fpInterfaces1
			use fpInterfaces2
			implicit None
			procedure (F_annihilation) :: F_ab_ann
			procedure (F_scattering) :: F_ab_sc
			procedure (collision_integrand_nue) :: f
			type(coll_args), intent(in) :: obj
		end function
	end interface
	interface
		pure real(dl) function collision_integrator_nunu(f, obj, F_nu_sc, F_nu_pa)
			use precision
			use variables
			use fpInterfaces1
			use fpInterfaces2
			implicit None
			procedure (Fnunu) :: F_nu_sc, F_nu_pa
			procedure (collision_integrand_nunu) :: f
			type(coll_args), intent(in) :: obj
		end function
	end interface
end module fpInterfaces3
