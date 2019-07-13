module Precision
implicit none

integer, parameter :: dl = KIND(1.d0)
integer, parameter :: sp = KIND(1.0)
end module Precision

module constants
	use precision
	implicit none

	real(dl), parameter :: four_thirds = 4.d0/3.d0
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

	real(dl), parameter :: zero = 0.0d0

	real(dl), parameter :: m_e = 0.5109989461*Mev2eV!eV
	real(dl), parameter :: m_e_cub = m_e**3
	real(dl), parameter :: m_mu = 105.6583745*Mev2eV!eV
	real(dl), parameter :: m_mu_o_m_e = m_mu/m_e
	real(dl), parameter :: x_muon_cut = 0.5d0!do not compute mu densities above this value, to avoid float overflow
	real(dl), parameter :: m_W = 80.385*Gev2eV!eV
	real(dl), parameter :: G_F = 1.1663787d-5/(Gev2eV*Gev2eV)
	real(dl), parameter :: G_Fsq = G_F * G_F
	real(dl), parameter :: sin2thW = 0.23129
	real(dl), parameter :: cos2thW = 1.d0-sin2thW
	real(dl), parameter :: alpha_fine = 1.d0/137.035999139d0
	real(dl), parameter :: electron_charge = sqrt(4*PI*alpha_fine)
	real(dl), parameter :: planck_mass = 1.220910e19*Gev2eV

	integer,  parameter :: maxFlavorNumber = 6
	integer,  parameter :: i_flavorNumber = 3
	!from PDG 2018: M. Tanabashi et al. (Particle Data Group), Phys.Rev.D, 98, 030001 (2018).
	real(dl), parameter :: i_theta12 = 0.297
	real(dl), parameter :: i_theta13 = 0.0215
	real(dl), parameter :: i_theta23 = 0.425
	real(dl), parameter :: i_dm21 = 7.37e-05
	real(dl), parameter :: i_dm31 = 0.00256

	real(dl), parameter :: leptDensFactor = -8*SQRT2*G_F*m_e**6/(3*m_W**2)
	real(dl), parameter :: collTermFactor = G_Fsq/(8.d0*PICub) * m_e_cub
	real(dl), parameter :: overallFactor = planck_mass / sqrt(PIx8D3)
	real(dl), parameter :: dampTermFactor = -7.d0*PISQ*PISQ/135.d0

	real(dl), parameter :: zid = (11.d0/4.d0)**(1.d0/3.d0)

	character(len=5), parameter :: dblfmt = "E17.9"
	character(len=10), parameter :: multidblfmt = "(*(E17.9))"
end module constants

module variables
	use precision
	implicit none

	logical :: timing_tests = .false.
	character(len=300) :: outputFolder
	logical :: firstWrite = .true.
	logical :: firstPoint = .false.
	logical :: checkpoint = .false.

	integer :: collision_offdiag
	logical :: damping_read_zero
	logical :: ftqed_temperature_corr
	logical :: ftqed_log_term
	logical :: ftqed_ord3
	logical :: save_fd, save_Neff, save_nuDens_evolution, save_w_evolution, save_z_evolution
	logical :: save_energy_entropy_evolution

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
		real(dl) :: y1, y2, y3, y4, x, z, dme2
		integer :: ix1, ix2, iy
	end type coll_args

	real(dl), dimension(:), allocatable :: nuMasses, nuFactor
	real(dl) :: tot_factor_active_nu, tot_factor_nu
	logical , dimension(:), allocatable :: sterile
	real(dl), dimension(:,:), allocatable :: mixMat, mixMatInv, nuMassesMat, leptonDensities
	type(cmplxMatNN) :: nuDensities
	real(dl), dimension(:,:), allocatable :: dampTermMatrixCoeff
	real(dl), dimension(:,:), allocatable :: GL_mat, GR_mat
	real(dl), dimension(:,:,:), allocatable :: GLR_vec
	real(dl), dimension(:,:), allocatable :: idMat
	real(dl), dimension(:), allocatable :: massSplittings
	real(dl), dimension(:,:), allocatable :: mixingAngles

	type(cmplxMatNN), dimension(:), allocatable :: nuDensMatVec, nuDensMatVecFD
	real(dl), dimension(:), allocatable :: nuDensVec
	integer :: ntot

	!technical settings
	integer :: verbose = 1
	real(dl) :: Nprintderivs = 100.d0
	logical :: use_gauss_laguerre
	integer :: Nx, Ny, Nylog
	real(dl) :: x_in, x_fin, y_min, y_max, y_cen, z_in, logx_in, logx_fin
	real(dl), dimension(:), allocatable :: x_arr, y_arr
	real(dl), dimension(:), allocatable :: y_gl, w_gl, w_gl_arr, w_gl_arr2
	real(dl), dimension(:), allocatable :: dy_arr, fy_arr
	integer :: maxiter
	real(dl) :: toler_jkyg, dlsoda_atol, dlsoda_rtol

	integer, parameter :: N_opt_xoz = 63
	real(dl), parameter :: opt_xoz_cut = 30.d0
	real(dl), dimension(:), allocatable :: opt_xoz, opt_xoz_w

	integer, parameter :: N_opt_y = 110
	real(dl), parameter :: opt_y_cut = 100.d0
	real(dl), dimension(:), allocatable :: opt_y, opt_y_w

	!used for 2d interpolation:
	integer, parameter :: interp_nx0 = 750, interp_nz0 = 250, interp_nxz0 = 1800
	integer :: interp_nx, interp_nz, interp_nxz
	integer, parameter :: interp_ny = 100
	real(dl), parameter :: interp_logy_min = -2.
	real(dl), parameter :: interp_logy_max = 1.5
	real(dl), parameter :: interp_zmin0 = 0.9d0, interp_zmax0 = 1.5d0
	real(dl), parameter :: very_early_x=0.1*0.5109989461/105.6583745!initial temperature is 10 times the muon mass
	real(dl), parameter :: interp_logx_in=log10(very_early_x)
	real(dl) :: interp_zmin, interp_zmax
	real(dl), dimension(:), allocatable :: interp_xvec
	real(dl), dimension(:), allocatable :: interp_yvec
	real(dl), dimension(:), allocatable :: interp_zvec
	real(dl), dimension(:), allocatable :: interp_xozvec

	contains

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
end module variables

module ndInterfaces1
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
		real(dl) function nuDensity_integrator(i1, i2, reim)
			use precision
			integer, intent(in) :: i1, i2
			logical, intent(in), optional :: reim
		end function
	end interface
end module ndInterfaces1

module ndInterfaces2
	interface
		pure real(dl) function collision_integrand(a, b, o, F_ab_ann, F_ab_sc)
			use precision
			use variables
			use ndInterfaces1
			procedure (F_annihilation) :: F_ab_ann
			procedure (F_scattering) :: F_ab_sc
			integer, intent(in) :: a
			real(dl), intent(in) :: b
			type(coll_args), intent(in) :: o
		end function
	end interface
end module ndInterfaces2

module ndInterfaces3
	interface
		pure real(dl) function collision_integrator(f, obj, F_ab_ann, F_ab_sc)
			use precision
			use variables
			use ndInterfaces1
			use ndInterfaces2
			implicit None
			procedure (F_annihilation) :: F_ab_ann
			procedure (F_scattering) :: F_ab_sc
			procedure (collision_integrand) :: f
			type(coll_args), intent(in) :: obj
		end function
	end interface
end module ndInterfaces3
