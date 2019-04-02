module Precision
implicit none

integer, parameter :: dl = KIND(1.d0)
integer, parameter :: sp = KIND(1.0)
end module Precision

module constants
	use precision
	implicit none

	real(dl), parameter :: SQRT2 = sqrt(2.d0)
	real(dl), parameter :: PI  =3.141592653589793238463d0
	real(dl), parameter :: PIx2=2.d0*PI
	real(dl), parameter :: PIx8D3=8.d0*PI/3.d0
	real(dl), parameter :: PID2=PI/2.d0
	real(dl), parameter :: PID3=PI/3.d0
	real(dl), parameter :: PISQ = PI*PI
	real(dl), parameter :: PISQD2 = PISQ/2.d0
	real(dl), parameter :: PISQD15=PISQ/15.d0
	real(dl), parameter :: PICub = PI*PI*PI
	real(dl), parameter :: Gev2eV = 1.d9, Mev2eV = 1.d6

	real(dl), parameter :: zero = 0.0d0

	real(dl), parameter :: m_e = 0.5109989461*Mev2eV!eV
	real(dl), parameter :: m_e_cub = m_e**3
	real(dl), parameter :: m_mu = 105.6583745*Mev2eV!eV
	real(dl), parameter :: m_mu_o_m_e = m_mu/m_e
	real(dl), parameter :: m_W = 80.385*Gev2eV!eV
	real(dl), parameter :: G_F = 1.1663787d-5/(Gev2eV*Gev2eV)
	real(dl), parameter :: G_Fsq = G_F * G_F
	real(dl), parameter :: sin2thW =  0.23129
	real(dl), parameter :: alpha_fine = 1.d0/137.035999139d0
	real(dl), parameter :: planck_mass = 1.220910e19*Gev2eV

	integer,  parameter :: maxflavorNumber = 4
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
	real(dl), parameter :: coeff_dw_dx = -30.d0/(7.d0*PISQ*PISQ)

	real(dl), parameter :: zid = (11.d0/4.d0)**(1.d0/3.d0)

	real(dl), parameter :: fe_l = 0.d0, fe_u = 100.d0

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
	logical :: dme2_temperature_corr
	logical :: save_fd, save_nuDens_evolution, save_w_evolution, save_z_evolution

	!variables that will be read from config file
	logical :: giveSinSq
	integer :: flavorNumber, flavNumSqu
	real(dl) :: theta12, dm21
	real(dl) :: theta13, theta23, dm31
	real(dl) :: theta14, theta24, theta34, dm41

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
	real(dl), dimension(:,:), allocatable :: dampTermMatrixCoeff
	real(dl), dimension(:,:), allocatable :: GL_mat, GR_mat
	real(dl), dimension(:,:,:), allocatable :: GLR_vec
	real(dl), dimension(:,:), allocatable :: idMat
	real(dl), dimension(:,:), allocatable :: xcutsCollInt

	type(cmplxMatNN), dimension(:), allocatable :: nuDensMatVec, nuDensMatVecFD
	real(dl), dimension(:), allocatable :: nuDensVec
	integer :: ntot

	!technical settings
	integer :: verbose = 1
	real(dl) :: Nprintderivs = 100.d0
	integer :: Nx, Ny, Nylog
	real(dl) :: x_in, x_fin, y_min, y_max, y_cen, z_in, logx_in, logx_fin
	real(dl), dimension(:), allocatable :: x_arr, y_arr
	real(dl), dimension(:), allocatable :: dy_arr, fy_arr
	integer :: maxiter
	real(dl) :: toler, toler_dme2, toler_ed, toler_jkyg, dlsoda_atol, dlsoda_rtol

	!used for 2d interpolation:
	integer, parameter :: interp_nx0 = 750, interp_nz0 = 250, interp_nxz0 = 1500
	integer :: interp_nx, interp_nz, interp_nxz
	integer, parameter :: interp_ny = 40!not used in real code, it's for interpolations of D and Pi functions
	real(dl), parameter :: interp_zmin0 = 0.9d0, interp_zmax0 = 1.5d0
	real(dl) :: interp_zmin, interp_zmax
	real(dl), dimension(:), allocatable :: interp_xvec
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
