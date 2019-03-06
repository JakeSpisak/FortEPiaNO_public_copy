module Precision
implicit none

integer, parameter :: dl = KIND(1.d0)
integer, parameter :: sp = KIND(1.0)
end module Precision

module constants
	use precision
	implicit none
	
	real(dl), parameter :: SQRT2 = sqrt(2.d0)
	real(dl), parameter :: SQRT3 = sqrt(3.d0)
	real(dl), parameter :: PI  =3.141592653589793238463d0
	real(dl), parameter :: PIx2=2.d0*PI
	real(dl), parameter :: PIx4=4.d0*PI
	real(dl), parameter :: PIx8=8.d0*PI
	real(dl), parameter :: PIx8D3=8.d0*PI/3.d0
	real(dl), parameter :: PID2=PI/2.d0
	real(dl), parameter :: PID3=PI/3.d0
	real(dl), parameter :: PID4=PI/4.d0
	real(dl), parameter :: SQRTPI = sqrt(PI)
	real(dl), parameter :: PISQ = PI*PI
	real(dl), parameter :: PISQD2 = PISQ/2.d0
	real(dl), parameter :: PISQD4 = PISQ/4.d0
	real(dl), parameter :: PISQD15=PISQ/15.d0
	real(dl), parameter :: PISQD30=PISQ/30.d0
	real(dl), parameter :: PICub = PI*PI*PI
	real(dl), parameter :: e_neper = 2.718281828459045235d0
	real(dl), parameter :: gamma_par = 0.577215664901532861d0
	real(dl), parameter :: zeta3  = 1.2020569031595942853997_dl
	real(dl), parameter :: zeta5  = 1.0369277551433699263313_dl
	real(dl), parameter :: zeta7  = 1.0083492773819228268397_dl
	real(dl), parameter :: Gev2eV = 1.d9, Mev2eV = 1.d6, kev2eV = 1.d3
	
	real(dl), parameter :: zero = 0.0d0
	real(dl), parameter :: largeNum = 1.e20
	
	real(dl), parameter :: c = 299792458!m/s
	real(dl), parameter :: hbar = 1.054571800d-34 !Js
	real(dl), parameter :: m_e = 0.5109989461*Mev2eV!eV
	real(dl), parameter :: m_e_sq = m_e**2
	real(dl), parameter :: m_e_cub = m_e**3
	real(dl), parameter :: m_W = 80.385*Gev2eV!eV
	real(dl), parameter :: G_F = 1.1663787d-5/(Gev2eV*Gev2eV)
	real(dl), parameter :: G_Fsq = G_F * G_F
	real(dl), parameter :: sin2thW =  0.23129
	real(dl), parameter :: alpha_fine = 1.d0/137.035999139d0
	real(dl), parameter :: planck_mass = 1.220910e19*Gev2eV
	
	integer,  parameter :: maxflavorNumber = 4
	integer,  parameter :: i_flavorNumber = 3
	!from PDG 2016: C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016).
	real(dl), parameter :: i_theta12 = 0.5840
	real(dl), parameter :: i_theta13 = 0.1485
	real(dl), parameter :: i_theta23 = 0.7954
	real(dl), parameter :: i_dm12 = 7.53e-05
	real(dl), parameter :: i_dm23 = 0.00244
	real(dl), parameter :: i_deltaCP13 = 0.
	
	real(dl), parameter :: ymed = 3.15137
	real(dl), parameter :: leptDensFactor = -8*SQRT2*G_F*m_e**6/(3*m_W**2)
	real(dl), parameter :: collTermFactor = G_Fsq/(8.d0*PICub) * m_e_cub
	real(dl), parameter :: overallFactor = planck_mass / sqrt(PIx8D3)
	real(dl), parameter :: dampTermFactor = -7.d0*PISQ*PISQ/135.d0

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
	
	!variables that will be read from config file
	logical :: massOrdering, giveSinSq
	integer :: flavorNumber, flavNumSqu
	real(dl) :: m_lightest
	real(dl) :: theta12, dm12
	real(dl) :: theta13, theta23, dm23, deltaCP13
	real(dl) :: theta14, theta24, theta34, dm14

	!complex matrix, it will host the neutrino density matrix.
	!Intended to have the relative shape correction with respect to the FD in the diagonal
	type cmplxMatNN
		real(dl), dimension(:,:), allocatable :: re, im
		real(dl) :: x, y, z
		logical :: a=.false. !allocated?
	end type cmplxMatNN
	
	type nuDensArgs
		real(dl) :: x,z
		integer iFl
	end type nuDensArgs
	
	type coll_args
		real(dl) :: y1, y2, y3, y4, x, z, dme2
		integer :: ix1, ix2, iy
	end type coll_args
	
	real(dl), dimension(:), allocatable :: nuMasses, nuFactor
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
	
	pure subroutine deallocateCmplxMat(m)
		type(cmplxMatNN), intent(inout) :: m
		
		if (m%a) then
			m%a=.false.
			deallocate(m%re, m%im)
		end if
	end subroutine deallocateCmplxMat
end module variables
