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
	real(dl), parameter :: PISQD15=PISQ/15.d0
	real(dl), parameter :: PISQD30=PISQ/30.d0
	real(dl), parameter :: PICub = PI*PI*PI
	real(dl), parameter :: e_neper = 2.718281828459045235d0
	real(dl), parameter :: gamma_par = 0.577215664901532861d0
	
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
	
	integer,  parameter :: i_flavorNumber = 3
	!from PDG 2016: C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016).
	real(dl), parameter :: i_theta12 = 0.5840
	real(dl), parameter :: i_theta13 = 0.1485
	real(dl), parameter :: i_theta23 = 0.7954
	real(dl), parameter :: i_dm12 = 7.53e-05
	real(dl), parameter :: i_dm23 = 0.00244
	real(dl), parameter :: i_deltaCP13 = 0.
	
	real(dl), parameter :: i_photonTempToday = 2.7255
	real(dl), parameter :: i_HubbleParam = 70.
	
	real(dl), parameter :: leptDensFactor = -8*SQRT2*G_F*m_e**6/(3*m_W**2)
	
	character(len=5), parameter :: dblfmt = "E15.7"
end module constants

module variables
	use precision
	implicit none
	
	character(len=300) :: outputFolder
	logical :: firstWrite = .true.
	logical :: firstPoint = .false.
	logical :: checkpoint = .false.
	
	logical :: coll_scatt_em, coll_scatt_ep, coll_annih_epem
	logical :: dme2_temperature_corr
	
	!variables that will be read from config file
	logical :: massOrdering, only_1a_1s
	integer :: flavorNumber, flavNumSqu
	real(dl) :: m_lightest
	real(dl) :: theta12, dm12
	real(dl) :: theta13, theta23, dm23, deltaCP13
	real(dl) :: theta14, theta24, theta34, dm14
	
	real(dl) :: photonTemperatureToday, hubbleParam

	!matrix
	type cmplxMatNN
		real(dl), dimension(:,:), allocatable :: re, im
		real(dl) :: x, y, z
		logical :: a=.false. !allocated?
	end type cmplxMatNN
	
	type coll_args
		type(cmplxMatNN) :: na, nb
		real(dl) :: y1, y2, y3, y4, x, z, dme2
		logical :: s1, s2, s3, s4
		integer :: ix1, ix2, a, b
		integer :: rI !return real or imaginary part
	end type coll_args
	
	real(dl), dimension(:), allocatable :: nuMasses, nuFactor
	logical , dimension(:), allocatable :: sterile
	real(dl), dimension(:,:), allocatable :: mixMat, mixMatInv, nuMassesMat, leptonDensities
	real(dl), dimension(:,:), allocatable :: GL_mat, GR_mat
	real(dl), dimension(:,:,:), allocatable :: GLR_vec
	real(dl), dimension(:,:), allocatable :: idMat
	
	type(cmplxMatNN), dimension(:), allocatable :: nuDensMatVec
	real(dl), dimension(:), allocatable :: nuDensVec
	integer :: ntot
	
	!technical settings
	integer :: verbose = 1
	integer :: Nx, Ny
	real(dl) :: x_in, x_fin, y_min, y_max, z_in, logx_in, logx_fin, logy_min, logy_max
	real(dl), dimension(:), allocatable :: x_arr, y_arr
	integer :: maxiter
	real(dl) :: toler, dlsoda_atol, dlsoda_rtol

	!used for 2d interpolation:
	integer, parameter :: interp_nx = 100, interp_nz = 50
	real(dl), parameter :: interp_zmin = 1.0d0, interp_zmax = 1.5d0
	real(dl), dimension(:), allocatable :: interp_xvec
	real(dl), dimension(:), allocatable :: interp_zvec
	real(dl), dimension(:), allocatable :: interp_xozvec
	
	contains
	
	subroutine allocateCmplxMat(m)
		type(cmplxMatNN) :: m
		
		if (.not. m%a) then
			m%a=.true.
			if (.not.allocated(m%re)) &
				allocate(m%re(flavorNumber,flavorNumber))
			if (.not.allocated(m%im)) &
				allocate(m%im(flavorNumber,flavorNumber))
		end if
	end subroutine allocateCmplxMat
	
	subroutine deallocateCmplxMat(m)
		type(cmplxMatNN) :: m
		
		if (.not. m%a) then
			m%a=.false.
			deallocate(m%re, m%im)
		end if
	end subroutine deallocateCmplxMat
end module variables
