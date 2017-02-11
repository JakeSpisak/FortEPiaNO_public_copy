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
	real(dl), parameter :: PID2=PI/2.d0
	real(dl), parameter :: PID4=PI/4.d0
	real(dl), parameter :: SQRTPI = sqrt(PI)
	real(dl), parameter :: PISQ = PI*PI
	real(dl), parameter :: PICub = PI*PI*PI
	real(dl), parameter :: e_neper = 2.718281828459045235d0
	real(dl), parameter :: gamma_par = 0.577215664901532861d0
	
	real(dl), parameter :: Gev2eV = 1.d9, Mev2eV = 1.d6, kev2eV = 1.d3
	
	real(dl), parameter :: zero = 0.0d0
	real(dl), parameter :: largeNum = 1.e20
	
	real(dl), parameter :: c = 299792458!m/s
	real(dl), parameter :: hbar = 1.054571800d-34 !Js
	real(dl), parameter :: m_e = 0.5109989461*Mev2eV!eV
	real(dl), parameter :: G_F = 1.1663787d-5/(Gev2eV*Gev2eV)
	real(dl), parameter :: G_Fsq = G_F * G_F
	real(dl), parameter :: sin2thW =  0.23129
	
	integer,  parameter :: i_flavorNumber = 3
	!from PDG 2016: C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016).
	real(dl), parameter :: i_theta12 = 0.5839958715755919
	real(dl), parameter :: i_theta13 = 0.1485320301705043
	real(dl), parameter :: i_theta23 = 0.7953988301841435
	real(dl), parameter :: i_dm12 = 7.53e-05
	real(dl), parameter :: i_dm23 = 0.00244
	real(dl), parameter :: i_deltaCP13 = 0.
	
	real(dl), parameter :: i_photonTempToday = 2.7255
	real(dl), parameter :: i_HubbleParam = 70.
	
	!variables that will be read from config file
	logical :: massOrdering
	integer :: flavorNumber
	real(dl) :: m_lightest
	real(dl) :: theta12, dm12
	real(dl) :: theta13, theta23, dm23, deltaCP13
	real(dl) :: theta14, theta24, theta34, dm14
	
	real(dl) :: photonTemperatureToday, hubbleParam

end module constants
