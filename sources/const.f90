module Precision
implicit none

integer, parameter :: dl = KIND(1.d0)
integer, parameter :: sp = KIND(1.0)
end module Precision

module constants
	use precision
	implicit none
	
	real(dl), parameter :: largeNum = 1.e20
	
	!from PDG 2016: C. Patrignani et al.(Particle Data Group), Chin. Phys. C, 40, 100001 (2016).
	real(dl), parameter :: nh_theta12 = 0.5839958715755919
	real(dl), parameter :: nh_theta13 = 0.1485320301705043
	real(dl), parameter :: nh_theta23 = 0.7953988301841435
	real(dl), parameter :: nh_dm12 = 7.53e-05
	real(dl), parameter :: nh_dm23 = 0.00244
	real(dl), parameter :: deltaCP = 0.


end module constants
