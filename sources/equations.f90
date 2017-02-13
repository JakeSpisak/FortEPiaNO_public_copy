module ndEquations
	use precision
	use constants
	use ndErrors
	use ndinteractions
	use ndcosmology
	implicit none
	
	real(dl), parameter :: upper = 1.d8, tol=1.d-3
	contains 
	
	function J_int(o, u)
		real(dl) :: J_int, esuo
		real(dl), intent(in) :: o, u
		
		esuo=exp(sqrt(u*u+o*o))
		J_int = u*u * esuo / ((1.d0+esuo)**2)
	end function J_int
	
	function J_func(o)
		real(dl) :: J_func
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
		
		J_func = rombint_obj(o, J_int,0,upper,tol)/PISQ
	end function J_func
	
	function Jprime_int(o, u)
		real(dl) :: Jprime_int
		real(dl), intent(in) :: o, u
		real(dl) :: uuoo, sqrtuuoo, expsqrtuuoo
		
		uuoo=u*u+o*o
		sqrtuuoo = sqrt(uuoo)
		expsqrtuuoo = exp(sqrtuuoo)
		
		Jprime_int = u*u * expsqrtuuoo * (1.d0-expsqrtuuoo) /(sqrtuuoo*(expsqrtuuoo+1)**3)
	end function Jprime_int
	
	function Jprime(o)
		real(dl) :: Jprime
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
	
		Jprime = rombint_obj(o, Jprime_int,0,upper,tol) * o / PISQ
	end function Jprime
	
	function K_int(o, u)
		real(dl) :: K_int,suo
		real(dl), intent(in) :: o, u
		
		suo=sqrt(u*u+o*o)
		K_int = u*u / (suo * (1.d0+exp(suo)))
	end function K_int
	
	function K_func(o)
		real(dl) :: K_func
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
		
		k_func = rombint_obj(o, k_int,0,upper,tol)/PISQ
	end function K_func
	
	function Kprime_int(o, u)
		real(dl) :: Kprime_int
		real(dl), intent(in) :: o, u
		real(dl) :: uuoo, sqrtuuoo, expsqrtuuoo
		
		uuoo=u*u+o*o
		sqrtuuoo = sqrt(uuoo)
		expsqrtuuoo = exp(sqrtuuoo)
		
		Kprime_int = - u*u / (uuoo*sqrtuuoo*(expsqrtuuoo+1)**2) * &
			(1.d0 + expsqrtuuoo*(sqrtuuoo+1.d0))
	end function Kprime_int
	
	function Kprime(o)
		real(dl) :: Kprime
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
	
		Kprime = rombint_obj(o, Kprime_int,0,upper,tol) * o / PISQ
	end function Kprime
	
	function Y_int(o, u)
		real(dl) :: Y_int, esuo
		real(dl), intent(in) :: o, u
		
		esuo=exp(sqrt(u*u+o*o))
		Y_int = u**4 * esuo / ((1.d0+esuo)**2)
	end function Y_int
	
	function Y_func(o)
		real(dl) :: Y_func
		real(dl), intent(in) :: o
		real(dl) :: rombint_obj
		external rombint_obj
		
		Y_func = rombint_obj(o, Y_int,0,upper,tol)/PISQ
	end function Y_func
	
	function G12_func(o)
		real(dl), dimension(2) :: G12_func
		real(dl), intent(in) :: o
		real(dl) :: ko, jo, kp, jp, tmp
		ko=k_func(o)
		jo=j_func(o)
		kp=kprime(o)
		jp=jprime(o)
		tmp = (kp/6.d0 - ko*kp + jp/6.d0 +jp*ko + jo*kp)
		
		G12_func(1) = PIx2*alpha_fine *(&
			(ko/3.d0 + 2*ko*ko - jo/6.d0 -ko*jo)/o + &
			tmp )
		G12_func(2) = PIx2*alpha_fine*( &
			o * tmp &
			- 4.d0*( (ko+jo)/6.d0 + ko*jo -ko*ko/2.d0) )
	end function G12_func
	
	subroutine init_interp_jkyg12
	
	end subroutine init_interp_jkyg12
	
	function dz_o_dx(x,z)
		real(dl) :: dz_o_dx
		real(dl), intent(in) :: x,z
		real(dl) :: x_o_z
		
		x_o_z = x/z
		
		
		
		
	end function dz_o_dx

end module ndEquations
