module ndInteractions
	use precision
	use constants
	use ndErrors
	implicit none
	
	contains
	
	function E_k_m(k,m)
		real(dl), intent(in) :: k,m
		real(dl) :: E_k_m
		E_k_m = sqrt(k*k+m*m)
	end function E_k_m
	
	function fermiDirac(k,m,T)
		real(dl) :: fermiDirac, tmp
		real(dl), intent(in) :: k,m,T
		
		tmp = exp(E_k_m(k,m)/T) + 1.d0
		fermiDirac = 1.d0/tmp
	end function fermiDirac
	
	function dme2_e_i1(vec, k)
		real(dl) :: dme2_e_i1
		real(dl), intent(in) :: k
		real(dl), dimension(3) :: vec
		real(dl) :: m,T,p
		m=vec(1)
		T=vec(2)
		dme2_e_i1 = k*k/E_k_m(k,m)*fermiDirac(k,m,T)
	end function dme2_e_i1
	
	function dme2_e_i2(vec, k)
		real(dl) :: dme2_e_i2
		real(dl), intent(in) :: k
		real(dl), dimension(3) :: vec
		real(dl) :: m,T,p
		m=vec(1)
		T=vec(2)
		p=vec(3)
		dme2_e_i2 = k/E_k_m(k,m)*fermiDirac(k,m,T)*log(abs((p+k)/(p-k)))
	end function dme2_e_i2
	
	function dme2_electron(p,T)
		real(dl) :: dme2_electron, tmp, rombint_obj
		real(dl), intent(in) :: p,T
		real(dl), dimension(3) :: vec
		external rombint_obj
		
		vec(1) = m_e
		vec(2) = T
		vec(3) = p
		
		dme2_electron = PIx2 * alpha_fine * T*T/3.d0
		tmp = rombint_obj(vec, dme2_e_i1, 0.d0, 1d8, 1d-3, 200, 5)
		dme2_electron = dme2_electron + tmp*alpha_fine/pid4
		tmp = rombint_obj(vec, dme2_e_i2, 0.d0, 1d8, 1d-3, 200, 5)
		dme2_electron = dme2_electron + tmp*alpha_fine*m_e*m_e/(pid2*p)
	
	end function dme2_electron

end module
