subroutine derivatives(n, x, rho, rhodot)
	use precision
	use variables
	use ndEquations
	implicit none
	real(dl), dimension(:), allocatable :: rhodot
	type(cmplxMatNN) :: mat
	integer :: n, i, j, k, m
	real(dl) :: x, z
	real(dl), dimension(:) :: rho
	
	allocate(rhodot(n))
	print *,"a",rho,"b",size(rho)
	z = rho(n)
	call vec_2_densMat(rho(1:n-1))
	k=1
	do m=1, Ny
		mat = drhoy_dx_fullMat(x,z,m)
		do i=1, flavorNumber
			do j=i, flavorNumber
				rhodot(k) = mat%re(i,j)
				k=k+1
			end do
			if (i.lt.flavorNumber) then
				do j=i+1, flavorNumber
					rhodot(k) = mat%im(i,j)
					k=k+1
				end do
			end if
		end do
	end do
	rhodot(ntot) = dz_o_dx(x,z)
end subroutine derivatives

subroutine jdum

end subroutine jdum
