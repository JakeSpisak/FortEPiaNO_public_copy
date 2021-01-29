program prepare_gl_nodes
	use precision
	use fpErrors
	use utilities
	implicit none

	integer :: ix, iy
	real(dl), dimension(:), allocatable :: tyv, twv
	integer :: alpha
	integer, parameter :: uid = 8324, uod = 8123
	logical :: exists
	character(len=300) :: tmpstr

	do alpha=2, 3
		do ix=1, maxNgaulag
			write(tmpstr, "('GL_nodes/N',I4.4,'_alpha',I2.2,'.dat')") ix, alpha
			inquire(file=trim(tmpstr), exist=exists)
			if (.not. exists) then
				print *, ix, alpha
				call gaulag(tyv, twv, ix, 1.d0*alpha)
				open(file=trim(tmpstr), unit=uid, status="unknown", form="unformatted")
				do iy=1, ix
					write(uid) tyv(iy), twv(iy)
				end do
				close(uid)
			end if
		end do
	end do
end program
