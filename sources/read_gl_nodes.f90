program read_gl_nodes
	implicit none

	integer, parameter :: dl = KIND(1.d0)
	integer :: ix, iy
	real(dl), dimension(:), allocatable :: tyv, twv
	integer :: alpha
	integer, parameter :: uid = 8324, uod = 8123
	character(len=300) :: tmpstr
	logical :: exists

	do alpha=2, 3
		write(tmpstr, "('GL_nodes/expected_gl_nodes_alpha',I2.2,'.dat')") alpha
		open(file=trim(tmpstr), unit=uod, status="unknown", form="unformatted")
		do ix=1, 1500
			write(tmpstr, "('GL_nodes/N',I4.4,'_alpha',I2.2,'.dat')") ix, alpha
			inquire(file=trim(tmpstr), exist=exists)
			if (exists) then
				if (allocated(tyv)) &
					deallocate(tyv)
				if (allocated(twv)) &
					deallocate(twv)
			allocate(tyv(ix), twv(ix))
			open(file=trim(tmpstr), unit=uid, form="unformatted")
			do iy=1, ix
				read(uid) tyv(iy), twv(iy)
			end do
			close(uid)
			write(uod) ix, tyv(1), twv(1), tyv(ix), twv(ix)
			print *,ix, tyv(1), twv(1), tyv(ix), twv(ix)
		end if
	end do
	close(uod)
	end do
end program
