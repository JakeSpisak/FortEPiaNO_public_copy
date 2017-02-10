module ndMatrices
	use precision
	implicit none
	
	contains
	
	subroutine createIdentityMat(mat, nc)
		integer, intent(in) :: nc
		real(dl), dimension(:,:), allocatable, intent(out) :: mat
		integer :: ix
		
		if (.not. allocated(mat)) then
			allocate(mat(nc, nc))
		end if
		
		mat = 0.d0
		do ix = 1, nc
			mat(ix, ix) = 1.d0
		end do
	end subroutine createIdentityMat
	
	subroutine createDiagMat(mat, nc, d)
		integer, intent(in) :: nc
		real(dl), dimension(:), intent(in) :: d
		real(dl), dimension(:,:), allocatable, intent(out) :: mat
		integer :: ix
		
		if (.not. allocated(mat)) then
			allocate(mat(nc, nc))
		end if
		
		mat=0.d0
		do ix=1, nc
			mat(ix,ix) = d(ix)
		end do
	end subroutine createDiagMat
	
	subroutine createRotMat(mat, nc, i1, i2, theta)
		integer, intent(in) :: nc, i1, i2
		real(dl), intent(in) :: theta
		real(dl), dimension(:,:), allocatable, intent(out) :: mat
		integer :: ix1, ix2
		
		if (.not. allocated(mat)) then
			allocate(mat(nc, nc))
		end if
		if (ix1.lt.ix2) then
			ix1=i1
			ix2=i2
		else
			ix2=i1
			ix1=i2
		end if
		
		call createIdentityMat(mat, nc)
		
		mat(ix1, ix1) = cos(theta)
		mat(ix2, ix2) = cos(theta)
		mat(ix1, ix2) = -sin(theta)
		mat(ix2, ix1) = sin(theta)
	end subroutine createRotMat
	
	function trackMat(inM)
		real(dl), dimension(:,:), allocatable, intent(in) :: inM
		real(dl) :: trackMat
		integer :: i
		integer, dimension(2) :: s
		
		s = shape(inM)
		
		if (s(1) .ne. s(2)) then
			write (*,*) "Error: non-squared matrices have no defined track"
			call exit()
		end if
		
		trackMat = 0.d0
		do i=1, s(1)
			trackMat = trackMat + inM(i,i)
		end do
	end function trackMat

	subroutine transposeMat(inM, outM)
		real(dl), dimension(:,:), allocatable, intent(in) :: inM
		real(dl), dimension(:,:), allocatable, intent(out) :: outM
		integer :: i,j
		integer, dimension(2) :: sizeM
		
		sizeM=shape(inM)

		if (.not. allocated(outM)) allocate(outM(sizeM(2), sizeM(1)))
		
		outM = 0.d0
		
		do i=1, sizeM(1)
			do j=1, sizeM(2)
				outM(j,i) = inM(i,j)
			end do
		end do
	end subroutine transposeMat

	subroutine minoreMat(inM, ix1, ix2, outM)
		integer, intent(in) :: ix1, ix2
		real(dl), dimension(:,:), allocatable, intent(in) :: inM
		real(dl), dimension(:,:), allocatable, intent(out) :: outM
		integer :: i,j,r,c
		integer, dimension(2) :: sizeM
		
		sizeM=shape(inM)

		if (.not. allocated(outM)) allocate(outM(sizeM(1)-1, sizeM(2)-1))
		
		outM = 0.d0
		r=0
		do i=1, sizeM(1)
			c=0
			if (i.ne.ix1) then
				r=r+1
				do j=1, sizeM(2)
					if (j .ne. ix2) then
						c=c+1
						outM(r,c) = inM(i,j)
					end if
				end do
			end if
		end do
	end subroutine minoreMat
	
	recursive function determMat(mat) result(determ)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat
		real(dl), dimension(:,:), allocatable :: mi
		real(dl) :: determ
		integer :: ix
		integer, dimension(2) :: s
		
		s=shape(mat)
		determ = 0.d0
		
		if (s(1) .ne. s(2)) then
			write (*,*) "Error: non-squared matrices have no det"
			call exit()
		end if
		if (s(1) .eq. 1) then
			determ = mat(1,1)
		else if (s(1) .eq. 2) then
			determ = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
		else
			do ix=1, s(2)
				call minoreMat(mat, 1, ix, mi)
				determ = determ + determMat(mi) * mat(1,ix) * (1 - 2 * modulo(1+ix, 2) )
			end do
		end if
	end function determMat

	function cofactorMat(mat, i, j)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat
		real(dl), dimension(:,:), allocatable :: mi
		integer, intent(in) :: i,j
		real(dl) :: cofactorMat
		
		call minoreMat(mat, i, j, mi)
		cofactorMat = determMat(mi) * (1 - 2 * modulo(i+j, 2) )! for (-1)**(i+j)
		
	end function cofactorMat
	
	subroutine inverseMat(inM, outM)
		real(dl), dimension(:,:), allocatable, intent(in) :: inM
		real(dl), dimension(:,:), allocatable, intent(out) :: outM
		integer :: i,j,n
		integer, dimension(2) :: s
		real(dl) :: det
		
		s=shape(inM)
		det = determMat(inM)
		if (det.eq.0) then
			write (*,*) "Error: the matrix is degenerate!"
			call exit()
		end if
		if (s(1) .ne. s(2)) then
			write (*,*) "Error: non-squared matrices have no inverse"
			call exit()
		else
			n=s(1)
		end if

		if (.not. allocated(outM)) allocate(outM(n,n))
		outM=0
		do i=1, n
			do j=1, n
				outM(i,j) = cofactorMat(inM, j, i)/det
			end do
		end do
	
	end subroutine inverseMat

	subroutine multMat(mat1, mat2, outMat)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat1, mat2
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat
		integer :: i,j,k
		integer, dimension(2) :: size1, size2
		
		size1=shape(mat1)
		size2=shape(mat2)
		
		if (size1(2) .ne. size2(1)) then
			write (*,*) "Error! cannot multiply matrices with these sizes:"
			write (*,*) size1, size2
			call exit()
		end if
		
		if (.not. allocated(outMat)) allocate(outMat(size1(1), size2(2)))
		outMat=0.0d0
		
		!multiply
		do i=1, size1(1)
			do j=1, size2(2)
				do k=1, size1(2)
					outMat(i,j) = outMat(i,j) + mat1(i,k)*mat2(k,j)
				end do
			end do
		end do
	end subroutine multMat
	
	subroutine printMat(m)
		real(dl), dimension(:,:), allocatable, intent(in) :: m
		integer, dimension(2) :: s
		integer :: i
		
		s = shape(m)
		do i=1, s(1)
			write(*,"(*(E14.6))") m(i,:)
		end do
	end subroutine printMat

end module ndMatrices
