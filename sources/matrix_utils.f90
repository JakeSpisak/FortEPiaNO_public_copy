module ndMatrices
	use precision
	use ndErrors
	implicit none
	
	contains
	
	function numRows(mat)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat
		integer :: numRows
		integer, dimension(2) :: s
		
		s = shape(mat)
		numRows = s(1)
	end function numRows
	
	function numColumns(mat)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat
		integer :: numColumns
		integer, dimension(2) :: s
		
		s = shape(mat)
		numColumns = s(2)
	end function numColumns
	
	subroutine createEmptyMat(mat, nc)
		integer, intent(in) :: nc
		real(dl), dimension(:,:), allocatable, intent(out) :: mat
		
		if (.not. allocated(mat)) then
			allocate(mat(nc, nc))
		end if
		
		mat = 0.d0
	end subroutine createEmptyMat
	
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
			call criticalError("non-squared matrices have no defined track!")
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
			call criticalError("non-squared matrices have no det!")
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
			call criticalError("the matrix is degenerate!")
		end if
		if (s(1) .ne. s(2)) then
			call criticalError("non-squared matrices have no inverse!")
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

	subroutine multiplyMat(mat1, mat2, outMat)
		real(dl), dimension(:,:), intent(in) :: mat1, mat2
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat
		integer :: i,j,k
		integer, dimension(2) :: size1, size2
		character(len=12) :: sizestr
		
		size1=shape(mat1)
		size2=shape(mat2)
		write (sizestr,"(2('(',I1,',',I1,') '))") size1(:), size2(:)
		if (size1(2) .ne. size2(1)) then
			call criticalError("cannot multiply matrices with these sizes: "//sizestr)
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
	end subroutine multiplyMat

	subroutine tripleProdMat(mat1, mat2, mat3, outMat)
		real(dl), dimension(:,:), intent(in) :: mat1, mat2, mat3
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat
		real(dl), dimension(:,:), allocatable :: m1
		
		call multiplyMat(mat1, mat2, m1)
		call multiplyMat(m1, mat3, outMat)
	end subroutine tripleProdMat

	subroutine quadrupleProdMat(mat1, mat2, mat3, mat4, outMat)
		real(dl), dimension(:,:), intent(in) :: mat1, mat2, mat3, mat4
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat
		real(dl), dimension(:,:), allocatable :: m1, m2
		
		call multiplyMat(mat1, mat2, m1)
		call multiplyMat(m1, mat3, m2)
		call multiplyMat(m2, mat4, outMat)
	end subroutine quadrupleProdMat
	
	subroutine Commutator(mat1, mat2, outmat)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat1, mat2
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat
		
		call CommAntiComm(mat1, mat2, outmat, -1)
	end subroutine Commutator
	
	subroutine AntiCommutator(mat1, mat2, outmat)
		real(dl), dimension(:,:), allocatable, intent(in) :: mat1, mat2
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat
		
		call CommAntiComm(mat1, mat2, outmat, 1)
	end subroutine antiCommutator
	
	subroutine CommAntiComm(mat1, mat2, outmat, s)
		!use the sign s to switch between commutator and anticommutator
		real(dl), dimension(:,:), allocatable, intent(in) :: mat1, mat2
		real(dl), dimension(:,:), allocatable, intent(out) :: outMat
		real(dl), dimension(:,:), allocatable :: m1, m2
		integer, intent(in) :: s
		integer, dimension(2) :: size1, size2
		integer :: d
		character(len=2) :: sstr
		
		write(sstr, "(I2)") s
		
		size1=shape(mat1)
		size2=shape(mat2)
		if (abs(s).ne.1) &
			call criticalError("wrong sign in CommAntiComm (it must be either -1 or 1): "//sstr)
		
		if ((size1(2) .ne. size1(2)) .or. (size1(2) .ne. size2(1)) .or. (size2(1) .ne. size2(2))) then
			call criticalError("cannot compute commutator of non-square matrices!")
		else
			d = size1(1)
		end if
		
		allocate(m1(d,d), m2(d,d), outMat(d,d))
		call multiplyMat(mat1, mat2, m1)
		call multiplyMat(mat2, mat1, m2)
		outMat = m1 + s * m2
	end subroutine
	
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
