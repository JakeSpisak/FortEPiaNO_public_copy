subroutine assert_double(testname, num1, num2, tol)
	use precision
	implicit none
	character(len=*) :: testname
	real(dl) :: num1, num2, tol

	if (abs(num1 - num2) .lt. tol) then
		write(*, fmt="(a)", advance="no") "."
	else
		print *, testname, " failed"
		call exit()
	end if
end subroutine assert_double

program tests
	use precision
	use ndConfig
	use ndErrors
	use ndEquations
	implicit none

	print *, "starting tests..."
	call assert_double("D1 test 1", D1_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.4d0, 1d-7)
	call assert_double("D1 test 2", D1_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.4d0, 1d-7)
	call assert_double("D2 test 1", D2_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), -0.00133333d0, 1d-7)
	call assert_double("D2 test 2", D2_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), -0.0213333d0, 1d-7)
	call assert_double("D3 test 1", D3_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.000192d0, 1d-7)
	call assert_double("D3 test 2", D3_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.000192d0, 1d-7)
	
	print *, "all tests were successfull!"
end program tests
