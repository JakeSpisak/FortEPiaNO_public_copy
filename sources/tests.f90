subroutine assert_double(testname, num1, num2, tol)
	use precision
	implicit none
	character(len=*) :: testname
	real(dl) :: num1, num2, tol

	if (abs(num1 - num2) .lt. tol) then
		write(*, fmt="(a)", advance="no") "."
	else
		print *, testname, " failed"
		write (*,"(*(E15.7))") num1, num2, num1 - num2, tol
		call exit()
	end if
end subroutine assert_double

program tests
	use precision
	use ndConfig
	use ndErrors
	use ndEquations
	implicit none

	real(dl), dimension(2) :: temp_v2

	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "Starting tests"

	write(*,*) ""
	write(*,"(a)") "D_i functions (12 tests)"
	call assert_double("D1 test 1", D1_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.4d0, 1d-7)
	call assert_double("D1 test 2", D1_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.4d0, 1d-7)
	call assert_double("D1 test 3", D1_full(0.01d0,5.d0,2.6d0,2.41d0), 0.04d0, 1d-7)
	call assert_double("D1 test 4", D1_full(10.03d0,5.4d0,8.8d0,6.63d0), 21.6d0, 1d-4)
	call assert_double("D2 test 1", D2_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), -0.00133333d0, 1d-7)
	call assert_double("D2 test 2", D2_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), -0.0213333d0, 1d-7)
	call assert_double("D2 test 3", D2_full(0.01d0,5.d0,2.6d0,2.41d0), -1.333333333d-6, 1d-12)
	call assert_double("D2 test 4", D2_full(10.03d0,5.4d0,8.8d0,6.63d0), -209.952d0, 1d-4)
	call assert_double("D3 test 1", D3_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.000192d0, 1d-7)
	call assert_double("D3 test 2", D3_full(0.4d0, 0.2d0, 0.3d0, 0.1d0), 0.000192d0, 1d-7)
	call assert_double("D3 test 3", D3_full(0.01d0,5.d0,2.6d0,2.41d0), 2.50454d-5, 1d-10)
	call assert_double("D3 test 4", D3_full(10.03d0,5.4d0,8.8d0,6.63d0), 22692.2d0, 3d-1)

	write(*,*) ""
	write(*,"(a)") "Pi(yi,yj) functions (12 tests)"
	call assert_double("Pi_1_12 test 1", PI1_12_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.00933333d0, 1d-7)
	call assert_double("Pi_1_12 test 2", PI1_12_full(0.4d0, 0.3d0, 0.2d0, 0.1d0), 0.0893333d0, 1d-7)
	call assert_double("Pi_1_13 test 1", PI1_13_full(0.1d0, 0.2d0, 0.3d0, 0.4d0), 0.0106667d0, 1d-7)
	call assert_double("Pi_1_13 test 2", PI1_13_full(0.4d0, 0.3d0, 0.2d0, 0.1d0), 0.0106667d0, 1d-7)
	temp_v2 = PI2_nn_f(0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.3d0, 0.4d0)
	call assert_double("Pi_2_14 test 1", temp_v2(1), 0.00267733d0, 1d-7)
	call assert_double("Pi_2_13 test 1", temp_v2(2), 0.000810667d0, 1d-7)
	temp_v2 = PI2_nn_f(0.4d0, 0.3d0, 0.2d0, 0.1d0, 0.2d0, 0.1d0)
	call assert_double("Pi_2_14 test 2", temp_v2(1), 0.00267733d0, 1d-7)
	call assert_double("Pi_2_13 test 2", temp_v2(2), 0.000810667d0, 1d-7)
	temp_v2 = PI2_ne_f(0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.2d0, 0.4d0)
	call assert_double("Pi_2_14 test 1", temp_v2(1), 0.00267733d0, 1d-7)
	call assert_double("Pi_2_12 test 1", temp_v2(2), 0.00427733d0, 1d-7)
	temp_v2 = PI2_ne_f(0.4d0, 0.3d0, 0.2d0, 0.1d0, 0.3d0, 0.1d0)
	call assert_double("Pi_2_14 test 2", temp_v2(1), 0.00267733d0, 1d-7)
	call assert_double("Pi_2_12 test 2", temp_v2(2), 0.00427733d0, 1d-7)



	write(*,*) ""
	write(*,*) ""
	write(*,"(a)") "all tests were successfull!"
end program tests
