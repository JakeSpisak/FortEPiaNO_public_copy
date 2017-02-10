module ndErrors
	use precision
	implicit none
	
	integer :: totErrors
	character(len=50), parameter :: logFile = "messages.log"
	integer, parameter :: lfu=2411
	
	contains
	
	subroutine openLogFile
		open(unit=lfu, file=trim(logFile))
		totErrors=0
	end subroutine openLogFile
	
	subroutine closeLogFile
		write(*,*)   "Total errors: ",totErrors
		write(lfu,*) "Total errors: ",totErrors
		close(lfu)
	end subroutine closeLogFile
	
	subroutine addToLog(message)
		character(len=*), intent(in) :: message
		
		write(*,*) message
		write(lfu,*) message
	end subroutine addToLog
	
	subroutine error(message)
		character(len=*), intent(in) :: message
		
		totErrors = totErrors+1
		write(*,*)   "[Error] "//message
		write(lfu,*) "[Error] "//message
	end subroutine error
	
	subroutine criticalError(message)
		character(len=*), intent(in) :: message
		
		totErrors = totErrors+1
		write(*,*)   "[Critical error] "//message
		write(lfu,*) "[Critical error] "//message
		
		call closeLogFile
		call exit()
	end subroutine criticalError

end module
