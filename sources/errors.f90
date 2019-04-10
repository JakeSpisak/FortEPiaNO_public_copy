module ndErrors
	use precision
	use variables
	implicit none

	integer :: totErrors
	character(len=50) :: logFile = "messages.log"
	integer, parameter :: lfu=2411

	contains

	subroutine openLogFile
		write(*,*) "[log] Writing log into: log/", trim(logFile)
		open(unit=lfu, file="log/"//trim(logFile))
		totErrors=0
		call timeToLog
	end subroutine openLogFile

	subroutine closeLogFile
		write(*,*)   "Total errors: ",totErrors
		write(lfu,*) "Total errors: ",totErrors
		call timeToLog
		close(lfu)
	end subroutine closeLogFile

	subroutine renameLogFile(newfilename)
		character(len=*) :: newfilename
		call rename("log/"//trim(logFile), newfilename)
	end subroutine renameLogFile

	subroutine timeToLog
		character(8)  :: date
		character(10) :: time
		integer,dimension(8) :: values
        call date_and_time(VALUES=values)
		write(*  , '("-- ",I0.2,"/",I0.2,"/",I4," - h",I2,":",I0.2,":",I0.2)') values(3), values(2), values(1), values(5),values(6),values(7)
		write(lfu, '("-- ",I0.2,"/",I0.2,"/",I4," - h",I2,":",I0.2,":",I0.2)') values(3), values(2), values(1), values(5),values(6),values(7)
	end subroutine timeToLog

	subroutine addToLog(message)
		character(len=*), intent(in) :: message

		write(*,*) message
		write(lfu,*) message
	end subroutine addToLog

	subroutine printVerbose(message,l)
		character(len=*), intent(in) :: message
		integer :: l

		if (l.le.verbose) &
			write(*,*) message
	end subroutine printVerbose

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
