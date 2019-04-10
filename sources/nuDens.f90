program nuDens
	use precision
	use ndConfig
	use ndErrors
	use ndEquations
	implicit none

	call initConfig

	call solver
	call finalresults

	call addToLog("Finished. closing log file.")
	call closeLogFile
end program nuDens
