program nuDens
	use precision
	use ndConfig
	use ndErrors
	use ndEquations
	implicit none

	call initConfig

	call init_interp_jkyg12
	call init_interp_dme2_e
	call init_interp_ElDensity

	call solver
	call finalresults

	call addToLog("Finished. closing log file.")
	call closeLogFile
	call renameLogFile(trim(outputFolder)//'/messages.log')
end program nuDens
