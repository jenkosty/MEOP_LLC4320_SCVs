# MEOP_LLC4320_SCVs

This repository contains all code necessary to reproduce the results of our study. 

Summary of matlab files:

formatingLLCdata.m - takes matfiles of LLC4320 output data and restructures the data in a more easily usable format. Also includes the snapshot-wide calculation of the Okubo-Weiss parameter. 

creatingLLCsealdata.m - creates a synthetic version of the MEOP dataset using one LLC4320 snapshot

applyingLillyDetectionAlgorithm.m - uses the jLab data analysis package (must be downloaded separately) to identify SCVs in the LLC4320 snapshots (using 2D data). Calls lillySCVdetections.m

SCVdetectionAlgorithmMain.m - Main script to run algorithm. Calls preppingTimeSeriesForDetectionAlgorithm.m. Creates separate files for the anticyclonic tests, cyclonic tests, and extra background tests. 

anticycloneOptimization.m - Takes the results from the Lilly detection algorithm and the anticyclonic+background tests and creates an optimal set of parameters. 

cycloneOptimization.m - Takes the results from the Lilly detection algorithm and the cyclonic+background tests and creates an optimal set of parameters. 




