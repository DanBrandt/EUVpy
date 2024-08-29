@echo off
REM ------------- Delete-HEUVAC-output-files.bat -------------------------
REM .. This batch file is used to safely delete the files created 
REM .. by a HEUVAC run.
REM ..

DEL *.txt 
Del *.dat 

REM fort.* files are sometimes created when HEUVAC is aborted during a run
Del fort.* 

PAUSE

EXIT