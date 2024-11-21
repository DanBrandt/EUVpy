@echo off
REM ..............................HEUVAC.bat......................
REM  This batch file runs the HEUVAC program to put solar EUV fluxes into 
REM  user specified bins
REM     Note that this file acts as both a batch file to run the model 
REM     and as the unit 5 (FOR005) input data file that the model reads 
REM     to change the model.

REM --------- Execute the run --------------
CALL debug\HEUVAC.exe HEUVAC.bat

echo ****


EXIT

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!... The information below here is read by HEUVAC to input control
!... parameters
 
$ASSIGN HEUVAC-scratch.TXT FOR005      !.. a temporary scratch file
 
!.. Unit 13 - fluxes summed into Torr et al. 37 bins
$ASSIGN U13-Torr-37-bins-F10-75-85.txt FOR013

!.. Unit 14 - fluxes summed user specified bins
$ASSIGN U14-flux-User-bins-10A-F10-75-85.txt FOR014

!.. Unit 15 - flux weighted cross section in user specified bins
$ASSIGN U15-XS-User-bins-10A-F10-75-85.txt FOR015

! &=START=& DO NOT DELETE or CHANGE THIS FLAG and all "!" in the following lines.                 
 75  ! F107:  Daily F10.7 index
 85  ! F107A: 81-day average F10.7 index
  10   ! BIN_SIZE: size of wavelength bins in Angstrom
 &! &=END=&  DO NOT DELETE or CHANGE THIS FLAG
2
