C.................... HEUVAC .......................
C.. This a test driver program for the high resolution 
C.. EUVAC model (HEUVAC). 
C.. You can specify the wavelength bins for the summing
C...... Written by Phil Richards, March 2004
      IMPLICIT NONE
	INTEGER I,J,K             !.. loop control variables
	INTEGER KMAX              !.. index of last wavelength bin
	INTEGER IDIM              !.. dimension for the bin arrays
	PARAMETER (IDIM=3333)
	REAL F107,F107A           !.. Solar activity indices
	REAL BIN_SIZE             !.. Width of wavelength bins
	REAL BINLAM(IDIM+1)             !.. Bin wavelengths
	REAL BIN_ANGS(IDIM+1)           !.. Array for the binned fluxes(photons/cm2/s)
      REAL WL_FAC                   !.. Unit Conversion factor
      REAL XS_OUT(15,IDIM)          !.. Weighted cross section
      REAL TORR_WL(37),TORR_FLX(37) !.. Wavelengths and EUV flux in Torr 37 bins

      WRITE(6,676)
 676  FORMAT(/6X,'... Test Driver for HEUVAC solar flux model......')

      !... opens files that are assigned in batch (.bat) file
      CALL OPEN_FILE()

 10   READ(5,*) F107             !.. daily F107 index
      READ(5,*) F107A            !.. 81 day average F107 index
      READ(5,*) BIN_SIZE         !.. Bin size in Angstrom

      !.. Set up the wavelength bins in Angstroms to sum fluxes
      BINLAM(1)=0
	DO K=1,IDIM
	  BINLAM(K+1)=BINLAM(K)+BIN_SIZE
	  IF(BINLAM(K+1).GT.2000) GOTO 15
	ENDDO
 15   CONTINUE
      KMAX=K-1  !.. the index of the last bin
      
      !.. Call the high resolution EUVAC model
      CALL HEUVAC(IDIM,F107,F107A,KMAX,BINLAM,BIN_ANGS,XS_OUT,
     >  TORR_WL,TORR_FLX)

	WRITE(14,*) ' EUV fluxes summed into user specified bins'
	WRITE(14,94) 
  94  FORMAT(' Bin     WL_LO    WL_HI    Ave_A    Ave_nm'
     >  ,2X,'Flux(ph/cm2/s)   Flux(W/m2/nm)')

      !... Output the solar EUV fluxes
	DO K=1,KMAX
	  !.. conversion factor for W/m2/del(lamba)
	  WL_FAC=5.03556E+11*(BINLAM(K)+BINLAM(K+1))/20.0
	  WRITE(14,'(I5,4F9.2,1P,9E15.3)') K,BINLAM(K),BINLAM(K+1),
     >    (BINLAM(K)+BINLAM(K+1))/2.0,
     >    (BINLAM(K)+BINLAM(K+1))/20.0,BIN_ANGS(K),
     >     BIN_ANGS(K)/WL_FAC/(BIN_SIZE/10)
	ENDDO

      !... Output the flux weighted cross photoionization sections
	WRITE(15,90)
 90   FORMAT(' Flux weighted photoionization cross sections. Based on'
     > ' Fennelly and Torr, Atomic Data and Nuclear Data Tables, 51,'
     > ,/2X,' 321-363, 1992. Units are 1.0E-18 cm-2. BINLO, BINHI ='
     > ' wavelength limits in Angstrom'
     > ,/3X,'BINLO   BINHI   OP4S   OP2D   OP2P   OP2E   OP4E'
     > '  O+_tot N2_AB   N2+    N2_N+  N2_ION O2_AB  O2+   O2_O+'
     > '  O2_ION   N+')

	DO K=1,KMAX-1
	  WRITE(15,'(2F8.2,22F7.2)') BINLAM(K),BINLAM(K+1),
     >     (XS_OUT(J,K),J=1,15)
      ENDDO

      !.. Write out the Torr et al. fluxes in 37 wavelength bins
	WRITE(13,*) '   #  Ave-WL   Flux(ph/cm2/s) in Torr 37 wavelength bins'
	DO K=1,37
	  WRITE(13,'(I6,F8.2,1P,E10.2)') K,TORR_WL(K),TORR_FLX(K)
      ENDDO

	STOP
	END
C:::::::::::::::::::::OPEN_FILE::::::::::::::::::::::::::::::::::::
C.... This routine is responsible for processing the control file(eg.HEUVAC.bat)
C.... which is specified in first argument.
C.... after opening the control file,this routine calls GET_INOUT_FILES()
C.... to extract all I/O files from the control file and open them with associated 
C.... unit numbers.
C.... Finally,this routine calls GET_CONTR_PARAMS() to extract some paramter 
C.... values from control file and write them to a file associated with unit 5.
C....

      SUBROUTINE OPEN_FILE() 
	    
C	USE DFLIB                     !use visual digital fortran library
	                                   
	INTEGER       INPUT           !Input unit for control file
	INTEGER       OUTPUT          !Output unit for log file
	INTEGER       DATAUNIT        !Unit number for data file
	CHARACTER*60  ARG_1           !First command line argument
	CHARACTER*60  FILENAME        !User input file name
	CHARACTER     ANSWER1         !User input for "y"or "n"
	CHARACTER*60  ARG_2           !Second command line argument
	INTEGER       FLAG            !0:polwinds;1:equwinds;2:dhdu
	
      !... the following three lines are not standard Fortran 
      INTEGER(2) n1, status         !n1:indicate which command line argument
                                    !status: needed by GETARG()
      !... get first command line argument
      n1 = 1                     
C      CALL GETARG(n1, ARG_1, status)
      ARG_1 = 'HEUVAC.bat'
	FILENAME=ARG_1 
	!... get second command line argument
	n2 = 2
C     CALL GETARG(n2,ARG_2, status)
        ARG_2 = 'NONE'
	FLAG=0
	IF(INDEX(ARG_2,'EQWIND').NE.0) THEN
	    FLAG=1
	ELSE IF (INDEX(ARG_2,'DHDU-NORTH').NE.0) THEN
	    FLAG=2
	ELSE IF (INDEX(ARG_2,'DHDU-SOUTH').NE.0) THEN
	    FLAG=3
	ENDIF  
		         	      
      !... open input control file
  	INPUT=49                      
      OPEN(INPUT,FILE=FILENAME,STATUS='OLD',ERR=350) 	
      !...extract I/O files and their unit numbers from .RUN file 
	
	CALL GET_INOUT_FILES(INPUT,FLAG)
	
      !... rewind .RUN file and extract parameters from it and  
	!... store these parameters into the file which is associated with unit 5
C	REWIND(INPUT)
C      CALL GET_CONTR_PARAMS(INPUT,5)
C	REWIND(5)
	GO TO 400
	    
 350  WRITE(*,*) 'FILE', FILENAME ,'Does NOT EXIST!'
 400 	CLOSE(INPUT)
      END


C:::::::::::::::::::::GET_INOUT_FILES::::::::::::::::::::::::::::::::::::
C.... This routine is used to extract input and output files from the
C.... control file (eg.HEUVAC.bat) and open them.
C.... This routine is called in the OPEN_FILE() routine. 
C.... ASSUME: each I/O file and its unit number is contained in line starting 
C.... with "$" and there is an "assign" or "ASSIGN" in front of the line.
C....
      SUBROUTINE GET_INOUT_FILES(INUNIT,FLAG)     
	                                                                              	                                       
      INTEGER        INUNIT                  !... Input unit in
	INTEGER        OUTUNIT                 !... Output unit for debugging
	INTEGER        FLAG                    !... 0:polwinds;1:equwinds;2:dhdu                     
      CHARACTER*200  CH_ARRAY                !... Character array
	CHARACTER*60   DIR                     !... Directory name 
      CHARACTER*60   FILENAME                !... File name extracted
	CHARACTER*10   UNIT                    !... Unit string extracted
	INTEGER        UNITNO                  !... Unit number
	INTEGER        LENGTH                  !... length of unit string 
	INTEGER        DIRFLAG                 !... flag to indicate directory
	INTEGER        DIRLEN                  !... directory length
	
        DATA OUTUNIT/0/         !..change to 100 for debug writes

        OPEN(5, FILE = 'HEUVAC-scratch.TXT')
        OPEN(13, FILE = 'Torr-37-bins.txt')
        OPEN(14, FILE = 'flux-User-bins-10A.txt')
        OPEN(15, FILE = 'XS-User-bins-10A.txt')
C
C
C
C        
C 50   READ(INUNIT,'(A200)',END=60) CH_ARRAY  !... read each line from input file
C      
C      DIRFLAG=0                              !... initial directory flag=0
C	
C      !... skip all starting blanks in each line and blank lines
C      DO I=1,100
C      IF (CH_ARRAY(I:I).NE.' ') GO TO 55
C      ENDDO
C	GO TO 50
C      !... identify the lines containing '$'
C 55   IF(CH_ARRAY(I:I).EQ.'$') THEN
C	   DO I1=I+1,100
C	      IF (CH_ARRAY(I1:I1).NE.' ') GO TO 57 
C	   ENDDO
C	
C	!... check the line contains "assign" or not
C 57      IF(CH_ARRAY(I1:I1+6).EQ.'ASSIGN'.
C     >      OR.CH_ARRAY(I1:I1+6).EQ.'assign') THEN
C               DO J=I1+7,100 
C	             IF (CH_ARRAY(J:J).NE.' ') GO TO 80
C               ENDDO
C	!... check the line contains "IF", "EQWIND","THEN" or not
C	   ELSE IF (FLAG.EQ.1.AND.INDEX(CH_ARRAY,'IF').NE.0.AND.
C     >            INDEX(CH_ARRAY,'EQWIND').NE.0.AND.
C     >	        INDEX(CH_ARRAY,'THEN').NE.0.AND.
C     >		    INDEX(CH_ARRAY,'ASSIGN').NE.0) THEN
C	            I1=INDEX(CH_ARRAY,'ASSIGN')
C
C	            DO J=I1+7,100 
C	                IF (CH_ARRAY(J:J).NE.' ') GO TO 80
C                  ENDDO
C               
C 80            DO K=J,100
C      !... use ":" to identify there is directory info before file name
C                   IF (CH_ARRAY(K:K).EQ.':')THEN
C		            DIRFLAG=1
C		            GO TO 86
C                   ENDIF
C               ENDDO
C	         K=J
C	!...  find out directory name and length
C
C 86            IF(K.GT.J) THEN
C                   DIR=CH_ARRAY(J:K-1)
C                   DIRLEN=K-J
C                   IF (OUTUNIT.EQ.100) WRITE(OUTUNIT,*)'DIR=',DIR
C	         ELSE
C	             DIR=''
C	             DIRLEN=0
C	         ENDIF
C               DO K1=K,100
C	             IF (CH_ARRAY(K1:K1).EQ.' ')GO TO 90
C               ENDDO
C
C	!...  add directory name to the file name    	      
C 90	         IF (DIRFLAG.EQ.1) THEN   
C		         FILENAME(1:DIRLEN)=DIR
C		         FILENAME(1+DIRLEN:1+DIRLEN)='\'
C		         FILENAME(1+DIRLEN+1:)=CH_ARRAY(K+1:K1-1)
C               ELSE
C		         FILENAME(1:)=CH_ARRAY(K:K1-1)
C	         ENDIF
C               IF (OUTUNIT.EQ.100) WRITE(OUTUNIT,*)'FILENAME=',FILENAME 
C	         DO L=K1,100
C	             IF (CH_ARRAY(L:L).NE.' ') GO TO 100
C	         ENDDO
C 100           DO M=L,100
C                   IF (CH_ARRAY(M:M).EQ.' ') GO TO 110
C               ENDDO
C	!... extract unit number
C 110           UNIT=CH_ARRAY(L+3:M-1)
C               LENGTH=M-1-L-3+1
C               IF (OUTUNIT.EQ.100) WRITE(OUTUNIT,*) 'UNIT=',UNIT
C
C	!... convert unit from string to integer
C	         CALL CHAR_TO_INTEGER(UNIT,LENGTH,UNITNO)        
C      !... open file with assigned unit number. Note unit 2 and 99
C	!... is assigned to binary files
C	         IF (UNITNO.EQ.2) THEN
C       	         OPEN(UNITNO,FILE=FILENAME,FORM='UNFORMATTED')
C               ELSE 
C	             IF (UNITNO.EQ.99) THEN
C	                 OPEN(UNITNO,FILE=FILENAME,STATUS='OLD',
C     >					 FORM='UNFORMATTED',ERR=300)
C	                 go to 50
C 300                   continue 
C                   ELSE
C	                 OPEN(UNITNO,FILE=FILENAME)
C	             ENDIF
C	         ENDIF
C       
C	         
C         ENDIF   !end of if with checking "assign"
C         	      
C      ENDIF      !end of if with checking "$"
C     
C      GOTO 50
 60   RETURN
      END

C
CC:::::::::::::::::::::CHAR_TO_INTEGER::::::::::::::::::::::::::::::::::::
CC.... This routine is used to convert string to integer.
CC.... It is called in GET_INOUT_FILES().
C
C
C      SUBROUTINE CHAR_TO_INTEGER(INCHAR, !... input string
C     >                           LENGTH, !... length of INCHAR
C     >                           OUTINT) !... integer out
C
C	CHARACTER*10 INCHAR    ! input string
C	INTEGER      OUTINT    ! converted integer
C      INTEGER      DIGIT     ! each digit in string
C	CHARACTER*10 NEWCHAR   ! temp variable
C	M=LENGTH
C	OUTINT=0
C	NEWCHAR=INCHAR
C	!... get each digit and convert to integer
C	DO I=1,LENGTH
C        DIGIT=ICHAR(NEWCHAR)-48
C	  M=M-1	  
C	  OUTINT=OUTINT+DIGIT*(10**M)        
C	  NEWCHAR=INCHAR(I+1:LENGTH)
C      ENDDO
C	
C	RETURN
C	END
C
CC:::::::::::::::::::::GET_CONTR_PARAMS::::::::::::::::::::::::::::::::::::
CC.... This routine is used to extract control parameter values from
CC.... the control file (HEUVAC.bat) and write them into files
CC.... instead of input from keyboard. 
CC.... Use &=START=& and &=END=& to identify the range of data.
CC.... Use ! to identify the data in a line
CC.... This routine is called in OPEN_FILE() routine.
C
C      SUBROUTINE GET_CONTR_PARAMS(INUNIT,OUTUNIT)
C	INTEGER INUNIT                 ! Input unit 
C	INTEGER OUTUNIT                ! Output unit 
C         
C      CHARACTER*200 LINE             ! Line buffer
C      
C	!... read each line from input file until find out flag
C  50	READ(INUNIT,'(A200)',END=100) LINE
C      DO I=1,200
C	  IF(INDEX(LINE,'&=START=&').NE.0) GO TO 60 
C	ENDDO
C	GO TO 50
C	!... parse the line and pick out the data
C  60  READ(INUNIT,'(A200)',END=100) LINE
C      
C	!... if meet end flag,return
C      IF(INDEX(LINE,'&=END=&').NE.0) GO TO 100
C	!... use "!" as a flag to identify data
C      I1=INDEX(LINE,'!')
C	!... output data into output file
C	IF(I1.NE.0) THEN
C	   WRITE(OUTUNIT,*) LINE(1:I1-1)
C	ENDIF    
C	GO TO 60  
C  100 RETURN 
C      END											      
	           
