
!+   Represent variation of solar output - its total & distribution.
! Subroutine Interface:
      SUBROUTINE SOLVAR (YEAR, SCS, NBANDS, FB, RA, OZWB1, RESTORE, &
     &     ErrorMessage, ErrorStatus)

      Use scvary_Mod
!      REAL SC                           !  Solar constant
!      PARAMETER ( SC = 1365. )
!20110812: define SC in coupling control
!
! Now get SC from namelist of coupling control
       use auscom_cpl_data_mod, only: SC,VOLCTS_val

      IMPLICIT NONE
! Description:
!  SOLVAR allows for solar variability: its effects on the total solar 
!   output, and on its distribution with wavelength - or at least some 
!   of those things that depend on the latter, namely spectrally 
!   averaged items in the shortwave spectral file.  It assumes that the 
!   empirical (but theoretically plausible) correlations between total 
!   solar output & its spectral distribution over recent solar cycles 
!   due to Lean et al apply across the board, even down to Maunder 
!   minimum, but the time variation follows Solanki & Krivova (2003).  
!   It is written to work only with SW spectral files which match 
!   HadCM3's / HadGEM1's in the relevant aspects.  Only the division 
!   of the total energy between the bands, ozone absorption in band 1 
!   & Rayleigh scattering are worth bothering with (& the last only 
!   because it's so easy).  The ozone absorption can be done just by 
!   tweaking the weights.  This code is for HadGAM1 at 5.5, based on 
!   the 4.4 code for HadAM3.
!
! Method:
!  First time in, save the "standard" values of the quantities 
!   concerned in the spectral file.  
!  Every time the year changes, the quantities to be altered in the 
!   spectral file are varied as functions of the solar "constant".
!  Every call re-set SCS, which is set fresh every timestep by SOLPOS.
!
! Current Code Owner: William Ingram
!                    
! History:
! Version  Date     Comment
! =======  ====     =======
! 5.5-mods 6/04     4.4 mods upgraded.  (William Ingram & Gareth Jones)
!          8/04     Debugged (Peter Stott)
!          9/04     Changes made to ensure that Solanki has a mean
!                   of 1365, the value of SC the solar "constant" 
!                   over all the data  Note that the original
!                   way of specifying Lean was to make sure that the 
!                   average over the Lean data was 1365 which means that Lean
!                   is slightly higher as Lean includes the Maunder min
!                   Nobody knows what the 
!                   absolute values of solar irradiance are anyway.
!                   (Peter Stott)
!    
! Code description:
!   FORTRAN 77 + common extensions also in fortran 90.
!   This code is written to UM programming standards version 6
! 
! System component covered:  not applicable
! System task:               not applicable
!
! Declarations: these are of the form:-
!     INTEGER      ExampleVariable  !Description of variable

! 
! Global variables (*CALLed common blocks etc.)


! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER YEAR,                					&
					!  Current year
     &        NBANDS            					
					!  Number of bands in spectral file

      LOGICAL RESTORE              !  Is this call to restore the 
!     original standard values (rather than to set time-varying ones) ? 

      CHARACTER*80 ErrorMessage    !
      CHARACTER*80 errmessage2
!   Scalar arguments with intent(InOut):

      REAL SCS                     
!     Solar constant scaling factor, taking into consideration the
!     changes in solar irradiance & the earth-sun distance.
!   Array  arguments with intent(InOut):

      REAL FB(6),							&
		                  ! Fraction of insolation in each band
     &     RA(6),                 					& 
				  ! Rayleigh coefficients for each band
     &     OZWB1(5)               					
				  ! Ozone esft weights in band 1
!   ErrorStatus:
      INTEGER      ErrorStatus     !+ve = fatal error    

! Local parameters:
      REAL FC2B1D,                 					&
				  ! Used to convert FCSC to BAND1D
     &     FBSEN1,                					& 
				  ! Sensitivity of FB(1) to BAND1D
     &     FBSEN2,                					& 
				  ! Sensitivity of FB(2) to FCSC
     &     RASEN1,              					& 
				  ! Sensitivity of RA(1) to BAND1D
     &     RASEN2              						 
			          ! Sensitivity of RA(2) to FCSC 

      PARAMETER ( FC2B1D = 213. )
      PARAMETER ( FBSEN1 = 5.212, FBSEN2 = 0.1458 )
      PARAMETER ( RASEN1 = 4.360, RASEN2 = 0.2185 )
! Local scalars:
      REAL SOV                     ! Solar "constant" for the year
      REAL MNFB(6), MNRA(2), MNOZ(5)
      ! Reference values for FB, RA and OZWB1.
      REAL FCSC,                  					& 
				  ! Fractional change in solar constant
     &     SCS_LASTIN,   						&
				  ! Value SCS was re-set away from
     &     BAND1D                  

!     Factor representing the band-1 quantities' quadratic dependence
!     on solar constant (in band 2 they vary linearly).

!      REAL SCVARY (601)
!      INTEGER I, YSOL
!     define solar constants for year 2301 and after to follow 11 years solar cycle
      REAL,PARAMETER :: sc_cycle(1:11) =(/1366.53954,1366.30442,1366.02501,   &
          1365.78996,1365.67390,1365.71370,1365.89665,1366.16474,&
          1366.43289,1366.61586,1366.65563/)

      INTEGER I
! DATA statements replaced with a read from a text file later
! CDJ. 4/6/09

      REAL OZSEN(5)                ! Sensitivities of OZWB1 to BAND1D
      DATA OZSEN / -3.65009, -3.84821, -2.57906, 2.44108, 12.7060 /

      LOGICAL FIRST                ! First time in ?
      DATA FIRST  / .TRUE. /
      INTEGER LASTYR               ! Year last time through
      DATA LASTYR / -99999 /

      SAVE FIRST, LASTYR, MNRA, MNFB, MNOZ, FCSC, SCS_LASTIN
!- End of header  

      IF ( FIRST .AND. RESTORE ) THEN
        ErrorMessage =							&
     &  ' SOLVAR: Asked to restore original values on 1st CALL'
        WRITE (6, *) ErrorMessage
        ErrorStatus = 7
        RETURN
      ENDIF

      IF ( FIRST ) THEN

!  Check only on the number of bands & that there seem to be 5 esft 
!   terms for the first gas in the first band - any more complex test on, 
!   say, the total size of the file, or band limits could be confused by
!   adding more aerosols or re-labelling the (pseudo-)band limits.

        IF ( SUM(OZWB1) .GT. 1.0001 .OR. SUM(OZWB1) .LE. 0.9999 	&
     &                               .OR. NBANDS .NE. 6 ) THEN
          ErrorMessage = ' SOLVAR: Inappropriate spectral file'
          WRITE (6, *) NBANDS, OZWB1
          WRITE (6, *) ErrorMessage
          ErrorStatus = 8
          RETURN
        ENDIF

!       Store reference values.
        MNFB = FB
        MNRA = RA(1:2)
        MNOZ = OZWB1
        FIRST = .FALSE.

      ENDIF

      IF ( YEAR .LT. 1700 .OR. YEAR .GT. 2300 ) THEN
ErrorMessage = ' SOLVAR: years 1700-2300 covered in file, SOLVAR should be '
errmessage2='turned off before 1700. it cycles from year 2301 for every 11 years'

        WRITE (6, *) YEAR
        WRITE (6, *) ErrorMessage//errmessage2, sc_cycle(mod(YEAR-2301,11)+1)

        !ErrorStatus = 9
        !RETURN
      ENDIF

      IF ( RESTORE ) THEN

         FB = MNFB 
         RA(1:2) = MNRA 
         OZWB1 = MNOZ
         SCS = SCS_LASTIN
! Do this, rather than just re-divide by 1+FCSC, in hope of keeping
!  bit-reproducibility.

       ELSE

         IF ( YEAR .NE. LASTYR ) THEN

!          Calculate this year's solar constant & thence FCSC and BAND1D
!          IF YOU WANT TO INCREASE THE SIZE OF THE SOLAR CONSTANT VARIATION
!          DO NOT DO IT HERE AS IT WILL GIVE SILLY SPECTRAL STUFF.
           IF ( YEAR .GT. 2300 ) THEN
             SOV = sc_cycle(mod(YEAR-2301,11)+1)
           ELSE
             SOV = SCVARY(YEAR-1699)
           ENDIF

!          Set SC as the mean over 1700-2003 (original Met Office code scaled
!          this to 1365).
           SC = sum(SCVARY(1:304))/304.
           FCSC = SOV / SC - 1.
           write(6,*) "Solar constant SC, SOV", sc, sov, fcsc
           BAND1D = FCSC * ( 1. + FC2B1D * FCSC )

!          Set FB(1,2) using BAND1D & FCSC and the rest by normalization.
           FB(1) = ( 1. + FBSEN1 * BAND1D ) * MNFB (1)
           FB(2) = ( 1. + FBSEN2 * FCSC )   * MNFB (2)
           FB(3:6) = MNFB(3:6)   *  ( 1. - FB(1) - FB(2) ) 		&
     &     / SUM(MNFB(3:6))

!          Set RA(1,2) using BAND1D & FCSC: the rest are insignificant.
           RA(1) = ( 1. + RASEN1 * BAND1D ) * MNRA (1)
           RA(2) = ( 1. + RASEN2 * FCSC )   * MNRA (2)

!          Set OZWB1 using BAND1D but must then re-normalize.
           OZWB1 = ( 1. + OZSEN * BAND1D ) * MNOZ
           OZWB1 = OZWB1 / SUM(OZWB1)

           LASTYR = YEAR
 
         ENDIF

         SCS_LASTIN = SCS
!        TO AMPLIFY THE SOLAR CONSTANT VARIATION PUT A FACTOR IN HERE:
         SCS = SCS * ( 1. + FCSC )
 
      ENDIF

      RETURN
      END

