#if defined(A05_4A) || defined(A05_5A)
!
!    MODEL
!    VERSION  DATE
!      5.3   09-10-01  Reduce the value of DELTHST to speed-up code.
!                                                          S. Cusack
!
      REAL DELTHST  !  DIFFERENCE IN POTENTIAL TEMPERATURE BETWEEN
                    !  LEVELS ABOVE WHICH THE ATMOSPHERE IF ASSUMED
                    !  TO BE TOO STABLE TO CONVECT (K)
      PARAMETER (DELTHST = 0.5)
#endif
