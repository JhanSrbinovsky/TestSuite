! SCLFNC3A defines types of scaling for absorber amounts in two-stream
! radiation code

      ! null scaling function
      INTEGER,PARAMETER:: IP_SCALE_FNC_NULL=0

      ! power law scaling function
      INTEGER,PARAMETER:: IP_SCALE_POWER_LAW=1

      ! power law for p; quadratic for t
      INTEGER,PARAMETER:: IP_SCALE_POWER_QUAD  =2

      ! power law for p; quadratic for t with implicit doppler
      ! correction
      INTEGER,PARAMETER:: IP_SCALE_DOPPLER_QUAD=3

      ! Wenyi scaling for pressure and temperature
      INTEGER, PARAMETER:: IP_SCALE_WENYI=4

      ! number of scaling variables (values set in SCLFND3A)
      INTEGER  N_SCALE_VARIABLE(0: NPD_SCALE_FNC)

! SCLFNC3A end
