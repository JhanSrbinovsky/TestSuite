!L  COMDECK DOCOBS
!L  --------------
!L
!L    The AC Observation files are read in at the start of an
!L    assimilation (in routine RDOBS2). The observations are then
!L    merged and any observations not required by the assimilation
!L    to follow will be compressed out. The required observations are
!L    then written to a temporary file attached to Unit 15.
!L
!L    On each timestep the observations are read in from Unit 15 and
!L    assimilated. This is done until the AC Observation files
!L    need to be read in again which will normally be every 6 (3) hours
!L    for a global (elf) assimilation. A new set of observations will
!L    then be written to unit 15.
!L
!L    COMOBS contains informtion on the observations in store.
!L
!L    Observations are assigned to types. Each type provides one sort
!L    of information and has its own detailed format. Note that one
!L    observation report may be in several types, for example, a sonde
!L    observation may report temperature, wind and Relative Humidity and
!L    the data will come under types 201, 301 and 401 respectively.
!L
!L    There are NOBTYP types which are listed in OBSTYP. The
!L    The observation types are listed in decreasing order of the
!L    number of data values in NDATAV. Each type OBSTYP(J) has :-
!L
!L        1. NOBS(J) observations stored for it.
!L        2. NDATAV(J) data values in each obs.
!L        3. A total of NDATAV(J)*NOBS(J) data values.
!L
!L    The array MDISPOBT gives the displacement to the first
!L    observation of each type.
!L
!L    ie. MDISPOBT(1)=0 and MDISPOBT(J+1)=MDISPOBT(J)+NOBS(J)
!L
!L    The Total number of observations and data values after merging
!L    the observation files are TNOBS and TNDV respectively.
!L
!L    Maximum allowed values of NOBTYP is NOBTYPMX
!L    Maximum allowed values of NDATAV is NDATAVMX
!L    Maximum allowed values of TNDV   is TNDVMAX
!L
!L  MAXNLEV1    - Maximum no of levels for any obs type in OBSTYP + 1
!L  MDISPOBT    - Offset to first obs for each type in OBS.
!L  MISSD       - Value used for missing data in obs files (=5555.0)
!L  NDATAV      - No of data values for each type in OBSTYP.
!L  NDVHDR      - First NDVHDR data values common to all observations.
!L  NERLEV1     - No of data value corresponding to first error ratio
!L              - for each type in OBSTYP.
!L  NOBLEV      - No of levels for each type in OBSTYP.
!L  NOBS        - No of obs    for each type in OBSTYP.
!L  NOBTYP      - No of observation types in OBSTYP.
!L  OBLAYERB    - Pressure Levels of layer boundaries for each type.
!L              - Set to MISSD for obs types not concerned.
!L  OBLEVELS    - Obs Pressure Levels for each type.
!L              - Set to MISSD for obs types not concerned.
!L  OBLEVTYP    - Level Type obs data in on for each type in OBSTYP.
!L  OBSTYP      - List of observation types in observation files.
!L  OBS_LAT_N   - N ) Boundaries of area for
!L  OBS_LAT_S   - S ) which obs are to fetched for.
!L  OBS_LONG_W  - W )
!L  OBS_LONG_E  - E )
!L  OBS_REF_YY  - Year  ) Reference
!L  OBS_REF_MM  - Month ) Time/Date for observations
!L  OBS_REF_DD  - Day   ) = Start of assimilation.
!L  OBS_REF_HH  - Hour  )
!L  OBS_REF_MIN - Mins  )
!L  OBS_NO_ST   - Offset to first observation number for each type
!L  TIMEINT     - Interval in minutes between reading obs files.
!L  TIMENEXT    - Time in minutes relative to start of assimilation
!L              - when obs files are next to be read.
!L ---------------------------------------------------------------------
