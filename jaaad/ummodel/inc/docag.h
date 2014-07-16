!-----------------------------------------------------------------------
!L  COMDECK DOCAG
!L  -------------
!L  COMAG contains variables to define the Analysis Grid.
!L
!L  AG = Analysis Grid
!L  MG = Model    Grid
!L
!L  DEF_AGRES_ROWS - Default Ratio of MG Rows to AG Rows for each group.
!L                 - Use AGRES_ROWS in namelist to change defaults.
!L  DEF_AGRES_PTS  - Default Ratio of MG Points to AG Points for each
!L                 - group. (Equatorwards from AGLATDEC).
!L                 - Use AGRES_ROWS and AGRES_PTS in namelist to
!L                 - change defaults.
!L
!L  COLATITUDES & LONGITUDES are in RADIANS, not DEGREES.
!L
!L  Variables defining Analysis Grid.
!L  ---------------------------------
!L
!L  NROWSAG        - No of rows in AG.
!L  NPTSAG         - No of points in each AG row.
!L  NPTS0AG        - Offset to first point in each row in AG.
!L  NPTSAGMX       - Maximum no of points in AG rows.
!L  NPTSAGMN       - Minimum no of points in AG rows.
!L  MIN_AGPTS      - Minimum no of points on any AG ROW.
!L  LAGNP          - Logical set to T if first row of AG = North Pole
!L  LAGSP          - Logical set to T if last  row of AG = South Pole
!L  STAGROW1       - Stagger of first row of MG to first row of AG,
!L                 - (expressed as fraction of AG row spacing).
!L  STAGPT1        - Stagger of first point of MG to first point of AG,
!L                 - (expressed as fraction of AG row spacing).
!L  ROW1AG         - Co-latitude of first row   of AG.
!L  PT1AG          - Longitude   of first point of AG.
!L  DLATAG         - Co-latitude of row spacing of AG.
!L  DLONGAG        - Longitude spacing for each AG row.
!L  AGLATDEC       - Co-latitude at which no of pts in AG rows start
!L                 - to decrease.
!L  AGROWLEN       - Length of row in radians of latitude.
!L  COSROWAG       - Array storing 1/(DLONGAG(x)*cos(lat of row x))
!L                 - where x is the AG row no.
!L
!L  Variables defining Model Grid.
!L  ------------------------------
!L
!L  ROW1MG         - Co-latitude of first row of MG in use.
!L  ROW1MGTH       - Co-latitude of first row of MG for Theta.
!L  ROW1MGUV       - Co-latitude of first row of MG for Winds.
!L  PT1MGTH        - Longitude of first point on MG for Theta.
!L  PT1MGUV        - Longitude of first point on MG for Winds.
!L  DLATMG         - Co-latitude of row spacing of MG (+ve for N->S)
!L  DLONGMG        - Longitude spacing of MG.
!-----------------------------------------------------------------------
