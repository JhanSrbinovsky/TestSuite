
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ List parameters and cross reference info for GRIB data

Module rcf_GRIB_lookups_Mod

! Description: This module stores the tables (in the form of derived
!              types) and parameters used to cross reference STASH IDs
!              with those used by other centers.
!
!              To add data for another center just add another column
!              and use 'Is_not_def' for any ID already in the table
!              that any center has no equivalent for.
!
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     12/06/02   Original code. Roddy Sharp (frtz)
!  5.5     07/02/03   De-magic-number-ed code for readability. T.White.
!  6.4     19/12/06   Add ozone. Remove empty field. C Mathison
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!=======================================================================
!  The 'Lists' - used for storage of GRIB header data
!=======================================================================
! These parameters define
! 1) How many linked lists are used to store the fields being read
!    (1 list per variable - multiple 'fields' or 'levels' per list)
! 2) Which list each variable is stored in.
! 3) The numerical order of the lists defines the order in which they
!    Are written to the intermediary UM formatted dump
!
!
! Use all of RCF_STASHCODES for stash magic numbers
Use Rcf_Stashcodes_Mod

! Use all of RCF_ECMWFCODES for ECMWF magic numbers
Use Rcf_ECMWFcodes_Mod

Implicit None
Integer , Parameter :: grib_max_fields    = 19 ! no. of lists

Integer , Parameter              :: grib_LandMask_field    =  1
Integer , Parameter              :: grib_Orog_field        =  2
Integer , Parameter              :: grib_Surf_Temp_field   =  3
Integer , Parameter              :: grib_Surf_Pres_field   =  4
Integer , Parameter              :: grib_U_field           =  5
Integer , Parameter              :: grib_V_field           =  6
Integer , Parameter              :: grib_W_field           =  7
Integer , Parameter              :: grib_Temp_field        =  8
Integer , Parameter              :: grib_Q_field           =  9
Integer , Parameter              :: grib_Rel_Hum_field     = 10
Integer , Parameter              :: grib_Exner_field       = 11
Integer , Parameter              :: grib_Soil_Temp_field   = 12
Integer , Parameter              :: grib_Soil_Moist_field  = 13
Integer , Parameter              :: grib_ozone_field       = 14
Integer , Parameter              :: grib_NOX_field         = 15
Integer , Parameter              :: grib_CH4_field         = 16
Integer , Parameter              :: grib_CO_field          = 17
Integer , Parameter              :: grib_HCHO_field        = 18
Integer , Parameter              :: grib_GO3_field         = 19

Integer , Parameter              :: grib_Misc_field        =  0
! The misc field is used to store header data not stored in the
! other lists

!=======================================================================
!  Parameters : The Originating center no.s
!=======================================================================
! These parameters correspond the the values of the Originating centers
! supplied by the GRIB record - Block 1 - Octet 1. They are used within
! rcf_grib_assign.F90 to select which column in table A (below) to use
! to find the variables STASH code from it's given GRIB parameter code
! (Table 2 or Block 1 - Octet 5)

Integer, Parameter  :: GrbOrigECMWF  = 98  ! ECMWF
Integer, Parameter  :: GrbOrigIntrnl = -1  ! Internal data (created by
                                           ! GRIB code not read)
!=======================================================================
!  Parameters : The Table no.s
!=======================================================================
! These are set so that entries can be distingushed by table number.

Integer, Parameter :: GrbTblECMWFstd  = 128 ! Standard table
Integer, Parameter :: GrbTblECMWFgems = 210 ! GEMS Chemistry table

!=======================================================================
!  Table A: The Parameter ID cross reference table
!  Part 1 : Defining the table
!=======================================================================

Integer, Parameter  :: p_Max_Cols    =  3  ! Number of columns in table
                                           ! 1 per Center (inc UK Met O)

Integer, Parameter  :: p_Max_Rows    = 27  ! Current no. of paramters
                                           ! listed. No real relation to
                                           ! no. of lists

Integer, Parameter  :: Is_not_def    = -1  ! Used when no valid ID is
                                           ! available at the center in
                                           ! question.

Integer , Parameter :: lenDesc       = 20  ! length of description

Type Cross_Ref_Table
  Character(len=lenDesc) :: DescText           ! Textual desc of param
  Integer           :: CrossRefIDs(p_Max_Cols) ! The columns of ID nos
  Integer           :: Table_No(p_Max_Cols)    ! Grib table number for field
  Integer           :: List_No           ! List Number for dynamic lists
End Type Cross_Ref_Table

Type (Cross_Ref_Table)           :: Param_ID_CrossRef(p_Max_Rows)

!=======================================================================
!  Parameters : The ID columns
!=======================================================================
! These parameters define which column (in Table A) to use for a given
! set of variable ID no.s - e.g. column 1 is used to hold STASH no.s
! and column 2 for ECMWF's version of Table 2 parameter IDs

Integer, Parameter  :: p_STASH_IDCol  =  1  ! Column used to store STASH
                                            ! ID no's
Integer, Parameter  :: p_ECMWF_IDCol  =  2  ! Column used to store ECMWF
                                            ! ID no's

! The last column is used to store locally defined ID no.s for data
! which will be generated internally by the GRIB code itself and not
! read in from the GRIB data.
Integer, Parameter  :: p_Intrnl_IDCol =  p_Max_Cols

!=======================================================================
!  Table A: The Parameter ID cross reference table
!  Part 2 : Populating the Table
!=======================================================================
! To add a definition for a new parameter, just append another block of
! data definitions.
! To add a new center to get data from, you must add the ID code to
! each (all) of the CrossRefIDs.

! Example
!Data  Param_ID_CrossRef(1) % CrossRefIDs / -STASH-,-ECMWF-,-INTERNAL-/
!Data  Param_ID_CrossRef(1) % DescText    /'Description here    '/
!Data  Param_ID_CrossRef(1) % List_No     /list no. parameter/

Data  Param_ID_CrossRef( 1) % CrossRefIDs / stashcode_u, &
                                            ECMWFcode_u, &
                                            Is_not_def/
Data  Param_ID_CrossRef( 1) % DescText    /'U Wind              '/
Data  Param_ID_CrossRef( 1) % List_No     /grib_U_field/
Data  Param_ID_CrossRef( 1) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef( 2) % CrossRefIDs / stashcode_v, &
                                            ECMWFcode_v, &
                                            Is_not_def/
Data  Param_ID_CrossRef( 2) % DescText    /'V Wind              '/
Data  Param_ID_CrossRef( 2) % List_No     /grib_V_field/
Data  Param_ID_CrossRef( 2) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef( 3) % CrossRefIDs / stashcode_theta, &
                                            ECMWFcode_T,     &
                                            Is_not_def/
Data  Param_ID_CrossRef( 3) % DescText    /'Temperature         '/
Data  Param_ID_CrossRef( 3) % List_No     /grib_Temp_field/
Data  Param_ID_CrossRef( 3) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef( 4) % CrossRefIDs / stashcode_q,        &
                                            ECMWFcode_spec_hum, &
                                            Is_not_def/
Data  Param_ID_CrossRef( 4) % DescText    /'Specific Humidity   '/
Data  Param_ID_CrossRef( 4) % List_No     /grib_Q_field/
Data  Param_ID_CrossRef( 4) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/


Data  Param_ID_CrossRef( 5) % CrossRefIDs / stashcode_soil_temp, &
                                            ECMWFcode_soil_temp_1, &
                                            Is_not_def/
Data  Param_ID_CrossRef( 5) % DescText    /'Soil Temp Lvl 1     '/
Data  Param_ID_CrossRef( 5) % List_No     /grib_Soil_Temp_field/
Data  Param_ID_CrossRef( 5) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef( 6) % CrossRefIDs / stashcode_soil_temp, &
                                            ECMWFcode_soil_temp_2, &
                                            Is_not_def/
Data  Param_ID_CrossRef( 6) % DescText    /'Soil Temp Lvl 2     '/
Data  Param_ID_CrossRef( 6) % List_No     /grib_Soil_Temp_field/
Data  Param_ID_CrossRef( 6) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef( 7) % CrossRefIDs / stashcode_soil_temp, &
                                            ECMWFcode_soil_temp_3, &
                                            Is_not_def/
Data  Param_ID_CrossRef( 7) % DescText    /'Soil Temp Lvl 3     '/
Data  Param_ID_CrossRef( 7) % List_No     /grib_Soil_Temp_field/
Data  Param_ID_CrossRef( 7) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef( 8) % CrossRefIDs / stashcode_soil_temp, &
                                            ECMWFcode_soil_temp_4, &
                                            Is_not_def/
Data  Param_ID_CrossRef( 8) % DescText    /'Soil Temp Lvl 4     '/
Data  Param_ID_CrossRef( 8) % List_No     /grib_Soil_Temp_field/
Data  Param_ID_CrossRef( 8) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef( 9) % CrossRefIDs / stashcode_tstar,     &
                                            ECMWFcode_skin_temp, &
                                            Is_not_def/
Data  Param_ID_CrossRef( 9) % DescText    /'Skin Temperature    '/
Data  Param_ID_CrossRef( 9) % List_No     /grib_Surf_Temp_field/
Data  Param_ID_CrossRef( 9) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(10) % CrossRefIDs / stashcode_tstar, &
                                            Is_not_def,      &
                                            Is_not_def/
Data  Param_ID_CrossRef(10) % DescText    /'Surface Temperature '/
Data  Param_ID_CrossRef(10) % List_No     /grib_Surf_Temp_field/
Data  Param_ID_CrossRef(10) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/


Data  Param_ID_CrossRef(11) % CrossRefIDs / stashcode_lsm, &
                                            ECMWFcode_lsm, &
                                            Is_not_def/
Data  Param_ID_CrossRef(11) % DescText    /'Land/Sea Mask       '/
Data  Param_ID_CrossRef(11) % List_No     /grib_LandMask_field/
Data  Param_ID_CrossRef(11) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(12) % CrossRefIDs / stashcode_orog,   &
                                            ECMWFcode_geopot, &
                                            Is_not_def/
Data  Param_ID_CrossRef(12) % DescText    /'Geopotential        '/
Data  Param_ID_CrossRef(12) % List_No     /grib_Orog_field/
Data  Param_ID_CrossRef(12) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(13) % CrossRefIDs / stashcode_w, &
                                            ECMWFcode_w, &
                                            Is_not_def/
Data  Param_ID_CrossRef(13) % DescText    /'W Wind              '/
Data  Param_ID_CrossRef(13) % List_No     /grib_W_field/
Data  Param_ID_CrossRef(13) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(14) % CrossRefIDs / stashcode_pstar, &
                                            ECMWFcode_pstar, &
                                            Is_not_def/
Data  Param_ID_CrossRef(14) % DescText    /'Surface Pressure    '/
Data  Param_ID_CrossRef(14) % List_No     /grib_Surf_Pres_field/
Data  Param_ID_CrossRef(14) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(15) % CrossRefIDs / stashcode_q, &
                                            Is_not_def,  &
                                            Is_not_def/
Data  Param_ID_CrossRef(15) % DescText    /'Relative Humidity   '/
Data  Param_ID_CrossRef(15) % List_No     /grib_Rel_Hum_field/
Data  Param_ID_CrossRef(15) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(16) % CrossRefIDs / stashcode_exner, &
                                            Is_not_def, 1/
Data  Param_ID_CrossRef(16) % DescText    /'Exner on P levels   '/
Data  Param_ID_CrossRef(16) % List_No     /grib_Exner_field/
Data  Param_ID_CrossRef(16) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(17) % CrossRefIDs / stashcode_pstar, &
                                            ECMWFcode_log_p, &
                                            Is_not_def/
Data  Param_ID_CrossRef(17) % DescText    /'Log Surface Pressure'/
Data  Param_ID_CrossRef(17) % List_No     /grib_Surf_Pres_field/
Data  Param_ID_CrossRef(17) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(18) % CrossRefIDs / stashcode_soil_moist, &
                                            ECMWFcode_soil_moist_1, &
                                            Is_not_def/
Data  Param_ID_CrossRef(18) % DescText    /'Vol Soil Moist Lvl 1'/
Data  Param_ID_CrossRef(18) % List_No     /grib_Soil_Moist_field/
Data  Param_ID_CrossRef(18) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(19) % CrossRefIDs / stashcode_soil_moist, &
                                            ECMWFcode_soil_moist_2, &
                                            Is_not_def/
Data  Param_ID_CrossRef(19) % DescText    /'Vol Soil Moist Lvl 2'/
Data  Param_ID_CrossRef(19) % List_No     /grib_Soil_Moist_field/
Data  Param_ID_CrossRef(19) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(20) % CrossRefIDs / stashcode_soil_moist, &
                                            ECMWFcode_soil_moist_3, &
                                            Is_not_def/
Data  Param_ID_CrossRef(20) % DescText    /'Vol Soil Moist Lvl 3'/
Data  Param_ID_CrossRef(20) % List_No     /grib_Soil_Moist_field/
Data  Param_ID_CrossRef(20) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(21) % CrossRefIDs / stashcode_soil_moist, &
                                            ECMWFcode_soil_moist_4, &
                                            Is_not_def/
Data  Param_ID_CrossRef(21) % DescText    /'Vol Soil Moist Lvl 4'/
Data  Param_ID_CrossRef(21) % List_No     /grib_Soil_Moist_field/
Data  Param_ID_CrossRef(21) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(22) % CrossRefIDs / stashcode_ozone, &
                                            ECMWFcode_ozone, &
                                            Is_not_def/
Data  Param_ID_CrossRef(22) % DescText    /'Ozone               '/
Data  Param_ID_CrossRef(22) % List_No     /grib_Ozone_field/
Data  Param_ID_CrossRef(22) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           is_not_def/

Data  Param_ID_CrossRef(23) % CrossRefIDs / stashcode_NO2, &
                                            ECMWFcode_NOX, &
                                            Is_not_def/
Data  Param_ID_CrossRef(23) % DescText    /'Nitrogen Oxides     '/
Data  Param_ID_CrossRef(23) % List_No     /grib_NOX_field/
Data  Param_ID_CrossRef(23) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           is_not_def/

Data  Param_ID_CrossRef(24) % CrossRefIDs / stashcode_CH4, &
                                            ECMWFcode_CH4, &
                                            Is_not_def/
Data  Param_ID_CrossRef(24) % DescText    /'Methane             '/
Data  Param_ID_CrossRef(24) % List_No     /grib_CH4_field/
Data  Param_ID_CrossRef(24) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           is_not_def/

Data  Param_ID_CrossRef(25) % CrossRefIDs / stashcode_CO, &
                                            ECMWFcode_CO, &
                                            Is_not_def/
Data  Param_ID_CrossRef(25) % DescText    /'Carbon monoxide     '/
Data  Param_ID_CrossRef(25) % List_No     /grib_CO_field/
Data  Param_ID_CrossRef(25) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           is_not_def/

Data  Param_ID_CrossRef(26) % CrossRefIDs / stashcode_HCHO, &
                                            ECMWFcode_HCHO, &
                                            Is_not_def/
Data  Param_ID_CrossRef(26) % DescText    /'Formaldehyde        '/
Data  Param_ID_CrossRef(26) % List_No     /grib_HCHO_field/
Data  Param_ID_CrossRef(26) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           is_not_def/

Data  Param_ID_CrossRef(27) % CrossRefIDs / stashcode_O3,  &
                                            ECMWFcode_GO3, &
                                            Is_not_def/
Data  Param_ID_CrossRef(27) % DescText    /'GEMS ozone          '/
Data  Param_ID_CrossRef(27) % List_No     /grib_GO3_field/
Data  Param_ID_CrossRef(27) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           is_not_def/

End Module rcf_GRIB_lookups_Mod
