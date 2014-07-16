!*L --------------------- Comdeck: CENVIR   ----------------------------
!
!   Purpose: COMDECK defining Character enviroment variables used
!            by portable IO to open and close files
!
!   Author : R A Stratton      Date : 22/10/92
!
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL 3.2     28/05/93  Add file BAS_IND at unit number 58. M.Carter.
!LL
!LL 3.1     15/01/93  Increase no. of unit nos. from 1-99  to 1-199
!LL                   Dummy names have been set up temporarily for
!LL                   files 104-119. R.Rawlins
!LL
!LL 3.3     09/03/94  Separate data statements into COMDECK
!LL                   CENVIRDT. Also includes mods originally
!LL                   in RB221193 : Add source terms at unit no.110
!LL                   P.Burton and R.T.H Barnes
!LL

!   Vn3.0  12/02/93 - Environment variables PERTURB and TRANSP put in
!                     positions 37 and 97 respectively in character
!                     array FT_ENVIRON, and the appropriate character
!                     lengths put in LEN_FT_ENVIR. C. S. Douglas
!
! Type declarations
!
      CHARACTER*8 FT_ENVIRON(199)  ! Array holding enviroment variables
!                                  for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
!


!
!Common Blocks for character and integer arrays
!
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
!
