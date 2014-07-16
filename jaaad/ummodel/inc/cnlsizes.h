!LL  COMDECK: CNLSIZES
!LL
!LL  This calls the :                        COMDECK
!LL                 type declarations        TYPSIZE
!LL                 COMMON blocks            TYPSIZE
!LL                 NAMELISTs                NAMSIZE
!LL  for the array size variables defined by the User Interface.
!LL
!LL  Additions need to be added in all 3 positions in TYPSIZE and
!LL  NAMSIZE COMDECKS for sizes obtained from the UI. Other sizes are
!LL  derived from the UI set and only need entries in type declarations
!LL  and common blocks in TYPSIZE.
!LL
!LL  DATA is split into sections denoting which elements belong to
!LL  the set provided in the two UI provided files: STASHC# and NAMLST#.
!LL
!LL  A Common Block is used to allow for data to be read from
!LL  a NAMELIST.
!LL
!LL  This COMDECK should only be called by UM_SHELL and READSIZE.
!LL  Portability pre-processing does not allow arrays to be dimensioned
!LL  via variables on a common block within the same routine, but sizes
!LL  can be passed down as arguments for dynamic allocation at a lower
!LL  level of subroutine.
!LL
!LL  THIS COMDECK CANNOT CALL THE TYPSIZE COMDECKS AS *DEF CALLS WOULD
!LL  CHANGE THE CONTENTS OF THE NAMELIST. ALL NAMELIST ELEMENTS ARE
!LL  PROVIDED BY THE UI.
!LL
!LL  Model            Modification history:
!LL version  Date
!LL  3.2   05/05/93   M.CARTER: DECK CNLSIZES on MOD SET MC050593
!LL                   creation: added for dynamic allocation
!LL
!----------------------------------------------------------------------
!
#include "typsize.h"
#include "namsize.h"
