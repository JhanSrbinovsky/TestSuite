#if defined(C70_1A) || defined(UTILIO) || defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Block Data Subprogram : BLKDATA
!  
!   Purpose : Holds DATA for any variables that are in common blocks,
!             so that they are initialised in only one place.
!  
!   Written by P.Burton
!  

      BLOCK DATA BLKDATA

#include "cenvir.h"
#include "cntl_io.h"
#include "parvars.h"
#include "decomptp.h"
#include "decompdb.h"
#include "decompdt.h"
#include "acparm.h"
#include "comobs.h"
#include "c_gwave.h"
#include "aercmp3a.h"
#include "entcnst.h"
#include "entcnstdt.h"
#include "c_lspdrp.h"
#include "c_lspdrpdt.h"
#if defined(PUMF) || defined(CUMF)
#include "sx_size.h"
#endif
#if defined(MAKEBC)
#include "typsize.h"
#endif

#include "cenvirdt.h"
#include "cntliodt.h"
#include "comobsdt.h"
#include "nstypes.h"
#include "cmaxsize.h"
#include "parcomdt.h"
#if defined(PUMF) || defined(CUMF) || defined(MAKEBC)
#include "sx_dt.h"
#endif
#if !defined(UTILIO) && !defined(FLDIO) && defined(SCMA)
#include "s_soilpr.h"
#include "s_soilpt.h"
#endif

      END BLOCKDATA BLKDATA
#endif
