#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer. 
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution. 
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Fast-j routine for calculating online photolysis rates
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Current Code Owner:       Colin Johsnon/Fiona O'Connor
!                            Oliver Wild
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
      REAL FUNCTION FASTJ_FLINT (tint,t1,t2,t3,f1,f2,f3)
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     Three-point linear interpolation function
!-----------------------------------------------------------------------

      REAL :: tint,t1,t2,t3,f1,f2,f3

      IF (tint <= t2)  THEN
        IF (tint <= t1)  THEN
          fastj_flint  = f1
        ELSE
          fastj_flint = f1 + (f2 - f1)*(tint -t1)/(t2 -t1)
        ENDIF
      ELSE
        IF (tint >= t3)  THEN
          fastj_flint  = f3
        ELSE
          fastj_flint = f2 + (f3 - f2)*(tint -t2)/(t3 -t2)
        ENDIF
      ENDIF

      RETURN
      END FUNCTION FASTJ_FLINT
#endif
