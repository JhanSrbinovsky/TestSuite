#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      MODULE HEIGHT_MOD
      USE IN_STOCHEM_GRD

      PRIVATE

!kk   select fixed ind_max(.true.) or automatic detection (.false.)

      LOGICAL,PARAMETER                      :: fixed_max=.false.

      INTEGER,PARAMETER                      :: ind_max_fixed=10000
      INTEGER,SAVE                           :: ind_max   !max index
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)  :: ind_array !index array

      REAL,SAVE                              :: scale_factor
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)     :: eta_theta_scale
      LOGICAL,SAVE                           :: first=.true.

      INTERFACE HEIGHT_ETA_THETA
         MODULE PROCEDURE HEIGHT_ETA_THETA
         MODULE PROCEDURE HEIGHT_ETA_THETA_TODo
      END INTERFACE ! HEIGHT_ETA_THETA

      INTERFACE HEIGHT_ETA_RHO
         MODULE PROCEDURE HEIGHT_ETA_RHO
      END INTERFACE ! HEIGHT_ETA_RHO

      INTERFACE HEIGHT_ETA_STOCHEM
         MODULE PROCEDURE HEIGHT_ETA_STOCHEM
      END INTERFACE ! HEIGHT_ETA_STOCHEM

      PUBLIC HEIGHT_INI
      PUBLIC HEIGHT_ETA_THETA, HEIGHT_ETA_RHO, HEIGHT_ETA_STOCHEM

      CONTAINS

      SUBROUTINE HEIGHT_INI
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  14/01/04   Created. K. Ketelsen
!  6.1  07/09/04   Added '&' for F90 compatability as needed.
!                  M.G. Sanderson
!
!kk   The HEIGHT_Eta_theta uses table lookup to reduce operation count
!kk   In the original version, the height index is computed in a search
!kk   loop. starting with index one, the loop was exit as soon as the
!kk   level was found. An avarage of nmetlev/2 stride through the
!kk   search loop was required. This code did not vectorize.
!kk
!kk   Using table lookup, an index array (ind_array) of size ind_max
!kk   has been created. The values of eta_theta are scaled to a
!kk   maximum of ind_max. ind_array is filled that the index matches
!kk   the relation to the eta_theta values.

!kk   ind_table=POS(3,i)*scale_factor gives the index in the lookup tabl
!kk   ind_array(ind_table) gives a first guess of the index in eta_theta
!kk
!kk   To find the exact ind_theta, POS(3,i) has to be compared with
!kk   eta_theta(ind) and eta_theta(ind+1).

!kk   ind_max has to be > 3*(1./min(eta_theta(i)-eta_theta(i-1))), i.e.
!kk   for every eta_theta value at least 3 ind_array entries are
!kk   required.

        IMPLICIT NONE
!
        INTEGER              :: i,j
        REAL                 :: xmin

        IF (first) THEN                ! call only once
          IF (fixed_max) THEN
             ind_max=ind_max_fixed
          ELSE
             xmin = eta_theta(nmetlev)
             DO i=1,nmetlev
               xmin = min(xmin,eta_theta(i)-eta_theta(i-1))
             END DO
             xmin = xmin/eta_theta(nmetlev)         !normalize to 1.0
             ind_max = 4*(1./xmin)                  !3* might be ok
          END IF
          ALLOCATE(ind_array(0:ind_max))
          ALLOCATE(eta_theta_scale(0:nmetlev))

          scale_factor = ind_max/eta_theta(nmetlev)

          eta_theta_scale=eta_theta*scale_factor

          WRITE(6,*) 'Initialize HEIGHT table lookup'
          WRITE(6,*) 'ind_max = ',ind_max,'  scale_factor = ',          &
     &      scale_factor

!         Create lookup table

          j = 0
          DO i=0,ind_max-1
            IF (REAL(i) >= eta_theta_scale(j)) THEN
              j = j+1
            END IF
            ind_array(i) = j-1
          END DO
          first = .false.
        END If
        ind_array(ind_max) = nmetlev-1

        RETURN
      END SUBROUTINE HEIGHT_INI

      SUBROUTINE HEIGHT_ETA_THETA(pos,ka)
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  14/01/04   Created. K. Ketelsen
!  6.1  07/09/04   Added '&' for F90 compatability as needed.
!                  M.G. Sanderson
!
        IMPLICIT NONE

        REAL,DIMENSION(:,:),            INTENT(IN)  :: pos
        INTEGER,DIMENSION(:),           INTENT(OUT) :: ka

!--     local variables

        INTEGER                   :: i,n1
        INTEGER                   :: ind_table,ind_theta

        DO i=1,size(ka)
          ind_table = pos(3,i)*scale_factor
          ind_theta = ind_array(ind_table)
          IF (pos(3,i) > eta_theta(ind_theta+1) .AND.                   &
     &      eta_theta(ind_theta) /= eta_theta(ind_theta+1)) THEN
            ka(i) = ind_theta+1
          ELSE
            ka(i) = ind_theta
          END IF
        END DO

        RETURN
      END SUBROUTINE HEIGHT_ETA_THETA

      SUBROUTINE HEIGHT_ETA_THETA_TODO(pos,todo,ka)
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  14/01/04   Created. K. Ketelsen
!  6.1  07/09/04   Added '&' for F90 compatability as needed.
!                  M.G. Sanderson
!
        IMPLICIT NONE

        REAL,DIMENSION(:,:),            INTENT(IN)  :: pos
        LOGICAL,DIMENSION(:),           INTENT(IN)  :: todo
        INTEGER,DIMENSION(:),           INTENT(OUT) :: ka

!--     local variables

        INTEGER                   :: i,n1
        INTEGER                   :: ind_table,ind_theta

        DO i=1,size(ka)
          IF (todo(i)) THEN
            ind_table = pos(3,i)*scale_factor
            ind_theta = ind_array(ind_table)
            IF (pos(3,i) > eta_theta(ind_theta+1) .AND.                 &
     &        eta_theta(ind_theta) /= eta_theta(ind_theta+1)) THEN
              ka(i) = ind_theta+1
            ELSE
              ka(i) = ind_theta
            END IF
          END IF
        END DO

        RETURN
      END SUBROUTINE HEIGHT_ETA_THETA_TODO

      SUBROUTINE HEIGHT_ETA_RHO(pos,ka)
        IMPLICIT NONE

        REAL,DIMENSION(:,:),            INTENT(IN)  :: pos
        INTEGER,DIMENSION(:),           INTENT(OUT) :: ka

        INTEGER                   :: i,n1

        DO i=1,size(ka)
          n1 = 1
          IF (pos(3,i) > eta_rho(nmetlev)) THEN
            ka(i)=nmetlev
          ELSE
            DO
              IF (pos(3,i)<=eta_rho(n1)) EXIT
              n1=n1+1
            END DO
            ka(i)=n1-1
          END IF
        END DO

        RETURN
      END SUBROUTINE HEIGHT_ETA_RHO


      SUBROUTINE HEIGHT_ETA_STOCHEM(pos,ka)
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  14/01/04   Created. K. Ketelsen
!  6.1  07/09/04   Added '&' for F90 compatability as needed.
!                  M.G. Sanderson
!
        IMPLICIT NONE

        REAL,DIMENSION(:,:),            INTENT(IN)  :: pos
        INTEGER,DIMENSION(:),           INTENT(OUT) :: ka

        INTEGER                   :: i,n1

        DO i=1,size(ka)
          n1 = 1
          DO
            IF (pos(3,i) <= eta_stochem(n1)) EXIT
            n1=n1+1
          END DO
          ka(i)=n1
        END DO

        RETURN

      END SUBROUTINE HEIGHT_ETA_STOCHEM

      END MODULE HEIGHT_MOD
#endif
