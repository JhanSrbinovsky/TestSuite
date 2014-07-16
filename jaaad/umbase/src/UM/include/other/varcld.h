! ----------------------- header file: VARCLD  -----------------------
!
! Description:
!
!   Interface block for Subroutine Gen_VarDiagCloud.
!
!
! Current Code Owner: Martin Sharpe.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.1   17/08/04   Original code. Martin C. Sharpe
!   6.4   13/12/06   Updated to match changes to subroutine for
!                    incrementing qcf. Martin C. Sharpe
!
! Code Description:
!   Language: FORTRAN 90

      INTERFACE
        SUBROUTINE Gen_VarDiagCloud (                                   &
     &                             field_size,                          &
                                                   ! in
     &                             p_theta_levels,                      &
                                                   ! in
     &                             RHc,                                 &
                                                   ! in
     &                             IncrementIce,                        &
                                                   ! in
     &                             qT,                                  &
                                                   ! inout
     &                             CMessage,                            &
                                                   ! inout
     &                             ICode,                               &
                                                   ! inout
     &                             CL,                                  &
                                                   ! out   (optional)
     &                             qCL,                                 &
                                                   ! out   (optional)
     &                             dCLdp,                               &
                                                   ! out   (optional)
     &                             dCLdqT,                              &
                                                   ! out   (optional)
     &                             dCLdT,                               &
                                                   ! out   (optional)
     &                             dqCLdp,                              &
                                                   ! out   (optional)
     &                             dqCLdqT,                             &
                                                   ! out   (optional)
     &                             dqCLdT,                              &
                                                   ! out   (optional)
     &                             CF,                                  &
                                                   ! out   (optional)
     &                             qCF,                                 &
                                                   ! out   (optional)
     &                             qCFmax,                              &
                                                   ! out   (optional)
     &                             dCFdp,                               &
                                                   ! out   (optional)
     &                             dCFdqcf,                             &
                                                   ! out   (optional)
     &                             dCFdT,                               &
                                                   ! out   (optional)
     &                             dqCFdqCL,                            &
                                                   ! out   (optional)
     &                             dqCFdT,                              &
                                                   ! out   (optional)
     &                             T,                                   &
                                                   ! in    (optional)
     &                             BGqcl,                               &
                                                   ! in    (optional)
     &                             BGqcf,                               &
                                                   ! in    (optional)
     &                             BGT,                                 &
                                                   ! in    (optional)
     &                              TL)
                                                   ! inout (optional)

           IMPLICIT NONE

           INTEGER,      INTENT(IN)              ::                     &
     &         field_size

           REAL,         INTENT(IN)              ::                     &
     &         p_theta_levels (field_size),                             &
     &         RHc            (field_size)

           LOGICAL,      INTENT(IN)              ::                     &
     &         IncrementIce

           REAL,         INTENT(IN),    OPTIONAL ::                     &
     &         T     (:),                                               &
     &         BGqcl (:),                                               &
     &         BGqcf (:),                                               &
     &         BGT   (:)

           INTEGER,      INTENT(INOUT)           ::                     &
     &         ICode

           REAL,         INTENT(INOUT)           ::                     &
     &         qT (field_size)

           REAL,         INTENT(INOUT), OPTIONAL ::                     &
     &         TL        (:)

          CHARACTER(80), INTENT(INOUT)           ::                     &
     &         CMessage

           REAL,         INTENT(OUT),   OPTIONAL ::                     &
     &         CL             (:),                                      &
     &         qCl            (:),                                      &
     &         dCLdp          (:),                                      &
     &         dCLdqT         (:),                                      &
     &         dCLdT          (:),                                      &
     &         dqCLdp         (:),                                      &
     &         dqCLdqT        (:),                                      &
     &         dqCLdT         (:),                                      &
     &         CF             (:),                                      &
     &         qCF            (:),                                      &
     &         qCFmax         (:),                                      &
     &         dCFdp          (:),                                      &
     &         dCFdqcf        (:),                                      &
     &         dCFdT          (:),                                      &
     &         dqCFdqCL       (:),                                      &
     &         dqCFdT         (:)

         END SUBROUTINE Gen_VarDiagCloud
      END INTERFACE
