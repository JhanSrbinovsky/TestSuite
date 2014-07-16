! COCTMOM start
      REAL :: ddxt(imt_iso,km_iso,nt,0:1)  ! d(tracer)/dx at j=j:j+1
      REAL :: ddyt(imt_iso,km_iso,nt,0:2)   ! d(tracer)/dy at j=j-1:j
      REAL :: ddzt(imt_iso,0:km_iso,nt,0:2) ! d(tracer)/dz at j=j:j+1
      REAL :: alphai(imt_iso,km_iso,0:2)   ! temp expansion coeff at j+1
      REAL :: betai(imt_iso,km_iso,0:2)    ! salt expansion coeff at j+1
      REAL :: adv_vntiso(imt_gmm,km_gmm,0:1) ! GM meridional velocity

      ! meridional isopyc tracer flx
      REAL :: diff_fn(imt_iso,km_iso,nt,0:1)
      REAL :: stn(imt_gmm,km_gmm,0:2)
      REAL :: sbn(imt_gmm,km_gmm,0:2)
      REAL :: ath_tn(imt_gmm,km_gmm,0:2)
      REAL :: ath_bn(imt_gmm,km_gmm,0:2)
      REAL :: Bi_ez(imt_iso,km_iso,0:1,0:1,0:2)
      REAL :: Bi_nz(imt_iso,km_iso,0:1,0:1,0:2)
      REAL :: Bi_bx(imt_iso,km_iso-1,0:1,0:1,0:2)
      REAL :: Bi_by(imt_iso,km_iso-1,0:1,0:1,0:2)
! COCTMOM end
