!     ------------------------------------------------------------------
!     Module setting gaussian points and weights: only even
!     orders are used at present.
!
!     First order:
      DATA GAUSS_POINT(1,1)/0.00000000000000E+00_real64/
      DATA GAUSS_WEIGHT(1,1)/2.00000000000000E+00_real64/
!
!     Second order:
      DATA GAUSS_POINT(1,2)/-5.77350269189626E-01_real64/
      DATA GAUSS_POINT(2,2)/5.77350269189626E-01_real64/
      DATA GAUSS_WEIGHT(1,2)/1.00000000000000E+00_real64/
      DATA GAUSS_WEIGHT(2,2)/1.00000000000000E+00_real64/
!
!     Third order:
      DATA GAUSS_POINT(1,3)/-7.74596669241484E-01_real64/
      DATA GAUSS_POINT(2,3)/0.00000000000000E+00_real64/
      DATA GAUSS_POINT(3,3)/7.74596669241484E-01_real64/
      DATA GAUSS_WEIGHT(1,3)/5.55555555555556E-01_real64/
      DATA GAUSS_WEIGHT(2,3)/8.88888888888889E-01_real64/
      DATA GAUSS_WEIGHT(3,3)/5.55555555555556E-01_real64/
!
!     Fourth order:
      DATA GAUSS_POINT(1,4)/-8.61136311594053E-01_real64/
      DATA GAUSS_POINT(2,4)/-3.39981043584856E-01_real64/
      DATA GAUSS_POINT(3,4)/3.39981043584856E-01_real64/
      DATA GAUSS_POINT(4,4)/8.61136311594053E-01_real64/
      DATA GAUSS_WEIGHT(1,4)/3.47854845137454E-01_real64/
      DATA GAUSS_WEIGHT(2,4)/6.52145154862546E-01_real64/
      DATA GAUSS_WEIGHT(3,4)/6.52145154862546E-01_real64/
      DATA GAUSS_WEIGHT(4,4)/3.47854845137454E-01_real64/
!
!     Fifth order:
      DATA GAUSS_POINT(1,5)/-9.06179845938664E-01_real64/
      DATA GAUSS_POINT(2,5)/-5.38469310105683E-01_real64/
      DATA GAUSS_POINT(3,5)/0.00000000000000E+00_real64/
      DATA GAUSS_POINT(4,5)/5.38469310105683E-01_real64/
      DATA GAUSS_POINT(5,5)/9.06179845938664E-01_real64/
      DATA GAUSS_WEIGHT(1,5)/2.36926885056189E-01_real64/
      DATA GAUSS_WEIGHT(2,5)/4.78628670499366E-01_real64/
      DATA GAUSS_WEIGHT(3,5)/4.67913934572691E-01_real64/
      DATA GAUSS_WEIGHT(4,5)/4.78628670499366E-01_real64/
      DATA GAUSS_WEIGHT(5,5)/2.36926885056189E-01_real64/
!
!     Sixth order:
      DATA GAUSS_POINT(1,6)/-9.32469514203152E-01_real64/
      DATA GAUSS_POINT(2,6)/-6.61209386466265E-01_real64/
      DATA GAUSS_POINT(3,6)/-2.38619186083197E-01_real64/
      DATA GAUSS_POINT(4,6)/2.38619186083197E-01_real64/
      DATA GAUSS_POINT(5,6)/6.61209386466265E-01_real64/
      DATA GAUSS_POINT(6,6)/9.32469514203152E-01_real64/
      DATA GAUSS_WEIGHT(1,6)/1.71324492379170E-01_real64/
      DATA GAUSS_WEIGHT(2,6)/3.60761573048139E-01_real64/
      DATA GAUSS_WEIGHT(3,6)/4.67913934572691E-01_real64/
      DATA GAUSS_WEIGHT(4,6)/4.67913934572691E-01_real64/
      DATA GAUSS_WEIGHT(5,6)/3.60761573048139E-01_real64/
      DATA GAUSS_WEIGHT(6,6)/1.71324492379170E-01_real64/
!
!     Seventh order:
      DATA GAUSS_POINT(1,7)/-9.49107912342759E-01_real64/
      DATA GAUSS_POINT(2,7)/-7.41531185599394E-01_real64/
      DATA GAUSS_POINT(3,7)/-4.05845151377397E-01_real64/
      DATA GAUSS_POINT(4,7)/0.00000000000000E+00_real64/
      DATA GAUSS_POINT(5,7)/4.05845151377397E-01_real64/
      DATA GAUSS_POINT(6,7)/7.41531185599394E-01_real64/
      DATA GAUSS_POINT(7,7)/9.49107912342759E-01_real64/
      DATA GAUSS_WEIGHT(1,7)/1.29484966168870E-01_real64/
      DATA GAUSS_WEIGHT(2,7)/2.79705391489277E-01_real64/
      DATA GAUSS_WEIGHT(3,7)/3.81830050505119E-01_real64/
      DATA GAUSS_WEIGHT(4,7)/4.17959183673469E-01_real64/
      DATA GAUSS_WEIGHT(5,7)/3.81830050505119E-01_real64/
      DATA GAUSS_WEIGHT(6,7)/2.79705391489277E-01_real64/
      DATA GAUSS_WEIGHT(7,7)/1.29484966168870E-01_real64/
!
!     Eighth order:
      DATA GAUSS_POINT(1,8)/-9.60289856497536E-01_real64/
      DATA GAUSS_POINT(2,8)/-7.96666477413627E-01_real64/
      DATA GAUSS_POINT(3,8)/-5.25532409916329E-01_real64/
      DATA GAUSS_POINT(4,8)/-1.83434642495650E-01_real64/
      DATA GAUSS_POINT(5,8)/1.83434642495650E-01_real64/
      DATA GAUSS_POINT(6,8)/5.25532409916329E-01_real64/
      DATA GAUSS_POINT(7,8)/7.96666477413627E-01_real64/
      DATA GAUSS_POINT(8,8)/9.60289856497536E-01_real64/
      DATA GAUSS_WEIGHT(1,8)/1.01228536290376E-01_real64/
      DATA GAUSS_WEIGHT(2,8)/2.22381034453374E-01_real64/
      DATA GAUSS_WEIGHT(3,8)/3.13706645877887E-01_real64/
      DATA GAUSS_WEIGHT(4,8)/3.62683783378362E-01_real64/
      DATA GAUSS_WEIGHT(5,8)/3.62683783378362E-01_real64/
      DATA GAUSS_WEIGHT(6,8)/3.13706645877887E-01_real64/
      DATA GAUSS_WEIGHT(7,8)/2.22381034453374E-01_real64/
      DATA GAUSS_WEIGHT(8,8)/1.01228536290376E-01_real64/
!
!     Ninth order:
      DATA GAUSS_POINT(1,9)/-9.68160239507626E-01_real64/
      DATA GAUSS_POINT(2,9)/-8.36031107326636E-01_real64/
      DATA GAUSS_POINT(3,9)/-6.13371432700590E-01_real64/
      DATA GAUSS_POINT(4,9)/-3.24253423403809E-01_real64/
      DATA GAUSS_POINT(5,9)/0.00000000000000E+00_real64/
      DATA GAUSS_POINT(6,9)/3.24253423403809E-01_real64/
      DATA GAUSS_POINT(7,9)/6.13371432700590E-01_real64/
      DATA GAUSS_POINT(8,9)/8.36031107326636E-01_real64/
      DATA GAUSS_POINT(9,9)/9.68160239507626E-01_real64/
      DATA GAUSS_WEIGHT(1,9)/8.1274388361574E-02_real64/
      DATA GAUSS_WEIGHT(2,9)/1.80648160694857E-01_real64/
      DATA GAUSS_WEIGHT(3,9)/2.60610696402935E-01_real64/
      DATA GAUSS_WEIGHT(4,9)/3.12347077040003E-01_real64/
      DATA GAUSS_WEIGHT(5,9)/3.30239355001260E-01_real64/
      DATA GAUSS_WEIGHT(6,9)/3.12347077040003E-01_real64/
      DATA GAUSS_WEIGHT(7,9)/2.60610696402935E-01_real64/
      DATA GAUSS_WEIGHT(8,9)/1.80648160694857E-01_real64/
      DATA GAUSS_WEIGHT(9,9)/8.1274388361574E-02_real64/
!
!     Tenth order:
      DATA GAUSS_POINT(1,10)/-9.73906528517172E-01_real64/
      DATA GAUSS_POINT(2,10)/-8.65063366688985E-01_real64/
      DATA GAUSS_POINT(3,10)/-6.79409568299024E-01_real64/
      DATA GAUSS_POINT(4,10)/-4.33395394129427E-01_real64/
      DATA GAUSS_POINT(5,10)/-1.48874338981631E-01_real64/
      DATA GAUSS_POINT(6,10)/1.48874338981631E-01_real64/
      DATA GAUSS_POINT(7,10)/4.33395394129427E-01_real64/
      DATA GAUSS_POINT(8,10)/6.79409568299024E-01_real64/
      DATA GAUSS_POINT(9,10)/8.65063366688985E-01_real64/
      DATA GAUSS_POINT(10,10)/9.73906528517172E-01_real64/
      DATA GAUSS_WEIGHT(1,10)/6.6671344308688E-02_real64/
      DATA GAUSS_WEIGHT(2,10)/1.49451349150581E-01_real64/
      DATA GAUSS_WEIGHT(3,10)/2.19086362515982E-01_real64/
      DATA GAUSS_WEIGHT(4,10)/2.69266719309996E-01_real64/
      DATA GAUSS_WEIGHT(5,10)/2.95524224714753E-01_real64/
      DATA GAUSS_WEIGHT(6,10)/2.95524224714753E-01_real64/
      DATA GAUSS_WEIGHT(7,10)/2.69266719309996E-01_real64/
      DATA GAUSS_WEIGHT(8,10)/2.19086362515982E-01_real64/
      DATA GAUSS_WEIGHT(9,10)/1.49451349150581E-01_real64/
      DATA GAUSS_WEIGHT(10,10)/6.6671344308688E-02_real64/
!
!     ------------------------------------------------------------------
