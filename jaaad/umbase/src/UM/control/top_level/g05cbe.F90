#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  Purpose of this deck : provide a number of routines equivalent
!LL  to the random number generation routines of the NAG library.
!LL  Where possible, the interface of the routines/function
!LL  matches its equivalent in the NAG library.
!LL
!LL  Written for T3E & HP platforms indifferently (then run at their
!LL  own default precision).
!LL
!LL  Model            Modification history :
!LL version  Date
!LL  4.5     04.10.98 Deck introduced as part of the Single Column
!LL                   Model
!LL
!LL  Programming standards : Unified Model documentation paper no. 4
!LL
!LL  Logical components :
!LL
!LL  System task :
!LL
!LL  Documentation : Unified Model documentation paper no ??
!LL
!
!------------------------------------------------------------------
!     Function inspired from the 'Numerical Recipies in fortran 77'
!     well known book. Generates a random number uniformly distributed
!     over [0,1].
!     It should not be called direct externally from outside this
!     deck. Instead it is accessed through the G05xx
!     routines/functions when they need a random number.
!
!     Arguments : nil
!     Returned value : a random number uniformly distributed on [0,1].
!     Side effects : the status of the /ran_jc/ common block
!                    is updated at each call.

!----------------------------------------------------------------
!
!     Initialize the random number generator RAN1_JC (and therefore
!     also the functions like g05xxx which use RAN_JC) to give
!     repeatable sequence
!
!     Argument : iseed input, not destroyed after the call
!     Side effect : initialize the /ran_jc/ common block
!
!     This routine *needs* to be called first when using the
!     random routines G05xxE
      subroutine G05CBE(iseed)
      Implicit none
!     Argument
      Integer iseed
!     Local variables
      Integer ntab
      Parameter (ntab=32)
      Integer idum, iv(ntab), iy ! This common block constitutes the
      Common /ran_jc/ idum,iv,iy ! 'memory' of the random timeserie.

      Integer j                 ! index


!     ----------------------
!     Zero the arrays of the ran1_jc function
      Do j = 1, ntab
        iv(j) = 0
      Enddo
      iy = 0
!     initialize the ran1_jc function with the user provided value:
      idum = - int(iseed)       ! idum must be negative when
                                ! initializing RAN1_JC.

      Return
      END SUBROUTINE G05CBE

!----------------------------------------------------------------
!
!     Save the state of the random number generator
!
!     Arguments : they are all output value. It is the
!     task of the user to store those values somewhere
!     until s/he wants to reinitialise to an identical
!     serie of random numbers with G05CGE, called
!     with the same arguments.
!     Side effect : nil


!----------------------------------------------------------------
!
!     Restore the state of the random number generator
!
!     Arguments : they are all input values.
!     They are the values which should be extracted
!     from the memory of the random routines by G05CFE
!     The arguments are not destroyed after the call.
!     Side effect : the /ran_jc/ common block is
!     updated with the values provided as argument.

!----------------------------------------------------------------
!
!     Pseudo-random real numbers, Gaussian distribution of
!     mean m, standard deviation s.
!
!     Also inspired from the 'Numerical Recipies in fortran 77' book.
!
!     Arguments : Input (m,s) ; are not destroyed after the call.
!     calls RAN_JC() in order to access a uniform random number
!     generator over [0,1].
!     Side effect : changes the values of the /ran_jc/ common
!     block
#endif
