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
      Subroutine G05CFE(idum_out,iv_out,iy_out)
      Implicit none

      Integer ntab
      Parameter (ntab=32)
!     Arguments to extract from the memory of the random routines
      Integer idum_out,iv_out(ntab),iy_out
!     Local variables
      Integer idum,iv(ntab),iy
      Common /ran_jc/ idum,iv,iy ! This common block constitutes the
                                ! 'memory' of the timeserie.
      Integer j
!

!     Copy the values straight from the common into the _out
!     variables.
      Do j = 1, ntab
        iv_out(j) = iv(j)
      Enddo
      idum_out = idum
      iy_out = iy

      Return
      END SUBROUTINE G05CFE


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
