! TYPDUMG start
! Description:
!    Generic form of dump component arguments, providing a sub-model
!    independent interface into STASH or other generic modules. Used
!    in conjunction with argument list ARGDUMG. Sub-model specific
!    argument lists in ARGDUMGA, ARGDUMGO provide the top level
!    interface. Definitions are equivalent to TYPDUMA,TYPDUMO originals.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 5.0    11/05/99   Original code. R. Rawlins
!
      INTEGER :: LEN_INTHD
      INTEGER :: LEN_REALHD
      INTEGER :: LEN1_LEVDEPC
      INTEGER :: LEN2_LEVDEPC
      INTEGER :: len_ixsts
      INTEGER :: len_spsts

      INTEGER :: FIXHD(LEN_FIXHD)
      INTEGER :: INTHD(LEN_INTHD)
      INTEGER :: LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)
      INTEGER :: ixsts(len_ixsts)
      INTEGER :: spsts(len_spsts)
      REAL    :: REALHD(LEN_REALHD)
      REAL    :: LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC+1)
! TYPDUMG end
