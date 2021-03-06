! ==============================================================================
! Subroutine: TINT11 (LIST,EB)
! 
! Purpose: Performs a table lookup for the internal energy for a mismatch
!          pair between two basepairs in a helix.
!
! Method: Uses the MFOLD 3.0 energy function table for RNA @ T=37.
!
! Arguments:
!
!          LIST - Array of length 4 containing the nucleotides in
!                 numerical code (A=1,C=2,G=3,U=4) for the
!                 following locations:
!
!                         (3)         
!                 5' (1) A . X (5) 3' 
!                 3' (2) U . Y (6) 5' 
!                         (4)         
!
!                 where LIST(1) = letter code for position 1 etc.
!
!            EB - (OUTPUT) MFOLD 3.0 internal loop energy of the sequence
!                 provided in LIST.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependencies:
!
! Modules -
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE TINT11 (LIST,EB)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: list(6)
        REAL, INTENT(INOUT) :: eb

        !=== VARIABLES ===!

        INTEGER :: i,j,i1,i2,i3,i4,i5,i6
        REAL :: au(4,24),cg(4,24),gc(4,24)
        REAL :: ua(4,24),gu(4,24),ug(4,24)

        DATA (au(1,i),i=1,24) / 1.70, 1.70, 1.70, 1.70, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /
        DATA (au(2,i),i=1,24) / 1.70, 1.70, 1.70, 1.70, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, & 
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /
        DATA (au(3,i),i=1,24) / 1.70, 1.70,-0.40, 1.70, 1.10, 1.10, &
                              &-1.00, 1.10, 1.10, 1.10,-1.00, 1.10, &
                              & 1.70, 1.70,-0.40, 1.70, 1.70, 1.70, &
                              &-0.40, 1.70, 1.70, 1.70,-0.40, 1.70 /
        DATA (au(4,i),i=1,24) / 1.70, 1.70, 1.70, 1.50, 1.10, 1.10, &
                              & 1.10, 1.00, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.20, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /

        DATA (cg(1,i),i=1,24) / 1.10, 1.10, 1.10, 1.10, 0.40,-0.40, &
                              & 0.40, 0.40, 1.10, 0.40, 0.40, 0.40, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 /
        DATA (cg(2,i),i=1,24) / 1.10, 1.10, 1.10, 1.10, 0.30, 0.50, &
                              & 0.40, 0.50, 0.40, 0.40, 0.40, 0.40, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 /
        DATA (cg(3,i),i=1,24) / 1.10, 1.10,-1.00, 1.10,-0.10, 0.40, &
                              &-1.70, 0.40, 0.40, 0.40,-1.40, 0.40, &
                              & 1.10, 1.10,-1.00, 1.10, 1.10, 1.10, &
                              &-1.00, 1.10, 1.10, 1.10,-1.00, 1.10 /
        DATA (cg(4,i),i=1,24) / 1.10, 1.10, 1.10, 1.10, 0.40, 0.00, &
                              & 0.40,-0.30, 0.40, 0.40, 0.40, 0.40, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 /

        DATA (gc(1,i),i=1,24) / 1.10, 1.10, 1.10, 1.10, 0.80, 0.40, &
                              & 0.40, 0.40, 0.40, 0.30,-0.10, 0.40, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 /
        DATA (gc(2,i),i=1,24) / 1.10, 1.10, 1.10, 1.10, 0.40, 0.40, &
                              & 0.40, 0.40,-0.40, 0.50, 0.40, 0.00, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 /
        DATA (gc(3,i),i=1,24) / 1.10, 1.10,-1.00, 1.10, 0.40, 0.40, &
                              &-2.10, 0.40, 0.40, 0.40,-1.70, 0.40, &
                              & 1.10, 1.10,-1.00, 1.10, 1.10, 1.10, &
                              &-1.00, 1.10, 1.10, 1.10,-1.00, 1.10 /
        DATA (gc(4,i),i=1,24) / 1.10, 1.10, 1.10, 1.10, 0.40, 0.40, &
                              & 0.40,-0.70, 0.40, 0.50, 0.40,-0.30, &
                              & 1.10, 1.10, 1.10, 1.00, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 /

        DATA (ua(1,i),i=1,24) / 1.70, 1.70, 1.70, 1.70, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /
        DATA (ua(2,i),i=1,24) / 1.70, 1.70, 1.70, 1.70, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /
        DATA (ua(3,i),i=1,24) / 1.70, 1.70,-0.40, 1.70, 1.10, 1.10, &
                              &-1.00, 1.10, 1.10, 1.10,-1.00, 1.10, &
                              & 1.70, 1.70,-0.40, 1.70, 1.70, 1.70, &
                              &-0.40, 1.70, 1.70, 1.70,-0.40, 1.70 /
        DATA (ua(4,i),i=1,24) / 1.70, 1.70, 1.70, 1.80, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.50, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /

        DATA (gu(1,i),i=1,24) / 1.70, 1.70, 1.70, 1.70, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /
        DATA (gu(2,i),i=1,24) / 1.70, 1.70, 1.70, 1.70, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /
        DATA (gu(3,i),i=1,24) / 1.70, 1.70,-0.40, 1.70, 1.10, 1.10, &
                              &-1.00, 1.10, 1.10, 1.10,-1.00, 1.10, &
                              & 1.70, 1.70,-0.40, 1.70, 1.70, 1.70, &
                              &-0.40, 1.70, 1.70, 1.70,-0.40, 1.70 /
        DATA (gu(4,i),i=1,24) / 1.70, 1.70, 1.70, 1.70, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /

        DATA (ug(1,i),i=1,24) / 1.70, 1.70, 1.70, 1.70, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /
        DATA (ug(2,i),i=1,24) / 1.70, 1.70, 1.70, 1.70, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /
        DATA (ug(3,i),i=1,24) / 1.70, 1.70,-0.40, 1.70, 1.10, 1.10, &
                              &-1.00, 1.10, 1.10, 1.10,-1.00, 1.10, &
                              & 1.70, 1.70,-0.40, 1.70, 1.70, 1.70, &
                              &-0.40, 1.70, 1.70, 1.70,-0.40, 1.70 /
        DATA (ug(4,i),i=1,24) / 1.70, 1.70, 1.70, 1.70, 1.10, 1.10, &
                              & 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, &
                              & 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 /


        !         (3)         !
        ! 5' (1) A . X (5) 3' !
        ! 3' (2) U . Y (6) 5' !
        !         (4)         !

        i1 = list(1)
        i2 = list(2)
        i3 = list(3)
        i4 = list(4)
        i5 = list(5)
        i6 = list(6)

        IF ( i5 == 1 .and. i6 == 4 ) j = 0 + i4
        IF ( i5 == 2 .and. i6 == 3 ) j = 4 + i4
        IF ( i5 == 3 .and. i6 == 2 ) j = 8 + i4
        IF ( i5 == 4 .and. i6 == 1 ) j =12 + i4
        IF ( i5 == 3 .and. i6 == 4 ) j =16 + i4
        IF ( i5 == 4 .and. i6 == 3 ) j =20 + i4

        IF ( i1 == 1 .and. i2 == 4 ) eb = eb + au(i3,j)
        IF ( i1 == 2 .and. i2 == 3 ) eb = eb + cg(i3,j)
        IF ( i1 == 3 .and. i2 == 2 ) eb = eb + gc(i3,j)
        IF ( i1 == 4 .and. i2 == 1 ) eb = eb + ua(i3,j)
        IF ( i1 == 3 .and. i2 == 4 ) eb = eb + gu(i3,j)
        IF ( i1 == 4 .and. i2 == 3 ) eb = eb + ug(i3,j)

        RETURN

      END SUBROUTINE TINT11
