! ==============================================================================
! SUBROUTINE: RUN
! 
! Description: A program for computing the folding kinetics of an RNA
!              sequence using Turner energies.
!
!              Please See/Cite:
!              Dykeman,E.C.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            12/05/2018   Original Code
!
! Author(s): Alex Reis
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE RUN(FASTA, STR0, STRF, NSIM, TMAX, TIMES, TRAJOUT, EOUT, STRFINAL, EFINAL, FPT)

        IMPLICIT NONE

        CHARACTER (LEN=10001), INTENT(IN) :: fasta
        CHARACTER (LEN=10000), INTENT(IN) :: str0
        CHARACTER (LEN=10000), INTENT(IN) :: strf
        DOUBLE PRECISION,      INTENT(IN) :: tmax
        INTEGER,               INTENT(IN) :: nsim
        DOUBLE PRECISION,      INTENT(IN) :: times(100)

        CHARACTER,        INTENT(OUT) :: trajout(nsim,100,10001)
        REAL,             INTENT(OUT) :: eout(nsim,100)
        DOUBLE PRECISION, INTENT(OUT) :: fpt(nsim)
        
        CHARACTER,        INTENT(OUT) :: strfinal(nsim,10001)
        REAL,             INTENT(OUT) :: efinal(nsim)

        CALL KFOLD(fasta,str0,strf,nsim,tmax,times,trajout,eout,strfinal,efinal,fpt)

        RETURN

      END SUBROUTINE RUN
