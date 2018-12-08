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

    SUBROUTINE RUN(FASTA, STR0, STRF, EF, NSIM, TMAX, TRAJOUT, EOUT, FPT, EFPT)

        CHARACTER (LEN=10001), INTENT(IN) :: fasta
        CHARACTER (LEN=10001), INTENT(IN) :: str0
        CHARACTER (LEN=10001), INTENT(IN) :: strf
        REAL,                  INTENT(IN) :: ef
        DOUBLE PRECISION,      INTENT(IN) :: tmax
        INTEGER,               INTENT(IN) :: nsim

        CHARACTER,        INTENT(OUT) :: trajout(nsim,100,10001)
        REAL,             INTENT(OUT) :: eout(nsim,100)
        DOUBLE PRECISION, INTENT(OUT) :: fpt(nsim)
        DOUBLE PRECISION, INTENT(OUT) :: efpt(nsim)

        CALL KFOLD(fasta,str0,strf,ef,nsim,tmax,trajout,eout,fpt,efpt)

        RETURN

    END SUBROUTINE RUN
