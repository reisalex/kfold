! =============================================================================
! Program: KFOLD
! 
! Description: A program for computing the folding kinetics of an RNA
!              sequence using Turner energies.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependancies:
!
! Modules - Class_RNAFold, RNAVar
! Functions -
! Subroutines - READDATA SSAREACTION V2CT
!
! Author(s): Eric Dykeman, Alex Reis
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! =============================================================================

      SUBROUTINE KFOLD(FASTA,STR0,STRF,EF,NSIM,TMAX,TRAJOUT,EOUT,FPT,EFPT)

        USE RNAVar, ONLY : mxnt

        USE Class_RNAFold

        IMPLICIT NONE

        !=== VARIABLES FOR F2PY INTERFACE ===!

        CHARACTER (LEN=mxnt+1), INTENT(IN) :: fasta
        CHARACTER (LEN=mxnt),   INTENT(IN) :: str0
        CHARACTER (LEN=mxnt),   INTENT(IN) :: strf
        REAL,                   INTENT(IN) :: ef
        INTEGER,                INTENT(IN) :: nsim
        DOUBLE PRECISION,       INTENT(IN) :: tmax

        CHARACTER,              INTENT(OUT) :: trajout(nsim,100,mxnt+1)
        REAL,                   INTENT(OUT) :: eout(nsim,100)
        DOUBLE PRECISION,       INTENT(OUT) :: fpt(nsim)
        DOUBLE PRECISION,       INTENT(OUT) :: efpt(nsim)

        !=== VARIABLES ===!

        TYPE(RNA_STRUC) :: rna

        CHARACTER :: seq(mxnt)
        CHARACTER :: fld(mxnt)

        INTEGER :: iseq(mxnt)
        INTEGER :: ibpi(mxnt)
        INTEGER :: ibpf(mxnt)

        INTEGER :: i,j,k,n,nn,is,io,oindex
        INTEGER :: isim,narg,iseed,seedgen

        DOUBLE PRECISION :: random,tstart,time
        DOUBLE PRECISION :: tout,dt

        LOGICAL :: fstop,estop

        REAL :: e

        !=== DEFAULT SETTINGS ===!

        !=== Max Time (microseconds) ===!

        tstart = 0.0d0

        efpt(:) = 1.0d12
        fpt(:)  = 1.0d12

        iseed = SEEDGEN()

        !=== Inital / Final Structure ===!

        ibpi(:) = 0
        ibpf(:) = 0

        nn = LEN_TRIM(fasta)

        IF ( nn > mxnt ) THEN
          WRITE(*,*)'ERROR: Maximum number of nt = ',mxnt
          WRITE(*,*)'Increase mxnt in rnavar.f90'
          STOP
        ENDIF

        !=== Convert fasta, str0, strf ===!

        READ(fasta,'(10000A1)')(seq(k),k=1,nn)

        READ(str0,'(10000A1)')(fld(k),k=1,nn)
        CALL V2CT (ibpi,fld,'C',nn)

        READ(strf,'(10000A1)')(fld(k),k=1,nn)
        CALL V2CT (ibpf,fld,'C',nn)

        IF ( ALL( (/1.0d0, 10.0d0, 100.0d0, 1000.0d0, 10000.0d0, 100000.0d0/)/=tmax ) )THEN
          WRITE(*,*)'ERROR: tmax must be an nth power of 10, i.e. 1, 10, 100'
          CALL EXIT(0)
        ENDIF

        !=== SECTION 1 - Read in data ===!

        CALL READDATA

        !=== Setup RNA ===!

        CALL CONVERT (seq,iseq,nn)
        CALL SETUPNUC (nn)

        rna% seq(:) = seq(:)
        rna% iseq(:) = iseq(:)
        rna% n = nn

        !=== SECTION 2 - Perform RNA Kinetics ===! 

        DO isim=1,nsim

          io = 1
          oindex = 1
          dt = 1.0d-2

          tout = dt
          time = tstart

          fstop = .TRUE.
          estop = .TRUE.

          ! Initial structure always provided
          rna% ibsp(:) = ibpi(:)
 
          CALL LOOP_INIT (rna)

          !=== STOCHASTIC SIMULATION ===!

          DO WHILE ( time < tmax )

            CALL SSAREACTION (rna,iseed,time,tout,fld,e)

            !=== Increment tout ===!

            IF ( time > tout ) THEN

              trajout(isim,oindex,1:nn) = fld(1:nn)
              eout(isim,oindex) = e

              oindex = oindex + 1

              tout = tout + dt

              io = io + 1

              IF ( io > 9 ) THEN
                io = 1
                dt = dt * 10.0d0
              ENDIF

            ENDIF

            !=== Check for stop structure ===!

            IF ( fstop ) THEN

              j = 0

              DO i=1,nn
              IF ( ibpf(i) == rna%ibsp(i) ) j = j + 1
              ENDDO

              IF ( j == nn ) THEN
                fpt(isim) = time
                fstop = .FALSE.
              ENDIF

            ENDIF

            !=== Check for stop energy threshold ===!

            IF ( estop ) THEN

              IF ( e <= ef ) THEN
                efpt(isim) = time
                estop = .FALSE.
              ENDIF

            ENDIF

          ENDDO

        ENDDO

        RETURN

      END SUBROUTINE KFOLD
