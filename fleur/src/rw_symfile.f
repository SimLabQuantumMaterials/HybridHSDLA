      MODULE m_rwsymfile
!----------------------------------------------------------------------!
!     writes spacegroup operations                                     ! 
!     and                                                              |  
!     rec. lattice vectors (for external k-point generator)            !
!----------------------------------------------------------------------!
      CONTAINS
      SUBROUTINE rw_symfile(
     >                      rw,symfh,symfn,nopd,bmat,
     X                      mrot,tau,nop,nop2,symor)

      IMPLICIT NONE

!===> Arguments
      CHARACTER(len=1), INTENT (IN)    :: rw
      CHARACTER(len=7), INTENT (IN)    :: symfn
      INTEGER,          INTENT (IN)    :: nopd,symfh
      REAL,             INTENT (IN)    :: bmat(3,3)
      INTEGER,          INTENT (INOUT) :: nop,nop2
      INTEGER,          INTENT (INOUT) :: mrot(3,3,nopd)
      REAL,             INTENT (INOUT) :: tau(3,nopd)
      LOGICAL,          INTENT (INOUT) :: symor

!===> Variables
      INTEGER i,j,n,ios,no2,no3,gen1
      REAL    t,d
      LOGICAL ex,op
      CHARACTER(len=3) :: type
      CHARACTER(len=7) :: sym2fn

      sym2fn = 'sym.out'

      IF ( SCAN(rw,'wW') > 0 ) THEN

!===> write symfile

        OPEN (symfh, file=sym2fn, status='unknown', err=911, iostat=ios)
        WRITE (symfh,*) nop,nop2,symor,'    ! nop,nop2,symor '
        DO n = 1, nop
           WRITE (symfh,'(a1,i3)') '!', n
           WRITE (symfh,'(3i5,5x,f10.5)') 
     &           ((mrot(i,j,n),j=1,3),tau(i,n),i=1,3)
        ENDDO
!        WRITE (symfh,*) '! reciprocal lattice vectors'
!        WRITE (symfh,'(3f25.15)') ((bmat(i,j),j=1,3),i=1,3)

      ELSEIF ( SCAN(rw,'rR') > 0 ) THEN

!===> read symfile
        OPEN (symfh, file=trim(symfn),status='old',err=911,iostat=ios)
        READ (symfh,*) nop,nop2,symor
        IF (symfn.EQ.'sym.out') THEN
          gen1 = 0
        ELSEIF (trim(symfn).EQ.'sym') THEN
          gen1 = 1
        ELSE
          STOP 'rw_symfile: symfn should be sym or sym.out'
        ENDIF
        DO n = 1 + gen1, nop + gen1
          READ (symfh,*)
          READ (symfh,*) 
     &         ((mrot(i,j,n),j=1,3),tau(i,n),i=1,3)
        ENDDO
        IF (symor) THEN
          DO n=1,nop
            t= tau(1,n)**2 + tau(2,n)**2 + tau(3,n)**2
            IF (t .GT. 1.e-8) STOP 'rw_symfile: not symmorphic'
          ENDDO
        ELSE
          DO n=1,nop
            DO i = 1,3
             IF (ABS(tau(i,n)-0.33333) < 0.00001) THEN
               tau(i,n) = 1./3.
             ENDIF
             IF (ABS(tau(i,n)+0.33333) < 0.00001) THEN
               tau(i,n) = -1./3.
             ENDIF
             IF (ABS(tau(i,n)-0.66667) < 0.00001) THEN
               tau(i,n) = 2./3.
             ENDIF
             IF (ABS(tau(i,n)+0.66667) < 0.00001) THEN
               tau(i,n) = -2./3.
             ENDIF
             IF (ABS(tau(i,n)) > 0.00001) THEN
             IF (ABS(ABS(tau(i,n))-0.5) > 0.00001) THEN
               STOP 'rw_symfile: complex phases not implemented!'
             ENDIF
             ENDIF
            ENDDO
          ENDDO
        ENDIF

        DO n = 1,nop
!
! Determine the kind of symmetry operation we have here
!
          d = mrot(1,1,n)*mrot(2,2,n)*mrot(3,3,n) +
     +        mrot(1,2,n)*mrot(2,3,n)*mrot(3,1,n) +
     +        mrot(2,1,n)*mrot(3,2,n)*mrot(1,3,n) -
     +        mrot(1,3,n)*mrot(2,2,n)*mrot(3,1,n) -
     +        mrot(2,3,n)*mrot(3,2,n)*mrot(1,1,n) -
     +        mrot(2,1,n)*mrot(1,2,n)*mrot(3,3,n)
          t =  mrot(1,1,n) + mrot(2,2,n) + mrot(3,3,n)

          IF (d.EQ.-1) THEN
            type = 'm  '
            IF (t.EQ.-3) type = 'I  '
          ELSEIF (d.EQ.1) THEN
            IF (t.EQ.-1) type = 'c_2'
            IF (t.EQ. 0) type = 'c_3'
            IF (t.EQ. 1) type = 'c_4'
            IF (t.EQ. 2) type = 'c_6'
            IF (t.EQ. 3) type = 'E  '
          ELSE
            STOP 'spg2gen: determinant =/= +/- 1'
          ENDIF
 
          WRITE (6,FMT=8020) n, type
 8020     FORMAT (/,1x,i3,' : ',a3)
          DO i = 1,3
             WRITE (6,FMT=8030) (mrot(i,j,n),j=1,3),tau(i,n)
          ENDDO
 8030     FORMAT (5x,3i3,3x,f4.1)
        ENDDO

      ELSE
        STOP 'ERROR! rw_symfile #1'
      ENDIF
      CLOSE (symfh)
      RETURN

! === errors
 911  CONTINUE
      WRITE(*,*) 'Error in inquire. IOS=',ios
 912  CONTINUE
      WRITE(*,*) 'Error in open. IOS=',ios
      STOP 'writesymfile: i/o ERROR'

      END SUBROUTINE rw_symfile
      END MODULE m_rwsymfile
