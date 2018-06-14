      MODULE m_jcoff
c-------------------------------------------------------------------
c     Determines the cone angles and phase shifts of the spin 
c     vectors on magnetic atoms for the calculation of the
c     interaction constants Jij from the Heisenberg model
c                                   M. Lezaic '04
c-------------------------------------------------------------------
      CONTAINS
      SUBROUTINE jcoff(
     >                 i_J,j_J,phn,irank,ntypd,natd,ntype,taual,neq,
     >                 qss,nmagn,l_magn,l_disp,thetaJ,M,
     X                 alph1,alph,beta,l_wr)

      USE m_constants, ONLY : pimach
      IMPLICIT NONE

c     .. Scalar arguments ..

      INTEGER, INTENT (IN)   :: ntype,nmagn
      INTEGER, INTENT (IN)   :: ntypd,natd
      INTEGER, INTENT (IN)   :: i_J,j_J,phn,irank
      REAL,    INTENT (IN)   :: thetaJ
      LOGICAL, INTENT (IN)   :: l_disp

c     .. Array arguments ..

      REAL,    INTENT (IN)   :: qss(3),taual(3,natd),M(ntypd)
      REAL,    INTENT (INOUT):: alph(ntypd),beta(ntypd)
      REAL,    INTENT (INOUT):: alph1(ntypd)
      INTEGER, INTENT (IN)   :: neq(ntypd)
      LOGICAL, INTENT (IN)   :: l_magn(ntypd)
      LOGICAL, INTENT (INOUT):: l_wr

c     .. Local scalars ..

      INTEGER n,itype,iatom
      REAL    tpi,pihalf,phi
c-------------------------------------------------------------------
      tpi = 2.0 * pimach()
      pihalf=0.5 * pimach()

      IF (.not.l_disp) THEN
        n = 0
        alph(:) = 0
        beta(:) = 0
c-mpi-
        IF (irank.eq.0) THEN
          IF (phn.eq.1) THEN
            WRITE(6,*) 'Calculating Jij(q) for the atom types:'
          ENDIF
        ENDIF
c-mpi-
        DO itype=1,ntype
          IF (l_magn(itype)) THEN
            n=n+1

            IF (n.eq.i_J) THEN
              beta(itype) = thetaJ
              IF (phn.eq.2) THEN
                alph(itype) = pihalf
              ENDIF
c-mpi-
              IF (irank.eq.0) THEN
                IF (phn.eq.1) THEN
                  WRITE(6,*)'atom type ',itype
                ENDIF
                WRITE(114,*) ' atom type ','     phase     ',
     &              '   cone angle ',' magnetic moment  '
                WRITE(114,5000)itype,alph(itype),beta(itype),M(itype)
              ENDIF
c-mpi-
            ENDIF ! n = i_J

            IF (n.eq.j_J) THEN
              beta(itype) = thetaJ 
c-mpi-
              IF (irank.eq.0) THEN
                IF (phn.eq.1) THEN
                  WRITE(6,*)'atom type ',itype
                ENDIF
                WRITE(114,5000)itype,alph(itype),beta(itype),M(itype)
              ENDIF
c-mpi-
            ENDIF  ! n = j_J

          ENDIF    ! l_magn
        ENDDO      ! itype
        
c-mpi-
        IF (irank.eq.0) THEN
          WRITE(6,*) 'The spin-spiral setup is:'
          WRITE(6,*) ' atom type ','     phase     ','   cone angle ',
     &               ' magnetic moment  '
          DO itype=1,ntype
            WRITE(6,5000)itype,alph(itype),beta(itype),M(itype)
          ENDDO
 5000     FORMAT(3x,i4,4x,f14.10,1x,f14.10,4x,f8.5)
        ENDIF
c-mpi-
      ELSE ! l_disp = T

        alph(:)=alph1(:)
        IF (l_wr) THEN
          IF (irank.eq.0) THEN   
            WRITE(6,*) 'The magnon spectrum will be calculated'
            WRITE(6,*) 'for the following setup:'
            DO itype=1,ntype
              WRITE(6,*)'atom type=',itype
              WRITE(6,*)'starting phase=',alph(itype)
              WRITE(6,*)'cone angle=',beta(itype)
              WRITE(6,*)
            ENDDO
          ENDIF
        ENDIF
        l_wr=.false.

      ENDIF ! l_disp

      iatom = 1
      DO itype = 1,ntype
        phi = tpi*(  qss(1)*taual(1,iatom)
     +             + qss(2)*taual(2,iatom)
     +             + qss(3)*taual(3,iatom) )
        alph(itype) = alph(itype) + phi
        iatom = iatom + neq(itype)
      ENDDO

      END SUBROUTINE jcoff
      END MODULE m_jcoff
