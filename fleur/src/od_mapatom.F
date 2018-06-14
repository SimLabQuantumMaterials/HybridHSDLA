      MODULE m_od_mapatom
      CONTAINS 
      SUBROUTINE od_mapatom(
     >        odd,natd,ntype,taual,neq,mrot,tau,amat,bmat,
     <        ngopr,invtab,multab,invsat,invsatnr)

c      written by Y.Mokrousov in order to generate the arrays connected
c      to the operations, transforming atoms into each other,
c      for more details look in mapatom.F.    year 2004 

      USE m_od_types, ONLY : od_dim
      USE m_dotir, ONLY : dotirl
      USE m_dotset

      IMPLICIT NONE

      TYPE (od_dim), INTENT (IN) :: odd
      INTEGER, INTENT (IN) :: natd,ntype
      REAL,    INTENT (IN) :: taual(3,natd)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3)
      INTEGER, INTENT (IN) :: neq(ntype)
      REAL,    INTENT (IN) :: mrot(3,3,odd%nop)
      REAL,    INTENT (IN) :: tau(3,odd%nop)

      INTEGER, INTENT (OUT) :: ngopr(natd)
      INTEGER, INTENT (OUT) :: invtab(odd%nop)
      INTEGER, INTENT (OUT) :: multab(odd%nop,odd%nop)
      INTEGER, INTENT (OUT) :: invsat(natd) 
      INTEGER, INTENT (OUT) :: invsatnr(natd)

      REAL ij,pps(3),norm,aamat(3,3),bbmat(3,3)
      INTEGER i,j,n1,k,n,n2,np1,na,ix,iy,iz,nat1,nat2,na2
      REAL mt(3,3),sum_tau_lat(3),sum_taual(3)
      REAL, PARAMETER :: del = 1.0e-4

      CALL dotset(
     >            amat,bmat,
     <            aamat,bbmat)

      n1 = 1
      DO 45 n = 1,ntype
        n2 = n1 + neq(n) - 1
        IF (neq(n).EQ.1) THEN
           ngopr(n2) = 1
           WRITE (6,FMT=8010) n2,n2,ngopr(n2)
           n1 = n1 + neq(n)
           GOTO 45
        END IF
        DO 55 na = n1,n2
           DO np1 = 1,odd%nop
                  pps(1) = mrot(1,1,np1)*taual(1,n1) +
     +                     mrot(1,2,np1)*taual(2,n1) +
     +                     mrot(1,3,np1)*taual(3,n1)
                  pps(2) = mrot(2,1,np1)*taual(1,n1) +
     +                     mrot(2,2,np1)*taual(2,n1) +
     +                     mrot(2,3,np1)*taual(3,n1)
                  pps(3) = mrot(3,1,np1)*taual(1,n1) +
     +                     mrot(3,2,np1)*taual(2,n1) +
     +                     mrot(3,3,np1)*taual(3,n1) +
     +                               tau(3,np1)/amat(3,3)
                  IF (abs(taual(1,na)-pps(1)).LE.1.e-4 .AND.
     &                abs(taual(2,na)-pps(2)).LE.1.e-4 .AND.
     &                abs(taual(3,na)-pps(3)).LE.1.e-4) THEN
                      ngopr(na) = np1
                      WRITE (6,FMT=8010) na,n1,ngopr(na)
 8010                 FORMAT (5x,'atom',i3,' can be mapped into atom',
     +                     i3,' through group  operation',i4)
                      GOTO 55
                   END IF
           END DO
 55     CONTINUE
        n1 = n1 + neq(n)
 45   CONTINUE

c---> defining inverse operations for the Hamiltonian and forces
c     where we do not need to consider the translational part

      DO 300 n1 = 1,odd%nop
         DO 290 n2 = 1,odd%nop
            DO 250 i = 1,3
               DO 240 j = 1,3
                  ij = 0.
                  DO 230 k = 1,3
                     ij = ij + mrot(i,k,n1)*mrot(k,j,n2)
  230             CONTINUE
                  mt(i,j) = ij
  240          CONTINUE
  250       CONTINUE
            DO 280 n = 1,odd%nop
               DO 270 i = 1,3
                  DO 260 j = 1,3
                     IF (abs(mt(i,j) - mrot(i,j,n)).LE.1.e-06)
     &                                         GO TO 260
                     GO TO 280
  260             CONTINUE
  270          CONTINUE
               multab(n1,n2) = n
               IF (n.EQ.1) invtab(n1) = n2
               GO TO 290
  280       CONTINUE
            WRITE (6,FMT=8050) n1,n2
 8050       FORMAT (' error - n1,n2=',2i3)
            STOP 'mult'
  290    CONTINUE
  300 CONTINUE

 8060 FORMAT (1x,24i3)

      WRITE (6,FMT='(//," inverse operations",//)')

      DO n1 = 1,odd%nop
         WRITE (6,FMT=8060) n1,invtab(n1)
      END DO

      DO na = 1,natd
         invsat(na) = 0
         invsatnr(na) = 0
      END DO
#if ( defined(CPP_INVERSION) )

      IF (odd%invs) THEN
         WRITE (6,FMT=*)
         nat1 = 1
         DO n = 1,ntype
            nat2 = nat1 + neq(n) - 1
            DO na = nat1,nat2 - 1
               IF (invsat(na).EQ.0) THEN
                  DO na2 = na + 1,nat2
                     DO i = 1,3
                        sum_taual(i) = taual(i,na) + taual(i,na2)
                     END DO
                     DO ix = -2,2
                       sum_tau_lat(1) = sum_taual(1) + real(ix)
                       DO iy = -2,2
                         sum_tau_lat(2) = sum_taual(2) + real(iy)
                         DO iz = -2,2
                           sum_tau_lat(3) = sum_taual(3) + real(iz)
                           norm = sqrt(dotirl(
     +                              sum_tau_lat,sum_tau_lat,aamat))
                           IF (norm.LT.del) THEN
                              invsat(na) = 1
                              invsat(na2) = 2
                              invsatnr(na)  = na2
                              invsatnr(na2) = na
                              WRITE (6,FMT=9000) n,na,na2
                              GOTO 310
                           END IF
                        END DO
                      END DO
                    END DO
 310                CONTINUE
                  END DO
               END IF
            END DO
            nat1 = nat1 + neq(n)
         END DO
      END IF
#endif
      WRITE (6,FMT=*) invsat
 9000 FORMAT ('atom type',i3,': atom',i3,' can be mapped into atom',i3,
     +       ' via 3d inversion')

      END SUBROUTINE od_mapatom
      END MODULE m_od_mapatom
