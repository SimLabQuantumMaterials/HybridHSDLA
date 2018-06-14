      MODULE m_mapatom
c*******************************************************************
c     determines the group operation which maps the representive
c     atom into its equivalent atoms     c.l.fu
c*******************************************************************
      CONTAINS
      SUBROUTINE mapatom(
     >                   nop,ntypd,natd,ntype,ntypsy,
     >                   amat,bmat,mrot,neq,tau,taual,l_f,n_u,
     >                   l_soc,theta,phi,
     <                   invs,ngopr,invsat,invsatnr,bbmat,
     <                   multab,invtab,invarop,invarind)
c
c     if (l_f) setup multab,invtab,invarop,invarind for force_a12 & 21
c***********************************************************************
c the contribution to the hamiltonian and to the overlap matrix in a
c system with inversion symmetry from one muffin tin is the complex
c conjugate of the contribution from the "invers" muffin tin. this fact
c can be exploited to save cpu-time in hssphn. Thus, it is nessessary to
c know whether an atom can be mapped onto an equivalent atom via 3d
c inversion. where both atoms have to belong to the same unit cell, i.e.
c the are not related to each other by a lattice translation. therefore,
c an array invsatom is set up.
c invsatom(natom) =
c 0 if the atom cannot be mapped onto an eqivalent atom via inversion
c 1 if the atom can be mapped onto an eqivalent atom via inversion, and
c   has a smaller atom index than the related atom
c 2 if the atom can be mapped onto an eqivalent atom via inversion, and
c   has a bigger atom index than the related atom
c p.kurz aug. 1996
c***********************************************************************
c
      USE m_dotset
      USE m_dotir, ONLY : dotirl
      USE m_socsym

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN)  :: nop,ntypd,natd
      INTEGER, INTENT (IN)  :: ntype,n_u
      REAL,    INTENT (IN)  :: theta,phi
      LOGICAL, INTENT (IN)  :: l_f,invs,l_soc
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: mrot(3,3,nop),neq(ntypd),ntypsy(natd)
      INTEGER, INTENT (OUT) :: ngopr(natd)
      INTEGER, INTENT (OUT) :: invarind(natd)
      INTEGER, INTENT (OUT) :: invarop(natd,nop)
      INTEGER, INTENT (OUT) :: invtab(nop)
      INTEGER, INTENT (OUT) :: multab(nop,nop)
      INTEGER, INTENT (OUT) :: invsat(natd)
      INTEGER, INTENT (OUT) :: invsatnr(natd)
      REAL,    INTENT (IN)  :: tau(3,nop),taual(3,natd)
      REAL,    INTENT (IN)  :: amat(3,3),bmat(3,3)
      REAL,    INTENT (OUT) :: bbmat(3,3)
C     ..
C     .. Local Scalars ..
      REAL s3,norm
      INTEGER i,icount,j,j1,j2,j3,jop,n,na,nat1,nat2,nb,na_r
      INTEGER k,ij,n1,n2,ix,iy,iz,na2
      REAL, PARAMETER :: del = 1.0e-4
C     ..
C     .. Local Arrays ..
      INTEGER mt(3,3),mp(3,3)
      REAL aamat(3,3),sum_tau_lat(3),sum_taual(3)
      REAL gam(3),gaminv(3),gamr(3),sr(3),ttau(3)
      LOGICAL error(nop)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real,sqrt
C     ..
      CALL dotset(
     >            amat,bmat,
     <            aamat,bbmat)
c
      IF (l_soc) THEN  ! check once more here...
        CALL soc_sym(
     >               nop,mrot,theta,phi,amat,
     <               error)
      ELSE
        error(:) = .false.
      ENDIF
                               
      WRITE (6,FMT=8000)
 8000 FORMAT (/,/,5x,'group operations on equivalent atoms:')
      nat1 = 1
      DO 40 n = 1,ntype
         nat2 = nat1 + neq(n) - 1
         ngopr(nat1) = 1
c+gu
         na_r = nat1
         DO 30 na = nat1,nat2
            IF (ntypsy(na).NE.ntypsy(na_r)) na_r = na
c-gu
            DO i = 1,3
               gam(i) = taual(i,na)
            END DO
            invarind(na) = 0
            icount = 0
            DO 10 jop = 1,nop
               DO i = 1,3
                  gamr(i) = 0.
                  DO j = 1,3
                     gamr(i) = gamr(i) + mrot(i,j,jop)*gam(j)
                  END DO
                  gamr(i) = gamr(i) + tau(i,jop)
               END DO
               DO i = 1,3
                  gaminv(i) = gamr(i) - taual(i,na)
                  gamr(i)   = gamr(i) - taual(i,nat1) ! cf local_sym
               END DO
               IF (icount.EQ.0) THEN
                  DO j3 = -2,2
                     sr(3) = gamr(3) + real(j3)
                     DO j2 = -2,2
                        sr(2) = gamr(2) + real(j2)
                        DO j1 = -2,2
                           sr(1) = gamr(1) + real(j1)
                           s3 = sqrt(dotirl(sr,sr,aamat))
                           IF ((s3.LT.del).AND.(.not.error(jop))) THEN
                              icount = icount + 1
                              ngopr(na) = jop
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
c
c search for operations which leave taual invariant
c
               IF (l_f.OR.(n_u.GT.0)) THEN 
                  DO j3 = -2,2
                     sr(3) = gaminv(3) + real(j3)
                     DO j2 = -2,2
                        sr(2) = gaminv(2) + real(j2)
                        DO j1 = -2,2
                           sr(1) = gaminv(1) + real(j1)
                           s3 = sqrt(dotirl(sr,sr,aamat))
                           IF (s3.LT.del) THEN
                              invarind(na) = invarind(na) + 1
                              invarop(na,invarind(na)) = jop
                           END IF
                        END DO
                     END DO
                  END DO
               ENDIF
c
c end of operations
   10       CONTINUE
            IF (icount.GT.0) GO TO 20
c
            STOP 'mapatom'
   20       CONTINUE
            WRITE (6,FMT=8010) nat1,na,ngopr(na)
 8010       FORMAT (5x,'atom',i3,' can be mapped into atom',i3,
     +             ' through group  operation',i4)
c
c end of equivalent atoms
   30    CONTINUE
c
         nat1 = nat1 + neq(n)
c
c end of different types of atoms
   40 CONTINUE

c------------------------- FORCE PART -------------------------------
c+gu this is the remainder of spgset necessary for force calculations
c
      IF (l_f.OR.(n_u.GT.0)) THEN

      WRITE (6,FMT=
     +  '(//,"list of operations which leave taual invariant",/)')
      DO na = 1,nat2
         WRITE (6,FMT='("atom nr.",i3,3x,(t14,"ops are:",24i3))') na,
     +     (invarop(na,nb),nb=1,invarind(na))
      END DO

      ENDIF
c------------------------- FORCE PART ENDS --------------------------
!
!     check closure  ; note that:  {R|t} tau = R^{-1} tau -  R^{-1} t
!
!--->    loop over all operations
c
      WRITE (6,FMT=8040)
 8040 FORMAT (/,/,' multiplication table',/,/)
      multab = 0
      DO j=1,nop

!--->    multiply {R_j|t_j}{R_i|t_i}
         DO i=1,nop
            mp = matmul( mrot(:,:,j) , mrot(:,:,i) )
            ttau = tau(:,j) + matmul( mrot(:,:,j) , tau(:,i) )
            ttau = ttau - anint( ttau - 1.e-7 )

!--->    determine which operation this is
            DO k=1,nop
              IF( all( mp(:,:) == mrot(:,:,k) ) .AND.
     &            ALL( abs( ttau(:)-tau(:,k) ) < 1.e-7 ) ) THEN
                 IF (multab(j,i) .EQ. 0 ) THEN
                    multab(j,i) = k
                    IF (k .EQ. 1) invtab(j)=i
                 ELSE
                    WRITE(6,'(" Symmetry error: multiple ops")')
                    STOP 'Multiple ops'
                 ENDIF
              ENDIF
            ENDDO

            IF (multab(j,i).EQ.0) THEN
               WRITE (6,'(" Group not closed")')
               WRITE (6,'("  j , i =",2i4)') j,i
               STOP 'mapatom: group not closed'
            ENDIF
         ENDDO
      ENDDO

      DO n1 = 1,nop
         WRITE (6,FMT=8060) (multab(n1,n2),n2=1,nop)
      END DO
 8060 FORMAT (1x,24i3)
      WRITE (6,FMT='(//," inverse operations",//)')
      DO n1 = 1,nop
         WRITE (6,FMT=8060) n1,invtab(n1)
      END DO

      DO na = 1,natd
         invsat(na) = 0
         invsatnr(na) = 0
      END DO
!#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) ) ! soc moved to main gb`02
                                                     ! before call eigenso 
#if ( defined(CPP_INVERSION) )                                                                         
      IF (invs) THEN
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
      WRITE (6,FMT=*) invsat
 9000 FORMAT ('atom type',i3,': atom',i3,' can be mapped into atom',i3,
     +       ' via 3d inversion')
#endif

      END  SUBROUTINE mapatom
      END  MODULE m_mapatom
