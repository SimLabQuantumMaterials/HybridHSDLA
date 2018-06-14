      MODULE m_brysh2
******************************************************
c     maps the density back from one single vector into the
c     proper component of interstitial, m.t. and vacuum density
******************************************************
      CONTAINS 
      SUBROUTINE brysh2(
     >                  jspins,n3d,ntypd,nlhd,ntypsd,jmtd,mmap,n2d,natd,
     >                  l_noco,nq3,ntype,nlh,ntypsy,jri,nmzxy,nmz,
     >                  film,invs,invs2,nvac,nq2,neq,nmap,s_in,n_u,
     <                  n_mmp,odi,qpw,rho,rht,rhtxy,cdom,cdomvz,cdomvxy)

      USE m_od_types, ONLY : od_inp
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins,n3d,ntypd,nlhd,ntypsd,jmtd,mmap,n_u
      INTEGER, INTENT (IN) :: nmzxy,nmz,n2d,natd,nmap,nq3,nq2,ntype,nvac
      LOGICAL, INTENT (IN) :: film,invs,invs2,l_noco
C     ..
C     .. Array Arguments ..
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c-odim
      INTEGER, INTENT (IN) :: jri(ntypd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),neq(ntypd)
      REAL,    INTENT (IN) :: s_in(mmap)
      REAL,    INTENT (OUT) :: rho(jmtd,0:nlhd,ntypd,jspins)
      REAL,    INTENT (OUT) :: rht(nmz,2,jspins)
      COMPLEX, INTENT (OUT) :: qpw(n3d,jspins),cdom(n3d),cdomvz(nmz,2)
      COMPLEX, INTENT (OUT) :: rhtxy(nmzxy,odi%n2d-1,2,jspins)
      COMPLEX, INTENT (OUT) :: cdomvxy(nmzxy,odi%n2d-1,2)
      COMPLEX, INTENT (OUT) :: n_mmp(-3:3,-3:3,n_u,jspins)
C     ..
C     .. Local Scalars ..
      INTEGER i,iv,j,js,k,l,n,na
C
      j=0
      DO 20 js = 1,jspins
         IF (invs) THEN
            DO i = 1,nq3
               j = j + 1
               qpw(i,js) = cmplx(s_in(j),0.0)
            END DO
         ELSE
            DO i = 1,nq3
               j = j + 1
               qpw(i,js) = cmplx(s_in(j),s_in(j+nq3))
            END DO
            j = j + nq3
         ENDIF
         na = 1
         DO n = 1,ntype
            DO l = 0,nlh(ntypsy(na))
               DO i = 1,jri(n)
                  j = j + 1
                  rho(i,l,n,js) = s_in(j)
               END DO
            END DO
            na = na + neq(n)
         END DO
         IF (film) THEN
            DO iv = 1,nvac
               DO k = 1,nmz
                  j = j + 1
                  rht(k,iv,js) = s_in(j)
               END DO
               DO k = 1,odi%nq2-1
                  DO i = 1,nmzxy
                     j = j + 1
                     rhtxy(i,k,iv,js) = cmplx(s_in(j),0.0)
                  END DO
               END DO
               IF (.not.invs2) THEN
                  DO k = 1,odi%nq2-1
                     DO i = 1,nmzxy
                        j = j + 1
                        rhtxy(i,k,iv,js) = rhtxy(i,k,iv,js) +
     +                                       cmplx(0.0,s_in(j))
                     END DO
                  END DO
               END IF
            END DO
         END IF
 20   CONTINUE   

      IF (l_noco) THEN
c--->    off-diagonal part of the density matrix
         DO i = 1,nq3
            j = j + 1
            cdom(i) = cmplx(s_in(j),0.0)
         END DO
         DO i = 1,nq3
            j = j + 1
            cdom(i) = cdom(i) + cmplx(0.0,s_in(j))
         END DO
         IF (film) THEN
            DO iv = 1,nvac
               DO k = 1,nmz
                  j = j + 1
                  cdomvz(k,iv) = cmplx(s_in(j),0.0)
               END DO
               DO k = 1,odi%nq2-1
                  DO i = 1,nmzxy
                     j = j + 1
                     cdomvxy(i,k,iv) = cmplx(s_in(j),0.0)
                  END DO
               END DO
            END DO
            DO iv = 1,nvac
               DO k = 1,nmz
                  j = j + 1
                  cdomvz(k,iv) = cdomvz(k,iv) + cmplx(0.0,s_in(j))
               END DO
               DO k = 1,odi%nq2-1
                  DO i = 1,nmzxy
                     j = j + 1
                     cdomvxy(i,k,iv) = cdomvxy(i,k,iv)
     +                               + cmplx(0.0,s_in(j))
                  END DO
               END DO
            END DO
         END IF
      ENDIF

      IF ( n_u > 0 ) THEN
        DO js = 1,jspins
          DO n = 1, n_u
            DO k = -3, 3
              DO i = -3, 3
                j = j + 1
                n_mmp(i,k,n,js) = cmplx(s_in(j),s_in(j+1))
                j = j + 1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF (j.NE.nmap) THEN
         WRITE (6,'(a,i5,a,i5)') 'j = ',j,' nmap = ',nmap
         STOP ' brysh2 : j =/= nmap '
      ENDIF

      END SUBROUTINE brysh2
      END MODULE m_brysh2
