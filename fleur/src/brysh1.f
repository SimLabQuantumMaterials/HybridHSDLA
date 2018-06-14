      MODULE m_brysh1
******************************************************
c      shifts the charge density of the interstitial, m.t
c      and vacuum part in one single vector
c      in the spin polarized case the arrays consist of 
c      spin up and spin down densities
******************************************************
      CONTAINS
      SUBROUTINE brysh1(
     >                  jspins,n3d,ntypd,nlhd,ntypsd,jmtd,mmap,n2d,natd,
     >                  l_noco,nq3,ntype,nlh,ntypsy,jri,nmzxy,nmz,
     >                  film,invs,invs2,nvac,nq2,neq,intfac,vacfac,
     >                  qpw,rho,rht,rhtxy,cdom,cdomvz,cdomvxy,n_u,
     >                  n_mmp,odi,
     <                  nmap,nmaph,mapmt,mapvac,mapvac2,sout) 

      USE m_od_types, ONLY : od_inp
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins,n3d,ntypd,nlhd,ntypsd,jmtd,mmap
      INTEGER, INTENT (IN) :: nmzxy,nmz,n2d,natd,nq3,nq2,ntype,nvac,n_u
      REAL,    INTENT (IN) :: intfac,vacfac
      LOGICAL, INTENT (IN) :: film,invs,invs2,l_noco
      INTEGER, INTENT (OUT):: mapmt,mapvac,mapvac2,nmap,nmaph
C     ..
C     .. Array Arguments ..
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim
      COMPLEX, INTENT (IN) :: qpw(n3d,jspins)
      COMPLEX, INTENT (IN) :: cdomvz(nmz,2),cdomvxy(nmzxy,odi%n2d-1,2)
      COMPLEX, INTENT (IN) :: cdom(n3d),rhtxy(nmzxy,odi%n2d-1,2,jspins)
      COMPLEX, INTENT (IN) :: n_mmp(-3:3,-3:3,n_u,jspins)
      INTEGER, INTENT (IN) :: jri(ntypd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),neq(ntypd)
      REAL,    INTENT (IN) :: rho(jmtd,0:nlhd,ntypd,jspins)
      REAL,    INTENT (IN) :: rht(nmz,2,jspins)
      REAL,    INTENT (OUT):: sout(mmap)
C     ..
C     .. Local Scalars ..
      INTEGER i,iv,j,js,k,l,n,nall,na,nvaccoeff,nvaccoeff2,mapmtd
c
c--->  put input into arrays sout 
c      in the spin polarized case the arrays consist of 
c      spin up and spin down densities

      j=0
      DO 10 js = 1,jspins
         DO i = 1,nq3
            j = j + 1
            sout(j) = real(qpw(i,js))
         END DO
         IF (.not.invs) THEN
            DO i = 1,nq3
               j = j + 1
               sout(j) = aimag(qpw(i,js))
            END DO
         ENDIF
         mapmt=0
         na = 1
         DO n = 1,ntype
            DO l = 0,nlh(ntypsy(na))
               DO i = 1,jri(n)
                  mapmt = mapmt +1
                  j = j + 1
                  sout(j) = rho(i,l,n,js)
               END DO
            END DO
            na = na + neq(n)
         END DO
         mapvac=0
         IF (film) THEN
            DO iv = 1,nvac
               DO k = 1,nmz
                  mapvac = mapvac + 1
                  j = j + 1
                  sout(j) = rht(k,iv,js)
               END DO
               DO k = 1,odi%nq2-1
                  DO i = 1,nmzxy
                     mapvac = mapvac + 1
                     j = j + 1
                     sout(j) =  real(rhtxy(i,k,iv,js))
                  END DO
               END DO
               IF (.not.invs2) THEN
                  DO k = 1,odi%nq2-1
                     DO i = 1,nmzxy
                        mapvac = mapvac + 1
                        j = j + 1
                        sout(j) =  aimag(rhtxy(i,k,iv,js))
                     END DO
                  END DO
               END IF
            END DO
         END IF
         IF (js .EQ. 1) nmaph = j
 10   CONTINUE   

      mapvac2=0
      IF (l_noco) THEN
c--->    off-diagonal part of the density matrix
         DO i = 1,nq3
            j = j + 1
            sout(j) = real(cdom(i))
         END DO
         DO i = 1,nq3
            j = j + 1
            sout(j) = aimag(cdom(i))
         END DO
         IF (film) THEN
            DO iv = 1,nvac
               DO k = 1,nmz
                  mapvac2 = mapvac2 + 1
                  j = j + 1
                  sout(j) = real(cdomvz(k,iv))
               END DO
               DO k = 1,odi%nq2-1
                  DO i = 1,nmzxy
                     mapvac2 = mapvac2 + 1
                     j = j + 1
                     sout(j) =  real(cdomvxy(i,k,iv))
                  END DO
               END DO
            END DO
            DO iv = 1,nvac
               DO k = 1,nmz
                  mapvac2 = mapvac2 + 1
                  j = j + 1
                  sout(j) = aimag(cdomvz(k,iv))
               END DO
               DO k = 1,odi%nq2-1
                  DO i = 1,nmzxy
                     mapvac2 = mapvac2 + 1
                     j = j + 1
                     sout(j) =  aimag(cdomvxy(i,k,iv))
                  END DO
               END DO
            END DO
            nvaccoeff2 = 2*nmzxy*(odi%nq2-1)*nvac + 2*nmz*nvac
            IF (mapvac2 .NE. nvaccoeff2) THEN
               WRITE (6,*)'The number of vaccum coefficients off the'
               WRITE (6,*)'off-diagonal part of the density matrix is'
               WRITE (6,*)'inconsitent:'
               WRITE (6,8000) mapvac2,nvaccoeff2
 8000          FORMAT ('mapvac2= ',i12,'nvaccoeff2= ',i12)
               STOP 'brysh1: # of vacuum coeff. inconsistent'
            ENDIF
         END IF
      ENDIF ! noco

      IF ( n_u > 0 ) THEN     ! lda+U
        DO js = 1,jspins
          DO n = 1, n_u
            DO k = -3, 3
              DO i = -3, 3
                j = j + 1 
                sout(j) = real(n_mmp(i,k,n,js))
                j = j + 1 
                sout(j) = aimag(n_mmp(i,k,n,js))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF (film) THEN
         nvaccoeff = vacfac*nmzxy*(odi%nq2-1)*nvac + nmz*nvac
         IF (mapvac .NE. nvaccoeff) THEN
            WRITE(6,*)'The number of vaccum coefficients is'
            WRITE(6,*)'inconsitent:'
            WRITE (6,8010) mapvac,nvaccoeff
 8010       FORMAT ('mapvac= ',i12,'nvaccoeff= ',i12)
            STOP 'brysh1: # of vacuum coeff. inconsistent'
         ENDIF
      ENDIF

      mapmtd = ntypd*(nlhd+1)*jmtd
      IF (mapmt .GT. mapmtd) THEN
         WRITE(6,*)'The number of mt coefficients is larger than the'
         WRITE(6,*)'dimensions:'
         WRITE (6,8040) mapmt,mapmtd
 8040    FORMAT ('mapmt= ',i12,' > mapmtd= ',i12)
         STOP 'brysh1: mapmt > mapmtd (dimensions)'
      ENDIF

      nmap = j
      nall = (intfac*nq3 + mapmt + mapvac + 49*2*n_u )*jspins
      IF (l_noco) nall = nall + 2*nq3 + mapvac2
      IF (nall.NE.nmap) THEN
         WRITE(6,*)'The total number of charge density coefficients is'
         WRITE(6,*)'inconsitent:'
         WRITE (6,8020) nall,nmap
 8020    FORMAT ('nall= ',i12,'not equal nmap= ',i12)
         WRITE (6,'(a,i5,a,i5)') 'nall = ',nall,' nmap = ',nmap
         STOP 'brysh1: total # of charge density coeff. inconsistent'
      ENDIF
      IF (nmap.GT.mmap) THEN 
         WRITE(6,*)'The total number of charge density coefficients is'
         WRITE(6,*)'larger than the dimensions:'
         WRITE (6,8030) nmap,mmap
 8030    FORMAT ('nmap= ',i12,' > mmap= ',i12)
         STOP 'brysh1: nmap > mmap (dimensions)'
      ENDIF

      END SUBROUTINE brysh1
      END MODULE m_brysh1
