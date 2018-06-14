      MODULE m_metric
c*********************************************************
c     multiplicates the vector s_in with the metric G
c     output vector sout
c********************************************************* 
      CONTAINS
      SUBROUTINE metric(
     >                  z1,jspd,ntypd,nmzd,ntypsd,jmtd,mmap,mmaph,nmaph,
     >                  mapmt,mapvac2,jspins,l_noco,s_in,ufft,kimax,
     >                  igfft,pgfft,nq3,omtil,k1d,k2d,k3d,n2d,n3d,natd,
     >                  ntype,neq,ntypsy,nlh,jri,rmsh,dx,invs,invs2,
     >                  film,nvac,nmz,nmzxy,nq2,area,delz,nstr2,n_u,odi,
     <                  sout,lpot) 

      USE m_metrz0
      USE m_convol
      USE m_set, ONLY : dset
      USE m_od_types, ONLY : od_inp
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,ntypd,nmzd,ntypsd,jmtd,mmap,mmaph
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,n2d,n3d,natd,mapmt,mapvac2
      INTEGER, INTENT (IN) :: nmaph,jspins,nmz,nmzxy,nq2,nq3,ntype
      INTEGER, INTENT (IN) :: nvac,kimax,n_u
      REAL,    INTENT (IN) :: area,delz,omtil,z1
      LOGICAL, INTENT (IN) :: film,invs,invs2,l_noco
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),neq(ntypd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: nstr2(n2d),ntypsy(natd)
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (IN) :: ufft(0:27*k1d*k2d*k3d-1)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd),s_in(mmap)
      REAL,    INTENT (OUT):: sout(mmap)
      LOGICAL, OPTIONAL,INTENT (IN) :: lpot !do we mix a potential??
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c-odim
C     ..
C     .. Local Scalars ..
      INTEGER :: imap,ivac,iz,j,js,k2,l,n,iv2c,iv2,na,ioff
      REAL    :: dvol,dxn,dxn2,dxn4,volnstr2
      LOGICAL :: l_pot
C     ..
C     .. Local Arrays ..
      REAL,    ALLOCATABLE :: g(:),wght(:)
      COMPLEX, ALLOCATABLE :: ag3(:),fg3(:)
c
c     calculate the coefficients of the g-matrix 
c     for the m.t. and vacuum region
c
      IF (PRESENT(lpot)) THEN
         l_pot=lpot
      ELSE
         l_pot=.false.
      ENDIF
      ALLOCATE (g(mmaph),wght(nmzd),ag3(n3d),fg3(n3d))
      CALL dset(nmaph,0.0,g,1)
      IF (invs) THEN
         imap=nq3
      ELSE
         imap=2*nq3
      ENDIF
      iv2=1
      IF (.not.invs2) iv2=2
c
c metric for MT is r^2 dr/di di = r^3 dx ; divided by r^4 to 
c compensate use of n(r) * r^2 in array rho
c simpson integration used, weights for first and last point: 1
c weights forthe rest alternating: 2 or 4    
c
      na = 1
      DO 10 n = 1,ntype
         dxn = neq(n)*dx(n)/3.0e0
         dxn2 =2.0e0 *dxn
         dxn4 =4.0e0 *dxn
c
         DO l = 0,nlh(ntypsy(na))
            imap = imap + 1
            g(imap) = dxn/rmsh(1,n)
            IF (.not.l_pot) THEN
               DO j = 2,jri(n)-1,2
                  imap = imap + 2
                  g(imap-1) = dxn4/rmsh(j,n) 
                  g(imap) = dxn2/rmsh(j+1,n) 
               ENDDO
C CHANGE JR 96/12/01
c take care when jri(n) is even
               imap=imap+1-mod(jri(n),2)
               g(imap) = dxn/rmsh(jri(n),n)
            ELSE
c
c             for the potential multiply by r^4
c
               DO j = 2,jri(n)-1,2
                  imap = imap + 2
                  g(imap-1) = dxn4*rmsh(j,n)**3 
                  g(imap) = dxn2*rmsh(j+1,n)**3
               ENDDO
               imap=imap+1-mod(jri(n),2)
               g(imap) = dxn*rmsh(jri(n),n)**3
            ENDIF
   
         ENDDO
         na = na + neq(n)
 10   CONTINUE
c
c vacuum contribution
c
      IF (film) THEN
         dvol = area*delz
c     nvac=1 if (zrfs.or.invs)
         IF ( nvac.eq.1 ) dvol = dvol + dvol
         IF (odi%d1) THEN
            dvol = area*delz
         END IF
         DO 20 ivac = 1,nvac
c                                      G||=0 components
c
c---> use 7-point simpson integration in accordance to intgz0.f
c     calculate weights for integration
c
            call metr_z0(nmz,wght)
            DO iz = 1,nmz
               imap = imap + 1
               IF (odi%d1) THEN
                  g(imap) = wght(iz)*dvol*(z1+(iz-1)*delz)
               ELSE
                  g(imap) = wght(iz)*dvol
               ENDIF
            ENDDO
c                                      G||.ne.0 components
c     calculate weights for integration
            call metr_z0(nmzxy,wght)
            

            DO iv2c=1,iv2
               DO k2 = 1,odi%nq2-1
                  IF (odi%d1) THEN
                     DO iz = 1,nmzxy
                        imap = imap + 1
                        g(imap) = wght(iz)*odi%nst2(k2)*
     *                         dvol*(z1+(iz-1)*delz)
                     ENDDO
                  ELSE
                     volnstr2= dvol*nstr2(k2)
                     DO iz = 1,nmzxy
                        imap = imap + 1
                        g(imap) = wght(iz)*volnstr2
                     ENDDO
                  END IF
               ENDDO
            ENDDO
 20      CONTINUE
      END IF
      IF ( imap.NE.nmaph ) THEN
CHANGE PODL
            WRITE (6,'(" IMAP, NMAPH, MMAPH",3i10)')
            WRITE (16,'(" IMAP, NMAPH, MMAPH",3i10)')
     +                   imap,nmaph,mmaph
C
                 STOP ' metric dim >< nmaph '
      END IF
c
c     multiplicate the metric with the vector
c
      DO 40 js = 1,jspins
c    map s_in on a complex help array ag3
         IF (invs) THEN
            DO imap = 1,nq3
               ag3(imap) = cmplx(s_in(imap+nmaph*(js-1)),0.0)
            ENDDO
         ELSE
            DO imap = 1,nq3
               ag3(imap) = cmplx(s_in(imap+nmaph*(js-1)),
     +                        s_in(imap+nq3+nmaph*(js-1)))
            ENDDO
         ENDIF
         CALL convol(
     >               k1d,k2d,k3d,n3d,
     <               fg3,
     >               ag3,nq3,
     =               kimax,igfft,pgfft,ufft)
         IF (invs) THEN
            DO imap = 1,nq3
               sout(imap+nmaph*(js-1)) = omtil*real(fg3(imap))
            ENDDO
            DO imap = nq3+1,nmaph
               sout(imap+nmaph*(js-1)) = g(imap)*s_in(imap+nmaph*(js-1))
            ENDDO
         ELSE
            DO imap = 1,nq3
               sout(imap+nmaph*(js-1)) = omtil*real(fg3(imap))
               sout(imap+nq3+nmaph*(js-1)) = omtil*aimag(fg3(imap))
            ENDDO
            DO imap = 2*nq3+1,nmaph
               sout(imap+nmaph*(js-1)) = g(imap)*s_in(imap+nmaph*(js-1))
            ENDDO
         ENDIF
 40   CONTINUE

      IF (l_noco) THEN
         DO imap = 1,nq3
            ag3(imap) = cmplx(s_in(2*nmaph + imap),
     +                        s_in(2*nmaph + nq3 + imap))
         ENDDO
         CALL convol(
     >        k1d,k2d,k3d,n3d,
     <        fg3,
     >        ag3,nq3,
     =        kimax,igfft,pgfft,ufft)
         DO imap = 1,nq3
            sout(2*nmaph + imap)       = omtil*real(fg3(imap))
            sout(2*nmaph + nq3 + imap) = omtil*aimag(fg3(imap))
         ENDDO
         IF (film) THEN
c--->    js runs over the real and imaginary part of the vacuum density
c--->    coefficients (not the spin).
            IF (invs) THEN
              ioff = nq3 + mapmt
            ELSE
              ioff = 2*nq3 + mapmt
            ENDIF
            DO js = 1,2
               DO imap = 1,mapvac2/2
                  sout(2*nmaph + 2*nq3 + mapvac2/2*(js-1) + imap) = 
     +               g(ioff + imap)
     +               *s_in(2*nmaph + 2*nq3 + mapvac2/2*(js-1) + imap)
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF ( n_u > 0 )  THEN
        j = jspins*nmaph
        IF (l_noco) THEN
          j = j +  2*nq3 
          IF (film) THEN
            j = j + mapvac2
          ENDIF
        ENDIF
        DO imap = j+1, j+49*2*jspins*n_u
          sout(imap) = s_in(imap)
        ENDDO
      ENDIF
      DEALLOCATE (g,wght,ag3,fg3)
         
      END SUBROUTINE metric
      END MODULE m_metric
