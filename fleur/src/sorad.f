      MODULE m_sorad
c*********************************************************************
c     1. generates radial spin-orbit matrix elements
c     based on m.weinert's radsra and radsrd subroutines
c*********************************************************************
      CONTAINS
      SUBROUTINE sorad(
     > jmtd,jspd,ntypd,nwdd,lmaxd,nlod,nlo,llo,l_dulo,
     > ulo_der,ello,jspins,ntyp,nw,vr,r0,dx,jri,lmax,el0,spav,
     < rsopp,rsopdpd,rsoppd,rsopdp,ddn,
     < rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     < us,dus,uds,duds,ulos,dulos,uulon,dulon)

      USE m_constants, ONLY : c_light
      USE m_intgr,     ONLY : intgr0
      USE m_sointg
      USE m_radsra
      USE m_radsrd
      USE m_radsrdn

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jmtd,jspd,ntypd,nwdd,lmaxd,nlod
      INTEGER, INTENT (IN) :: jri,jspins,ntyp,nw,lmax,nlo
      REAL,    INTENT (IN) :: dx,r0
      LOGICAL, INTENT (IN) :: spav ! if T, spin-averaged pot is used
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: llo(nlod),ulo_der(nlod)
      REAL,    INTENT (IN) :: vr(jmtd,jspd),ello(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: el0(0:lmaxd,ntypd,jspd,nwdd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod)
      REAL,    INTENT (OUT) :: rsopp  (ntypd,lmaxd,2,2)
      REAL,    INTENT (OUT) :: rsoppd (ntypd,lmaxd,2,2)
      REAL,    INTENT (OUT) :: rsopdp (ntypd,lmaxd,2,2)
      REAL,    INTENT (OUT) :: rsopdpd(ntypd,lmaxd,2,2)
      REAL,    INTENT (OUT) :: rsoplop (ntypd,nlod,2,2)
      REAL,    INTENT (OUT) :: rsoplopd(ntypd,nlod,2,2)
      REAL,    INTENT (OUT) :: rsopdplo(ntypd,nlod,2,2)
      REAL,    INTENT (OUT) :: rsopplo (ntypd,nlod,2,2)
      REAL,    INTENT (OUT) :: rsoploplop(ntypd,nlod,nlod,2,2)
      REAL,    INTENT (OUT) ::   us(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (OUT) ::  dus(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (OUT) ::  uds(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (OUT) :: duds(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (OUT) ::  ddn(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (OUT) ::  ulos(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: dulos(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: uulon(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: dulon(nlod,ntypd,jspd)
C     ..
C     .. Local Scalars ..
      REAL ddn1,e,c,ulops,dulops,duds1
      INTEGER i,j,ir,jspin,l,noded,nodeu,ilo,ilop
C     ..
C     .. Local Arrays ..
      REAL, ALLOCATABLE :: p(:,:),pd(:,:),q(:,:),qd(:,:),plo(:,:)
      REAL, ALLOCATABLE :: plop(:,:),glo(:,:),fint(:),pqlo(:,:)
      REAL, ALLOCATABLE :: filo(:,:)
      REAL, ALLOCATABLE :: v0(:),vso(:,:),qlo(:,:)
C     ..
      c = c_light(1.0)

      IF (jri.GT.jmtd) STOP 'sorad: jri.GT.jmtd'
      ALLOCATE ( p(jmtd,2),pd(jmtd,2),q(jmtd,2),plo(jmtd,2),fint(jmtd),
     +   qlo(jmtd,2),plop(jmtd,2),qd(jmtd,2),v0(jmtd),vso(jmtd,2) )
c
      DO l = 0,lmax 

         DO jspin = 1,jspins
!
!--->    calculate normalized function at e: p and q 
!
             e = el0(l,ntyp,jspin,nw)
             CALL radsra(
     >                   e,l,vr(1,jspin),r0,dx,jri,jmtd,c,
     <                   us(l,ntyp,jspin),dus(l,ntyp,jspin),
     <                   nodeu,p(1,jspin),q(1,jspin))
!                     
!--->    calculate orthogonal energy derivative at e : pd and qd
!
             CALL radsrd(
     >                   e,l,vr(1,jspin),r0,dx,jri,jmtd,c,
     <                   uds(l,ntyp,jspin),duds(l,ntyp,jspin),
     <           ddn(l,ntyp,jspin),noded,pd(1,jspin),qd(1,jspin),
     >                   p(1,jspin),q(1,jspin),dus(l,ntyp,jspin))

         END DO     ! end of spin loop
!
!---> in case of jspins=1
!

         IF (jspins.EQ.1) THEN
           DO i = 1,jri
               p(i,2) =  p(i,1)
              pd(i,2) = pd(i,1)
           ENDDO
         ENDIF
!
!---> common spin-orbit integrant V   (average spin directions)
!                                  SO
         v0(:) = 0.0
         IF (jspins.EQ.1) THEN
           v0(1:jri) = vr(1:jri,1)
           e = el0(l,ntyp,1,nw)
         ELSE
           DO i = 1,jri
             v0(i) = (vr(i,1)+vr(i,jspins))/2.
           END DO
           e = (el0(l,ntyp,1,nw)+el0(l,ntyp,jspins,nw))/2.
         END IF

         CALL sointg(
     >               e,vr,v0,r0,dx,jri,jmtd,c,jspins,
     <               vso)
         IF (spav) THEN
           DO i= 1,jmtd
             vso(i,1)= (vso(i,1)+vso(i,2))/2.
             vso(i,2)= vso(i,1)
           ENDDO
         ENDIF

!                        s       s'            .s       s'
!-->  radial integrals <u  |V  |u  > = rsopp, <u  |V  |u  > = rsopdp etc.
!                            SO                     SO

        IF (l.GT.0) THEN ! there is no spin-orbit for s-states
        DO i = 1, 2
        DO j = 1, 2
           rsopp(ntyp,l,i,j) = radso( p(1,i), p(1,j),vso(1,i),jri,r0,dx)
          rsopdp(ntyp,l,i,j) = radso(pd(1,i), p(1,j),vso(1,i),jri,r0,dx)
          rsoppd(ntyp,l,i,j) = radso( p(1,i),pd(1,j),vso(1,i),jri,r0,dx)
         rsopdpd(ntyp,l,i,j) = radso(pd(1,i),pd(1,j),vso(1,i),jri,r0,dx)
        ENDDO
        ENDDO
        ENDIF ! l>0
!
!--->  Check for local orbitals with same l
!
       DO ilo = 1, nlo
         IF (llo(ilo).EQ.l) THEN

           DO jspin = 1,jspins
             e = ello(ilo,ntyp,jspin)
             CALL radsra(
     >                   e,l,vr(1,jspin),r0,dx,jri,jmtd,c,
     <                   ulos(ilo,ntyp,jspin),dulos(ilo,ntyp,jspin),
     <                   nodeu,plo(1,jspin),qlo(1,jspin))

c+apw+lo
             IF (l_dulo(ilo).or.ulo_der(ilo).ge.1) THEN !  calculate energy derivative (of order ulo_der) at e
               ALLOCATE (glo(jmtd,2),pqlo(jmtd,2),filo(jmtd,2))
               pqlo(1:jri,1)=plo(1:jri,jspin)
               pqlo(1:jri,2)=qlo(1:jri,jspin)
               i = ulo_der(ilo)
               IF(l_dulo(ilo)) i=1
               CALL radsrdn(
     >                  e,l,vr(1,jspin),r0,dx,jri,jmtd,c,
     <                  ulos(ilo,ntyp,jspin),duds1,ddn1,noded,glo,filo,!filo is a dummy array
     >                  pqlo,dulos(ilo,ntyp,jspin),i)
               ddn1 = sqrt(ddn1)
               IF(l_dulo(ilo)) ddn1=1.0
               plo(1:jri,jspin) = glo(1:jri,1)/ddn1
               qlo(1:jri,jspin) = glo(1:jri,2)/ddn1
               dulos(ilo,ntyp,jspin) = duds1/ddn1
               DEALLOCATE (glo,pqlo,filo)
             ENDIF
c-apw+lo
           ENDDO

           IF (jspins.EQ.1) THEN
             plo(1:jri,2) = plo(1:jri,1)
             e = (ello(ilo,ntyp,1) + el0(l,ntyp,1,nw) )/2
           ELSE
             e = (ello(ilo,ntyp,1) +  ello(ilo,ntyp,jspins) +
     +            el0(l,ntyp,1,nw) + el0(l,ntyp,jspins,nw) )/4
           END IF
           CALL sointg(
     >                 e,vr,v0,r0,dx,jri,jmtd,c,jspins,
     <                 vso)
           IF (spav) THEN
             DO i= 1,jmtd
               vso(i,1)= (vso(i,1)+vso(i,2))/2.
               vso(i,2)= vso(i,1)
             ENDDO
           ENDIF
        
           DO i = 1, 2
             DO j = 1, 2
               rsoplop (ntyp,ilo,i,j) = radso(plo(1,i),p (1,j),vso(1,i),
     +                                                        jri,r0,dx)
               rsoplopd(ntyp,ilo,i,j) = radso(plo(1,i),pd(1,j),vso(1,i),
     +                                                        jri,r0,dx)
               rsopplo (ntyp,ilo,i,j) = radso(p (1,i),plo(1,j),vso(1,i),
     +                                                        jri,r0,dx)
               rsopdplo(ntyp,ilo,i,j) = radso(pd(1,i),plo(1,j),vso(1,i),
     +                                                        jri,r0,dx)
             ENDDO
           ENDDO

           DO i = 1,jspins
             fint(:) = plo(:,i) *  p(:,i) + qlo(:,i) *  q(:,i)
             CALL intgr0(fint,r0,dx,jri,uulon(ilo,ntyp,i))
             fint(:) = plo(:,i) * pd(:,i) + qlo(:,i) * qd(:,i)
             CALL intgr0(fint,r0,dx,jri,dulon(ilo,ntyp,i))
           ENDDO

           DO ilop = 1, nlo
             IF (llo(ilop).EQ.l) THEN

               DO jspin = 1,jspins
                 e = ello(ilop,ntyp,jspin)
                 CALL radsra(
     >                       e,l,vr(1,jspin),r0,dx,jri,jmtd,c,
     <                       ulops,dulops,nodeu,plop(1,jspin),q)
c+apw+lo
                 IF (l_dulo(ilo).or.ulo_der(ilo).ge.1) THEN ! calculate orthogonal energy derivative at e
                   ALLOCATE (glo(jmtd,2),pqlo(jmtd,2),filo(jmtd,2))
                   pqlo(1:jri,1)=plop(1:jri,jspin)
                   pqlo(1:jri,2)=q(1:jri,1)
                   i = ulo_der(ilo)
                   IF(l_dulo(ilo)) i=1
                   CALL radsrdn(
     >                         e,l,vr(1,jspin),r0,dx,jri,jmtd,c,
     <                         ulops,duds1,ddn1,noded,glo,filo,!filo is a dummy array
     >                         pqlo,dulops,i)
                   plop(1:jri,jspin) = glo(1:jri,1)
                   DEALLOCATE (glo,pqlo,filo)
                 ENDIF
c-apw+lo
               ENDDO
         
               IF (jspins.EQ.1) THEN
                 plop(1:jri,2) = plop(1:jri,1)
                 e = (ello(ilo,ntyp,1) + ello(ilop,ntyp,1) )/2
               ELSE
                 e = (ello(ilo,ntyp,1) +  ello(ilo,ntyp,jspins) +  
     +                ello(ilop,ntyp,1) + ello(ilop,ntyp,jspins) )/4
               END IF
               CALL sointg(
     >                     e,vr,v0,r0,dx,jri,jmtd,c,jspins,
     <                     vso)
               IF (spav) THEN
                 DO i= 1,jmtd
                   vso(i,1)= (vso(i,1)+vso(i,2))/2.
                   vso(i,2)= vso(i,1)
                 ENDDO
               ENDIF

               DO i = 1, 2
                 DO j = 1, 2
                   rsoploplop(ntyp,ilo,ilop,i,j) =
     +                      radso(plo(1,i),plop(1,j),vso(1,i),jri,r0,dx)
                 ENDDO
               ENDDO

             ENDIF
           ENDDO

         ENDIF
       ENDDO ! end of lo-loop

      ENDDO ! end of l-loop

      DEALLOCATE ( p,pd,q,qd,plo,plop,qlo,fint,v0,vso )
!      rsoplop (:,:,:,:) = 0.0
!      rsoplopd(:,:,:,:) = 0.0
!      rsopplo (:,:,:,:) = 0.0
!      rsopdplo(:,:,:,:) = 0.0
!      rsoploplop(:,:,:,:,:) = 0.0

      END SUBROUTINE sorad
!--------------------------------------------------------------------
      REAL FUNCTION radso(
     >                    a,b,vso,jri,r0,dx)
!
!     compute radial spin-orbit integrals
!
      USE m_intgr, ONLY : intgr0
      IMPLICIT NONE
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jri
      REAL,    INTENT (IN) :: dx,r0
!     ..
!     .. Array Arguments ..
      REAL,    INTENT (IN) :: a(jri),b(jri),vso(jri)
!     ..
!     .. Local Arrays ..
      REAL q(jri)
!     ..
      q(1:jri) = a(1:jri)*b(1:jri)*vso(1:jri)
      CALL intgr0(
     >            q,r0,dx,jri,
     <            radso)

      RETURN
      END FUNCTION radso
!--------------------------------------------------------------------
      END MODULE m_sorad
