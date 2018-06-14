      MODULE m_vacfun
      CONTAINS
      SUBROUTINE vacfun(
     >                  nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,jspd,
     >                  jsp,jspins,l_noco,qss,ipot,
     >                  nmzxy,nmz,invs2,delz,ig2,ig,rgphs,
     >                  bbmat,ivac,evac,bkpt,
     >                  vxy,vz,kvac1,kvac2,nv2,
     <                  tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv,wronk)
c*********************************************************************
c     determines the necessary values and derivatives on the vacuum
c     boundary (ivac=1 upper vacuum; ivac=2, lower) for energy
c     parameter evac.  also sets up the 2d hamiltonian matrices
c     necessary to update the full hamiltonian matrix.
c               m. weinert
c*********************************************************************

      USE m_intgr, ONLY : intgz0
      USE m_dotir, ONLY : dotirp
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nmzxyd,nmzd,nv2d,k1d,k2d,k3d,jspd
      INTEGER, INTENT (IN) :: jsp,jspins,ivac,nmzxy,nmz,n2d,n3d,ipot
      LOGICAL, INTENT (IN) :: invs2,l_noco
      REAL,    INTENT (IN) :: delz
      REAL,    INTENT (OUT) :: wronk
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nv2(jspd)
      INTEGER, INTENT (IN) :: kvac1(nv2d,jspd),kvac2(nv2d,jspd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      COMPLEX, INTENT (IN) :: vxy(nmzxyd,n2d-1)
      COMPLEX, INTENT (OUT):: tddv(nv2d,nv2d),tduv(nv2d,nv2d)
      COMPLEX, INTENT (OUT):: tudv(nv2d,nv2d),tuuv(nv2d,nv2d)
      REAL,    INTENT (IN) :: vz(nmzd,2,4),bbmat(3,3),evac(2,jspd)
      REAL,    INTENT (IN) :: qss(3)
      REAL,    INTENT (IN) :: bkpt(3),rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (OUT):: udz(nv2d,jspd),uz(nv2d,jspd)
      REAL,    INTENT (OUT):: dudz(nv2d,jspd),duz(nv2d,jspd)
      REAL,    INTENT (OUT):: ddnv(nv2d,jspd)
C     ..
C     .. Local Scalars ..
      REAL ev,phase,scale,xv,yv,vzero
      INTEGER i,i1,i2,i3,ik,ind2,ind3,jk,np1,jspin,jsp1,jsp2
      LOGICAL tail
C     ..
C     .. Local Arrays ..
      REAL u(nmzd,nv2d,jspd),ud(nmzd,nv2d,jspd),v(3),x(nmzd),
     +     qssbti(2,2)
C     ..
C     .. External Subroutines ..
      EXTERNAL vacudz,vacuz
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC aimag,cmplx,conjg,real
C     ..
      tail = .true.
      np1 = nmzxy + 1
c--->    wronksian for the schrodinger equation given by an identity
      wronk = 2.0
c---> setup the spin-spiral q-vector
      qssbti(1,1) = - qss(1)/2
      qssbti(2,1) = - qss(2)/2
      qssbti(1,2) = + qss(1)/2
      qssbti(2,2) = + qss(2)/2
c--->    generate basis functions for each 2-d k+g
      DO jspin = 1,jspins
         DO 20 ik = 1,nv2(jspin)
            v(1) = bkpt(1) + kvac1(ik,jspin) + qssbti(1,jspin)
            v(2) = bkpt(2) + kvac2(ik,jspin) + qssbti(2,jspin)
            v(3) = 0.0
            ev = evac(ivac,jspin) - 0.5*dotirp(v,v,bbmat)
            vzero = vz(nmzd,ivac,jspin)
            CALL vacuz(ev,vz(1,ivac,jspin),vzero,nmz,delz,
     +                uz(ik,jspin),duz(ik,jspin),u(1,ik,jspin))
            CALL vacudz(ev,vz(1,ivac,jspin),vzero,nmz,delz,
     +                  udz(ik,jspin),dudz(ik,jspin),ddnv(ik,jspin),
     +                  ud(1,ik,jspin),duz(ik,jspin),u(1,ik,jspin))
c--->       make sure the solutions satisfy the wronksian
            scale = wronk/ (udz(ik,jspin)*duz(ik,jspin)-
     -                         dudz(ik,jspin)*uz(ik,jspin))
            udz(ik,jspin) = scale*udz(ik,jspin)
            dudz(ik,jspin) = scale*dudz(ik,jspin)
            ddnv(ik,jspin) = scale*ddnv(ik,jspin)
            DO 10 i = 1,nmz
               ud(i,ik,jspin) = scale*ud(i,ik,jspin)
 10         CONTINUE
 20      CONTINUE
      ENDDO
c--->    set up the tuuv, etc. matrices
      IF (l_noco) THEN
         IF (ipot.EQ.1) THEN
            jsp1 = 1
            jsp2 = 1
         ELSEIF (ipot.EQ.2) THEN
            jsp1 = 2
            jsp2 = 2
         ELSEIF (ipot.EQ.3) THEN
            jsp1 = 2
            jsp2 = 1
         ENDIF
      ELSE
         jsp1 = jsp
         jsp2 = jsp
      ENDIF
      DO 120 ik = 1,nv2(jsp1)
         DO 110 jk = 1,nv2(jsp2)

c--->     determine the warping component of the potential
          i1 = kvac1(ik,jsp1) - kvac1(jk,jsp2)
          i2 = kvac2(ik,jsp1) - kvac2(jk,jsp2)
          i3 = 0
          ind3 = ig(i1,i2,i3)
          IF (ind3.EQ.0) GO TO 110
          phase = rgphs(i1,i2,i3)
          ind2 = ig2(ind3)
          IF (ind2.EQ.0) THEN
             WRITE (16,FMT=8000) ik,jk
             WRITE (6,FMT=8000) ik,jk
 8000        FORMAT (' **** error in map2 for 2-d stars',2i5)
             STOP 'vacfun'
          END IF
c--->     get the proper warping index (vxy starts with the 2nd star)
          ind2 = ind2 - 1
          IF (ind2.NE.0) THEN
c--->       only the warping part, 1st star (G=0) is done later

c--->       obtain the warping matrix elements
c--->       note that the tail correction (tail=.true.) is included for
c--->       the integrals, i.e. the integrand is from infinity inward

c--->       tuuv
            DO 30 i = 1,nmzxy
               x(np1-i) = u(i,ik,jsp1)*u(i,jk,jsp2)*real(vxy(i,ind2))
 30         CONTINUE
            CALL intgz0(x,delz,nmzxy,xv,tail)
            DO 40 i = 1,nmzxy
               x(np1-i) = u(i,ik,jsp1)*u(i,jk,jsp2)*aimag(vxy(i,ind2))
 40         CONTINUE
            CALL intgz0(x,delz,nmzxy,yv,tail)
            tuuv(ik,jk) = phase*cmplx(xv,yv)

c--->       tddv
            DO 50 i = 1,nmzxy
               x(np1-i) = ud(i,ik,jsp1)*ud(i,jk,jsp2)*real(vxy(i,ind2))
 50         CONTINUE
            CALL intgz0(x,delz,nmzxy,xv,tail)
            DO 60 i = 1,nmzxy
               x(np1-i) =ud(i,ik,jsp1)*ud(i,jk,jsp2)*aimag(vxy(i,ind2))
 60         CONTINUE
            CALL intgz0(x,delz,nmzxy,yv,tail)
            tddv(ik,jk) = phase*cmplx(xv,yv)

c--->       tudv
            DO 70 i = 1,nmzxy
               x(np1-i) = u(i,ik,jsp1)*ud(i,jk,jsp2)*real(vxy(i,ind2))
 70         CONTINUE
            CALL intgz0(x,delz,nmzxy,xv,tail)
            DO 80 i = 1,nmzxy
               x(np1-i) = u(i,ik,jsp1)*ud(i,jk,jsp2)*aimag(vxy(i,ind2))
 80         CONTINUE
            CALL intgz0(x,delz,nmzxy,yv,tail)
            tudv(ik,jk) = phase*cmplx(xv,yv)

c--->       tduv
            DO 90 i = 1,nmzxy
               x(np1-i) = ud(i,ik,jsp1)*u(i,jk,jsp2)*real(vxy(i,ind2))
 90         CONTINUE
            CALL intgz0(x,delz,nmzxy,xv,tail)
            DO 100 i = 1,nmzxy
               x(np1-i) = ud(i,ik,jsp1)*u(i,jk,jsp2)*aimag(vxy(i,ind2))
 100        CONTINUE
            CALL intgz0(x,delz,nmzxy,yv,tail)
            tduv(ik,jk) = phase*cmplx(xv,yv)

          ELSE

c--->       diagonal (film muffin-tin) terms
            IF ((ipot.EQ.1) .OR. (ipot.EQ.2)) THEN
               tuuv(ik,ik) = cmplx(evac(ivac,jsp1),0.0)
               tddv(ik,ik) = cmplx(evac(ivac,jsp1)*ddnv(ik,jsp1),0.0)
               tudv(ik,ik) = cmplx(0.5,0.0)
               tduv(ik,ik) = cmplx(0.5,0.0)
            ELSE

c--->          tuuv
               DO i = 1,nmz
                  x(nmz+1-i) = u(i,ik,jsp1)*u(i,jk,jsp2)*vz(i,ivac,3)
               ENDDO
               CALL intgz0(x,delz,nmz,xv,tail)
               DO i = 1,nmz
                  x(nmz+1-i) = u(i,ik,jsp1)*u(i,jk,jsp2)*vz(i,ivac,4)
               ENDDO
               CALL intgz0(x,delz,nmz,yv,tail)
               tuuv(ik,jk) = cmplx(xv,yv)
               
c--->          tddv
               DO i = 1,nmz
                  x(nmz+1-i) = ud(i,ik,jsp1)*ud(i,jk,jsp2)*vz(i,ivac,3)
               ENDDO
               CALL intgz0(x,delz,nmz,xv,tail)
               DO i = 1,nmz
                  x(nmz+1-i) = ud(i,ik,jsp1)*ud(i,jk,jsp2)*vz(i,ivac,4)
               ENDDO
               CALL intgz0(x,delz,nmz,yv,tail)
               tddv(ik,jk) = cmplx(xv,yv)

c--->          tudv
               DO i = 1,nmz
                  x(nmz+1-i) = u(i,ik,jsp1)*ud(i,jk,jsp2)*vz(i,ivac,3)
               ENDDO
               CALL intgz0(x,delz,nmz,xv,tail)
               DO i = 1,nmz
                  x(nmz+1-i) = u(i,ik,jsp1)*ud(i,jk,jsp2)*vz(i,ivac,4)
               ENDDO
               CALL intgz0(x,delz,nmz,yv,tail)
               tudv(ik,jk) = cmplx(xv,yv)
               
c--->          tduv
               DO i = 1,nmz
                  x(nmz+1-i) = ud(i,ik,jsp1)*u(i,jk,jsp2)*vz(i,ivac,3)
               ENDDO
               CALL intgz0(x,delz,nmz,xv,tail)
               DO i = 1,nmz
                  x(nmz+1-i) = ud(i,ik,jsp1)*u(i,jk,jsp2)*vz(i,ivac,4)
               ENDDO
               CALL intgz0(x,delz,nmz,yv,tail)
               tduv(ik,jk) = cmplx(xv,yv)
            ENDIF

          ENDIF
 110    CONTINUE
 120  CONTINUE

      END SUBROUTINE vacfun
      END MODULE m_vacfun

