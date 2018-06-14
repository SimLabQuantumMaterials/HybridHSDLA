      MODULE m_od_vacfun
      CONTAINS
      SUBROUTINE od_vacfun(
     >     m_cyl,z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,jspd,
     >     jsp,jspins,l_noco,qss,ipot,ig,ig1,tpi,
     >     nmzxy,nmz,delz,ig2,n2d_1,
     >     bbmat,ivac,evac,bkpt,MM,vM,
     >     vxy,vz,kvac3,nv2,
     <     tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv)
c*********************************************************************
c     determines the necessary values and derivatives on the cylindrical
c     vacuum boundary for energy parameter evac. also sets up the
c     hamiltonian matrices necessary to update the full hamiltonian
c     matrix.
c     Y. Mokrousov, June 2002
c*********************************************************************
      
      USE m_intgr, ONLY : intgz0
      USE m_dotir, ONLY : dotirp
      USE m_constants, ONLY : pimach
      
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
    
      INTEGER, INTENT (IN) :: nmzxyd,nmzd,nv2d,k1d,k2d,k3d,jspd,n3d
      INTEGER, INTENT (IN) :: jsp,jspins,nmzxy,nmz,ipot,MM,n2d,vM
      INTEGER, INTENT (IN) :: ivac,n2d_1,m_cyl
      REAL,    INTENT (IN) :: delz,z1,tpi
      LOGICAL, INTENT (IN) :: l_noco
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nv2(jspd)
      INTEGER, INTENT (IN) :: kvac3(nv2d,jspd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      INTEGER, INTENT (IN) :: ig1(-k3d:k3d,-MM:MM),ig2(n3d)
      COMPLEX, INTENT (IN) :: vxy(nmzxyd,n2d_1-1)
      COMPLEX, INTENT (OUT):: tddv(-vM:vM,-vM:vM,nv2d,nv2d)
      COMPLEX, INTENT (OUT):: tduv(-vM:vM,-vM:vM,nv2d,nv2d)
      COMPLEX, INTENT (OUT):: tudv(-vM:vM,-vM:vM,nv2d,nv2d)
      COMPLEX, INTENT (OUT):: tuuv(-vM:vM,-vM:vM,nv2d,nv2d)
      REAL,    INTENT (IN) :: vz(nmzd,2,4),bbmat(3,3),evac(2,jspd)
      REAL,    INTENT (IN) :: bkpt(3),qss(3)
      REAL,    INTENT (OUT):: udz(-vM:vM,nv2d,jspd),uz(-vM:vM,nv2d,jspd)
      REAL,    INTENT (OUT):: dudz(-vM:vM,nv2d,jspd)
      REAL,    INTENT (OUT):: duz(-vM:vM,nv2d,jspd)
      REAL,    INTENT (OUT):: ddnv(-vM:vM,nv2d,jspd)
C     ..
C     .. Local Scalars ..
      REAL ev,scale,xv,yv,vzero,v1,wronk
      INTEGER i,ik,jk,np1,jspin,jsp1,jsp2,m,l
      INTEGER i1,i2,i3,ind1,ind3
      LOGICAL tail
C     ..
C     .. Local Arrays ..
      REAL wdz(-vM:vM,nv2d,jspd),wz(-vM:vM,nv2d,jspd)
      REAL dwdz(-vM:vM,nv2d,jspd),dwz(-vM:vM,nv2d,jspd)
      REAL u(nmzd,-vM:vM,nv2d,jspd),ud(nmzd,-vM:vM,nv2d,jspd)
      REAL v(3),x(nmzd)
      REAL vr0(nmzd,2,4)
      REAL w(nmzd,-vM:vM,nv2d,jspd),wd(nmzd,-vM:vM,nv2d,jspd)
      REAL qssbti(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL vacudz,vacuz
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC aimag,cmplx,conjg,real,sqrt
C     ..

      tail = .true.
      np1 = nmzxy + 1
      
c     wronksian for the schrodinger equation given by an identity

      wronk = 2.0

      qssbti(1) = - qss(3)/2.
      qssbti(2) = + qss(3)/2.

      tuuv(:,:,:,:) = cmplx(0.,0.)
      tudv(:,:,:,:) = cmplx(0.,0.)
      tduv(:,:,:,:) = cmplx(0.,0.)
      tddv(:,:,:,:) = cmplx(0.,0.)

c     generate basis functions for each 1-d k_z+g_z and m

      DO jspin = 1,jspins
         DO 20 ik = 1,nv2(jspin)
            DO 25 m = 0,vM
               v(1) = 0.0
               v(2) = 0.0
               v(3) = bkpt(3) + kvac3(ik,jspin) + qssbti(jspin)
               ev = evac(ivac,jspin) - 0.5*dotirp(v,v,bbmat)
c     constructing of the 'pseudopotential'
               DO 55 i=1,nmzd
                  v1 = 1./(8.*((z1+(i-1)*delz)**2))
     -                 -(m*m)/(2.*((z1+(i-1)*delz)**2))
                  vr0(i,ivac,jspin) = vz(i,ivac,jspin)-v1
 55            CONTINUE
               vzero = vr0(nmzd,ivac,jspin)
c     obtaining solutions with the 'pseudopotential'
               
               CALL vacuz(ev,vr0(1,ivac,jspin),vzero,nmz,delz,
     +              wz(m,ik,jspin),dwz(m,ik,jspin),w(1,m,ik,jspin))
               CALL vacudz(ev,vr0(1,ivac,jspin),vzero,nmz,delz,
     +              wdz(m,ik,jspin),dwdz(m,ik,jspin),ddnv(m,ik,jspin),
     +              wd(1,m,ik,jspin),dwz(m,ik,jspin),w(1,m,ik,jspin))
c     make sure the solutions satisfy the wronksian
               scale = wronk/ (wdz(m,ik,jspin)*dwz(m,ik,jspin)-
     -              dwdz(m,ik,jspin)*wz(m,ik,jspin))
               wdz(m,ik,jspin) = scale*wdz(m,ik,jspin)
               dwdz(m,ik,jspin) = scale*dwdz(m,ik,jspin)
               ddnv(m,ik,jspin) = scale*ddnv(m,ik,jspin)
               IF (m.GT.0) THEN
                  wdz(-m,ik,jspin) = wdz(m,ik,jspin)
                  dwdz(-m,ik,jspin) = dwdz(m,ik,jspin)
                  ddnv(-m,ik,jspin) = ddnv(m,ik,jspin)
               END IF
               DO 10 i = 1,nmz
                  wd(i,m,ik,jspin) = scale*wd(i,m,ik,jspin)
                  w(i,m,ik,jspin) = scale*w(i,m,ik,jspin)
                  IF (m.GT.0) THEN
                     wd(i,-m,ik,jspin) = wd(i,m,ik,jspin)
                     w(i,-m,ik,jspin) = w(i,m,ik,jspin)
                  END IF   
 10            CONTINUE
c     constructing 'real' solutions
               DO 65 i=1,nmz
                  u(i,m,ik,jspin)=w(i,m,ik,jspin)/sqrt(z1+(i-1)*delz)
                  ud(i,m,ik,jspin)=wd(i,m,ik,jspin)/sqrt(z1+(i-1)*delz)
                  IF (m.GT.0) THEN
                     u(i,-m,ik,jspin) = u(i,m,ik,jspin)
                     ud(i,-m,ik,jspin) = ud(i,m,ik,jspin)
                  END IF
 65            CONTINUE    
               duz(m,ik,jspin)=(-dwz(m,ik,jspin))/sqrt(z1)-
     -              wz(m,ik,jspin)/(2.0*((z1)**(1.5)))
               uz(m,ik,jspin)=wz(m,ik,jspin)/sqrt(z1)
               dudz(m,ik,jspin)=(-dwdz(m,ik,jspin))/sqrt(z1)-
     -              wdz(m,ik,jspin)/(2.0*((z1)**(1.5))) 
               udz(m,ik,jspin)=wdz(m,ik,jspin)/sqrt(z1)
               IF (m.GT.0) THEN
                  duz(-m,ik,jspin) = duz(m,ik,jspin)
                  uz(-m,ik,jspin) = uz(m,ik,jspin)
                  dudz(-m,ik,jspin) = dudz(m,ik,jspin)
                  udz(-m,ik,jspin) = udz(m,ik,jspin)
               END IF
 25         CONTINUE 
 20      CONTINUE
      ENDDO
      
c     set up the tuuv, etc. matrices

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
         DO 115 jk = 1,nv2(jsp2)
            i1 = 0
            i2 = 0
            i3 = kvac3(ik,jsp1) - kvac3(jk,jsp2)
            ind3 = ig(i1,i2,i3)
            IF (ind3.EQ.0) GO TO 115
            DO 110 m = -vM,vM
               DO 105 l = -vM,vM
                IF (l.EQ.m .OR. (iabs(m).LE.m_cyl
     &                             .AND. iabs(l).LE.m_cyl)) THEN
c     determine the warping component of the potential
                  ind1 = ig1(i3,m-l)
                  IF (ind1.NE.0) THEN
                     IF(ind1.NE.1) THEN
                        ind1 = ind1 - 1
c     only the warping part
c--->             tuuv
                        DO i = 1,nmzxy
                           x(np1-i) = w(i,m,ik,jsp1)*w(i,l,jk,jsp2)
     *                          *real(vxy(i,ind1))
                        ENDDO
                        CALL intgz0(x,delz,nmzxy,xv,tail)
                        DO i = 1,nmzxy
                           x(np1-i) = w(i,m,ik,jsp1)*w(i,l,jk,jsp2)
     *                          *aimag(vxy(i,ind1))
                        ENDDO
                        CALL intgz0(x,delz,nmzxy,yv,tail)
                        tuuv(m,l,ik,jk) = cmplx(xv,yv)
c--->             tddv
                        DO i = 1,nmzxy
                           x(np1-i) = wd(i,m,ik,jsp1)*wd(i,l,jk,jsp2)
     *                          *real(vxy(i,ind1))
                        ENDDO
                        CALL intgz0(x,delz,nmzxy,xv,tail)
                        DO i = 1,nmzxy
                           x(np1-i) = wd(i,m,ik,jsp1)*wd(i,l,jk,jsp2)
     *                          *aimag(vxy(i,ind1))
                        ENDDO
                        CALL intgz0(x,delz,nmzxy,yv,tail)
                        tddv(m,l,ik,jk) = cmplx(xv,yv)
c--->             tudv
                        DO i = 1,nmzxy
                           x(np1-i) = w(i,m,ik,jsp1)*wd(i,l,jk,jsp2)
     *                          *real(vxy(i,ind1))
                        ENDDO
                        CALL intgz0(x,delz,nmzxy,xv,tail)
                        DO i = 1,nmzxy
                           x(np1-i) = w(i,m,ik,jsp1)*wd(i,l,jk,jsp2)
     *                          *aimag(vxy(i,ind1))
                        ENDDO
                        CALL intgz0(x,delz,nmzxy,yv,tail)
                        tudv(m,l,ik,jk) = cmplx(xv,yv)
c--->             tduv
                        DO i = 1,nmzxy
                           x(np1-i) = wd(i,m,ik,jsp1)*w(i,l,jk,jsp2)
     *                          *real(vxy(i,ind1))
                        ENDDO
                        CALL intgz0(x,delz,nmzxy,xv,tail)
                        DO i = 1,nmzxy
                           x(np1-i) = wd(i,m,ik,jsp1)*w(i,l,jk,jsp2)
     *                          *aimag(vxy(i,ind1))
                        ENDDO
                        CALL intgz0(x,delz,nmzxy,yv,tail)
                        tduv(m,l,ik,jk) = cmplx(xv,yv)
                     ELSE
c--->          diagonal terms
                      IF ((ipot.EQ.1).OR.(ipot.EQ.2)) THEN
                        tuuv(m,m,ik,ik) = cmplx(evac(ivac,jsp1),0.0)
                        tddv(m,m,ik,ik) = cmplx(evac(ivac,jsp1)*
     *                       ddnv(m,ik,jsp1),0.0)
                        tudv(m,m,ik,ik) = cmplx(0.5,0.0)
                        tduv(m,m,ik,ik) = cmplx(0.5,0.0)
                      ELSE
c--->             tuuv
                        DO i = 1,nmz
                           x(nmz+1-i) = w(i,m,ik,jsp1)*w(i,l,jk,jsp2)
     *                          *vz(i,ivac,3)
                        ENDDO
                        CALL intgz0(x,delz,nmz,xv,tail)
                        DO i = 1,nmz
                           x(nmz+1-i) = w(i,m,ik,jsp1)*w(i,l,jk,jsp2)
     *                          *vz(i,ivac,4)
                        ENDDO
                        CALL intgz0(x,delz,nmz,yv,tail)
                        tuuv(m,l,ik,jk) = cmplx(xv,yv)
c--->             tddv
                        DO i = 1,nmz
                           x(nmz+1-i) = wd(i,m,ik,jsp1)*wd(i,l,jk,jsp2)
     *                          *vz(i,ivac,3)
                        ENDDO
                        CALL intgz0(x,delz,nmz,xv,tail)
                        DO i = 1,nmz
                           x(nmz+1-i) = wd(i,m,ik,jsp1)*wd(i,l,jk,jsp2)
     *                          *vz(i,ivac,4)
                        ENDDO
                        CALL intgz0(x,delz,nmz,yv,tail)
                        tddv(m,l,ik,jk) = cmplx(xv,yv)
c--->             tudv
                        DO i = 1,nmz
                           x(nmz+1-i) = w(i,m,ik,jsp1)*wd(i,l,jk,jsp2)
     *                          *vz(i,ivac,3)
                        ENDDO
                        CALL intgz0(x,delz,nmz,xv,tail)
                        DO i = 1,nmz
                           x(nmz+1-i) = w(i,m,ik,jsp1)*wd(i,l,jk,jsp2)
     *                          *vz(i,ivac,4)
                        ENDDO
                        CALL intgz0(x,delz,nmz,yv,tail)
                        tudv(m,l,ik,jk) = cmplx(xv,yv)
c--->             tduv
                        DO i = 1,nmz
                           x(nmz+1-i) = wd(i,m,ik,jsp1)*w(i,l,jk,jsp2)
     *                          *vz(i,ivac,3)
                        ENDDO
                        CALL intgz0(x,delz,nmz,xv,tail)
                        DO i = 1,nmz
                           x(nmz+1-i) = wd(i,m,ik,jsp1)*w(i,l,jk,jsp2)
     *                          *vz(i,ivac,4)
                        ENDDO
                        CALL intgz0(x,delz,nmz,yv,tail)
                        tduv(m,l,ik,jk) = cmplx(xv,yv)

                      ENDIF !ipot
                     END IF !ind1 ne 1
                  ENDIF    ! ind1 ne 0
                END IF
 105           CONTINUE 
 110        CONTINUE
 115     CONTINUE 
 120  CONTINUE

       
      RETURN
      END SUBROUTINE od_vacfun
      END MODULE m_od_vacfun




      
