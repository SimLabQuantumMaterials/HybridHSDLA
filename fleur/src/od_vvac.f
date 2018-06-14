      MODULE m_od_vvac
      CONTAINS
      SUBROUTINE od_vvac(
     >     n3d,k1d,k2d,k3d,n2d,jspd,nvac,
     >     ig,ig2,sk2,phi2,nq3,nq2,
     >     z1,nmz,nmzd,delz,
     >     psq,rht,
     <     vz)

c     subroutine which calculates the non warped part of the
c     vacuum potential (m=0,gz=0)
c                               Y. Mokrousov
c     the potential in this subroutine can be defined in two
c     equivalent ways, which nevertheless give a bit defferent 
c     results, 2nd one seems to be more precise 

      USE m_constants, only : pimach
      USE m_qsf
      USE m_od_cylbes

      IMPLICIT NONE


      INTEGER, INTENT (IN) :: n3d,k1d,k2d,k3d,n2d,nmz,nmzd
      INTEGER, INTENT (IN) :: nq3,nq2
      INTEGER, INTENT (IN) :: jspd,nvac
      REAL,    INTENT (IN) :: z1,delz
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      INTEGER, INTENT (IN) :: ig2(n3d)
      REAL,    INTENT (IN) :: sk2(n2d),phi2(n2d)
      COMPLEX, INTENT (IN) :: psq(n3d)
      REAL,    INTENT (IN) :: rht(nmzd,2,jspd)
      REAL,    INTENT (OUT) :: vz(nmzd,2,jspd)

      COMPLEX  rhobar
      INTEGER  k1,k2,irec3,irec2,i,j,ivac,imz,imz1
      REAL     g2,phi,fpi,tpi,a(nmzd)
      REAL     fJ,z,zp
      REAL     rht1(nmzd)
      REAL     f2(nmzd),f22(nmzd)

      INTRINSIC cmplx

      fpi = 4*pimach()
      tpi = 2*pimach()
       
      DO i = 1,nmz
         f2(i) = 0.
         f22(i) = 0.
         DO ivac = 1,nvac
            vz(i,ivac,1) = 0.
         END DO
      END DO
        
         
      rhobar = -psq(1)
            
      DO 10 k1 = -k1d,k1d
         DO 20 k2 = -k2d,k2d
            irec3 = ig(k1,k2,0)
            IF (irec3.NE.0) THEN
               irec2 = ig2(irec3)
               IF (irec2.NE.1) THEN
                  g2 = sk2(irec2)
                  phi = phi2(irec2)
                  CALL od_cylbes(1,z1*g2,fJ)
                  rhobar = rhobar - 2.*psq(irec3)*
     *                 cmplx(fJ/(g2*z1),0.0)
                  
               END IF
            END IF
 20      CONTINUE
 10   CONTINUE

c----> 1st equivalent way      
      
      DO 30 i=1,nmz
         rht1(i) = fpi*(z1+(i-1)*delz)*rht(i,1,1)
 30   CONTINUE
      
      CALL qsf(delz,rht1(1),f2(1),nmz,1)
      
      DO 40 i = 1,nmz
         f2(i) = tpi*z1*z1*rhobar-f2(i)
 40   CONTINUE
      
      DO 50 i = 1,nmz
         DO 60 j = 1,nmz
            IF (j.lt.i) THEN
               f22(j) = 0.0
            ELSE
               f22(j) = f2(j)/(z1+delz*(j-1))
            END IF
 60      CONTINUE
         CALL qsf(delz,f22(1),a,nmz,0)
         DO 70 ivac =1,nvac
            vz(i,ivac,1) = -a(1)
 70      CONTINUE
 50   CONTINUE

c----> 2nd equivalent way (via the Green function)

      DO imz = 1,nmz
         z = z1 + (imz-1)*delz
         DO imz1 = 1,nmz
            zp = z1 +  (imz1-1)*delz
            IF (imz1.LE.imz) THEN
               rht1(imz1) = fpi*log(z)*zp*rht(imz1,1,1)
            ELSE
               rht1(imz1) = fpi*log(zp)*zp*rht(imz1,1,1)
            END IF
            
         END DO
         CALL qsf(delz,rht1,a,nmz,0)
         vz(imz,1,1) = tpi*log(z)*(z1*z1)*rhobar - a(1)
      END DO

      RETURN
      END SUBROUTINE od_vvac
      END MODULE m_od_vvac


      
   
   






      



      

