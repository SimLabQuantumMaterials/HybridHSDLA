      MODULE m_od_vvacis
      CONTAINS
      SUBROUTINE od_vvacis(
     >     n2d_1,jspd,nmzd,nmzxyd,nmzxy,nq2_1,
     >     kv2_1,nvac,z1,delz,MM,nq3,nstr2,
     >     k3d,bmat,sk3,amat,nq2,nstr2_1,
     >     k1d,k2d,sk2,phi2,ig1,ig,ig2,n2d,n3d,
     >     igfft2,pgfft2,kimax2,kv2,kv3,omtil,
     >     rht,rhtxy,psq,vz,zrfs,invs,
     >     kimax,nstr,igfft,pgfft,
     <     vxy,vpw)

c     generates m^2+gz^2.neq.0 coefficients of the vacuum potential
c     and fourier coefficients of the interstitial potential
c     in the case of 1-dimensional calculations
c     based on the using of the green's function for the 
c     radial equation for the potential:
c     [1/r(d/dr(r*d/dr))-m^2/r^2-gz^2]V(r) = 4pi*rho(r)
c     when gz.neq.0 green's function is represented in
c     the following way: 
c          gf = 4pi*I_m(gz*r<)*K_m(gz*r>)
c     where I and K are modified bessel funcions
c     in the case of gz.eq.0 the green's function is:
c            gf = 2pi*((r</r>)^m)/m
c                 Y.Mokrousov, autumn 2002

c     Fully symmetrized version, Mai 2003
      
      USE m_constants, ONLY : pimach
      USE m_od_cylbes
      USE m_modcyli
      USE m_modcylk
      USE m_vacp5_0
      USE m_vacp5_z
      USE m_visp5_0
      USE m_visp5_z
      USE m_angle
      USE m_qsf
      USE m_fft2d
      USE m_set,   ONLY : dset
      
      IMPLICIT NONE
      
      LOGICAL, INTENT (IN) :: zrfs,invs
      INTEGER, INTENT (IN) :: n2d_1,jspd,nmzxyd,MM,nmzd
      INTEGER, INTENT (IN) :: k1d,k2d,n2d,n3d,nq3
      INTEGER, INTENT (IN) :: nvac,nmzxy,nq2_1,k3d,nq2
      INTEGER, INTENT (IN) :: kimax2,kimax
      REAL,    INTENT (IN) :: z1,delz,omtil
      INTEGER, INTENT (IN) :: nstr2(n2d),nstr(n3d)
      INTEGER, INTENT (IN) :: nstr2_1(n2d_1)
      INTEGER, INTENT (IN) :: kv2_1(2,n2d_1),kv2(2,n2d)
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,2)
      INTEGER, INTENT (IN) :: kv3(3,n3d)
      INTEGER, INTENT (IN) :: igfft2(0:(2*k1d+1)*(2*k2d+1)-1,2)
      REAL,    INTENT (IN) :: pgfft2(0:(2*k1d+1)*(2*k2d+1)-1)
      REAL,    INTENT (IN) :: bmat(3,3),amat(3,3)
      REAL,    INTENT (IN) :: sk2(n2d),phi2(n2d),sk3(n3d)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      INTEGER, INTENT (IN) :: ig2(n3d),ig1(-k3d:k3d,-MM:MM)
      COMPLEX, INTENT (INOUT) :: psq(n3d)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
      REAL,    INTENT (IN) :: vz(nmzd,2,jspd) 
      REAL,    INTENT (IN) :: rht(nmzd,2,jspd)
      COMPLEX, INTENT (IN) :: rhtxy(nmzxyd,n2d_1-1,2,jspd)
      COMPLEX, INTENT (OUT):: vxy(nmzxyd,n2d_1-1,2,jspd)
      COMPLEX, INTENT (OUT):: vpw(n3d,jspd)
      
c     local

c------> reciprocal vectors staff

      INTEGER irec2,irec3,k1,k2,gzi,gzi1,k3,gzmin
      REAL    gz,gx,gy,g,g2,phi

c------> different garbage

      INTEGER m,imz,imz1,i,ivac,iirec1,m1,j,irc1,im,irc
      INTEGER ix,iy,iz,s,n
      REAL    x,y,r,z,b,q,zf,tol_21
      REAL    tpi,fpi,mult
      REAL    ani1,ani2,rhti
      INTEGER ivfft1,ivfft2,ivfft2d
      COMPLEX ic
     
c------> different staff with special functions

      REAL  IR(1:k3d,0:MM)
      REAL, ALLOCATABLE :: KK(:),II(:)
      REAL  fJ,fJ2,fJ1,iJ2,IIIR,fJ3
      REAL, ALLOCATABLE :: III(:),IIII(:,:,:)
      REAL, ALLOCATABLE :: fJJ(:,:),iJJ(:,:)

c------> used for the optimization 

      INTEGER irec1(4),l,lmin,lmax,l1
       
c------> some factors used for the construction of the ch. density

      REAL, ALLOCATABLE :: fact(:)
      COMPLEX aa,a

c------> values of the potential on the vacuum boundary

      COMPLEX val_help
      COMPLEX, ALLOCATABLE :: val(:),val_m(:,:)
      
c------> Charge density

      COMPLEX, ALLOCATABLE :: rxy(:)
     
c---------------> POTENTIALS

c-> interstitial potential by its gz,m - components and total on the 
c-> real grid 

      COMPLEX, ALLOCATABLE :: vis(:,:,:),vis_tot(:,:)

c-> potential in the vacuum caused by the vacuum charge density      

      COMPLEX, ALLOCATABLE :: pvac(:)

c-> potential in the vacuum caused by the interst. density

      COMPLEX, ALLOCATABLE :: pint(:)
 
c-> radial components of the potential caused by the vacuum  
c-> density on the real grid  

      COMPLEX, ALLOCATABLE :: vis_help(:,:),vpw_help(:)

c-> real grid vis_z,vis_0 for the fft2d

      REAL,ALLOCATABLE :: af2(:),bf2(:)

c-> radial grids

      INTEGER, ALLOCATABLE :: rmap(:,:)
      REAL   , ALLOCATABLE :: rr(:)
      INTEGER nrd

c--> time

      REAL    gxy0,fxy0
      COMPLEX gxy(n2d-1)
      COMPLEX fxy(n2d-1)

      INTRINSIC real,aimag

c--------- preparations ---------->

      ALLOCATE ( KK(nmzxyd),II(nmzxyd),III(9*k1d*k2d),
     &     IIII(9*k1d*k3d,1:k3d,0:MM),
     &     rmap(0:3*k1d-1,0:3*k2d-1),rr(1:9*k1d*k2d),
     &     fact(nmzxyd),val(n2d_1),fJJ(0:MM+1,n2d),
     &     val_m(-k3d:k3d,-MM:MM),rxy(nmzxyd),iJJ(0:MM+1,1:k3d),
     &     vis(0:3*k1d-1,0:3*k2d-1,n2d_1),
     &     vis_tot(0:3*k1d-1,0:3*k2d-1),pvac(nmzxyd),
     &     pint(nmzxyd),vis_help(0:3*k1d-1,0:3*k2d-1),
     &     af2(0:9*k1d*k2d-1),bf2(0:9*k1d*k2d-1),vpw_help(n3d) )

      fpi = 4.*pimach()
      tpi = 2.*pimach()
      ivfft2d = 9*k1d*k2d
      ic = cmplx(0.,1.)
      ivfft1 = 3*k1d
      ivfft2 = 3*k2d
      ani1 = 1./real(ivfft1)
      ani2 = 1./real(ivfft2)
      DATA   tol_21 /1.0e-21/

c--------- initializations -------->

c----> vpw in the '1st aproximation' (V - tilde)

      vpw(1,1) = cmplx(0.,0.)

      DO 140 irec3 = 2,nq3
                  
         g = sk3(irec3)

         vpw(irec3,1) = fpi*psq(irec3)/(g*g)
         
 140  CONTINUE

      DO irc1 = 2,nq2_1
         DO i = 1,nmzxy
            vxy(i,irc1-1,1,1) = cmplx(0.,0.)
         END DO
      END DO

c----> values of the potential in the 1st approximation on the boundary
c----> if nstr2.ne.1 then it should be changed!!!

      DO m = 0,MM+1
         DO k3 = 1,k3d
            CALL modcyli(m,bmat(3,3)*k3*z1,iJJ(m,k3))
         END DO
         DO irec2 = 1,nq2
            g2 = sk2(irec2)
            IF (irec2.NE.0) THEN
               CALL od_cylbes(m,g2*z1,fJJ(m,irec2))
            END IF
         END DO
      END DO

      DO irc1 = 1,nq2_1
         val(irc1) = cmplx(0.,0.)
         m = kv2_1(2,irc1)
         IF (m.LT.0) THEN
            mult = real((-1)**m)
         ELSE
            mult = 1.
         END IF
         k3 = kv2_1(1,irc1)
         DO irec2 = 1,nq2
            phi = phi2(irec2)
            irec3 = ig(kv2(1,irec2),kv2(2,irec2),k3)
            IF (irec3.NE.0) THEN
               val(irc1) = val(irc1) +
     +              (ic**m)*vpw(irec3,1)*exp(-ic*
     *              m*phi)*fJJ(iabs(m),irec2)*
     *              nstr2(irec2)*mult
            END IF
         END DO
      END DO

c-----> preparing arrays for radial grids:selecting from
c-----> x & y only the radius

      nrd = 0

      DO ix = 0,ivfft1 - 1
         DO 7 iy = 0,ivfft2 - 1
            x = ix*ani1
            IF (x.GT.0.5) x = x - 1.
            y = iy*ani2
            IF (y.GT.0.5) y = y - 1.
            r = sqrt((x*amat(1,1) + y*amat(1,2))**2 +
     +               (x*amat(2,1) + y*amat(2,2))**2)
            DO i = 1,nrd
               IF (abs(r-rr(i)).LE.1.e-6) THEN
                  rmap(ix,iy) = i
                  GOTO 7
               END IF
            END DO
            nrd = nrd + 1
            rmap(ix,iy) = nrd
            rr(nrd) = r
 7       CONTINUE
      END DO

      DO gzi = -k3d,k3d
         DO m = -MM,MM
            IF (m.LT.0) THEN
               mult = real((-1)**m)
            ELSE
               mult = 1.
            END IF
            val_m(gzi,m) = cmplx(0.,0.)
            irc1 = ig1(gzi,m)
            IF (irc1.NE.0) THEN
               val_m(gzi,m) = val(irc1)
            ELSE
               DO irec2 = 1,nq2
                  phi = phi2(irec2)
                  irec3 = ig(kv2(1,irec2),kv2(2,irec2),gzi)
                  IF (irec3.NE.0) THEN
                     val_m(gzi,m) = val_m(gzi,m) +
     +                    (ic**m)*vpw(irec3,1)*exp(-ic*
     *                    m*phi)*fJJ(iabs(m),irec2)*
     *                    nstr2(irec2)*mult
                  END IF
               END DO
            END IF
         END DO
      END DO
        
c-------  During the following staff we miss the m=0,gz=0 -> irec1 = 1
c-------  component of the potential which should be also added
c-------  in order to get the total potential for the FFT's

      DO ix = 0,ivfft1 - 1
         DO iy = 0,ivfft2 - 1
            irc = rmap(ix,iy)
            r = rr(irc)
            IF (r.GT.z1) THEN
               zf = (r-z1)/delz + 1.0
               im = zf
               q = zf - im
               vis(ix,iy,1) = 0.5*(q-1.)*
     *              (q-2.)*vz(im,1,1) -
     +              q*(q-2.)*vz(im+1,1,1) +
     +              0.5*q*(q-1.)*vz(im+2,1,1)
            ELSE
               vis(ix,iy,1) = 
     +              vz(1,1,1) - val(1) + tpi*
     *              psq(1)*(z1*z1 - r*r)/2.
            END IF
            DO irc1 = 2,nq2_1
               vis(ix,iy,irc1) = cmplx(0.,0.)
            END DO
         END DO
      END DO

      DO m = 0,MM
         DO k3 = 1,k3d
            CALL modcyli(m,bmat(3,3)*k3*z1,IR(k3,m))
         END DO
      END DO

      DO i = 1,nrd
         DO m = 0,MM
            DO k3 = 1,k3d
               CALL modcyli(m,bmat(3,3)*k3*rr(i),IIII(i,k3,m))
            END DO
         END DO
      END DO

c------- cycle by positive gz---------->
      
      DO 100 gzi = 0,k3d                        ! gz

c------- cycle by positive m ---------->
 
      DO 110 m = 0,MM
        
c-------------------------------------->

      IF (m.NE.0 .OR. gzi.NE.0) THEN ! m^2 + gz^2.ne.0
            
         irec1(1) = ig1(gzi,m)
         irec1(2) = ig1(gzi,-m)
         irec1(3) = ig1(-gzi,m)
         irec1(4) = ig1(-gzi,-m)
      
      DO l = 1,3
         DO l1 = l+1,4
            IF (irec1(l).EQ.irec1(l1)) irec1(l1) = 0
         END DO
      END DO

c---> if all the irec1 not equal to zero
     
      s = 0
      
      DO l = 1,4
         s = s + irec1(l)*irec1(l)
      END DO   

c---> preparing special functions which depend only on the 
c---> absolute value of m and gz
      
      IF (s.NE.0 .AND. gzi.NE.0) THEN
         gz = bmat(3,3)*gzi
         DO  imz = 1,nmzxy
            z = z1 + delz*(imz-1)  
            CALL modcylk(m,gz*z,KK(imz))
            CALL modcyli(m,gz*z,II(imz))
         END DO
         IIIR = II(1)
         DO irc = 1,nrd
            III(irc) = IIII(irc,gzi,m)
         END DO
      ELSEIF (s.EQ.0) THEN
         GO TO 110
      ELSEIF (s.NE.0 .AND. gzi.EQ.0) THEN  
         gz = 0.
      END IF

c---> now we start the cycle by +-m,gz
      
      DO 777 l = 1,4

       IF (irec1(l).NE.0) THEN

c--------------------------------------->   

       DO ix = 0,ivfft1-1
          DO iy = 0,ivfft2-1
             vis_help(ix,iy) = cmplx(0.,0.)
          END DO
       END DO
       
       DO i = 1,nmzxy
          fact(i) = 0.
          rxy(i) = cmplx(0.,0.)
          pvac(i) = cmplx(0.,0.)
          pint(i) = cmplx(0.,0.)
       END DO
       
c-------------- gz = 0 ------------------------------------->
c----------------------------------------------------------->
      
      IF (gzi.EQ.0) THEN

c----- this form of the density is just more easy to use

         DO imz = 1,nmzxy
            rxy(imz) = rhtxy(imz,irec1(l)-1,1,1)
         END DO
         
c----- vacuum potential caused by the vacuum density

         CALL vacp5_0(
     >        nmzxyd,nmzxy,z1,tpi,rxy,m,delz,
     <        pvac,fact)


c----- vacuum potential caused by the interstitial density
         
         aa = cmplx(0.,0.)
         m1 = kv2_1(2,irec1(l))
         gzi1 = kv2_1(1,irec1(l))

         IF (m1.LT.0) THEN
            mult = real((-1)**m)
         ELSE
            mult = 1.
         END IF
               
         DO irec2 = 1,nq2
            irec3 = ig(kv2(1,irec2),kv2(2,irec2),0)
            IF (irec3.NE.0 .AND. irec2.NE.1) THEN
               phi = phi2(irec2)
               g2 = sk2(irec2)
               aa = aa +    
     +              (ic**m1)*psq(irec3)*
     *              exp(cmplx(0.,-m1*phi))*
     *              cmplx(mult*fJJ(m+1,irec2)/g2,0.)*nstr2(irec2)
            END IF
         END DO
               
c----- total vacuum potential
               
         DO 160 imz = 1,nmzxy  
                  
            pint(imz) =  fact(imz)*aa 
                  
            vxy(imz,irec1(l)-1,1,1) = pvac(imz) + pint(imz)
                  
 160     END DO
               
c----- array val further is a boundary values of the
c----- potential V- \tilde \tilde which is created to compensate 
c----- the influence of the outer 'noice' charge density - which 
c----- is just a periodical continuation of the interstitial charge
c----- density, V - \tilde and V - \tilde\tilde are then added in
c----- order to obtain the real interstitial potential           

         val_help = vxy(1,irec1(l)-1,1,1) - val(irec1(l))
               
c----- potential \tilde\tilde{V} is a solution of the Laplase equation
c----- in the interstitial with the boundary conditions val_0 and val_z
c----- further, it is generated on the uniform grid in the unit cell
c----- \tilde{\Omega}, in the space between the cylindrical 
c----- interstitial boundary and the squre boundaries it is put to
c----- the vacuum potential
          
         CALL visp5_0(
     >        nmzxyd,nmzxy,delz,m,ivfft1,ivfft2,l,
     >        rxy,ani1,ani2,z1,amat,
     >        pvac,pint,tpi,jspd,val_help,
     =        vis_help)

         DO ix = 0,ivfft1 - 1
            DO iy = 0,ivfft2 - 1
               vis(ix,iy,irec1(l)) = vis_help(ix,iy)
            END DO
         END DO

         
c------- gz.NEQ.0--------------------------------------------->
      ELSE   
c------------------------------------------------------------->
         
         DO  imz = 1,nmzxy
            rxy(imz) = rhtxy(imz,irec1(l)-1,1,1)
         END DO   
         
c----- vacuum potential caused by the vacuum density        

         CALL vacp5_z(
     >        nmzxyd,nmzxy,z1,delz,fpi,II,KK,rxy,m,
     <        pvac)
         
c----- vacuum potential caused by the intst. density

         a = cmplx(0.,0.)
         m1 = kv2_1(2,irec1(l))
         gzi1 = kv2_1(1,irec1(l))
         
         IF (m1.LT.0) THEN
            mult = real((-1)**m)
         ELSE
            mult = 1.
         END IF  
 
         DO irec2 = 1,nq2
            irec3 = ig(kv2(1,irec2),kv2(2,irec2),gzi1)
            IF (irec3.ne.0) THEN
               g = sk3(irec3)
               g2 = sk2(irec2)
               phi = phi2(irec2)
               b = z1*( gz*iJJ(m+1,gzi)*fJJ(m,irec2) + 
     +              g2*II(1)*fJJ(m+1,irec2))/(g*g)
               a = a +    
     +              (ic**m1)*mult*exp(-ic*m1*phi)
     *              *psq(irec3)*b*nstr2(irec2)
            END IF
         END DO

c----- total vacuum potential ---------------
            
         DO imz = 1,nmzxy
            pint(imz) = fpi*a*KK(imz) 
            vxy(imz,irec1(l)-1,1,1) =  pint(imz) + pvac(imz)  
         END DO
            
         val_help = vxy(1,irec1(l)-1,1,1) - val(irec1(l))

         CALL visp5_z(
     >        nmzxyd,nmzxy,delz,m,ivfft1,ivfft2,IIIR,
     >        rxy,ani1,ani2,z1,amat,pvac,pint,tpi,l,
     >        fpi,val_help,III,m1,gz,rmap,rr,
     <        vis_help)
         
         DO ix = 0,ivfft1 - 1
            DO iy = 0,ivfft2 - 1
               vis(ix,iy,irec1(l)) = vis_help(ix,iy)
            END DO
         END DO
         
c---  end of vacuum contibution to vis ---------------------->
         
      END IF                    ! gz.neq.0
      
c---- > finishing the cycle by +-m,gz
      
      END IF
      
 777  END DO
      
      END IF                    ! m^2+gz^2.neq.0
      
 110  CONTINUE                  ! m

 100  CONTINUE                  ! gz 

c*************************************************************
c-------> finding the Fourier coefficients of the potential
c-------  in the interstitial
c-------  the scheme (making the potential continuous) is as follows:
c-------  1. transforming the \tilde{V} from pw-representation
c-------     to the real grid, putting to zero outside the
c-------     cylindrical interstitial
c-------  2. Adding the \tilde\tilde{V} on the real grid 
c-------  3. Transforming back to the pw-representation 
c-------  Everything is done in \tilda {\Omega}

      gzmin = -k3d

      IF (zrfs.OR.invs) gzmin = 0
      
      DO k3 = gzmin,k3d          ! collect the Fourier components

         fxy0 = 0.

         rhti = 0.
         
         DO irec2 = 1,nq2 - 1
            fxy(irec2) = cmplx(0.,0.)
         END DO

         DO irec2 = 1,nq2
            irec3 = ig(kv2(1,irec2),kv2(2,irec2),k3)
            IF (irec3.NE.0) THEN
               IF (irec2.EQ.1) THEN
                  fxy0 = real(vpw(irec3,1))
                  rhti = aimag(vpw(irec3,1))
               ELSE
                  fxy(irec2-1) = vpw(irec3,1)
               END IF
            END IF
         END DO

         CALL dset(ivfft2d,0.0,af2,1)
         CALL dset(ivfft2d,0.0,bf2,1)

         CALL fft2d(
     >        k1d,k2d,n2d,
     <        af2(0),bf2,
     >        fxy0,rhti,fxy,
     >        1,nq2,kimax2,+1,
     >        igfft2,pgfft2,nstr2)
         

c--> although some harmonics of the vacuum potential are equal to zero
c--> the same harmonics of the interstitial potential \Tilde{V} are not,
c--> because this potential is determined from pw coefficients of the
c--> charge density, which has the periodicity of the rectangular 
c--> unit cell, hence, these harmonics should be extracted from the 
c--> final interstitial potential

         i = 0
         
         DO iy = 0,3*k2d-1
            DO ix = 0,3*k1d-1
               x = ix*ani1
               IF (x.GT.0.5) x = x - 1.
               y = iy*ani2
               IF (y.GT.0.5) y = y - 1.
               r = sqrt((x*amat(1,1) + y*amat(1,2))**2 +
     +                  (x*amat(2,1) + y*amat(2,2))**2)               
               phi = angle(x*amat(1,1) + y*amat(1,2),
     >                     x*amat(2,1) + y*amat(2,2))
               vis_tot(ix,iy) = cmplx(0.,0.)
               j = rmap(ix,iy)
               DO m = -MM,MM
                  irc1 = ig1(k3,m)
                  IF (irc1.NE.0) THEN
                     vis_tot(ix,iy) = vis_tot(ix,iy) +
     +                    vis(ix,iy,irc1)
                  ELSE
                     IF (k3.EQ.0) THEN
                        IF (r.LE.z1) THEN
                           vis_tot(ix,iy) = vis_tot(ix,iy) - 
     +                          val_m(k3,m)*
     *                          exp(cmplx(0.,m*phi))*
     *                          (r**(iabs(m)))/(z1**(iabs(m)))
                        END IF
                     ELSE
                      IF (r.LE.z1) THEN
                        vis_tot(ix,iy) = vis_tot(ix,iy) -
     +                    val_m(k3,m)*IIII(j,iabs(k3),iabs(m))*
     *                    exp(cmplx(0.,m*phi))/IR(iabs(k3),iabs(m))
                      END IF
                     END IF
                  END IF
               END DO
               
               IF (r.LE.z1) THEN
                  af2(i) = af2(i) + real(vis_tot(ix,iy))
                  bf2(i) = bf2(i) + aimag(vis_tot(ix,iy))
c$$$                 af2(i) = real(vis_tot(ix,iy))
c$$$                 bf2(i) = aimag(vis_tot(ix,iy))
               ELSE
                  af2(i) = real(vis_tot(ix,iy))
                  bf2(i) = aimag(vis_tot(ix,iy))
               END IF
               i = i+1
            END DO
         END DO     

         rhti = 0.
               
         CALL fft2d(
     >        k1d,k2d,n2d,
     >        af2(0),bf2,
     <        gxy0,rhti,gxy,
     >        1,nq2,kimax2,-1,
     >        igfft2,pgfft2,nstr2) 

         DO irec2 = 1,nq2
            irec3 = ig(kv2(1,irec2),kv2(2,irec2),k3)
            IF (irec3.NE.0) THEN
               IF (irec2.EQ.1) THEN
                  vpw_help(irec3) = cmplx(gxy0,rhti)
               ELSE
                  vpw_help(irec3) = gxy(irec2-1)
               END IF
            END IF
         END DO 
      
      END DO                    ! gz -> Vpw(.,.,gz)

      DO irec3 = 1,nq3
         vpw(irec3,1) = vpw_help(irec3)
c$$$         vpw(irec3,1) = vpw(irec3,1) + vpw_help(irec3)
      END DO


      DO irc1 = 2,nq2_1
        DO imz = 1,nmzxy
           IF (abs(vxy(imz,irc1-1,1,1)).LE.tol_21)
     &          vxy(imz,irc1-1,1,1) = cmplx(0.,0.)
        END DO
      END DO


      
      DEALLOCATE ( KK,II,III,IIII,fact,val,val_m,rxy,
     &     vis,vis_tot,pvac,pint,vis_help,rmap,rr,
     &     af2,bf2,vpw_help,fJJ,iJJ )

      RETURN
      END SUBROUTINE od_vvacis
      END MODULE m_od_vvacis












