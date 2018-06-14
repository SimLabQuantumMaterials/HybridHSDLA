      MODULE m_od_abvac
      CONTAINS
      SUBROUTINE od_abvac(
     >     z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,
     >     ig,ig1,tpi,qssbti,
     >     nmzxy,nmz,delz,ig2,n2d_1,
     >     bbmat,wronk,evac,bkpt,MM,vM,
     >     vz,kvac3,nv2,
     <     uz,duz,u,udz,dudz,ddnv,ud)
c**************************************************************
c      determines the nesessary values and derivatives on the 
c      vacuum cylindrical boundary for finding a and b coefficients
c      for the construcing vacuum charge density in vacden.F
c                          Y.Mokrousov, 7th of october 2002
c*************************************************************** 
      USE m_constants, ONLY : pimach
      USE m_dotir, ONLY : dotirp
      IMPLICIT NONE

c     .. scalar Arguments..
      real wronk
      integer, intent (in) :: nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n3d
      integer, intent (in) :: nmzxy,nmz,MM,n2d,vM
      integer, intent (in) :: n2d_1,nv2
      real,    intent (in) :: delz,z1,tpi,evac
c     ..array arguments..
      
      integer, intent (in) :: kvac3(nv2d)
      real,    intent (in) :: bkpt(3),qssbti 
      integer, intent (in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: ig1(-k3d:k3d,-MM:MM),ig2(n3d)
      real,    intent (in) :: vz(nmzd),bbmat(3,3)
      real,    intent (out):: udz(nv2d,-vM:vM)
      real,    intent (out):: uz(nv2d,-vM:vM)
      real,    intent (out):: dudz(nv2d,-vM:vM)
      real,    intent (out):: duz(nv2d,-vM:vM)
      real,    intent (out):: u(nmzd,nv2d,-vM:vM)
      real,    intent (out):: ud(nmzd,nv2d,-vM:vM)
      real,    intent (out):: ddnv(nv2d,-vM:vM)
c     ..local scalars..
      real ev,scale,xv,yv,vzero,v1
      integer i,ik,jk,jspin,jsp1,jsp2,m,l
      integer i1,i2,i3,ind1,ind3
c     .. local arrays..
      real wdz(nv2d,-vM:vM),wz(nv2d,-vM:vM)
      real dwdz(nv2d,-vM:vM),dwz(nv2d,-vM:vM)
      real v(3),x(nmzd)
      real vr0(nmzd)
      real w(nmzd,nv2d,-vM:vM),wd(nmzd,nv2d,-vM:vM)
c     ..external subroutines..
      external vacudz,vacuz
c     ..intrinsic functions..
      intrinsic aimag,cmplx,conjg,real,sqrt

c     wronksian for the schrodinger equation given by an identity

      wronk = 2.0

      do 20 ik = 1,nv2
         do 25 m = 0,vM
            v(1) = 0.0
            v(2) = 0.0
            v(3) = bkpt(3) + kvac3(ik) + qssbti
            ev = evac - 0.5*dotirp(v,v,bbmat)
            
c     constructing of the 'pseudopotential'
            
            do 55 i=1,nmz
               v1 = 1./(8.*((z1+(i-1)*delz)**2))
     -              -(m*m)/(2.*((z1+(i-1)*delz)**2))
               vr0(i) = vz(i)-v1
 55         continue
            
            vzero = vr0(nmz)
            
c     obtaining solutions with the 'pseudopotential'
            
            CALL vacuz(ev,vr0(1),vzero,nmz,delz,
     +           wz(ik,m),dwz(ik,m),w(1,ik,m))
            CALL vacudz(ev,vr0(1),vzero,nmz,delz,
     +           wdz(ik,m),dwdz(ik,m),ddnv(ik,m),
     +           wd(1,ik,m),dwz(ik,m),w(1,ik,m))
            
            scale = wronk/(wdz(ik,m)*dwz(ik,m)-
     -           dwdz(ik,m)*wz(ik,m))
            wdz(ik,m) = scale*wdz(ik,m)
            dwdz(ik,m) = scale*dwdz(ik,m)
            ddnv(ik,m) = scale*ddnv(ik,m)
            IF (m.GT.0) THEN
               wdz(ik,-m) = wdz(ik,m)
               dwdz(ik,-m) = dwdz(ik,m)
               ddnv(ik,-m) = ddnv(ik,m)
            END IF
            do 10 i = 1,nmz
               wd(i,ik,m) = scale*wd(i,ik,m)
               w(i,ik,m) = scale*w(i,ik,m)
               IF (m.GT.0) THEN
                  wd(i,ik,-m) = wd(i,ik,m)
                  w(i,ik,-m) = w(i,ik,m)
               END IF   
 10         continue
            
c     constructing 'real' solutions
            
            do 65 i=1,nmz
               u(i,ik,m)=w(i,ik,m)/sqrt(z1+(i-1)*delz)
               ud(i,ik,m)=wd(i,ik,m)/sqrt(z1+(i-1)*delz)
               IF (m.GT.0) THEN
                  u(i,ik,-m) = u(i,ik,m)
                  ud(i,ik,-m) = ud(i,ik,m)
               END IF
 65         continue
            duz(ik,m)=(-dwz(ik,m))/sqrt(z1)-
     -           wz(ik,m)/(2.0*((z1)**(1.5)))
            uz(ik,m)=wz(ik,m)/sqrt(z1)
            dudz(ik,m)=(-dwdz(ik,m))/sqrt(z1)-
     -           wdz(ik,m)/(2.0*((z1)**(1.5)))
            udz(ik,m)=wdz(ik,m)/sqrt(z1)
            IF (m.GT.0) THEN
               duz(ik,-m) = duz(ik,m)
               uz(ik,-m) = uz(ik,m)
               dudz(ik,-m) = dudz(ik,m)
               udz(ik,-m) = udz(ik,m)
            END IF

 25      continue 
 20   continue
      
      return
      END SUBROUTINE od_abvac
      END MODULE m_od_abvac
