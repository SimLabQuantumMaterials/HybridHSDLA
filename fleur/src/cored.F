      SUBROUTINE cored(
     >                 jspins,jspin,ntype,neq,zatom,rmsh,dx,jri,film,
     >                 frcor,ncst,rmt,rho,jmtd,jspd,msh,
     >                 nlhd,nstd,ntypd,
     >                 vr,gw,
     <                 qint,rhc,seig)

c     *******************************************************
c     *****   set up the core densities for compounds.  *****
c     *****                      d.d.koelling           *****
c     *******************************************************

      USE m_intgr, ONLY : intgr3,intgr0,intgr1
      USE m_constants, ONLY : c_light
      USE m_setcor
      USE m_differ

      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jmtd,jspd,msh,nlhd,nstd,ntypd
      INTEGER, INTENT (IN) :: jspin,jspins,ntype,gw
      REAL,    INTENT (OUT) :: seig
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: ncst(ntypd),jri(ntypd),neq(ntypd)
      REAL   , INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL   , INTENT (IN) :: rmt(ntypd),zatom(ntypd)
      REAL   , INTENT (IN) :: vr(jmtd,ntypd)
      REAL,    INTENT (INOUT) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (OUT) :: rhc(msh,ntypd),qint(ntypd,jspd)
      LOGICAL, INTENT (IN) ::  film,frcor
C     ..
C     .. Local Scalars ..
      REAL e,fj,fl,fn,q,rad,rhos,rhs,sea,sume,t2,tec,weight
      REAL sfp,c,d,dxx,rn,rnot,z,t1,rr,r,lambda,bmu
      INTEGER i,j,jatom,jm,korb,n,ncmsh,nm,nm1,nst,kk,l,ierr
C     ..
C     .. Local Arrays ..
      REAL rhcs(msh),rhoc(msh),rhoss(msh),vrd(msh),f(0:3)
      REAL occ(nstd),a(msh),b(msh),ain(msh),ahelp(msh)
      REAL occ_h(nstd,2)
      INTEGER kappa(nstd),nprnc(nstd)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,iabs,isign
C     ..
      c = c_light(1.0)
      sfp = 4. * sqrt( atan (1.) )
      seig = 0.
      IF (jspin.EQ.1) THEN
         OPEN (17,file='cdnc',form='unformatted',status='unknown')
      ENDIF
c
      IF (frcor) THEN
         IF (jspin.EQ.1) REWIND 17
         DO 50 n = 1,ntype
            jm = jri(n)
            rnot = rmsh(1,n) ; dxx = dx(n)
            ncmsh = nint( log( (rmt(n)+10.0)/rnot ) / dxx + 1 )
            ncmsh = min( ncmsh, msh )
c     --->    read in core density
            READ (17) (rhc(i,n),i=1,ncmsh)
c     --->    update spherical charge density
            DO 30 i = 1,jri(n)
               rhoc(i) = rhc(i,n)
               rho(i,0,n,jspin) = rho(i,0,n,jspin) + rhoc(i)/sfp
   30       CONTINUE
c     --->    read in kinetic enrgy of the core
            READ (17) tec
c     ---> for total energy calculations, determine the sum of the
c     ---> eigenvalues by requiring that the core kinetic energy
c     ---> remains constant.
            DO 40 i = 1,jri(n)
               rhoc(i) = rhoc(i)*vr(i,n)/rmsh(i,n)
   40       CONTINUE
            nm = jri(n)
            CALL intgr3(rhoc,rmsh(1,n),dx(n),nm,rhos)
            sea = tec + rhos
            WRITE (16,FMT=8030) n,jspin,tec,sea
            WRITE (6,FMT=8030) n,jspin,tec,sea
            seig = seig + neq(n)*sea
   50    CONTINUE
c     --->    read in qint
         READ (17) (qint(n,jspin),n=1,ntype)
         RETURN
      END IF

!#ifdef CPP_CORE
!      IF (jspin.EQ.1) THEN
!        OPEN (45,file='slaterf',form='formatted',status='unknown')
!      ENDIF
!#endif
c     ---> set up densities
      DO 150 jatom = 1,ntype
         sume = 0.
         z = zatom(jatom)
c         rn = rmt(jatom)
         dxx = dx(jatom)
         bmu = 0.0
         CALL setcor(
     >               z,nstd,jspd,jspins,bmu,
     <               nst,kappa,nprnc,occ_h)
         IF ((bmu > 99.)) THEN
              occ(1:nst) = jspins *  occ_h(1:nst,jspin)
         ELSE
              occ(1:nst) = occ_h(1:nst,1) 
         ENDIF
         rnot = rmsh(1,jatom)
         d = exp(dx(jatom))
         ncmsh = nint( log( (rmt(jatom)+10.0)/rnot ) / dxx + 1 )
         ncmsh = min( ncmsh, msh )
         rn = rnot* (d** (ncmsh-1))
         WRITE (6,FMT=8000) z,rnot,dxx,jri(jatom)
         WRITE (16,FMT=8000) z,rnot,dxx,jri(jatom)
         DO 60 j = 1,jri(jatom)
            rhoss(j) = 0.
            vrd(j) = vr(j,jatom)
   60    CONTINUE
c
#ifdef CPP_CORE
c--->    linear extension of the potential with slope t1 / a.u.
         t1=0.125
         t1 = max( (vrd(jri(jatom)) - vrd(jri(jatom)-1)*d)*
     +              d / (rmt(jatom)**2 * (d-1) ) , t1)
         t2=vrd(jri(jatom))/rmt(jatom)-rmt(jatom)*t1
         rr = rmt(jatom)
#else
         t2 = vrd(jri(jatom)) / ( jri(jatom) - ncmsh )
#endif
         IF ( jri(jatom) .LT. ncmsh) THEN
            DO 70 i = jri(jatom) + 1,ncmsh
               rhoss(i) = 0.
#ifdef CPP_CORE
               rr = d*rr
               vrd(i) = rr*( t2 + rr*t1 )
c               vrd(i) = 2*vrd(jri(jatom)) - rr*( t2 + rr*t1 )
#else
               vrd(i) = vrd(jri(jatom)) + t2* (i-jri(jatom))
#endif
c
   70       CONTINUE
         END IF

!#ifndef CPP_CORE
         nst = ncst(jatom)        ! for lda+U
!#endif
         IF (gw.eq.1) WRITE(15) nst,rmsh(1:jri(jatom),jatom)
         DO 90 korb = 1,nst
!#ifndef CPP_CORE
            IF (occ(korb).EQ.0) GOTO 90
!#endif
            fn = nprnc(korb)
            fj = iabs(kappa(korb)) - .5e0
            weight = 2*fj + 1.e0
            IF (bmu > 99.) weight = occ(korb)
            fl = fj + (.5e0)*isign(1,kappa(korb))
            e = -2* (z/ (fn+fl))**2
            CALL differ(
     >                  fn,fl,fj,c,z,dxx,rnot,rn,d,ncmsh,vrd,
     X                  e,
     <                  a,b,ierr)
            WRITE (6,FMT=8010) fn,fl,fj,e,weight
            WRITE (16,FMT=8010) fn,fl,fj,e,weight
            IF (ierr.NE.0) STOP 'error in core-level routine'
            IF (gw.eq.1) WRITE (15) NINT(fl),weight,e,a(1:jri(jatom)),
     &                                                b(1:jri(jatom))
!_ldau
!#ifdef CPP_CORE
!            IF (korb.GT.ncst(jatom)) THEN
!               l = int(fl)
!               DO kk = 0, 2*l, 2                      ! F0 for s, F0 + F2 for p etc.
!                 lambda = 3.5449*sqrt((l+1.0)/100)     ! screening (TF) sqrt(4pi N(ef))
!                 IF (kk.GT.0) lambda = 2*lambda
!                 r = rnot
!                 DO i = 1, ncmsh
!                   ain(i) = a(i)**2 * r**(-kk-1)      ! prepare inner integrand
!!                   IF (kk.EQ.0) THEN
!                     ain(i) = ain(i) * exp(-r*lambda)
!!                   ENDIF
!                   r = r * d
!                 ENDDO
!                 CALL intgr1(ain,rnot,dxx,ncmsh,           ! integrate
!     <                                  ahelp)
!                 r = rnot
!                 DO i = 1, ncmsh-1
!                   ain(i) = a(i)**2 * r**kk * (ahelp(ncmsh) - ahelp(i))
!!                   IF (kk.EQ.0) THEN
!                     ain(i) = ain(i) * exp(r*lambda)
!!                   ENDIF
!                   r = r * d
!                 ENDDO
!                 CALL intgr0(ain,rnot,dxx,ncmsh-1,           ! integrate 2nd r
!     <                                   f(kk/2))
!               ENDDO
!!               write(*,*) (27.21*2*f(kk),kk=0,l)
!!               WRITE(45,'(i3,4f20.10)') l,(27.21*2*f(kk),kk=0,l)
!            ELSE 
!            IF (occ(korb).EQ.0) GOTO 90
!#endif
!-ldau
            sume = sume + weight*e/jspins
            DO j = 1,ncmsh
               rhcs(j) = weight* (a(j)**2+b(j)**2)
               rhoss(j) = rhoss(j) + rhcs(j)
            ENDDO
!#ifdef CPP_CORE
!            ENDIF
!#endif
   90    CONTINUE
          
c     ---->update spherical charge density rho with the core density.
c     ---->for spin-polarized (jspins=2), take only half the density
         nm = jri(jatom)
         DO 100 j = 1,nm
            rhoc(j) = rhoss(j)/jspins
            rho(j,0,jatom,jspin) = rho(j,0,jatom,jspin) + rhoc(j)/sfp
  100    CONTINUE

         rhc(1:ncmsh,jatom)   = rhoss(1:ncmsh) / jspins
         rhc(ncmsh+1:msh,jatom) = 0.0

         seig = seig + neq(jatom)*sume
!         WRITE (17) (rhoc(i),i=1,nm)
         WRITE (17) (rhc(i,jatom),i=1,ncmsh)
         DO 120 i = 1,nm
            rhoc(i) = rhoc(i)*vr(i,jatom)/rmsh(i,jatom)
  120    CONTINUE
         CALL intgr3(rhoc,rmsh(1,jatom),dx(jatom),nm,rhs)
         tec = sume - rhs
         WRITE (6,FMT=8030) jatom,jspin,tec,sume
         WRITE (16,FMT=8030) jatom,jspin,tec,sume
         WRITE (17) tec
  130    CONTINUE
c     ---> simpson integration
         rad = rmt(jatom)
         q = rad*rhoss(nm)/2.
         DO 140 nm1 = nm + 1,ncmsh - 1,2
            rad = d*rad
            q = q + 2*rad*rhoss(nm1)
            rad = d*rad
            q = q + rad*rhoss(nm1+1)
  140    CONTINUE
         q = 2*q*dxx/3
c+sb
         WRITE (6,FMT=8020) q/jspins
         WRITE (16,FMT=8020) q/jspins
c-sb
         qint(jatom,jspin) = q*neq(jatom)
c   ---neu---
c      qsave(jatom)=q
c      nm=jri(jatom)
c      call intgr1(rhoss,rnot(jatom),dx(jatom),nm,rhcs)
c      corr=qsave(jatom)/rhcs(nm)
c      do 1070 j=1,nm
c      rho(j,jatom)=rho(j,jatom)+rhoss(j)*corr
c 1070 continue
c      write(16,'(a,i2,a,f10.5)') 'core charge of atom', jatom
c     + 'multiplied with',corr
  150 CONTINUE

#ifdef CPP_CORE
      IF (jspin.EQ.jspins) THEN
        CLOSE (45)
      ENDIF
#endif

c      qint=0.
      WRITE (17) (qint(n,jspin),n=1,ntype)
c
      IF (jspin.EQ.jspins) CLOSE (17)
      RETURN

 8000 FORMAT (/,/,10x,'z=',f4.0,5x,'r(1)=',e14.6,5x,'dx=',f8.6,5x,
     +       'm.t.index=',i4,/,15x,'n',4x,'l',5x,'j',4x,'energy',7x,
     +       'weight')
 8010 FORMAT (12x,2f5.0,f6.1,f10.4,f10.0)
 8020 FORMAT (f20.8,'  electrons lost from core.')
 8030 FORMAT (10x,'atom type',i3,'  (spin',i2,') ',/,10x,
     +       'kinetic energy=',e20.12,5x,'sum of the eigenvalues=',
     +       e20.12)
      END
