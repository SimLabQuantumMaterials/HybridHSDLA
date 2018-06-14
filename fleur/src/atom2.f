      MODULE m_atom2
c     *************************************************************
c     fully relativistic atomic program based on the subroutines
c     differ, outint and inwint by d.d.koelling
c     erich wimmer     august 1981
c     modified for use in start-up generator.  m.w. april 1982
c     modified by adding stabilizing well. m. weinert may 1990
c     *************************************************************
      CONTAINS
      SUBROUTINE atom2(
     >                 msh,mshd,nstd,zatom,rnot1,dx,icorr,total,krla,
     >                 igrd,ndvgrd,idsprs,isprsv,sprsv,
     >                 jspins,jspd,bmu,qdel,jri,
     <                 rhoss,nst,lnum,eig,vbar)

      USE m_intgr, ONLY : intgr1,intgr0
      USE m_constants, ONLY : c_light, pimach
      USE m_xcall, ONLY : vxcall
      USE m_potl0
      USE m_setcor
      USE m_differ

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT (IN)  :: msh,mshd,nstd,jspins,jspd,jri
      INTEGER,INTENT (IN)  :: igrd,ndvgrd,idsprs,isprsv,icorr,krla
      REAL,   INTENT (IN)  :: dx,rnot1,zatom,sprsv,bmu,qdel
      REAL,   INTENT (OUT) :: rhoss(mshd,jspd),eig(nstd,jspd),vbar(jspd)
      INTEGER,INTENT (OUT) :: nst,lnum(nstd)
      LOGICAL,INTENT (IN)  :: total
C     ..
C     .. Local Scalars ..
      REAL c,d,delrv,dist,distol,e,fisr,fj,fl,fn,fpi,h,
     +     p,p1,pmax,pmin,r,r3,rn,rnot,z,zero,bmu_l
      INTEGER i,inr0,it,itmax,k,l,n,ispin,kk,ierr,msh_l
      LOGICAL conv,lastit
C     ..
C     .. Local Arrays ..
      REAL a(msh),b(msh),dens(msh),occ(nstd,jspd)
      REAL rad(msh),rev(nstd,jspd),ahelp(msh),ain(msh),
     +     rh(msh),vr(msh),f(0:3),
     +     vr1(msh,jspd),vr2(msh,jspd),vxc(mshd,jspd)
      INTEGER kappa(nstd),nprnc(nstd)
C     ..
C     .. External Subroutines ..
      EXTERNAL stpot1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,iabs,isign,log,max,min,sqrt
C     ..
C     .. Data statements ..
c---->     distol set from 1.0e-6 to 1.0e-3
      DATA zero,distol/0.0e0,1.0e-3/
C     ..
      c = c_light(1.0)
      fpi = 4.0 * pimach()
c
      WRITE (6,FMT=8000)
 8000 FORMAT (' subroutine atom2 entered')
      z = zatom
      n = msh
      rnot = rnot1
      itmax = 100
      pmin = 0.01e0
      pmax = 0.2e0
      h = dx
      d = exp(h)
      r = rnot
      DO 10 i = 1,n
         rad(i) = r
         r = r*d
   10 CONTINUE
      rn = rad(n)
      bmu_l = bmu
      CALL setcor(
     >            z,nstd,jspd,jspins,bmu_l,
     <            nst,kappa,nprnc,occ)
c
c--->   for electric field case (sigma.ne.0), add the extra charge
c--->   to the uppermost level; ignore the possible problem that
c--->   the occupations may not be between 0 and 2
      IF (jspins.EQ.1) THEN
        occ(nst,1) = occ(nst,1) + qdel
      ELSE
        occ(nst,1) = occ(nst,1) + qdel/2.
        occ(nst,jspins) = occ(nst,jspins) + qdel/2.
      ENDIF
c
      CALL stpot1(
     >           msh,n,z,rad,
     <           vr1)
      DO i = 1,n
         vr1(i,jspins) = vr1(i,1)
      ENDDO
c
c     start iterating
c
      lastit = .false.
      conv = .true.
      delrv = 0.100
      inr0 = log(5.0/rnot)/h + 1.5
      DO 180 it = 1,itmax
         DO ispin = 1,jspins
c
c---->     load potential
           DO i = 1,n
              vr(i) = vr1(i,ispin)
           ENDDO
c----> adding stabilizing well: m. weinert
           DO i = inr0,n
              vr(i) = vr(i) + rad(i)*delrv* (rad(i)-rad(inr0))
           ENDDO
c---->     note that vr contains r*v(r) in hartree units
           DO i = 1,n
             rhoss(i,ispin) = zero
           ENDDO
           DO 90 k = 1,nst
              fn = nprnc(k)
              fj = iabs(kappa(k)) - 0.5e0
              fl = fj + 0.5e0*isign(1,kappa(k))
              e = -2* (z/ (fn+fl))**2
              ierr = -1 
              msh_l = msh 
              DO WHILE (ierr.NE.0) 
                CALL differ(
     >                      fn,fl,fj,c,z,h,rnot,rn,d,msh_l,vr,
     X                      e,
     <                      a,b,ierr)
                msh_l = msh_l - 1
                IF (msh-msh_l > 100) STOP 'atom2'
              ENDDO
              DO i = msh_l+1, msh
                a(i) = a(msh_l) 
                b(i) = b(msh_l) 
              ENDDO  
              DO i = 1,n
                 rh(i) = occ(k,ispin)* (a(i)**2+b(i)**2)
              ENDDO
c+ldau
              IF (lastit) THEN                         ! calculate slater interals
                 l = int(fl)
!                 write(*,*) nprnc(k),l 
                 DO kk = 0, 2*l, 2                      ! F0 for s, F0 + F2 for p etc.
                   r = rnot
                   DO i = 1, n
                     ain(i) = a(i)**2 * r**(-kk-1)      ! prepare inner integrand
                     r = r * d
                   ENDDO
                   CALL intgr1(ain,rnot,h,n,           ! integrate
     <                         ahelp)
                   r = rnot
                   DO i = 1, n-1
                     ain(i) =  a(i)**2 * r**kk * (ahelp(n) - ahelp(i))
                     r = r * d
                   ENDDO
                   CALL intgr0(ain,rnot,h,n-1,           ! integrate 2nd r
     <                         f(kk/2))

                 ENDDO
!                 write(*,*) (27.21*2*f(kk),kk=0,l)
              ENDIF
c-ldau
              eig(k,ispin) = e
c---->       calculate <r>
              DO i = 1,n
                 a(i) = (a(i)**2+b(i)**2)*rad(i)
              ENDDO
              CALL intgr1(a,rnot,h,n,b)
              rev(k,ispin) = b(n)
              DO i = 1,n
                 rhoss(i,ispin) = rhoss(i,ispin) + rh(i)
              ENDDO
   90      ENDDO
         ENDDO
c
c     solve poisson's equation
c
         DO i = 1,n
            dens(i) = rhoss(i,1)
         ENDDO
         IF (jspins.EQ.2) THEN
           DO i = 1,n
              dens(i) = dens(i) + rhoss(i,jspins)
           ENDDO
         ENDIF
         CALL intgr1(dens,rnot,h,n,a)
         DO 110 i = 1,n
            rh(i) = dens(i)/rad(i)
  110    CONTINUE
         CALL intgr1(rh,rnot,h,n,b)
         fisr = b(n)
         DO 120 i = 1,n
            vr(i) = (a(i)+rad(i)* (fisr-b(i))-z)
  120    CONTINUE
c+ta
         DO ispin = 1, jspins
           DO i = 1,n
             rhoss(i,ispin) = rhoss(i,ispin) / (fpi*rad(i)**2)
           ENDDO
         ENDDO
         IF ((igrd.EQ.0).AND.(icorr.NE.-1)) THEN
c
           CALL  vxcall(6,icorr,krla,jspins,
     >                   mshd,msh,rhoss,
     <                   vxc)
c
         ELSEIF ((igrd.GT.0).OR.(icorr.EQ.-1)) THEN
c
           CALL potl0(
     >                mshd,jspd,jspins,icorr,n,dx,rad,rhoss,
     X                vxc,
     >                ndvgrd,idsprs,isprsv,sprsv)
c
         ENDIF
         DO ispin = 1, jspins
           DO i = 1,n
             vr2(i,ispin) = vr(i) + vxc(i,ispin)*rad(i)
           ENDDO
         ENDDO
c-ta
c        determine distance of potentials
c
         r3 = rn**3
         dist = 0.0
         DO ispin = 1, jspins
           DO i = 1,n
              a(i) = (vr2(i,ispin)-vr1(i,ispin))**2
           ENDDO
           CALL intgr1(a,rnot,h,n,b)
           dist = dist + sqrt((3.0e0/r3)*b(n))
         ENDDO
         IF (lastit) GO TO 190
         IF (dist.LT.distol) lastit = .true.
c     mix new input potential
         p1 = 1.0e0/dist
         p = min(pmax,p1)
         p = max(p,pmin)
         WRITE (6,FMT=8060) it,dist,p
         p1 = 1.0e0 - p
         DO ispin = 1, jspins
           DO i = 1,n
             vr1(i,ispin) = p1*vr1(i,ispin) + p*vr2(i,ispin)
           ENDDO
         ENDDO
  180 CONTINUE
c
c output 
c
      WRITE (6,FMT=8030) dist
      conv = .false.
c     list eigenvalues
  190 IF (conv) WRITE (6,FMT=8040) it,dist
      DO ispin = 1,jspins
        WRITE (6,'(a8,i2)') 'spin No.',ispin 
        DO k = 1,nst
           fj = iabs(kappa(k)) - 0.5e0
           l = fj + 0.5e0*isign(1,kappa(k)) + 0.01e0
           lnum(k) = l
           WRITE (6,FMT=8050) nprnc(k),kappa(k),l,fj,
     +                        occ(k,ispin),eig(k,ispin),rev(k,ispin)
        ENDDO  
!
!--->   guess enpara if it doesn't exist, using floating energy parameters
!
        i = jri - (log(4.0)/dx+1.51)
        vbar(ispin) = vr1(i,ispin)/( rnot*exp(dx*(i-1)) )
        WRITE ( 6,'(/,'' reference energy = '',2f12.6,/)') vbar(ispin)
      ENDDO

 8030 FORMAT (/,/,/,' $$$ error: not converged, dist=',f10.6,/)
 8040 FORMAT (/,/,3x,'converged in',i4,' iterations to a distance of',
     +       e12.5,' har',/,/,3x,'n  kappa  l    j  ',5x,
     +       'occ.   eigenvalue (har)  <r>  ',/)
 8050 FORMAT (3x,i1,i5,i5,f6.1,2 (3x,f7.2,1x,2f12.6))
 8060 FORMAT ('it,dist,p=',i4,2f12.5)

      END SUBROUTINE atom2
      END MODULE m_atom2
