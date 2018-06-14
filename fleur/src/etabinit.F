      MODULE m_etabinit
c     *******************************************************
c     *****   set up etab via old core program          *****
c     *******************************************************
!     modified to run with core-levels as provided by setcor
!     ntab & ltab transport this info to core.F        gb`02
!------------------------------------------------------------
      CONTAINS
      SUBROUTINE etabinit(
     >                    jmtd,msh,ntypd,nstd,ntype,ncst,
     >                    zatom,rmt,dx,rmsh,jri,vr,
     <                    etab,ntab,ltab,nkmust)

      USE m_constants, ONLY : c_light
      USE m_setcor
      USE m_differ
      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jmtd,msh,ntypd,nstd,ntype
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),ncst(ntypd)
      REAL   , INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd),rmt(ntypd)
      REAL   , INTENT (IN) :: vr(jmtd,ntypd),zatom(ntypd)
      REAL   , INTENT (OUT):: etab(100,ntypd)
      INTEGER, INTENT (OUT):: ntab(100,ntype),ltab(100,ntype)
      INTEGER, INTENT (OUT):: nkmust(ntypd)
C     ..
C     .. Local Scalars ..
      REAL c,d,dxx,e,fj,fl,fn,rn,rnot,t2,weight,z,t1,rr,bmu
      INTEGER i,ic,iksh,ilshell,j,jatom,korb,l,
     +        nst,ncmsh,nsh,nshell,ipos,ierr
C     ..
C     .. Local Arrays ..
      INTEGER kappa(nstd),nprnc(nstd)
      REAL eig(nstd),occ(nstd),vrd(msh),a(msh),b(msh)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,iabs,isign,max
C
      c = c_light(1.0)
c
      WRITE (6,FMT=8020)
c
      ncmsh = msh
c     ---> set up densities
      DO 40 jatom = 1,ntype
         z = zatom(jatom)
         rn = rmt(jatom)
         dxx = dx(jatom)
         bmu = 0.0
         CALL setcor(
     >               z,nstd,1,1,bmu,
     <               nst,kappa,nprnc,occ)
         rnot = rmsh(1,jatom)
         d = exp(dx(jatom))
         rn = rnot* (d** (ncmsh-1))
         WRITE (6,FMT=8000) z,rnot,dxx,jri(jatom)
         WRITE (16,FMT=8000) z,rnot,dxx,jri(jatom)
         DO j = 1,jri(jatom)
            vrd(j) = vr(j,jatom)
         ENDDO
#ifdef CPP_CORE
c--->    linear extension of the potential with slope t1 / a.u.
         t1=0.125
         t2  = vrd(jri(jatom))/rmt(jatom)-rmt(jatom)*t1
         rr = rmt(jatom)
         d = exp(dx(jatom))
#else
         t2 = vrd(jri(jatom))/ (jri(jatom)-msh)
#endif
         IF (jri(jatom).LT.msh) THEN
            DO i = jri(jatom) + 1,msh
#ifdef CPP_CORE
               rr = d*rr
               vrd(i) = rr*( t2 + rr*t1 )
#else
               vrd(i) = vrd(jri(jatom)) + t2* (i-jri(jatom))
#endif
            ENDDO
         END IF

         nst = ncst(jatom)
         DO 30 korb = 1,nst
            fn = nprnc(korb)
            fj = iabs(kappa(korb)) - .5e0
            weight = 2*fj + 1.e0
            fl = fj + (.5e0)*isign(1,kappa(korb))
            e = -2* (z/ (fn+fl))**2
            CALL differ(
     >                  fn,fl,fj,c,z,dxx,rnot,rn,d,msh,vrd,
     X                  e,
     <                  a,b,ierr)
            IF (ierr.NE.0) STOP 'etabinit: error in core-levels'
            WRITE (6,FMT=8010) fn,fl,fj,e,weight
            WRITE (16,FMT=8010) fn,fl,fj,e,weight
            eig(korb) = e
   30    CONTINUE

         ic = 0
         DO korb = 1,nst
            fn = nprnc(korb)
            fj = iabs(kappa(korb)) - .5e0
            weight = 2*fj + 1.e0
            fl = fj + (.5e0)*isign(1,kappa(korb))
            DO i = 1, int(weight)
              ic = ic + 1
              IF (kappa(korb).GT.0) THEN
                ipos = ic + 1 + i 
              ELSEIF (kappa(korb).LT.-1) THEN
                ipos = ic - 2*(iabs(kappa(korb))-1) + max(i-2,0)
              ELSE
                ipos = ic
              ENDIF
              etab(ipos,jatom) = eig(korb)
              ntab(ipos,jatom) = nint(fn)
              ltab(ipos,jatom) = nint(fl)
            ENDDO
         ENDDO
         nkmust(jatom) = ic

         DO i=1,nkmust(jatom)
            write(6,'(f12.6,2i3)') etab(i,jatom),ntab(i,jatom),
     +                                           ltab(i,jatom)
         ENDDO

   40 CONTINUE
 8000 FORMAT (/,/,10x,'z=',f4.0,5x,'r(1)=',e14.6,5x,'dx=',f8.6,5x,
     +       'm.t.index=',i4,/,15x,'n',4x,'l',5x,'j',4x,'energy',7x,
     +       'weight')
 8010 FORMAT (12x,2f5.0,f6.1,f10.4,f10.0)
 8020 FORMAT (/,/,12x,'core e.v. initialization')

      END SUBROUTINE etabinit
      END MODULE m_etabinit
