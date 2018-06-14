      SUBROUTINE coredr(jspins,ntype,neq,zatom,rmsh,dx,jri,seig,film,
     =                  ncst,rmt,rho,jmtd,jspd,msh,nlhd,nstd,ntypd,
     >                  vrs,
     <                  qints,rhc)
c     *******************************************************
c     *****   set up the core densities for compounds   *****
c     *****   for relativistic core                     *****
c     *******************************************************

      USE m_etabinit
      USE m_spratm
      USE m_ccdnup

      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      REAL seig
      INTEGER jspins,ntype
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jmtd,jspd,msh,nlhd,nstd,ntypd
      INTEGER, INTENT (IN) :: ncst(ntypd),jri(ntypd),neq(ntypd)
      REAL   , INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL   , INTENT (IN) :: rmt(ntypd),zatom(ntypd)
      REAL   , INTENT (IN) :: vrs(jmtd,ntypd,jspd)
      REAL,    INTENT (INOUT) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (OUT) :: rhc(msh,ntypd,jspd),qints(ntypd,jspd)
      LOGICAL, INTENT (IN) ::  film
C     ..
C     .. Local Scalars ..
      REAL dxx,rnot,sume,t2,t2b,z,t1,rr,d,v1,v2
      INTEGER i,j,jatom,jspin,k,n,ncmsh
      LOGICAL exetab
C     ..
C     .. Local Arrays ..
      REAL br(jmtd,ntypd),brd(msh),etab(100,ntypd),
     +     rhcs(jmtd,ntypd,jspd),rhochr(msh),rhospn(msh),
     +     tecs(ntypd,jspd),vr(jmtd,ntypd),vrd(msh)
      INTEGER nkmust(ntypd),ntab(100,ntype),ltab(100,ntype)

C     ..
      ntab(:,:) = -1 ; ltab(:,:) = -1 ; etab(:,:) = 0.0
c
c setup potential and field
c
      IF (jspins.EQ.1) THEN
         DO n = 1,ntype
            DO j = 1,jmtd
               vr(j,n) = vrs(j,n,1)
               br(j,n) = 0.0
            END DO
         END DO
      ELSE
         DO n = 1,ntype
            DO j = 1,jmtd
               vr(j,n) = (vrs(j,n,1)+vrs(j,n,jspins))/2.
               br(j,n) = (vrs(j,n,jspins)-vrs(j,n,1))/2.
            END DO
         END DO
      END IF
c
c setup eigenvalues
      exetab = .false.
      INQUIRE (file='core.dat',exist=exetab)
      IF (exetab) THEN
         OPEN (58,file='core.dat',form='formatted',status='old')
         REWIND 58
         DO n = 1,ntype
            READ (58,FMT=*) nkmust(n)
            DO k = 1,nkmust(n)
               READ (58,FMT='(f12.6,2i3)') etab(k,n),ntab(k,n),
     +                                               ltab(k,n)

            END DO
         END DO
      ELSE
         OPEN (58,file='core.dat',form='formatted',status='new')
         CALL etabinit(
     >                 jmtd,msh,ntypd,nstd,ntype,ncst,
     >                 zatom,rmt,dx,rmsh,jri,vr,
     <                 etab,ntab,ltab,nkmust)
      END IF
c
      ncmsh = msh
      seig = 0.
c ---> set up densities
      DO jatom = 1,ntype
c
         DO j = 1,jri(jatom)
            vrd(j) = vr(j,jatom)
            brd(j) = br(j,jatom)
         END DO

#ifdef CPP_CORE
c--->    linear extension of the potential with slope t1 / a.u.
         rr = rmt(jatom)
         d = exp(dx(jatom))
         t1=0.125
!         t2  = vrd(jri(jatom))/rr - rr*t1
!         t2b = brd(jri(jatom))/rr - rr*t1
         t2  = vrs(jri(jatom),jatom,1)     /rr - rr*t1
         t2b = vrs(jri(jatom),jatom,jspins)/rr - rr*t1
#else
         t2 = vrd(jri(jatom))/ (jri(jatom)-ncmsh)
         t2b = brd(jri(jatom))/ (jri(jatom)-ncmsh)
#endif
         IF (jri(jatom).LT.ncmsh) THEN
            DO i = jri(jatom) + 1,ncmsh
#ifdef CPP_CORE
               rr = d*rr
               v1 = rr*( t2  + rr*t1 )
               v2 = rr*( t2b + rr*t1 )
               vrd(i) = 0.5*(v2 + v1)
               brd(i) = 0.5*(v2 - v1)
#else
               vrd(i) = vrd(jri(jatom)) + t2* (i-jri(jatom))
               brd(i) = brd(jri(jatom)) + t2b* (i-jri(jatom))
#endif
            END DO
         END IF

!        rr = rmsh(1,jatom)
!        do i =1, ncmsh
!          rr = d*rr
!         write(*,'(3f20.10)') rr,vrd(i),brd(i)
!        enddo
        
c
         rnot = rmsh(1,jatom)
         z = zatom(jatom)
         dxx = dx(jatom)

         CALL spratm(
     >               msh,vrd,brd,z,rnot,dxx,ncmsh,
     >               etab(1,jatom),ntab(1,jatom),ltab(1,jatom),
     <               sume,rhochr,rhospn)

         seig = seig + neq(jatom)*sume
c
c     rho_up=2(ir) = (rhochr(ir)  + rhospn(ir))*0.5
c     rho_dw=1(ir) = (rhochr(ir)  - rhospn(ir))*0.5
c
         IF (jspins.EQ.2) THEN
            DO j = 1,jri(jatom)
               rhcs(j,jatom,jspins) = (rhochr(j)+rhospn(j))*0.5
               rhcs(j,jatom,1) = (rhochr(j)-rhospn(j))*0.5
            END DO
         ELSE
            DO j = 1,jri(jatom)
               rhcs(j,jatom,1) = rhochr(j)
            END DO
         END IF
         IF (jspins.EQ.2) THEN
           DO j = 1,msh
             rhc(j,jatom,jspins) = (rhochr(j)+rhospn(j))*0.5
             rhc(j,jatom,1) = (rhochr(j)-rhospn(j))*0.5
           ENDDO
         ELSE
            DO j = 1,msh
             rhc(j,jatom,1) = rhochr(j)
            END DO
         END IF
c
c store atomic eigenvalues to file.58
         IF (jatom.EQ.1) REWIND 58
         WRITE (58,FMT=*) nkmust(jatom)
         DO k = 1,nkmust(jatom)
            WRITE (58,FMT='(f12.6,2i3)') etab(k,jatom),ntab(k,jatom),
     +                                                 ltab(k,jatom)
         END DO
c---->update spherical charge density rho with the core density.
         CALL ccdnup(
     >               jmtd,jspd,msh,nlhd,ntypd,jspins,jatom,
     >               jri,neq,dx,rho,rmsh,rmt,
     >               sume,vrs,rhochr,rhospn,
     =               tecs,qints)

      END DO ! loop over atoms (jatom)
c
c---->store core charge densities to file.17
      OPEN (17,file='cdnc',form='unformatted',status='unknown')
      REWIND 17
      DO jspin = 1,jspins
         DO jatom = 1,ntype
            WRITE (17) (rhcs(j,jatom,jspin),j=1,jri(jatom))
            WRITE (17) tecs(jatom,jspin)
         END DO
         WRITE (17) (qints(jatom,jspin),jatom=1,ntype)
      END DO
      CLOSE (17)
c
      RETURN
      END
