      MODULE m_radflo
      CONTAINS
      SUBROUTINE radflo(
     >                  ntypd,nlod,jspd,jmtd,lmaxd,
     >                  ntyp,jsp,ello,vr,jri,r0,dx,
     >                  f,g,llo,nlo,l_dulo,irank,ulo_der,
     <                  ulos,dulos,uulon,dulon,uloulopn,
     <                  uuilon,duilon,ulouilopn,flo)
c
c***********************************************************************
c generates the scalar relativistic wavefunctions (flo) needed for the
c local orbitals at atom type n for angular momentum l.
c the values of the function and the radial derivative on the sphere
c (ulos,dulos) boundaries, the overlap with the other radial functions
c (uulon,dulon) and between the different local orbitals (uloulopn) are
c also calculated.
c
c ulos    : the value of the radial function of a local orbital
c           at the muffin tin radius
c dulos   : the value of the radial derivative of the radial function
c           function of a local orbital at the muffin tin radius
c uulon   : overlap integral between the radial functions of a local
c           obital and the flapw radial function with the same l
c dulon   : overlap integral between the radial functions of a local
c           obital and the energy derivative of the flapw radial
c           function with the same l
c uloulopn: overlap integral between the radial functions of two
c           different local orbitals
c           (only needed if they have the same l)
c l_dulo  : whether we use a dot u as local orbital (see setlomap.F)
c
c p.kurz jul. 1996 gb 2001
c
c ulo_der : specifies the order of the energy derivative to be used
c           (0, 1, 2, ... for primitive, first, second derivatives etc.)
c uuilon  : overlap integral between the radial functions of the
c           integral (multiplied by ulo_der) of a local orbital and the
c           flapw radial function with the same l
c duilon  : overlap integral between the radial functions of the
c           integral of a local orbital and the energy derivative of the
c           flapw radial function with the same l
c ulouilopn: overlap integral between the radial functions of the
c           integral of a local orbital and another local orbital with
c           the same l.
c
c C. Friedrich Feb. 2005
c***********************************************************************
c
      USE m_intgr, ONLY : intgr0
      USE m_constants, ONLY : c_light
      USE m_radsra
      USE m_radsrdn
      USE m_differ
#include "cpp_double.h"
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,nlod,jspd,jmtd,lmaxd
      INTEGER, INTENT (IN) :: jri,jsp,ntyp,irank
      REAL,    INTENT (IN) :: dx,r0
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),nlo(ntypd)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
      REAL,    INTENT (IN) :: ello(nlod,ntypd),vr(jmtd)
      REAL,    INTENT (IN) :: f(jmtd,2,0:lmaxd),g(jmtd,2,0:lmaxd)
      REAL,    INTENT (OUT):: ulos(nlod,ntypd),uulon(nlod,ntypd)
      REAL,    INTENT (OUT):: dulon(nlod,ntypd),dulos(nlod,ntypd)
      REAL,    INTENT (OUT):: uloulopn(nlod,nlod,ntypd)
      REAL,    INTENT (OUT):: uuilon(nlod,ntypd),duilon(nlod,ntypd)
      REAL,    INTENT (OUT):: ulouilopn(nlod,nlod,ntypd)
      REAL,    INTENT (OUT):: flo(jmtd,2,nlod)
      LOGICAL, INTENT (IN) :: l_dulo(nlod)
C     ..
C     .. Local Scalars ..
      INTEGER i,j,k,l,ilo,jlo,nodelo,noded,msh,ierr
      REAL c,uds,duds,ddn,dus,rn,t1,t2,d,rmt,rr,fn,fl,fj,e
      LOGICAL ofdiag
C     ..
C     .. Local Arrays ..
      REAL dulo(jmtd),ulo(jmtd),glo(jmtd,2),filo(jmtd,2,nlod)
      REAL, ALLOCATABLE :: f_rel(:,:),vrd(:)
      REAL help(nlod+2,nlod+2)
C     ..
      c = c_light(1.0)
c
      IF (irank.EQ.0) WRITE (6,FMT=8000)
      ofdiag = .false.
c---> calculate the radial wavefunction with the appropriate
c---> energy parameter (ello)
      DO ilo = 1,nlo(ntyp)
         CALL radsra(
     >             ello(ilo,ntyp),llo(ilo,ntyp),vr,r0,dx,jri,jmtd,c,
     <             ulos(ilo,ntyp),dulos(ilo,ntyp),nodelo,flo(1,1,ilo),
     <             flo(1,2,ilo))
c
c+apw+lo
         IF (l_dulo(ilo).or.ulo_der(ilo,ntyp).ge.1) THEN
           IF (ulo_der(ilo,ntyp).le.8) THEN
c--->    calculate orthogonal energy derivative at e
           i = ulo_der(ilo,ntyp)
           IF(l_dulo(ilo)) i=1
           CALL radsrdn(
     >              ello(ilo,ntyp),llo(ilo,ntyp),vr,r0,dx,jri,jmtd,c,
     <              uds,duds,ddn,noded,glo,filo(1,1,ilo),
     >              flo(1,1,ilo),dulos(ilo,ntyp),i)
           DO i=1,jri
             flo(i,1,ilo) = glo(i,1)
             flo(i,2,ilo) = glo(i,2)
           ENDDO
           nodelo = noded
           ddn    = sqrt(ddn)
           IF(l_dulo(ilo)) ddn=1.0
           flo (:,:,ilo)   = flo (:,:,ilo)/ddn ! Normalize ulo (flo) if APW+lo is not used
           filo(:,:,ilo)   = filo(:,:,ilo)/ddn ! and scale its integral (filo) accordingly
           dulos(ilo,ntyp) = duds/ddn          !   (setabc1lo and slomat assume <flo|flo>=1)
           ulos (ilo,ntyp) =  uds/ddn          ! 
           ELSE
!
!          test:
!
! set up core-mesh
            d = exp(dx)
            rmt = r0
            DO i = 1, jri - 1
             rmt = rmt * d
            ENDDO 
            rn = rmt
            msh = jri
            DO WHILE (rn < rmt + 20.0)
               msh = msh + 1
               rn = rn*d
            ENDDO
            rn = r0*( d**(msh-1) )
            ALLOCATE ( f_rel(msh,2),vrd(msh) )

! extend core potential (linear with slope t1 / a.u.)

            DO j = 1, jri
              vrd(j) = vr(j)
            ENDDO
            t1=0.125
            t2 = vrd(jri)/rmt - rmt*t1
            rr = rmt
            DO j = jri + 1, msh
               rr = d*rr
               vrd(j) = rr*( t2 + rr*t1 )
            ENDDO
            e = ello(ilo,ntyp)
            fn = 6.0 ; fl = 1.0 ; fj = 0.5
            CALL differ(
     >                  fn,fl,fj,c,82.0,dx,r0,
     >                  rn,d,msh,vrd,
     X                  e,
     <                  f_rel(1,1),f_rel(1,2),ierr)

            f_rel(:,1) = 2 * f_rel(:,1)
            f_rel(:,2) = 2 * f_rel(:,2)
            rn = r0
            DO i = 1, jri 
             rn = rn * d
             write(123,'(5f20.15)') rn,f_rel(i,1),f_rel(i,2),
     >                              f(i,1,1),f(i,2,1)
            ENDDO

            stop
            ENDIF
           
         ENDIF
c-apw+lo
c
c--->    calculate the overlap between these fcn. and the radial functions
c--->    of the flapw basis with the same l
         DO i = 1,jri
            ulo(i) = f(i,1,llo(ilo,ntyp))*flo(i,1,ilo) +
     +               f(i,2,llo(ilo,ntyp))*flo(i,2,ilo)
            dulo(i) = g(i,1,llo(ilo,ntyp))*flo(i,1,ilo) +
     +                g(i,2,llo(ilo,ntyp))*flo(i,2,ilo)
         END DO
         CALL intgr0(ulo, r0,dx,jri,uulon(ilo,ntyp))
         CALL intgr0(dulo,r0,dx,jri,dulon(ilo,ntyp))
         IF (l_dulo(ilo)) dulon(ilo,ntyp) = 0.0
         IF (irank.EQ.0) 
     +     WRITE (6,FMT=8010) ilo,llo(ilo,ntyp),ello(ilo,ntyp),
     +     ulos(ilo,ntyp),dulos(ilo,ntyp),nodelo,uulon(ilo,ntyp),
     +     dulon(ilo,ntyp)
c
c--->   case LO = energy derivative (ulo_der>=1):
c--->   calculate the overlap between the LO-integral (filo) and the radial functions
         IF(ulo_der(ilo,ntyp).ge.1) then
           DO i=1,jri
            ulo(i) = f(i,1,llo(ilo,ntyp))*filo(i,1,ilo) +
     +               f(i,2,llo(ilo,ntyp))*filo(i,2,ilo)
            dulo(i) = g(i,1,llo(ilo,ntyp))*filo(i,1,ilo) +
     +                g(i,2,llo(ilo,ntyp))*filo(i,2,ilo)
           ENDDO
           CALL intgr0(ulo, r0,dx,jri,uuilon(ilo,ntyp))
           CALL intgr0(dulo,r0,dx,jri,duilon(ilo,ntyp))
         ELSE
           uuilon(ilo,ntyp)=0
           duilon(ilo,ntyp)=0
         ENDIF
c--->   calculate overlap between radial fcn. of different local
c--->   orbitals (only if both have the same l)
c         uloulopn(ilo,ilo,ntyp) = 1.0
c         DO jlo = 1, (ilo-1)
         DO jlo = 1, ilo
            IF (llo(ilo,ntyp).EQ.llo(jlo,ntyp)) THEN
               DO i = 1,jri
                  ulo(i) = flo(i,1,ilo)*flo(i,1,jlo) +
     +                     flo(i,2,ilo)*flo(i,2,jlo)
               END DO
               CALL intgr0(ulo,r0,dx,jri,uloulopn(ilo,jlo,ntyp))
               uloulopn(jlo,ilo,ntyp) = uloulopn(ilo,jlo,ntyp)
               ofdiag = .true.
            ELSE
               uloulopn(ilo,jlo,ntyp) = 0.0
               uloulopn(jlo,ilo,ntyp) = 0.0
            END IF
         END DO
      END DO
c
c---> case: one of LOs = energy derivative (ulo_der>=1):
c---> calculate overlap between LOs and integrals of LOs
      DO ilo = 1,nlo(ntyp)
        DO jlo = 1,nlo(ntyp)
          IF(ulo_der(jlo,ntyp).ge.1.and.llo(ilo,ntyp).eq.llo(jlo,ntyp))
     &        THEN
            DO i = 1,jri
              ulo(i) = flo(i,1,ilo)*filo(i,1,jlo) +
     &                 flo(i,2,ilo)*filo(i,2,jlo)
            ENDDO
            CALL intgr0(ulo,r0,dx,jri,ulouilopn(ilo,jlo,ntyp))
          ELSE
            ulouilopn(ilo,jlo,ntyp)=0.0
          ENDIF
        ENDDO
      ENDDO
c 
      IF ( (ofdiag).AND.(irank.EQ.0) ) THEN
         WRITE (6,FMT=*)
         WRITE (6,FMT=*)
     +     'overlap matrix between different local orbitals'
         WRITE (6,FMT=8020) (i,i=1,nlo(ntyp))
         DO ilo = 1,nlo(ntyp)
            WRITE (6,FMT='(i3,20e15.6)') ilo,
     +        (uloulopn(ilo,jlo,ntyp),jlo=1,nlo(ntyp))
         END DO
      END IF
c
c
c     Diagonalize overlap matrix of normalized MT functions
c       to check linear dependencies
c       help is overlap matrix of all MT functions (flapw & LO)
      IF (.false.) THEN
      DO l=0,maxval(llo(1:nlo(ntyp),ntyp))
        IF(all(llo(1:nlo(ntyp),ntyp).ne.l)) cycle
        help(1,1)=1.0
        help(1,2)=0.0
        help(2,1)=0.0
        help(2,2)=1.0
        DO i=1,jri
          ulo(i)=g(i,1,l)**2+g(i,2,l)**2
        ENDDO
        CALL intgr0(ulo,r0,dx,jri,ddn)
        ddn=sqrt(ddn)
        j=2
        k=2
        DO ilo=1,nlo(ntyp)
          IF(llo(ilo,ntyp).eq.l) THEN
            j=j+1
            help(1,j)=uulon(ilo,ntyp)
            help(2,j)=dulon(ilo,ntyp)/ddn
            k=2
            DO jlo=1,ilo
              IF(llo(jlo,ntyp).eq.l) THEN
                k=k+1
                help(k,j)=uloulopn(ilo,jlo,ntyp)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        WRITE(6,'(/A,I2)')
     &    'Overlap matrix of normalized MT functions '//
     &    'and its eigenvalues for l =',l
        DO i=1,j
          WRITE(6,'(30F11.7)') (help(k,i),k=1,i)
        ENDDO
        CALL CPP_LAPACK_ssyev('N','U',j,help,nlod+2,ulo,dulo,
     &    jmtd,i)
        IF(i.ne.0) STOP 'radflo: ssyev failed.'
        WRITE(6,'(/30F11.7)') (ulo(i),i=1,j)
      ENDDO
      ENDIF ! .false.
c

 8000 FORMAT (/,t20,'radial function for local orbitals',/,t2,'lo',t6,
     +       'l',t11,'energy',t29,'value',t42,'derivative',t56,'nodes',
     +       t63,'ovlp with u',t78,'ovlp with udot')
 8010 FORMAT (i3,i3,f10.5,5x,1p,2e16.7,i5,2e16.7)
 8020 FORMAT (' lo',i9,19i15)
      RETURN
      END SUBROUTINE radflo
      END MODULE m_radflo
