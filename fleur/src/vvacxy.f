      MODULE m_vvacxy
c     **********************************************************
c     g.ne.0 coefficient of vacuum coulomb potential           *
c     due to warped vacuum charged density                     *
c                                 c.l.fu, r.podloucky          *
c     **********************************************************
!     modified for thick films to avoid underflows gb`06
!---------------------------------------------------------------
      CONTAINS
      SUBROUTINE vvacxy(
     >                  n2d,jspd,nmzxyd,nmzxy,ng2,nvac,z1,delz,invs,sk2,
     >                  rhtxy,
     X                  vxy,
     <                  alphm)

      USE m_intgr, ONLY : intgz1
      USE m_constants, ONLY : pimach
      USE m_qsf
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n2d,jspd,nmzxyd
      INTEGER, INTENT (IN) :: nvac,nmzxy,ng2
      REAL,    INTENT (IN) :: z1,delz
      LOGICAL, INTENT (IN) :: invs
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN)    :: sk2(n2d)
      COMPLEX, INTENT (IN)    :: rhtxy(nmzxyd,n2d-1,2,jspd)
      COMPLEX, INTENT (INOUT) :: vxy(nmzxyd,n2d-1,2,jspd)
      COMPLEX, INTENT (OUT)   :: alphm(n2d,2)
C     ..
C     .. Local Scalars ..
      COMPLEX alph0,alph1,alph2,alphaz,betaz,test
      REAL g,vcons,z,tpi,maxexp,minexp,e_m,e_p
      INTEGER imz,ip,irec2,ivac
      LOGICAL tail
C     ..
C     .. Local Arrays ..
      REAL fra(nmzxyd),frb(nmzxyd),fia(nmzxyd),fib(nmzxyd)
      REAL alpha(nmzxyd,2,2),beta(nmzxyd,2,2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC aimag,cmplx,conjg,exp,real
C     ..
      tpi = 2 * pimach()
      maxexp = log(2.0)*MAXEXPONENT(tpi)
      minexp = log(2.0)*MINEXPONENT(tpi)
C     ..
c     2-dim star loop g.ne.0
      tail = .true.
      ip = nmzxy + 1
      DO 50 irec2 = 2,ng2
         g = sk2(irec2)
         vcons = tpi/g
         DO 20 ivac = 1,nvac
            z = z1
            DO 10 imz = 1,nmzxy
               IF (-g*z >= minexp ) THEN
                 e_m = exp( -g*z )
               ELSE
                 e_m = exp( minexp )
               ENDIF
               IF ( g*z <= maxexp ) THEN
                 e_p = exp( g*z )
               ELSE
                 e_p = exp( maxexp )
               ENDIF
               fra(ip-imz) = real(rhtxy(imz,irec2-1,ivac,1))* e_m
               fia(ip-imz) = aimag(rhtxy(imz,irec2-1,ivac,1))*e_m
               frb(imz) = real(rhtxy(imz,irec2-1,ivac,1))* e_p
               fib(imz) = aimag(rhtxy(imz,irec2-1,ivac,1))*e_p
               z = z + delz
   10       CONTINUE
            CALL intgz1(fra,delz,nmzxy,alpha(1,ivac,1),tail)
            CALL intgz1(fia,delz,nmzxy,alpha(1,ivac,2),tail)
            CALL qsf(delz,frb,beta(1,ivac,1),nmzxy,1)
            CALL qsf(delz,fib,beta(1,ivac,2),nmzxy,1)
   20    CONTINUE
         alph1 = cmplx(alpha(nmzxy,1,1),alpha(nmzxy,1,2))
         IF (nvac.EQ.1) THEN
            IF (invs) THEN
               alph2 = conjg(alph1)
            ELSE
               alph2 = alph1
            END IF
         ELSE
            alph2 = cmplx(alpha(nmzxy,2,1),alpha(nmzxy,2,2))
         END IF
         DO 40 ivac = 1,nvac
            z = z1
            IF (ivac.EQ.1) alph0 = alph2
            IF (ivac.EQ.2) alph0 = alph1
            DO 30 imz = 1,nmzxy
               betaz = cmplx(beta(imz,ivac,1),beta(imz,ivac,2))
               alphaz = cmplx(alpha(ip-imz,ivac,1),alpha(ip-imz,ivac,2))
               IF (-g*z >= minexp ) THEN
                 e_m = exp( -g*z  )
               ELSE
                 e_m = exp( minexp )
               ENDIF
               IF ( g*z  <= maxexp ) THEN
                 e_p = exp( g*z  )
               ELSE
                 e_p = exp( maxexp )
               ENDIF
               test = e_m*(alph0+betaz) + e_p*alphaz
               IF ( 2.0 * test == test ) test = cmplx(0.0,0.0)
               vxy(imz,irec2-1,ivac,1) = vxy(imz,irec2-1,ivac,1) +
     +                    vcons * test
               z = z + delz
   30       CONTINUE
   40    CONTINUE
         alphm(irec2-1,1) = alph1
         alphm(irec2-1,2) = alph2
   50 CONTINUE

      END SUBROUTINE vvacxy
      END MODULE m_vvacxy
