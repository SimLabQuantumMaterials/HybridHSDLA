      MODULE m_cdntot
c     ********************************************************
c     calculate the total charge density in the interstial.,
c     vacuum, and mt regions      c.l.fu
c     ********************************************************
      CONTAINS
      SUBROUTINE cdntot(
     >                  k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,
     >                  nlhd,nmzd,nmz,jspins,film,nvac,area,nq3,
     >                  ntype,neq,volmts,taual,z1,vol,volint,
     >                  symor,tau,mrot,rmt,sk3,bmat,ig2,ig,
     >                  nstr,kv3,delz,jri,dx,rmsh,invtab,
     >                  qpw,rho,rht,odi,
     <                  qtot,qistot)

      USE m_intgr, ONLY : intgr3
      USE m_constants, ONLY : pimach
      USE m_qsf
      USE m_pwint
      USE m_od_types, ONLY : od_inp

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd
      INTEGER, INTENT (IN) :: nvac,nq3,ntype,nmz,jspins,nlhd,nmzd
      REAL,    INTENT (IN) :: z1,vol,volint,area,delz
      LOGICAL, INTENT (IN) :: symor,film
      REAL,    INTENT (OUT):: qtot,qistot
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: qpw(n3d,jspd)
      INTEGER, INTENT (IN) :: kv3(3,n3d),jri(ntypd),invtab(nop)
      INTEGER, INTENT (IN) :: mrot(3,3,nop),neq(ntypd),nstr(n3d)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      REAL,    INTENT (IN) :: taual(3,natd),volmts(ntypd),tau(3,nop)
      REAL,    INTENT (IN) :: bmat(3,3),rmt(ntypd),sk3(n3d)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),dx(ntypd)
      REAL,    INTENT (IN) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (IN) :: rht(nmzd,2,jspd)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim
C     ..
C     .. Local Scalars ..
      COMPLEX x
      REAL q,qis,w,tpi,sfp
      INTEGER i,ivac,j,jspin,n,nz
C     ..
C     .. Local Arrays ..
      REAL qmt(ntypd),qvac(2),q2(nmz),rht1(nmzd,2,jspd)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real
C     ..
      tpi = 2 * pimach()
      sfp = sqrt( 2 * tpi )
c
      qtot = 0.e0
      qistot = 0.e0
      DO 40 jspin = 1,jspins
         q = 0.e0
c     -----mt charge
         DO 10 n = 1,ntype
            CALL intgr3(rho(1,0,n,jspin),rmsh(1,n),dx(n),jri(n),w)
            qmt(n) = w*sfp
            q = q + neq(n)*qmt(n)
   10    CONTINUE
c     -----vacuum region
         IF (film) THEN
            DO 20 ivac = 1,nvac
               DO nz = 1,nmz
                  IF (odi%d1) THEN
                     rht1(nz,ivac,jspin) = (z1+(nz-1)*delz)*
     *                    rht(nz,ivac,jspin)
                  ELSE
                     rht1(nz,ivac,jspin) =  rht(nz,ivac,jspin)
                  END IF
               END DO
               CALL qsf(delz,rht1(1,ivac,jspin),q2,nmz,0)
               qvac(ivac) = q2(1)*area
               IF (.NOT.odi%d1) THEN
                  q = q + qvac(ivac)*2./real(nvac)
               ELSE
                  q = q + area*q2(1)
               END IF
   20       CONTINUE
         END IF
c     -----is region
         qis = 0.
         DO 30 j = 1,nq3
            CALL pwint(
     >                 k1d,k2d,k3d,n3d,ntypd,natd,nop,invtab,odi,
     >                 ntype,neq,volmts,taual,z1,vol,volint,
     >                 symor,tpi,tau,mrot,rmt,sk3,bmat,ig2,ig,
     >                 kv3(1,j),
     <                 x)
            qis = qis + qpw(j,jspin)*x*nstr(j)
   30    CONTINUE
         qistot = qistot + qis
         q = q + qis
         WRITE (6,FMT=8000) jspin,q,qis, (qmt(n),n=1,ntype)
         IF (film) WRITE (6,FMT=8010) (i,qvac(i),i=1,nvac)
         WRITE (16,FMT=8000) jspin,q,qis, (qmt(n),n=1,ntype)
         IF (film) WRITE (16,FMT=8010) (i,qvac(i),i=1,nvac)
         qtot = qtot + q
   40 CONTINUE
      WRITE (6,FMT=8020) qtot
      WRITE (16,FMT=8020) qtot
 8000 FORMAT (/,10x,'total charge for spin',i3,'=',f12.6,/,10x,
     +       'interst. charge =   ',f12.6,/,
     +       (10x,'mt charge=          ',4f12.6,/))
 8010 FORMAT (10x,'vacuum ',i2,'  charge=  ',f12.6)
 8020 FORMAT (/,10x,'total charge  =',f12.6)

      END SUBROUTINE cdntot
      END MODULE m_cdntot
