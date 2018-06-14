      MODULE m_qfix
c     *******************************************************
c     check total charge and renormalize        c,l.fu
c     *******************************************************
      CONTAINS
      SUBROUTINE qfix(
     >                k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,nmzxyd,
     >                nlhd,nmzd,nmz,jspins,film,nvac,area,nq3,nmzxy,n2d,
     >                ntype,neq,volmts,taual,z1,vol,volint,nq2,invtab,
     >                symor,tau,mrot,rmt,sk3,bmat,ig2,ig,nlh,ntypsd,
     >                nstr,kv3,delz,jri,dx,rmsh,zatom,ntypsy,sigma,
     X                qpw,rhtxy,rho,rht,odi,
     <                fix)

      USE m_od_types, ONLY : od_inp
      USE m_cdntot
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd
      INTEGER, INTENT (IN) :: nvac,nq3,ntype,nmz,jspins,nlhd,nmzd
      INTEGER, INTENT (IN) :: nmzxyd,nmzxy,n2d,nq2,ntypsd
      REAL,    INTENT (IN) :: z1,vol,volint,area,delz,sigma
      LOGICAL, INTENT (IN) :: symor,film
      REAL,    INTENT (OUT) :: fix
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: kv3(3,n3d),jri(ntypd),invtab(nop)
      INTEGER, INTENT (IN) :: nlh(ntypsd),ntypsy(natd)
      INTEGER, INTENT (IN) :: mrot(3,3,nop),neq(ntypd),nstr(n3d)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      REAL,    INTENT (IN) :: taual(3,natd),volmts(ntypd),tau(3,nop)
      REAL,    INTENT (IN) :: bmat(3,3),rmt(ntypd),sk3(n3d)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),dx(ntypd),zatom(ntypd)
      COMPLEX,INTENT (INOUT) :: qpw(n3d,jspd)
      COMPLEX,INTENT (INOUT) :: rhtxy(nmzxyd,odi%n2d-1,2,jspd)
      REAL,   INTENT (INOUT) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL,   INTENT (INOUT) :: rht(nmzd,2,jspd)
C     ..
C     .. Local Scalars ..
      LOGICAL fixtot,l99
      REAL qtot,qis,zc
      INTEGER ivac,j,jm,jspin,k,lh,n,nl,ns,na
C     ..
      CALL cdntot(
     >            k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,
     >            nlhd,nmzd,nmz,jspins,film,nvac,area,nq3,
     >            ntype,neq,volmts,taual,z1,vol,volint,
     >            symor,tau,mrot,rmt,sk3,bmat,ig2,ig,
     >            nstr,kv3,delz,jri,dx,rmsh,invtab,
     >            qpw,rho,rht,odi,
     <            qtot,qis)
      zc = 0.
      DO 10 n = 1,ntype
         zc = zc + neq(n)*zatom(n)
   10 CONTINUE

c+roa (check via qfix file if total charge or only interstitial to fix)
      fixtot=.TRUE.
      INQUIRE(file='qfix',exist=l99)
      IF (l99) then
        OPEN (99,file='qfix',status='old',form='formatted')
        READ (99,'(1x,l1)',end=1199) fixtot
        IF (.NOT.fixtot ) THEN
           REWIND (99)
           WRITE (99,'(1x,l1,70a)') .TRUE.,
     +      ' (1x,l1) F..fix interstitial T..fix total charge '
        ENDIF
 1199   CLOSE (99)
      ENDIF
      zc = zc + 2*sigma
      IF ( fixtot ) THEN
c-roa
         fix = zc/qtot
         DO 100 jspin = 1,jspins
            na = 1
            DO 40 n = 1,ntype
               ns = ntypsy(na)
               lh = nlh(ns)
               jm = jri(n)
               DO 30 nl = 0,lh
                  DO 20 j = 1,jm
                     rho(j,nl,n,jspin) = fix*rho(j,nl,n,jspin)
   20             CONTINUE
   30          CONTINUE
               na = na + neq(n)
   40       CONTINUE
            DO 50 k = 1,nq3
               qpw(k,jspin) = fix*qpw(k,jspin)
   50       CONTINUE
            IF (film) THEN
               DO 90 ivac = 1,nvac
                  DO 60 n = 1,nmz
                     rht(n,ivac,jspin) = fix*rht(n,ivac,jspin)
   60             CONTINUE
                  DO 80 n = 1,nmzxy
                     DO 70 k = 2,odi%nq2
                        rhtxy(n,k-1,ivac,jspin) = fix*
     +                    rhtxy(n,k-1,ivac,jspin)
   70                CONTINUE
   80             CONTINUE
   90          CONTINUE
            END IF
  100    CONTINUE
         WRITE (6,FMT=8000) zc,fix
         CALL cdntot(
     >               k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,
     >               nlhd,nmzd,nmz,jspins,film,nvac,area,nq3,
     >               ntype,neq,volmts,taual,z1,vol,volint,
     >               symor,tau,mrot,rmt,sk3,bmat,ig2,ig,
     >               nstr,kv3,delz,jri,dx,rmsh,invtab,
     >               qpw,rho,rht,odi,
     <               qtot,qis)
c+roa 
      ELSE
         fix = (zc - qtot) / qis + 1.
         DO jspin = 1,jspins
            DO k = 1,nq3
               qpw(k,jspin) = fix*qpw(k,jspin)
            ENDDO
         ENDDO
         WRITE (6,FMT=8001) zc,fix
         CALL cdntot(
     >               k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,
     >               nlhd,nmzd,nmz,jspins,film,nvac,area,nq3,
     >               ntype,neq,volmts,taual,z1,vol,volint,
     >               symor,tau,mrot,rmt,sk3,bmat,ig2,ig,
     >               nstr,kv3,delz,jri,dx,rmsh,invtab,
     >               qpw,rho,rht,odi,
     <               qtot,qis)

      ENDIF
 8000 FORMAT (/,10x,'zc= ',f12.6,5x,'qfix=  ',f10.6)
 8001 FORMAT (/,' > broy only qis: ','zc= ',f12.6,5x,'qfix=  ',f10.6)
c-roa

      END SUBROUTINE qfix
      END MODULE m_qfix
