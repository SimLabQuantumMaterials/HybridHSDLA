      MODULE m_pwintsl
      CONTAINS
      SUBROUTINE pwint_sl(
     >                    k1d,k2d,k3d,n3d,ntypd,natd,nop,invtab,
     >                    ntype,neq,volmts,taual,zsl1,zsl2,
     >                    volsl,volintsl,symor,tpi,tau,mrot,
     >                    rmt,sk3,bmat,ig2,ig,nmtsl1,
     >                    kv,
     <                    x)
c     ******************************************************************
c     calculate the integral of a star function over the layer 
c     interstial region of a film                Yury Koroteev  
c                                   from  pwint.F  by  c.l.fu              
c     ******************************************************************
      USE m_spgrot
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,n3d,ntypd,natd,nop
      INTEGER, INTENT (IN) :: ntype
      REAL,    INTENT (IN) :: zsl1,zsl2,volsl,volintsl,tpi
      LOGICAL, INTENT (IN) :: symor
      COMPLEX, INTENT (OUT):: x
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: kv(3),mrot(3,3,nop),neq(ntypd)
      INTEGER, INTENT (IN) :: nmtsl1(ntypd),invtab(nop)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      REAL,    INTENT (IN) :: taual(3,natd),volmts(ntypd),tau(3,nop)
      REAL,    INTENT (IN) :: bmat(3,3),rmt(ntypd),sk3(n3d)
C     ..
C     .. Local Scalars ..
      COMPLEX s1,sfs
      REAL arg,g,s,srmt,gm,gp,zslm,zslp
      INTEGER ig2d,ig3d,n,nn,nat
C     ..
C     .. Local Arrays ..
      REAL ph(nop)
      INTEGER kr(3,nop)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,cos,sin
C     ..
      ig3d = ig(kv(1),kv(2),kv(3))
c
c     -----> interstitial contributions
c
      IF (ig3d.EQ.0) THEN
         x = (0.,0.)
         RETURN
      END IF
      IF (ig3d.EQ.1) THEN
         x = cmplx(volintsl,0.0)
         RETURN
      ELSE
         ig2d = ig2(ig3d)
         IF (ig2d.EQ.1) THEN
            zslm = 0.5*(zsl2 - zsl1) 
	    zslp = 0.5*(zsl2 + zsl1)
            g = kv(3)*bmat(3,3)		   
            gm = g*zslm
            gp = g*zslp
!            x = cmplx(volsl*cos(gp)*sin(gm)/gm,sin(gp))
            x = volsl*sin(gm)/gm*cmplx(cos(gp),sin(gp))
         ELSE
            x = (0.0,0.0)
         END IF
      END IF
c
c     -----> sphere contributions
c
      s = sk3(ig3d)
      nat = 1
      DO 20 n = 1,ntype
         srmt = s*rmt(n)
         CALL spgrot(
     >               nop,symor,tpi,mrot,tau,invtab,
     >               kv,
     <               kr,ph)
         sfs = (0.0,0.0)
         DO 10 nn = 1,nop
            arg = tpi* (kr(1,nn)*taual(1,nat)+kr(2,nn)*taual(2,nat)+
     +            kr(3,nn)*taual(3,nat))
            sfs = sfs + cmplx(cos(arg),sin(arg))*ph(nn)
   10    CONTINUE
         sfs = sfs/nop
c
c     -----3*ji(gr)/gr term
c
         s1 = 3.* (sin(srmt)/srmt-cos(srmt))/ (srmt*srmt)
         x = x - volmts(n)*nmtsl1(n)*s1*sfs
         nat = nat + neq(n)
   20 CONTINUE
c
      END SUBROUTINE pwint_sl
      END MODULE m_pwintsl
