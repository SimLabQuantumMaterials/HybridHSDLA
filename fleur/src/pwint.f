      MODULE m_pwint
c     ******************************************************************
c     calculate the integral of a star function over the interstial    *
c     region              c.l.fu                                       *
c     ******************************************************************
      CONTAINS
      SUBROUTINE pwint(
     >                 k1d,k2d,k3d,n3d,ntypd,natd,nop,invtab,odi,
     >                 ntype,neq,volmts,taual,z1,vol,volint,
     >                 symor,tpi,tau,mrot,rmt,sk3,bmat,ig2,ig,
     >                 kv,
     <                 x)

      USE m_spgrot
      USE m_od_cylbes
      USE m_od_types , ONLY : od_inp
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,n3d,ntypd,natd,nop
      INTEGER, INTENT (IN) :: ntype
      REAL,    INTENT (IN) :: z1,vol,volint,tpi
      LOGICAL, INTENT (IN) :: symor
      COMPLEX, INTENT (OUT):: x
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: kv(3),mrot(3,3,nop),neq(ntypd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (IN) :: invtab(nop)
      REAL,    INTENT (IN) :: taual(3,natd),volmts(ntypd),tau(3,nop)
      REAL,    INTENT (IN) :: bmat(3,3),rmt(ntypd),sk3(n3d)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim
C     ..
C     .. Local Scalars ..
      COMPLEX s1,sfs
      REAL arg,g,s,srmt,gr,fJ
      INTEGER ig2d,ig3d,n,nn,nat,ii
C     ..
C     .. Local Arrays ..
      REAL ph(nop)
      INTEGER kr(3,nop)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,cos,exp,sin
C     ..
      ig3d = ig(kv(1),kv(2),kv(3))
      IF (ig3d.EQ.0) THEN
         x = (0.,0.)
         RETURN
      END IF
      IF (ig3d.EQ.1) THEN
         x = cmplx(volint,0.0)
         RETURN
      ELSE
         IF (odi%d1) THEN
            IF (kv(3).EQ.0) THEN
               g = (kv(1)*bmat(1,1) + kv(2)*bmat(2,1))**2 +
     +             (kv(1)*bmat(1,2) + kv(2)*bmat(2,2))**2
               gr = sqrt(g)
               g  = gr*z1
               CALL od_cylbes(1,g,fJ)
               x = cmplx(2*vol*fJ/g,0.0)
            ELSE
               x = (0.0,0.0)
            END IF
          ELSE
             ig2d = ig2(ig3d)
             IF (ig2d.EQ.1) THEN
                g = kv(3)*bmat(3,3)*z1
                x = cmplx(vol*sin(g)/g,0.0)
             ELSE
                x = (0.0,0.0)
             END IF
          END IF
      END IF
c     -----> sphere contributions
      s = sk3(ig3d)
      nat = 1
      IF (.NOT.odi%d1) THEN
         DO 20 n = 1,ntype
            srmt = s*rmt(n)
            CALL spgrot(
     >           nop,symor,tpi,mrot,tau,invtab,
     >           kv,
     <           kr,ph)
            sfs = (0.0,0.0)
            DO 10 nn = 1,nop
               arg = tpi* (kr(1,nn)*taual(1,nat)+kr(2,nn)*taual(2,nat)+
     +              kr(3,nn)*taual(3,nat))
               sfs = sfs + exp(cmplx(0.0,arg))*ph(nn)
 10         CONTINUE
            sfs = sfs/nop
c     -----3*ji(gr)/gr term
            s1 = 3.* (sin(srmt)/srmt-cos(srmt))/ (srmt*srmt)
            x = x - volmts(n)*neq(n)*s1*sfs
            nat = nat + neq(n)
 20      CONTINUE
      ELSE
c-odim
         DO 21 n = 1,ntype
            DO ii = 1,neq(n)
               srmt = s*rmt(n)
               CALL spgrot(
     >              nop,symor,tpi,mrot,tau,invtab,
     >              kv,
     <              kr,ph)
               sfs = (0.0,0.0)
               DO 11 nn = 1,nop
                  arg = tpi* (kr(1,nn)*taual(1,nat)+
     +                 kr(2,nn)*taual(2,nat)+
     +                 kr(3,nn)*taual(3,nat))
                  sfs = sfs + exp(cmplx(0.0,arg))*ph(nn)
 11            CONTINUE
               sfs = sfs/nop
c     -----3*ji(gr)/gr term
               s1 = 3.* (sin(srmt)/srmt-cos(srmt))/ (srmt*srmt)
               x = x - volmts(n)*s1*sfs
               nat = nat + 1
            END DO
 21      CONTINUE
c+odim
      ENDIF

      END SUBROUTINE pwint
      END MODULE m_pwint
