      MODULE m_phasy1
c     ********************************************************************
c     calculate 4pi*i**l/nop(3)*sum(R){exp(iRG(taual-taur)*conjg(ylm(RG)) }
c     e. wimmer   oct.1984
c     ********************************************************************
      CONTAINS
      SUBROUTINE phasy1(
     >                  ntypd,n3d,natd,nop,lmaxd,ntype,neq,lmax,
     >                  fpi,taual,bmat,kv3,tau,mrot,symor,k,invtab,
     <                  pylm)

      USE m_ylm
      USE m_spgrot

      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,n3d,natd,nop,lmaxd
      INTEGER, INTENT (IN) :: ntype,k
      REAL,    INTENT (IN) :: fpi
      LOGICAL, INTENT (IN) :: symor
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd),kv3(3,n3d)
      INTEGER, INTENT (IN) :: mrot(3,3,nop),invtab(nop)
      REAL,    INTENT (IN) :: bmat(3,3),tau(3,nop),taual(3,natd)
      COMPLEX, INTENT (OUT):: pylm( (lmaxd+1)**2, ntypd )
C     ..
C     .. Local Scalars ..
      COMPLEX sf,ci,csf
      REAL x,tpi
      INTEGER j,l,m,n,na,lm,ll1
C     ..
C     .. Local Arrays ..
      COMPLEX ciall(0:lmaxd)
      REAL phas(nop),rg(3)
      INTEGER kr(3,nop)
      COMPLEX, ALLOCATABLE :: ylm(:,:)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,conjg,cos,sin
C     ..
      tpi = fpi / 2.0
      ci = cmplx(0.0,1.0)
      ciall(0) = fpi/nop
      DO l = 1,lmaxd
         ciall(l) = ciall(0)*ci**l
      ENDDO

      CALL spgrot(
     >           nop,symor,tpi,mrot,tau,invtab,
     >           kv3(1,k),
     <           kr,phas)
      ALLOCATE ( ylm( (lmaxd+1)**2, nop ) )
      DO j = 1,nop
          rg(1) = kr(1,j)*bmat(1,1) + kr(2,j)*bmat(2,1) +
     +            kr(3,j)*bmat(3,1)
          rg(2) = kr(1,j)*bmat(1,2) + kr(2,j)*bmat(2,2) +
     +            kr(3,j)*bmat(3,2)
          rg(3) = kr(1,j)*bmat(1,3) + kr(2,j)*bmat(2,3) +
     +            kr(3,j)*bmat(3,3)
          CALL ylm4(
     >              lmaxd,rg,
     <              ylm(1,j) )
      ENDDO
      ylm = conjg( ylm )

      na = 1
      DO n = 1,ntype
         DO lm = 1, (lmax(n)+1)**2
               pylm(lm,n) = cmplx(0.,0.)
         ENDDO
         DO j = 1,nop
            x = tpi* (kr(1,j)*taual(1,na)+kr(2,j)*taual(2,na)+
     +          kr(3,j)*taual(3,na))
            sf = cmplx(cos(x),sin(x))*phas(j)
            DO l = 0,lmax(n)
               ll1 = l*(l+1) + 1
               csf = ciall(l)*sf
               DO m = -l,l
                  lm = ll1 + m
                  pylm(lm,n) = pylm(lm,n) + csf*ylm(lm,j)
               ENDDO
            ENDDO
         ENDDO
         na = na + neq(n)
      ENDDO
      DEALLOCATE ( ylm )

      END SUBROUTINE phasy1
      END MODULE m_phasy1
