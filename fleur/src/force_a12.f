      MODULE m_forcea12
c ************************************************************
c Pulay 1st term  force contribution a la Rici et al., eq. A12
c
c ************************************************************
c
      CONTAINS
      SUBROUTINE force_a12(
     >                     lmaxd,ntypd,nobd,natd,nop,
     >                     jspd,mrot,ngopr,lmax,invarop,invarind,
     >                     ntype,neq,rmt,invtab,multab,amat,bmat,
     >                     we,jsp,ne,us,uds,acof,bcof,e1cof,e2cof,
     >                     lmd,nlod,llod,nlo,llo,acoflo,bcoflo,l_geo,
     >                     odi,ods,
     X                     force,f_a12)

      USE m_cotra, ONLY : cotra0,cotra1
      USE m_od_types, ONLY : od_inp, od_sym
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,ntypd,nobd,natd,nop,nlod,llod
      INTEGER, INTENT (IN) :: jspd,ne,ntype,jsp,lmd
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd)
      INTEGER, INTENT (IN) :: mrot(3,3,nop),ngopr(natd),lmax(ntypd)
      INTEGER, INTENT (IN) :: invarop(natd,nop),invarind(natd)
      INTEGER, INTENT (IN) :: multab(nop,nop),invtab(nop)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (IN) :: we(nobd),rmt(ntypd)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3)
      REAL,    INTENT (IN) ::  us(0:lmaxd,ntypd),uds(0:lmaxd,ntypd)
      COMPLEX, INTENT (IN) ::  acof(nobd,0:lmaxd*(lmaxd+2),natd)
      COMPLEX, INTENT (IN) ::  bcof(nobd,0:lmaxd*(lmaxd+2),natd)
      COMPLEX, INTENT (IN) :: e1cof(nobd,0:lmaxd*(lmaxd+2),natd)
      COMPLEX, INTENT (IN) :: e2cof(nobd,0:lmaxd*(lmaxd+2),natd)
      COMPLEX, INTENT (IN) :: acoflo(-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (IN) :: bcoflo(-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (INOUT) :: f_a12(3,ntypd)
      REAL,    INTENT (INOUT) :: force(3,ntypd,jspd)
      LOGICAL, INTENT (IN) :: l_geo(ntypd)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim

C     ..
C     .. Local Scalars ..
      COMPLEX a12,cil1,cil2,czero,ci
      REAL zero
      INTEGER i,ie,irinv,is,isinv,it,j,l,l1,l2,lm1,lm2,m,m1,m2,n,natom,
     +        natrun,ilo
C     ..
C     .. Local Arrays ..
      COMPLEX forc_a12(3),gv(3),acof_flapw(nobd,0:lmd),
     +        bcof_flapw(nobd,0:lmd)
      REAL aaa(2),bbb(2),ccc(2),ddd(2),eee(2),fff(2),gvint(3),
     +     starsum(3),vec(3),vecsum(3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,conjg,max,min,real,sign,sqrt,cmplx
C     ..
C     .. Statement Functions ..
      REAL alpha,beta,delta,epslon,gamma,phi
      INTEGER krondel
c
C     .. Data statements ..
c
      DATA czero/ (0.000,0.000)/,zero/0.000/
C     ..
C     .. Statement Function definitions ..
c  inline functions:
c
c Kronecker delta for arguments >=0 AND <0
c
c
      krondel(i,j) = min(abs(i)+1,abs(j)+1)/max(abs(i)+1,abs(j)+1)*
     +               (1+sign(1,i)*sign(1,j))/2
      alpha(l,m) = (l+1)*0.5e0*sqrt(real((l-m)* (l-m-1))/
     +             real((2*l-1)* (2*l+1)))
      beta(l,m) = l*0.5e0*sqrt(real((l+m+2)* (l+m+1))/
     +            real((2*l+1)* (2*l+3)))
      gamma(l,m) = (l+1)*0.5e0*sqrt(real((l+m)* (l+m-1))/
     +             real((2*l-1)* (2*l+1)))
      delta(l,m) = l*0.5e0*sqrt(real((l-m+2)* (l-m+1))/
     +             real((2*l+1)* (2*l+3)))
      epslon(l,m) = (l+1)*sqrt(real((l-m)* (l+m))/
     +               real((2*l-1)* (2*l+1)))
      phi(l,m) = l*sqrt(real((l-m+1)* (l+m+1))/real((2*l+1)* (2*l+3)))
C     ..
c
      ci = cmplx(0.0,1.0)
c
      natom = 1
      DO 10 n = 1,ntype
         IF (l_geo(n)) THEN
         DO i = 1,3
            forc_a12(i) = czero
         END DO
c
         DO natrun = natom,natom + neq(n) - 1
            DO i = 1,3
               gv(i) = czero
            END DO
c
c--->       the local orbitals do not contribute to
c--->       the term a12, because they vanish at the
c--->       mt-boundary. Therefore, the LO-contribution
c--->       to the a and b coefficients has to be
c--->       substracted before calculation a12.
c
            DO l1 = 0,lmax(n)
               DO m1 = -l1,l1
                  lm1 = l1* (l1+1) + m1
                  DO ie = 1,ne
                     acof_flapw(ie,lm1) = acof(ie,lm1,natrun)
                     bcof_flapw(ie,lm1) = bcof(ie,lm1,natrun)
                  ENDDO
               ENDDO
            ENDDO
            DO ilo = 1,nlo(n)
               l1 = llo(ilo,n)
               DO m1 = -l1,l1
                  lm1 = l1* (l1+1) + m1
                  DO ie = 1,ne
                     acof_flapw(ie,lm1) = acof_flapw(ie,lm1) -
     +                                    acoflo(m1,ie,ilo,natrun)
                     bcof_flapw(ie,lm1) = bcof_flapw(ie,lm1) -
     +                                    bcoflo(m1,ie,ilo,natrun)
                  ENDDO
               ENDDO
            ENDDO
c
            DO l1 = 0,lmax(n)
               cil1 = ci**l1
               DO m1 = -l1,l1
                  lm1 = l1* (l1+1) + m1
                  DO l2 = 0,lmax(n)
                     cil2 = ci**l2
                     DO m2 = -l2,l2
                        lm2 = l2* (l2+1) + m2
c
                        a12 = czero
                        DO ie = 1,ne
c
                           a12 = a12 + conjg(cil1*
     +                           ( acof_flapw(ie,lm1)*us(l1,n) +
     +                             bcof_flapw(ie,lm1)*uds(l1,n) ))*cil2*
     +                          ( e1cof(ie,lm2,natrun)*us(l2,n)+
     +                            e2cof(ie,lm2,natrun)*uds(l2,n))*we(ie)
                                        
                        END DO
                        aaa(1) = alpha(l1,m1)*krondel(l2,l1-1)*
     +                           krondel(m2,m1+1)
                        aaa(2) = alpha(l2,m2)*krondel(l1,l2-1)*
     +                           krondel(m1,m2+1)
                        bbb(1) = beta(l1,m1)*krondel(l2,l1+1)*
     +                           krondel(m2,m1+1)
                        bbb(2) = beta(l2,m2)*krondel(l1,l2+1)*
     +                           krondel(m1,m2+1)
                        ccc(1) = gamma(l1,m1)*krondel(l2,l1-1)*
     +                           krondel(m2,m1-1)
                        ccc(2) = gamma(l2,m2)*krondel(l1,l2-1)*
     +                           krondel(m1,m2-1)
                        ddd(1) = delta(l1,m1)*krondel(l2,l1+1)*
     +                           krondel(m2,m1-1)
                        ddd(2) = delta(l2,m2)*krondel(l1,l2+1)*
     +                           krondel(m1,m2-1)
                        eee(1) = epslon(l1,m1)*krondel(l2,l1-1)*
     +                           krondel(m2,m1)
                        eee(2) = epslon(l2,m2)*krondel(l1,l2-1)*
     +                           krondel(m1,m2)
                        fff(1) = phi(l1,m1)*krondel(l2,l1+1)*
     +                           krondel(m2,m1)
                        fff(2) = phi(l2,m2)*krondel(l1,l2+1)*
     +                           krondel(m1,m2)
c
                        gv(1) = gv(1) + (aaa(1)+bbb(1)-ccc(1)-ddd(1)+
     +                          aaa(2)+bbb(2)-ccc(2)-ddd(2))*0.5*
     +                          rmt(n)**2*a12
c
                        gv(2) = gv(2) + ci* (aaa(1)+bbb(1)+ccc(1)+
     +                          ddd(1)-aaa(2)-bbb(2)-ccc(2)-ddd(2))*0.5*
     +                          rmt(n)**2*a12
c
                        gv(3) = gv(3) + (eee(1)+eee(2)-fff(1)-fff(2))*
     +                          0.5*rmt(n)**2*a12
c
c  m1,m2 loops end
                     END DO
                  END DO
c  l1,l2 loops end
               END DO
            END DO
c
c  to complete summation over stars of k now sum
c  over all operations which leave (k+G)*R(natrun)*taual(natrun)
c  invariant. We sum over ALL these operations and not only
c  the ones needed for the actual star of k. Should be
c  ok if we divide properly by the number of operations
c  First, we find operation S where RS=T. T -like R- leaves
c  the above scalar product invariant (if S=1 then R=T).
c  R is the operation which generates position of equivalent atom
c  out of position of representative
c  S=R^(-1) T
c  number of ops which leave (k+G)*op*taual invariant: invarind
c  index of inverse operation of R: irinv
c  index of operation T: invarop
c  now, we calculate index of operation S: is
c
c  transform vector gv into internal coordinates
            DO i = 1,3
               vec(i) = real(gv(i)) /neq(n)
            END DO
            CALL cotra1(vec,gvint,bmat)
c
            DO i = 1,3
               vecsum(i) = zero
            END DO
!-gb2002
!            irinv = invtab(ngopr(natrun))
!            DO it = 1,invarind(natrun)
!               is = multab(irinv,invarop(natrun,it))
!c  note, actually we need the inverse of S but -in principle
!c  because {S} is agroup and we sum over all S- S should also
!c  work; to be lucid we take the inverse:
!                isinv = invtab(is)
!!               isinv = is
! Rotation is alreadt done in to_pulay, here we work only in the
! coordinate system of the representative atom (natom):
!!        
            DO it = 1,invarind(natom)
               is =invarop(natom,it)
               isinv = invtab(is)
               IF (odi%d1) isinv = ods%ngopr(natom)
!-gb 2002
c  now we have the wanted index of operation with which we have
c  to rotate gv. Note gv is given in cart. coordinates but
c  mrot acts on internal ones
               DO i = 1,3
                  vec(i) = zero
                  DO j = 1,3
                     IF (.NOT.odi%d1) THEN
                        vec(i) = vec(i) + mrot(i,j,isinv)*gvint(j)
                     ELSE
                        vec(i) = vec(i) + ods%mrot(i,j,isinv)*gvint(j)
                     END IF
                  END DO
               END DO
               DO i = 1,3
                  vecsum(i) = vecsum(i) + vec(i)
               END DO
c   end operator loop
            END DO
c
c   transform from internal to cart. coordinates
            CALL cotra0(vecsum,starsum,amat)
            DO i = 1,3
               forc_a12(i) = forc_a12(i) + starsum(i)/invarind(natrun)
            END DO
c
c  natrun loop end
         END DO
c
c
c
c     sum to existing forces
c
C  NOTE: force() is real and therefore takes only the
C  real part of forc_a12(). in general, force must be
C  real after k-star summation. Now, we put the proper
C  operations into real space. Problem: what happens
C  if in real space there is no inversion any more?
C  But we have inversion in k-space due to time reversal
C  symmetry, E(k)=E(-k)
C  We argue that k-space inversion is automatically taken
C  into account if force = (1/2)(forc_a12+conjg(forc_a12))
C  because time reversal symmetry means that conjg(PSI)
C  is also a solution of Schr. equ. if psi is one.
         DO i = 1,3
            force(i,n,jsp) = force(i,n,jsp) + real(forc_a12(i))
            f_a12(i,n)     = f_a12(i,n)     + forc_a12(i)
         END DO
c
c     write result moved to force_a8
c
         ENDIF
         natom = natom + neq(n)
   10 CONTINUE
c
      END SUBROUTINE force_a12
      END MODULE m_forcea12
