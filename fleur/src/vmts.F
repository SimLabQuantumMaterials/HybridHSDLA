      MODULE m_vmts
c     *******************************************************************
c     this subroutine calculates the lattice harmonics expansion coeffi-*
c     cients of the coulomb potential for all atom types                *
c                                c.l.fu, r.podloucky                    *
c     *******************************************************************
      CONTAINS
      SUBROUTINE vmts(irank,isize,
     >                n3d,jspd,memd,nlhd,ntypsd,jmtd,ntypd,lmaxd,
     >                ntype,nq3,zatom,dx,rmsh,jri,rmt,sk3,lmax,
     >                clnu,mlh,nmem,llh,nstr,ntypsy,nlh,invtab,
     >                nop,natd,neq,kv3,mrot,bmat,tau,taual,symor,
     >                vpw,rho,odi,ods,
     <                vr)

#include"cpp_double.h"
      USE m_constants, ONLY : pimach
      USE m_od_types, ONLY : od_inp, od_sym
      USE m_intgr, ONLY : intgr2
      USE m_phasy1
      USE m_sphbes
      USE m_od_phasy

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ...
      INTEGER, INTENT (IN) :: n3d,jspd,memd,nlhd,ntypsd,jmtd,ntypd,lmaxd
      INTEGER, INTENT (IN) :: ntype,nq3,nop,natd,irank,isize
      LOGICAL, INTENT (IN) :: symor
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),kv3(3,n3d),mrot(3,3,nop)
      REAL,    INTENT (IN) :: bmat(3,3),tau(3,nop),taual(3,natd)
      COMPLEX, INTENT (IN) :: vpw(n3d,jspd)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: jri(ntypd),ntypsy(natd),nstr(n3d)
      INTEGER, INTENT (IN) :: lmax(ntypd),invtab(nop)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      REAL,    INTENT (IN) :: zatom(ntypd),dx(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: rmt(ntypd),sk3(n3d)
      REAL,    INTENT (IN) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (OUT):: vr(jmtd,0:nlhd,ntypd,jspd)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      COMPLEX cp,sm
      REAL fpi,rmt2l,rmtl,ror,rr,rrlr,sfp,fpl21
      INTEGER i,jm,k,l,l21,lh,m,n,nd,nat,lm,n1
C     ..
C     .. Local Arrays ..
      COMPLEX vtl(0:nlhd,ntypd)
      COMPLEX pylm( (lmaxd+1)**2, ntypd )
      REAL    f1r(jmtd),f2r(jmtd),x1r(jmtd),x2r(jmtd)
      REAL    sbf(0:lmaxd),rrl(jmtd),rrl1(jmtd)
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER ierr(3)
      COMPLEX, ALLOCATABLE :: c_b(:)
C ..
C ..  External Subroutines
      EXTERNAL CPP_BLAS_ccopy,MPI_REDUCE
#endif
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg
C     ..
      fpi = 4 * pimach()
      sfp = sqrt( fpi )
C     ..
c     ----> calculate lattice harmonics expansion coefficients of the
c     ----> interstitial coulomb potential on the sphere boundaries
      nat = 1
      DO n = 1,ntype
         DO lh = 0,nlh(ntypsy(nat))
            vtl(lh,n) = 0.e0
         ENDDO
         nat = nat + neq(n) 
      ENDDO
#ifdef CPP_MPI
      vtl(:,:) = cmplx(0.0,0.0)
      CALL MPI_BCAST(vpw,n3d*jspd,CPP_MPI_COMPLEX,0,
     +                          MPI_COMM_WORLD,ierr)
#endif

c           ----> g=0 component
      IF (irank == 0) THEN
        DO n = 1,ntype
          vtl(0,n) = sfp*vpw(1,1)
        ENDDO
      ENDIF
c           ----> g.ne.0 components
      DO k = irank+2, nq3, isize
         cp = vpw(k,1)*nstr(k)
         IF (.NOT.odi%d1) THEN
            CALL phasy1(
     >                  ntypd,n3d,natd,nop,lmaxd,ntype,neq,lmax,
     >                  fpi,taual,bmat,kv3,tau,mrot,symor,k,invtab,
     <                  pylm)
         ELSE
c-odim
            CALL od_phasy(
     >           ntypd,n3d,natd,lmaxd,ntype,neq,lmax,
     >           fpi,taual,bmat,kv3,k,odi,ods,
     <           pylm)
c+odim
         END IF 
c
         nat = 1
         DO n = 1,ntype
            CALL sphbes(lmax(n),sk3(k)*rmt(n),sbf)
            nd = ntypsy(nat)
            DO lh = 0,nlh(nd)
               l = llh(lh,nd)
               sm = (0.,0.)
               DO jm = 1,nmem(lh,nd)
                  m = mlh(jm,lh,nd)
                  lm = l*(l+1) + m + 1 
                  sm = sm + conjg(clnu(jm,lh,nd))*pylm(lm,n)
               ENDDO
               vtl(lh,n) = vtl(lh,n) + cp*sbf(l)*sm
            ENDDO
            nat = nat + neq(n)
         ENDDO
      ENDDO

#ifdef CPP_MPI
        n1 = (nlhd+1)*ntypd
        ALLOCATE(c_b(n1))
        CALL MPI_REDUCE(vtl,c_b,n1,CPP_MPI_COMPLEX,MPI_SUM,0,
     +                                   MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n1, c_b, 1, vtl, 1)
        ENDIF
        DEALLOCATE (c_b)
#endif

c     ----> solution of the poisson's equation
      nat = 1
      DO 110 n = 1,ntype
         nd = ntypsy(nat)
         DO 100 lh = 0,nlh(nd)
            l = llh(lh,nd)
            l21 = 2*l + 1
            fpl21 = fpi/l21
            DO i = 1,jri(n)
               rrl(i) = rmsh(i,n)**l
               rrl1(i) = 1./( rrl(i) * rmsh(i,n) )
               x1r(i) = rrl(i)*rho(i,lh,n,1)
               x2r(i) = rrl1(i)*rho(i,lh,n,1)
            ENDDO
            CALL intgr2(x1r,rmsh(1,n),dx(n),jri(n),f1r)
            CALL intgr2(x2r,rmsh(1,n),dx(n),jri(n),f2r)
            rmtl = 1./rmt(n)**l
            rmt2l = 1./rmt(n)**l21
            DO 90 i = 1,jri(n)
               rrlr = rrl(i)*rmt2l
               ror = rrl(i)*rmtl
               vr(i,lh,n,1) = fpl21 * (rrl1(i)*f1r(i)-rrlr*f1r(jri(n))+
     +                   rrl(i) * (f2r(jri(n))-f2r(i))) + ror*vtl(lh,n)
   90       CONTINUE
  100    CONTINUE
         nat = nat + neq(n)
  110 CONTINUE
      DO 130 n = 1,ntype
         DO 120 i = 1,jri(n)
            rr = rmsh(i,n)/rmt(n)
            vr(i,0,n,1) = vr(i,0,n,1) - sfp* (1.-rr)/rmsh(i,n)*zatom(n)
  120    CONTINUE
  130 CONTINUE

      END SUBROUTINE vmts
      END MODULE m_vmts
