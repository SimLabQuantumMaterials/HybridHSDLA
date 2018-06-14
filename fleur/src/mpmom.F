      MODULE m_mpmom
c     ***********************************************************
c     determine the multipole moments of (original charge minus
c     plane wave charge) for each atom type
c                                     c.l.fu
c cf. M.Weinert J.Math.Phys. 22(11) (1981) p.2434 eq. (10)-(15)
c
c     qlmo(m,l,n) : mult.mom. of the mufftn-tin charge density
c     qlmp(m,l,n) : mult.mom. of the plane-wave charge density
c     qlm (m,l,n) : (output) difference of the former quantities
c     
c     ***********************************************************
      CONTAINS
      SUBROUTINE mpmom(irank,isize,
     >                 ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,memd,
     >                 ntype,nq3,jri,ntypsy,lmax,jspd,ncvd,
     >                 rmt,rmsh,dx,zatom,volmts,fpi,sfp,sk3,nstr,
     >                 nop,natd,neq,kv3,mrot,bmat,tau,taual,
     >                 clnu,mlh,nmem,llh,nlh,symor,invtab,
     >                 qpw,rho,odi,ods,
     <                 qlm)

#include"cpp_double.h"
      USE m_intgr, ONLY : intgr3
      USE m_phasy1
      USE m_sphbes
      USE m_od_types, ONLY : od_inp, od_sym
      USE m_od_phasy

      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd
      INTEGER, INTENT (IN) :: ntype,nq3,ncvd,nop,natd,irank,isize
      REAL,    INTENT (IN) :: sfp,fpi
      LOGICAL, INTENT (IN) :: symor
C
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),kv3(3,n3d),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: lmax(ntypd),ntypsy(natd),jri(ntypd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),nlh(ntypsd),nstr(n3d)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: invtab(nop)
      REAL,    INTENT (IN) :: bmat(3,3),tau(3,nop),taual(3,natd)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd),rmt(ntypd)
      REAL,    INTENT (IN) :: sk3(n3d),volmts(ntypd),zatom(ntypd)
      REAL,    INTENT (IN) :: rho(jmtd,0:nlhd,ntypd,jspd)
      COMPLEX, INTENT (IN) :: qpw(n3d,jspd),clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (OUT):: qlm(-lmaxd:lmaxd,0:lmaxd,ntypd)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      REAL fint,rl3,sk3r
      INTEGER j,jm,k,l,lh,m,mb,mem,mems,n,nd,nl,ns,nat,lm,ll1
      COMPLEX cil,nqpw,sk3i
C     ..
C     .. Local Arrays ..
      COMPLEX pylm( (lmaxd+1)**2, ntypd )
      COMPLEX qlmo(-lmaxd:lmaxd,0:lmaxd,ntypd)
      COMPLEX qlmp(-lmaxd:lmaxd,0:lmaxd,ntypd)
      REAL aj(0:lmaxd+ncvd+1)
      REAL f(jmtd)

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
      INTRINSIC cmplx
C     ..
C Initialise Arrays
C
      DO n = 1,ntype
         DO l = 0,lmax(n)
            DO m = -l,l
               qlm(m,l,n) = cmplx(0.e0,0.e0)
               qlmp(m,l,n) = cmplx(0.e0,0.e0)
               qlmo(m,l,n) = cmplx(0.e0,0.e0)
            ENDDO
         ENDDO
      ENDDO
C
C     multipole moments of original charge (q_{lm}^i)
C
      IF (irank == 0) THEN
      nat = 1
      DO n = 1, ntype
         ns = ntypsy(nat)
         jm = jri(n)
         DO nl = 0, nlh(ns)
            l = llh(nl,ns)
            DO j = 1, jm
               f(j) = (rmsh(j,n)**l)*rho(j,nl,n,1)
            ENDDO
            CALL intgr3(f,rmsh(1,n),dx(n),jm,fint)
            DO mb = 1, nmem(nl,ns)
               m = mlh(mb,nl,ns)
               qlmo(m,l,n) = qlmo(m,l,n) + clnu(mb,nl,ns)*fint
            ENDDO
         ENDDO
         qlmo(0,0,n) = qlmo(0,0,n) - zatom(n)/sfp
         nat = nat + neq(n)
      ENDDO
C
C     multipole moments of plane wave charge inside the spheres (q_{lm}^{Ii})
C
C     g eq 0 term : \sqrt{4 \pi}/3 R_i^3 \rho_I(0) \delta_{l,0}
C
      DO n = 1,ntype
         qlmp(0,0,n) = qpw(1,1)*nstr(1)*volmts(n)/sfp
      ENDDO
C
C     g ne 0 terms : \sum_{K \= 0} 4 \pi i^l \rho_I(K) R_i^{l+3} \times
C                    j_{l+1} (KR_i) / KR_i \exp{iK\xi_i} Y^*_{lm} (K)
C
      ENDIF ! irank == 0
#ifdef CPP_MPI
      CALL MPI_BCAST(qpw,n3d*jspd,CPP_MPI_COMPLEX,0,
     &                          MPI_COMM_WORLD,ierr)
#endif
      DO k = irank+2, nq3, isize
        IF (.NOT.odi%d1) THEN 
           CALL phasy1(
     >                 ntypd,n3d,natd,nop,lmaxd,ntype,neq,lmax,
     >                 fpi,taual,bmat,kv3,tau,mrot,symor,k,invtab,
     <                 pylm)
        ELSE
c-odim
            CALL od_phasy(
     >           ntypd,n3d,natd,lmaxd,ntype,neq,lmax,
     >           fpi,taual,bmat,kv3,k,odi,ods,
     <           pylm)
c+odim
        END IF
c
         nqpw = qpw(k,1)*nstr(k)
         DO n = 1,ntype
            sk3r = sk3(k)*rmt(n)
            CALL sphbes(lmax(n)+1,sk3r,aj)
            rl3 = rmt(n)**3
            sk3i = nqpw/sk3r
            DO l = 0,lmax(n)
               cil = aj(l+1)*sk3i*rl3
               ll1 = l*(l+1) + 1
               DO m = -l,l
                  lm = ll1 + m 
                  qlmp(m,l,n) = qlmp(m,l,n) + cil*pylm(lm,n)
               ENDDO
               rl3 = rl3*rmt(n)
            ENDDO                 ! l = 0, lmax(n)
         ENDDO                    ! n = 1, ntype
      ENDDO                       ! k = 2, nq3

#ifdef CPP_MPI
        n = (2*lmaxd+1)*(lmaxd+1)*ntypd
        ALLOCATE(c_b(n))
        CALL MPI_REDUCE(qlmp,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0,
     +                                   MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, qlmp, 1)
        ENDIF
        DEALLOCATE (c_b)
#endif

C
C eq.(15): \tilde q_(lm}^i = q_{lm}^i - q_{lm}^{Ii}
C
      IF (irank == 0) THEN
      DO n = 1,ntype
         DO l = 0,lmax(n)
            DO m = -l,l
               qlm(m,l,n) = qlmo(m,l,n) - qlmp(m,l,n)
            ENDDO
         ENDDO                 ! l = 0, lmax(n)
      ENDDO                    ! n = 1, ntype
C
C Output section
C
      nat = 1
      DO 180 n = 1,ntype
         WRITE (6,FMT=8000) n
         nd = ntypsy(nat)
         DO 170 lh = 0,nlh(nd)
            l = llh(lh,nd)
            mems = nmem(lh,nd)
            DO 160 mem = 1,mems
               m = mlh(mem,lh,nd)
               WRITE (6,FMT=8010) l,m,qlmo(m,l,n),qlmp(m,l,n)
c     write(16,1002) l,m,qlmo(m,l,n),qlmp(m,l,n)
  160       CONTINUE
  170    CONTINUE
         nat = nat + neq(n)
  180 CONTINUE
c
 8000 FORMAT (/,10x,'multipole moments for atom type=',i5,/,/,t3,'l',t7,
     +       'm',t27,'original',t57,'plane wave')
 8010 FORMAT (1x,i2,2x,i2,2x,2 (5x,2e15.5))
c
      ENDIF ! irank == 0

      END SUBROUTINE mpmom
      END MODULE m_mpmom
