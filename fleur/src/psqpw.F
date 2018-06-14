      MODULE m_psqpw
c     ***********************************************************
c     generates the fourier coefficients of pseudo charge density
c                                   c.l.fu
c         corrected april 1990   m.w.
c
c cf. M.Weinert J.Math.Phys. 22(11) (1981) p.2434 eq. (10)-(15)
c
!
!     parallelized 04/08 gb
c     ***********************************************************
      CONTAINS
      SUBROUTINE psqpw(irank,isize,
     >                 ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,memd,nmzd,
     >                 ntype,nq3,jri,ntypsy,lmax,jspd,ncvd,
     >                 rmt,rmsh,dx,zatom,volmts,sfp,sk3,nstr,
     >                 omtil,dvac,z1,area,film,nvac,neq,nmz,delz,
     >                 clnu,mlh,nmem,llh,nlh,bmat,kv3,ig2,
     >                 nop,natd,symor,mrot,tau,taual,invtab,
     >                 ncv,qpw,rho,rht,sigma,vol,odi,ods,l_xyav,
     <                 psq)
c
#include"cpp_double.h"
      USE m_constants, ONLY : pimach
      USE m_phasy1
      USE m_mpmom 
      USE m_sphbes
      USE m_qsf
      USE m_od_phasy
      USE m_od_cylbes
      USE m_od_types, ONLY : od_inp, od_sym

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: irank,isize
      INTEGER, INTENT (IN) :: ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,jspd,memd
      INTEGER, INTENT (IN) :: ntype,nq3,ncvd,nmzd,nmz,nvac,nop,natd
      REAL,    INTENT (IN) :: omtil,dvac,z1,area,delz,sigma,vol
      LOGICAL, INTENT (IN) :: film,symor,l_xyav
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: mrot(3,3,nop),invtab(nop)
      COMPLEX, INTENT (IN) :: qpw(n3d,jspd),clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: lmax(ntypd),ntypsy(natd),jri(ntypd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),nlh(ntypsd),nstr(n3d)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: kv3(3,n3d),ig2(n3d),neq(ntypd),ncv(ntypd)
      REAL,    INTENT (IN) :: tau(3,nop),taual(3,natd)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd),rmt(ntypd)
      REAL,    INTENT (IN) :: sk3(n3d),volmts(ntypd),zatom(ntypd)
      REAL,    INTENT (IN) :: rho(jmtd,0:nlhd,ntypd,jspd),bmat(3,3)
      REAL,    INTENT (IN) :: rht(nmzd,2,jspd)
      COMPLEX, INTENT (OUT):: psq(n3d)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      COMPLEX psint,sa,sl,sm
      REAL f,fact,fpo,gz,p,qvac,rmtl,s,sfp,fpi,fJ,gr,g
      INTEGER ivac,k,l,m,n,n1,nc,ncvn,lm,ll1,nz,nd
C     ..
C     .. Local Arrays ..
      COMPLEX pylm( (lmaxd+1)**2, ntypd )
      COMPLEX  qlm(-lmaxd:lmaxd,0:lmaxd,ntypd)
      REAL q2(nmzd),pn(0:lmaxd,ntypd),aj(0:lmaxd+ncvd+1)
      REAL rht1(nmz)
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
      INTRINSIC abs,conjg,real,sin
C     ..
      fpi = 4 * pimach()
      sfp = sqrt( fpi )
C
C Calculate multipole moments
C
      CALL mpmom(irank,isize,
     >           ntypd,ntypsd,nlhd,n3d,jmtd,lmaxd,memd,
     >           ntype,nq3,jri,ntypsy,lmax,jspd,ncvd,
     >           rmt,rmsh,dx,zatom,volmts,fpi,sfp,sk3,nstr,
     >           nop,natd,neq,kv3,mrot,bmat,tau,taual,
     >           clnu,mlh,nmem,llh,nlh,symor,invtab,
     >           qpw,rho,odi,ods,
     <           qlm)
#ifdef CPP_MPI
      psq(:) = cmplx(0.0,0.0)
      CALL MPI_BCAST(qpw,n3d*jspd,CPP_MPI_COMPLEX,0,
     +                          MPI_COMM_WORLD,ierr)
      nd = (2*lmaxd+1)*(lmaxd+1)*ntypd
      CALL MPI_BCAST(qlm,nd,CPP_MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
#endif
C
C pn(l,n) = (2l + 2nc(n) + 3)!! / (2l + 1)!! R^l  ;   ncv(n)=n+l in paper
C 
      DO n = 1,ntype
         rmtl = 1.0
         DO l = 0,lmax(n)
            IF (l.GE.ncv(n)) THEN
               pn(l,n) = 0.0
            ELSE
               p = 1.
               DO nc = l,ncv(n)
                  p = p* (2*nc+3)
               ENDDO
               pn(l,n) = p/rmtl
            END IF
            rmtl = rmtl*rmt(n)
         ENDDO
      ENDDO
C
C G eq 0 term (eq.29) : \tilde \rho_s (0) = \sqrt{4 pi} / \Omega \sum_i \tilde q_{00}^i
C
      s = 0.
      DO n = 1,ntype
         s = s + neq(n)*real(qlm(0,0,n))
      ENDDO
      IF (irank == 0) THEN
        psq(1) = qpw(1,1) + (sfp/omtil)*s
      ENDIF
C
C G ne 0 term (eq.28) : \tilde \rho_s (K) = 4 pi / \Omega \sum_{lmi} (-i)^l \exp{-iK\xi_i}
C                    (2n+3)!!/(2l+1)!! * 1/R_i^l * j_{n+1}(KR_i)/(KR_i)^{n-l+1} Y_{lm} (K)
C
      fpo = 1./omtil
      DO k = irank+2, nq3, isize
       IF (.NOT.odi%d1) THEN
         CALL phasy1(
     >               ntypd,n3d,natd,nop,lmaxd,ntype,neq,lmax,
     >               fpi,taual,bmat,kv3,tau,mrot,symor,k,invtab,
     <               pylm)
       ELSE
c-odim
         CALL od_phasy(
     >        ntypd,n3d,natd,lmaxd,ntype,neq,lmax,
     >        fpi,taual,bmat,kv3,k,odi,ods,
     <        pylm)
c+odim
       END IF
c
         sa = 0.
         DO n = 1,ntype
            ncvn = ncv(n)
            CALL sphbes(ncvn+1,sk3(k)*rmt(n),aj)
            sl = 0.
            DO l = 0,lmax(n)
               IF (l.GE.ncvn) GO TO 60
               n1 = ncvn - l + 1
               ll1 = l*(l+1) + 1
               sm = 0.
               DO m = -l,l
                  lm = ll1 + m 
                  sm = sm + qlm(m,l,n)*conjg(pylm(lm,n))
               ENDDO
   60          sl = sl + pn(l,n)/ ((sk3(k)*rmt(n))**n1)*aj(ncvn+1)*sm
            ENDDO
            sa = sa + neq(n)*sl
         ENDDO
         psq(k) = qpw(k,1) + fpo*sa
      ENDDO

#ifdef CPP_MPI
        ALLOCATE(c_b(n3d))
        CALL MPI_REDUCE(psq,c_b,n3d,CPP_MPI_COMPLEX,MPI_SUM,0,
     +                                    MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n3d, c_b, 1, psq, 1)
        ENDIF
        DEALLOCATE (c_b)
#endif
      IF (irank == 0) THEN
C
C Check: integral of the pseudo charge density within the slab
C
      IF (film .AND. .NOT.odi%d1) THEN
         psint = psq(1)*nstr(1)*dvac
         DO k = 2,nq3
            IF (ig2(k).EQ.1) THEN
               gz = kv3(3,k)*bmat(3,3)
               f = 2.*sin(gz*z1)/gz
               psint = psint + nstr(k)*psq(k)*f
            END IF
         ENDDO 
         psint = area*psint
      ELSEIF (film .AND. odi%d1) THEN
c-odim
         psint = (0.0,0.0)
         DO 95 k = 2,nq3
            IF (kv3(3,k).EQ.0) THEN
               g = (kv3(1,k)*bmat(1,1) + kv3(2,k)*bmat(2,1))**2 +
     +             (kv3(1,k)*bmat(1,2) + kv3(2,k)*bmat(2,2))**2
               gr = sqrt(g)
               CALL od_cylbes(1,gr*z1,fJ)
               f = 2*vol*fJ/(gr*z1)
               psint = psint + nstr(k)*psq(k)*f
            END IF
 95      CONTINUE
         psint = psint + psq(1)*nstr(1)*vol
c+odim
      ELSEIF (.NOT.film) THEN
         psint = psq(1)*nstr(1)*omtil
      ENDIF
      WRITE (6,FMT=8000) psint
      WRITE (16,FMT=8000) psint
 8000 FORMAT (/,10x,'integral of pseudo charge density inside the slab='
     +       ,5x,2f11.6)
      IF (.NOT.film) RETURN
C
C Normalized pseudo density
C
      IF (.NOT.odi%d1) THEN
         qvac = 0.0
         DO ivac = 1,nvac
            CALL qsf(delz,rht(1,ivac,1),q2,nmz,0)
            q2(1) = q2(1)*area
            qvac = qvac + q2(1)*2./real(nvac)
         ENDDO
         qvac = qvac - 2*sigma
      ELSE
c-odim
         qvac = 0.0
         DO 105 nz = 1,nmz
            rht1(nz) = (z1+(nz-1)*delz)*
     *           rht(nz,nvac,1)
 105     CONTINUE
         CALL qsf(delz,rht1(1),q2,nmz,0)
         qvac = area*q2(1)
c+odim
      END IF 
!      fact = abs(qvac/psint)
!      DO k = 1,nq3
!         psq(k) = fact*psq(k)
!      ENDDO
      IF (l_xyav) RETURN
      fact = (qvac + psint)/(nstr(1)*vol)
      psq(1) = psq(1) - fact
c-gu
      WRITE (6,FMT=8010) fact*1000
      WRITE (16,FMT=8010) fact*1000
 8010 FORMAT (/,10x,'                     1000 * normalization const. ='
     +       ,5x,2f11.6)
c
      ENDIF ! irank == 0 
      END SUBROUTINE psqpw
      END MODULE m_psqpw
