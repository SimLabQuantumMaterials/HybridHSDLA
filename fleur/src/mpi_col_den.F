      MODULE m_mpi_col_den
c
c collect all data calculated in cdnval on different pe's on pe 0
c
      CONTAINS
      SUBROUTINE mpi_col_den(
     >                    irank,isize,ntypsd,ntype,nlod,natd,llod,odi,
     >                       lmaxd,ntypd,nlhd,n3d,n2d,nmzd,nmzxyd,
     >                       jspins,nq2,nq3,nmz,nmzxy,lmax,nlh,ntypsy,
     >                       film,l_f,l_soc,l_mperp,l_fmpl,
     >                       jspin,nw,nlo,neq,n_u,llpd,
     X                       rhtxy,rht,qpw,ener,sqal,force,svac,pvac,
     X                       uu,dd,du,uunmt,ddnmt,udnmt,dunmt,sqlo,
     X                       aclo,bclo,cclo,acnmt,bcnmt,ccnmt,enerlo,
     X                       orb,orbl,orblo,mt21,lo21,uloulop21,
     X                       uunmt21,ddnmt21,udnmt21,dunmt21,
     X                       cdom,cdomvz,cdomvxy,l_noco,n_mmp)
c
#include"cpp_double.h"
      USE m_types, ONLY : t_orb,t_orbl,t_orblo,t_mt21,t_lo21
      USE m_od_types, ONLY : od_inp

      IMPLICIT NONE
      INCLUDE 'mpif.h'
C ..
C ..  Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,ntypd,n3d,n2d,nmzd,nmzxyd,nlhd
      INTEGER, INTENT (IN) :: nq3,nmz,nmzxy,irank,isize,jspins,nq2,nlod
      INTEGER, INTENT (IN) :: ntypsd,ntype,jspin,nw,natd,n_u,llod,llpd
      LOGICAL, INTENT (IN) :: film,l_f,l_soc,l_noco,l_mperp,l_fmpl
C ..
C ..  Array Arguments ..
      TYPE (od_inp), INTENT (IN) :: odi
      INTEGER, INTENT (IN) :: lmax(ntypd),nlh(ntypsd),ntypsy(natd)
      INTEGER, INTENT (IN) :: nlo(ntypd),neq(ntypd)
      COMPLEX, INTENT (INOUT) :: qpw(n3d)
      COMPLEX, INTENT (INOUT) :: rhtxy(nmzxyd,odi%n2d-1,2)
      REAL,    INTENT (INOUT) :: rht(nmzd,2),force(3,ntypd)
      REAL,    INTENT (INOUT) :: ener(0:3,ntypd),sqal(0:3,ntypd)
      REAL,    INTENT (INOUT) :: svac(2),pvac(2)
      REAL,  INTENT (INOUT) :: dd(0:lmaxd,ntypd)
      REAL,  INTENT (INOUT) :: du(0:lmaxd,ntypd)
      REAL,  INTENT (INOUT) :: uu(0:lmaxd,ntypd)
      REAL,  INTENT (INOUT) :: ddnmt(0:llpd,nlhd,ntypd)
      REAL,  INTENT (INOUT) :: dunmt(0:llpd,nlhd,ntypd)
      REAL,  INTENT (INOUT) :: udnmt(0:llpd,nlhd,ntypd)
      REAL,  INTENT (INOUT) :: uunmt(0:llpd,nlhd,ntypd)
      REAL,  INTENT (INOUT) :: sqlo(nlod,ntypd),enerlo(nlod,ntypd)
      REAL,  INTENT (INOUT) :: aclo(nlod,ntypd),bclo(nlod,ntypd)
      REAL,  INTENT (INOUT) :: cclo(nlod,nlod,ntypd)
      REAL,  INTENT (INOUT) :: acnmt(0:lmaxd,nlod,nlhd,ntypd)
      REAL,  INTENT (INOUT) :: bcnmt(0:lmaxd,nlod,nlhd,ntypd)
      REAL,  INTENT (INOUT) :: ccnmt(nlod,nlod,nlhd,ntypd)
      COMPLEX,INTENT(INOUT) :: ddnmt21((lmaxd+1)**2,nlhd,ntypd)
      COMPLEX,INTENT(INOUT) :: dunmt21((lmaxd+1)**2,nlhd,ntypd)
      COMPLEX,INTENT(INOUT) :: udnmt21((lmaxd+1)**2,nlhd,ntypd)
      COMPLEX,INTENT(INOUT) :: uunmt21((lmaxd+1)**2,nlhd,ntypd)
      COMPLEX,INTENT(INOUT) :: uloulop21(nlod,nlod,ntypd)
      COMPLEX,INTENT(INOUT) :: n_mmp(-3:3,-3:3,n_u),cdomvz(nmzd,2)
      COMPLEX,INTENT(INOUT) :: cdom(n3d),cdomvxy(nmzxyd,odi%n2d-1,2)
      TYPE (t_orb),  INTENT (INOUT) :: orb(0:lmaxd,-lmaxd:lmaxd,ntypd)
      TYPE (t_orbl), INTENT (INOUT) :: orbl(nlod,-llod:llod,ntypd)
      TYPE (t_orblo),INTENT (INOUT) :: orblo(nlod,nlod,-llod:llod,ntypd)
      TYPE (t_mt21), INTENT (INOUT) :: mt21(0:lmaxd,ntypd)
      TYPE (t_lo21), INTENT (INOUT) :: lo21(nlod,ntypd)
C ..
C ..  Local Scalars ..
      INTEGER :: n
C ..
C ..  Local Arrays ..
      INTEGER :: ierr(3)
      COMPLEX, ALLOCATABLE :: c_b(:)
      REAL,    ALLOCATABLE :: r_b(:)
C ..
C ..  External Subroutines
      EXTERNAL CPP_BLAS_scopy,CPP_BLAS_ccopy,MPI_REDUCE

c
c -> Collect qpw()
c
      n = n3d
      ALLOCATE(c_b(n))
      CALL MPI_REDUCE(qpw,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0,
     +                              MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_ccopy(n, c_b, 1, qpw, 1)
      ENDIF
      DEALLOCATE (c_b)
c
c -> Collect rhtxy()
c
      IF (film) THEN

      n = nmzxyd*(odi%n2d-1)*2
      ALLOCATE(c_b(n))
      CALL MPI_REDUCE(rhtxy,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0,
     +                                MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_ccopy(n, c_b, 1, rhtxy, 1)
      ENDIF
      DEALLOCATE (c_b)
c
c -> Collect rht()
c
      n = nmzd*2
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(rht,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                       MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, rht, 1)
      ENDIF
      DEALLOCATE (r_b)

      ENDIF
c
c -> Collect uu(),ud() and dd()
c
      n = (lmaxd+1)*ntypd
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(uu,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                      MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, uu, 1)
      ENDIF
      CALL MPI_REDUCE(du,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                      MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, du, 1)
      ENDIF
      CALL MPI_REDUCE(dd,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                      MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, dd, 1)
      ENDIF
      DEALLOCATE (r_b)
c
c--> Collect uunmt,udnmt,dunmt,ddnmt
c
      n = (llpd+1)*nlhd*ntypd
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(uunmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                         MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, uunmt, 1)
      ENDIF
      CALL MPI_REDUCE(udnmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                         MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, udnmt, 1)
      ENDIF
      CALL MPI_REDUCE(dunmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                         MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, dunmt, 1)
      ENDIF
      CALL MPI_REDUCE(ddnmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                         MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, ddnmt, 1)
      ENDIF
      DEALLOCATE (r_b)
c
c--> ener & sqal
c
      n=4*ntypd
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(ener,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                        MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, ener, 1)
      ENDIF
      CALL MPI_REDUCE(sqal,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                        MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, sqal, 1)
      ENDIF
      DEALLOCATE (r_b)
c
c--> svac & pvac
c
      IF ( film ) THEN

      n=2
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(svac,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                        MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, svac, 1)
      ENDIF
      CALL MPI_REDUCE(pvac,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                        MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, pvac, 1)
      ENDIF
      DEALLOCATE (r_b)

      ENDIF
c
c -> Collect force
c   
      IF (l_f) THEN

      n=3*ntypd
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(force,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                         MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, force, 1)
      ENDIF
      DEALLOCATE (r_b)

      ENDIF
c
c -> Optional the LO-coefficients: aclo,bclo,enerlo,cclo,acnmt,bcnmt,ccnmt
c
      IF (nlod.GE.1) THEN

      n=nlod*ntypd 
      ALLOCATE (r_b(n))
      CALL MPI_REDUCE(aclo,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                        MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, aclo, 1)
      ENDIF
      CALL MPI_REDUCE(bclo,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                        MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, bclo, 1)
      ENDIF
      CALL MPI_REDUCE(enerlo,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                          MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, enerlo, 1)
      ENDIF
      CALL MPI_REDUCE(sqlo,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                        MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, sqlo, 1)
      ENDIF
      DEALLOCATE (r_b)

      n = nlod * nlod * ntypd
      ALLOCATE (r_b(n))
      CALL MPI_REDUCE(cclo,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                        MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, cclo, 1)
      ENDIF
      DEALLOCATE (r_b)

      n = (lmaxd+1) * ntypd * nlod * nlhd
      ALLOCATE (r_b(n))
      CALL MPI_REDUCE(acnmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                         MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, acnmt, 1)
      ENDIF
      CALL MPI_REDUCE(bcnmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                         MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, bcnmt, 1)
      ENDIF
      DEALLOCATE (r_b)

      n = ntypd * nlhd * nlod**2
      ALLOCATE (r_b(n))
      CALL MPI_REDUCE(ccnmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                         MPI_COMM_WORLD,ierr)
      IF (irank.EQ.0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, ccnmt, 1)
      ENDIF
      DEALLOCATE (r_b)

      ENDIF
c
c ->  Now the SOC - stuff: orb, orblo and orblo
c
      IF (l_soc) THEN
!
! orb
        n=(lmaxd+1)*(2*lmaxd+1)*ntypd
        ALLOCATE (r_b(n))
        CALL MPI_REDUCE(orb(:,:,:)%uu,r_b,n,CPP_MPI_REAL,
     +                               MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, orb(:,:,:)%uu, 1)
        ENDIF
        CALL MPI_REDUCE(orb(:,:,:)%dd,r_b,n,CPP_MPI_REAL,
     +                               MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, orb(:,:,:)%dd, 1)
        ENDIF
        DEALLOCATE (r_b)

        ALLOCATE (c_b(n))
        CALL MPI_REDUCE(orb(:,:,:)%uup,c_b,n,CPP_MPI_COMPLEX,
     +                                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb(:,:,:)%uup, 1)
        ENDIF
        CALL MPI_REDUCE(orb(:,:,:)%ddp,c_b,n,CPP_MPI_COMPLEX,
     +                                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb(:,:,:)%ddp, 1)
        ENDIF
        CALL MPI_REDUCE(orb(:,:,:)%uum,c_b,n,CPP_MPI_COMPLEX,
     +                                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb(:,:,:)%uum, 1)
        ENDIF
        CALL MPI_REDUCE(orb(:,:,:)%ddm,c_b,n,CPP_MPI_COMPLEX,
     +                                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb(:,:,:)%ddm, 1)
        ENDIF
        DEALLOCATE (c_b)
! orbl
!
        n = nlod * (2*llod+1) * ntypd
        ALLOCATE (r_b(n))
       CALL MPI_REDUCE(orbl(:,:,:)%uulo,r_b,n,CPP_MPI_REAL,
     +                               MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, orbl(:,:,:)%uulo, 1)
        ENDIF
       CALL MPI_REDUCE(orbl(:,:,:)%dulo,r_b,n,CPP_MPI_REAL,
     +                               MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, orbl(:,:,:)%dulo, 1)
        ENDIF
        DEALLOCATE (r_b)

        ALLOCATE (c_b(n))
        CALL MPI_REDUCE(orbl(:,:,:)%uulop,c_b,n,CPP_MPI_COMPLEX,
     +                                 MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orbl(:,:,:)%uulop, 1)
        ENDIF
        CALL MPI_REDUCE(orbl(:,:,:)%dulop,c_b,n,CPP_MPI_COMPLEX,
     +                                 MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orbl(:,:,:)%dulop, 1)
        ENDIF
        CALL MPI_REDUCE(orbl(:,:,:)%uulom,c_b,n,CPP_MPI_COMPLEX,
     +                                 MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orbl(:,:,:)%uulom, 1)
        ENDIF
        CALL MPI_REDUCE(orbl(:,:,:)%dulom,c_b,n,CPP_MPI_COMPLEX,
     +                                 MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orbl(:,:,:)%dulom, 1)
        ENDIF
        DEALLOCATE (c_b)
! orblo
!
        n = nlod * nlod * (2*llod+1) * ntypd
        ALLOCATE (r_b(n))
       CALL MPI_REDUCE(orblo(:,:,:,:)%z,r_b,n,CPP_MPI_REAL,
     +                               MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, orblo(:,:,:,:)%z, 1)
        ENDIF
        DEALLOCATE (r_b)

        ALLOCATE (c_b(n))
        CALL MPI_REDUCE(orblo(:,:,:,:)%p,c_b,n,CPP_MPI_COMPLEX,
     +                                   MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orblo(:,:,:,:)%p, 1)
        ENDIF
        CALL MPI_REDUCE(orblo(:,:,:,:)%m,c_b,n,CPP_MPI_COMPLEX,
     +                                   MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orblo(:,:,:,:)%m, 1)
        ENDIF
        DEALLOCATE (c_b)

      ENDIF

c
c -> Collect the noco staff: 
c
      IF ( l_noco .AND. jspin.EQ.1 ) THEN

        n = n3d
        ALLOCATE(c_b(n))
        CALL MPI_REDUCE(cdom,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0,
     +                                 MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, cdom, 1)
        ENDIF
        DEALLOCATE (c_b)
c
        IF (film) THEN

          n = nmzxyd*(odi%n2d-1)*2
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(cdomvxy,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0,
     +                                      MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, cdomvxy, 1)
          ENDIF
          DEALLOCATE (c_b)
c
          n = nmzd*2
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(cdomvz,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0,
     +                                     MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, cdomvz, 1)
          ENDIF
          DEALLOCATE (c_b)

        ENDIF ! film
 

        IF (l_mperp) THEN
!
! -->     for (spin)-off diagonal part of muffin-tin
!
          n = (lmaxd+1) * ntypd
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(mt21(:,:)%uu,c_b,n,CPP_MPI_COMPLEX,
     +                                MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, mt21(:,:)%uu, 1)
          ENDIF
          CALL MPI_REDUCE(mt21(:,:)%ud,c_b,n,CPP_MPI_COMPLEX,
     +                                MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, mt21(:,:)%ud, 1)
          ENDIF
          CALL MPI_REDUCE(mt21(:,:)%du,c_b,n,CPP_MPI_COMPLEX,
     +                                MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, mt21(:,:)%du, 1)
          ENDIF
          CALL MPI_REDUCE(mt21(:,:)%dd,c_b,n,CPP_MPI_COMPLEX,
     +                                MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, mt21(:,:)%dd, 1)
          ENDIF
          DEALLOCATE (c_b)
!
! -->     lo,u coeff's:
!
          n = nlod * ntypd
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(lo21(:,:)%uulo,c_b,n,CPP_MPI_COMPLEX,
     +                                MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, lo21(:,:)%uulo, 1)
          ENDIF
          CALL MPI_REDUCE(lo21(:,:)%ulou,c_b,n,CPP_MPI_COMPLEX,
     +                                MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, lo21(:,:)%ulou, 1)
          ENDIF
          CALL MPI_REDUCE(lo21(:,:)%dulo,c_b,n,CPP_MPI_COMPLEX,
     +                                MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, lo21(:,:)%dulo, 1)
          ENDIF
          CALL MPI_REDUCE(lo21(:,:)%ulod,c_b,n,CPP_MPI_COMPLEX,
     +                                MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, lo21(:,:)%ulod, 1)
          ENDIF
          DEALLOCATE (c_b)
!
! -->     lo,lo' coeff's:
!
          n = nlod*nlod*ntypd
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(uloulop21,c_b,n,CPP_MPI_COMPLEX,
     +                                MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (irank.EQ.0) THEN
            CALL CPP_BLAS_ccopy(n, c_b, 1, uloulop21, 1)
          ENDIF
          DEALLOCATE (c_b)

          IF (l_fmpl) THEN
c
c-->        Full magnetization plots: Collect uunmt21, etc.
c
            n = (lmaxd+1)**2 *nlhd*ntypd
            ALLOCATE(c_b(n))
            CALL MPI_REDUCE(uunmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0
     +                                        ,MPI_COMM_WORLD,ierr)
            IF (irank.EQ.0) THEN
              CALL CPP_BLAS_ccopy(n, c_b, 1, uunmt21, 1)
            ENDIF
            CALL MPI_REDUCE(udnmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0
     +                                        ,MPI_COMM_WORLD,ierr)
            IF (irank.EQ.0) THEN
              CALL CPP_BLAS_ccopy(n, c_b, 1, udnmt21, 1)
            ENDIF
            CALL MPI_REDUCE(dunmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0
     +                                        ,MPI_COMM_WORLD,ierr)
            IF (irank.EQ.0) THEN
              CALL CPP_BLAS_ccopy(n, c_b, 1, dunmt21, 1)
            ENDIF
            CALL MPI_REDUCE(ddnmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0
     +                                        ,MPI_COMM_WORLD,ierr)
            IF (irank.EQ.0) THEN
              CALL CPP_BLAS_ccopy(n, c_b, 1, ddnmt21, 1)
            ENDIF
            DEALLOCATE (c_b)

          ENDIF ! fmpl
        ENDIF  ! mperp
      ENDIF   ! noco

c+lda+U
      IF ( n_u.GT.0 ) THEN
        n = 49*n_u 
        ALLOCATE(c_b(n))
        CALL MPI_REDUCE(n_mmp,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0,
     +                                  MPI_COMM_WORLD,ierr)
        IF (irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, n_mmp, 1)
        ENDIF
        DEALLOCATE (c_b)
      ENDIF
c-lda+U

      RETURN
      END SUBROUTINE mpi_col_den
      END MODULE m_mpi_col_den
