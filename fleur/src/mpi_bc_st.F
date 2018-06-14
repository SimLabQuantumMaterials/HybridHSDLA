      MODULE m_mpi_bc_st
c**********************************************************************
!     mpi_bc_st :  broadcast all information for qpw_to_nmt
!     mpi_col_st:  collect the density from pe's 
c**********************************************************************
      CONTAINS
      SUBROUTINE mpi_bc_st(
     >                     irank,memd,nlhd,ntypsd,jmtd,ntypd,n3d,
     X                     jspin,l_cutoff,ng3,ntype,nop,natd,symor,
     X                     clnu,qpwc,lmax,ntypsy,jri,nmem,nlh,mlh,nstr,
     X                     neq,kv3,mrot,invtab,llh,bmat,tau,taual,dx,
     X                     rmsh,sk3)
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: irank
      INTEGER, INTENT (IN) :: memd,nlhd,ntypsd,jmtd,ntypd,n3d,natd
      INTEGER :: jspin,l_cutoff,ng3,ntype,nop
      LOGICAL :: symor
!      INTEGER, INTENT (INOUT) :: jspin,l_cutoff,ng3,ntype,nop
!      LOGICAL, INTENT (INOUT) :: symor
C     ..
C     .. Array Arguments ..
      COMPLEX :: clnu(memd,0:nlhd,ntypsd),qpwc(n3d)
      INTEGER :: lmax(ntypd),ntypsy(natd),jri(ntypd)
      INTEGER :: nmem(0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER :: mlh(memd,0:nlhd,ntypsd),nstr(n3d)
      INTEGER :: neq(ntypd),kv3(3,n3d),mrot(3,3,nop)
      INTEGER :: invtab(nop),llh(0:nlhd,ntypsd)
      REAL    :: bmat(3,3),tau(3,nop),taual(3,natd)
      REAL    :: dx(ntypd),rmsh(jmtd,ntypd),sk3(n3d)
!      COMPLEX, INTENT (INOUT) :: clnu(memd,0:nlhd,ntypsd),qpwc(n3d)
!      INTEGER, INTENT (INOUT) :: lmax(ntypd),ntypsy(natd),jri(ntypd)
!      INTEGER, INTENT (INOUT) :: nmem(0:nlhd,ntypsd),nlh(ntypsd)
!      INTEGER, INTENT (INOUT) :: mlh(memd,0:nlhd,ntypsd),nstr(n3d)
!      INTEGER, INTENT (INOUT) :: neq(ntypd),kv3(3,n3d),mrot(3,3,nop)
!      INTEGER, INTENT (INOUT) :: invtab(nop),llh(0:nlhd,ntypsd)
!      REAL,    INTENT (INOUT) :: bmat(3,3),tau(3,nop),taual(3,natd)
!      REAL,    INTENT (INOUT) :: dx(ntypd),rmsh(jmtd,ntypd),sk3(n3d)
C     ..
C     .. Local Scalars ..
      INTEGER n,TYP_REAL,TYP_COMPLEX
C     ..
C     .. Local Arrays ..
      INTEGER ierr(3)
C     ..
C     .. External Subroutines.. 
      EXTERNAL MPI_BCAST,MPI_TYPE_CONTIGUOUS,MPI_TYPE_COMMIT
C     ..
      INCLUDE 'mpif.h'
c
#ifdef CPP_T90
      CALL MPI_TYPE_CONTIGUOUS(1,MPI_REAL,TYP_REAL,ierr)
      CALL MPI_TYPE_CONTIGUOUS(1,MPI_COMPLEX,TYP_COMPLEX,ierr)
#else
      CALL MPI_TYPE_CONTIGUOUS(1,MPI_DOUBLE_PRECISION,TYP_REAL,ierr)
      CALL MPI_TYPE_CONTIGUOUS(1,MPI_DOUBLE_COMPLEX,TYP_COMPLEX,ierr)
#endif
      CALL MPI_TYPE_COMMIT(TYP_REAL,ierr)
      CALL MPI_TYPE_COMMIT(TYP_COMPLEX,ierr)

      CALL MPI_BCAST(jspin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ntype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(l_cutoff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ng3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(symor,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c
c -> Broadcast the arrays:
 
      CALL MPI_BCAST(qpwc,n3d,TYP_COMPLEX,0,MPI_COMM_WORLD,ierr)
      n = memd*(nlhd+1)*ntypsd
      CALL MPI_BCAST(clnu,n,TYP_COMPLEX,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mlh,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nlh,ntypsd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      n = (nlhd+1)*ntypsd
      CALL MPI_BCAST(nmem,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(llh,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(jri,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ntypsy,natd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(neq,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmax,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mrot,9*nop,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dx,ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sk3,n3d,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bmat,9,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(taual,3*natd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tau,3*nop,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rmsh,jmtd*ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(invtab,nop,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kv3,3*n3d,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nstr,n3d,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      END SUBROUTINE mpi_bc_st
c*********************************************************************
      SUBROUTINE mpi_col_st(
     >                      jmtd,nlhd,ntypd,irank,
     X                      rho)
c
#include"cpp_double.h"
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      EXTERNAL MPI_REDUCE,CPP_BLAS_scopy
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jmtd,nlhd,ntypd,irank
      REAL, INTENT (INOUT) :: rho(jmtd,0:nlhd,ntypd)

      INTEGER n
      INTEGER ierr(3)
      REAL, ALLOCATABLE :: r_b(:)

      n = jmtd*(nlhd+1)*ntypd
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(rho,r_b,n,CPP_MPI_REAL,MPI_SUM,0,
     +                                       MPI_COMM_WORLD,ierr)
      IF (irank == 0) THEN
        CALL CPP_BLAS_scopy(n, r_b, 1, rho, 1)
      ENDIF
      DEALLOCATE(r_b) 

      END SUBROUTINE mpi_col_st
c*********************************************************************
      END MODULE m_mpi_bc_st
