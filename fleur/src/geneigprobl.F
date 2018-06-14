      MODULE m_geneigprobl
c**********************************************************
c     Solve the generalized eigenvalue problem
c     Frank Freimuth, November 2006
c**********************************************************
      CONTAINS
      SUBROUTINE geneigprobl(
     >                       nbasfcn,
     >                       nsize,neigd,lb,ub,l_J,a,b,
     <                       z,eig,ne)
#include"cpp_double.h"
      IMPLICIT NONE

! ... Arguments ...

      INTEGER, INTENT (IN) :: nbasfcn
      INTEGER, INTENT (IN) :: neigd
      INTEGER, INTENT (IN) :: nsize
      REAL,    INTENT (IN) :: lb,ub
      LOGICAL, INTENT (IN) :: l_J

      REAL,    INTENT(OUT) :: eig(:)
      INTEGER, INTENT(OUT) :: ne
#ifdef CPP_F90

#ifdef CPP_INVERSION
      REAL,  INTENT (INOUT) :: a(:),b(:)
      REAL,  INTENT (INOUT) :: z(:,:)
#else
      COMPLEX, INTENT (INOUT)::a(:),b(:)
      COMPLEX, INTENT (INOUT) :: z(:,:)
#endif

#else

#ifdef CPP_INVERSION
      REAL, ALLOCATABLE, INTENT (INOUT) :: a(:),b(:)
      REAL, ALLOCATABLE, INTENT (INOUT) :: z(:,:)
#else
      COMPLEX, ALLOCATABLE, INTENT (INOUT) :: a(:),b(:)
      COMPLEX, ALLOCATABLE, INTENT (INOUT) :: z(:,:)
#endif

#endif

! ... Local Variables ..

      INTEGER iind,ind1,ind2,info,lwork,liwork,lrwork,err
      REAL toler

      REAL,    ALLOCATABLE :: work(:)
      INTEGER, ALLOCATABLE :: iwork(:),isuppz(:)
#ifdef CPP_INVERSION
      REAL, ALLOCATABLE :: largea(:,:),largeb(:,:)
#else
      COMPLEX, ALLOCATABLE :: largea(:,:),largeb(:,:)
      COMPLEX,ALLOCATABLE :: cwork(:)
#endif

!**********************************
!expand from packed to full storage: full storage lapack-routines
!are faster than the packed lapack-routines.
!**********************************
!hamiltonian
      ALLOCATE ( largea(nsize,nsize), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating largea'
      iind = 0
      DO ind1 = 1, nsize
        DO ind2 = 1, ind1
           iind = iind+1
           largea(ind2,ind1) = a(iind)
        ENDDO
      ENDDO
!save some storage by deallocation of unused array
#ifndef CPP_F90
      DEALLOCATE (a)
#endif
!metric
      ALLOCATE ( largeb(nsize,nsize), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating largeb'
      iind=0
      DO ind1 = 1, nsize
        DO ind2 = 1, ind1
           iind = iind+1
           largeb(ind2,ind1) = b(iind)
        ENDDO
      ENDDO
!save some storage by deallocation of unused array
#ifndef CPP_F90
      DEALLOCATE (b)
#endif



#ifdef CPP_INVERSION
      CALL CPP_LAPACK_spotrf('U',nsize,largeb,nsize,info)
      IF (info/=0) STOP 'geneigprobl: error in spotrf'

      CALL CPP_LAPACK_ssygst(1,'U',nsize,largea,nsize,largeb,nsize,info)
      IF (info/=0) STOP 'geneigprobl: error in ssygst'

      toler = 2.0*tiny(toler)
      liwork = 10*nsize
      ALLOCATE ( iwork(liwork), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating iwork'

      lwork = 26*nsize
      ALLOCATE ( work(lwork), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating work'
      ALLOCATE ( isuppz(2*nsize), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating isuppz'
#ifndef CPP_F90
      ALLOCATE ( z(nbasfcn,neigd), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating z'
#endif

#ifndef CPP_F90
      IF (l_J) THEN
            CALL CPP_LAPACK_ssyevr('N','I','U',nsize,largea,
     &        nsize,lb,ub,1,neigd,toler,ne,eig,z,
     &        size(z,1),isuppz,work,lwork,iwork,liwork,info)
      ELSE
            CALL CPP_LAPACK_ssyevr('V','I','U',nsize,largea,
     &        nsize,lb,ub,1,neigd,toler,ne,eig,z,
     &        size(z,1),isuppz,work,lwork,iwork,liwork,info)
      ENDIF
#else
      eig = 0.0
#endif
      IF (info/=0) STOP 'geneigprobl: error in ssyevr'
      DEALLOCATE (isuppz,work,iwork)

      CALL CPP_LAPACK_strtrs('U','N','N',nsize,ne,largeb,
     &                         nsize,z,size(z,1),info)
      IF (info/=0) STOP 'geneigprobl: error in strtrs'
#else

      CALL CPP_LAPACK_cpotrf('U',nsize,largeb,nsize,info)
      IF (info/=0) STOP 'geneigprobl: error in cpotrf'

      CALL CPP_LAPACK_chegst(1,'U',nsize,largea,nsize,largeb,nsize,info)
      IF (info/=0) STOP 'geneigprobl: error in chegst'

      toler = 2.0*tiny(toler)
      liwork = 50*nsize
      ALLOCATE ( iwork(liwork), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating iwork'

      lwork = 20*nsize
      ALLOCATE( cwork(lwork), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating cwork'
      ALLOCATE( isuppz(10*nsize), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating isuppz'

      lrwork = 84*nsize
      ALLOCATE (work(lrwork), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating work'
#ifndef CPP_F90
      ALLOCATE (z(nbasfcn,neigd), stat=err )
      IF (err.NE.0) STOP 'geneigprobl: error allocating z'
#endif
#ifndef CPP_F90
      IF (l_J) THEN
         CALL CPP_LAPACK_cheevr('N','I','U',nsize,largea,
     &      nsize,lb,ub,1,neigd,toler,ne,eig,z,
     &      size(z,1),isuppz,cwork,lwork,work,lrwork,iwork,liwork,info)
      ELSE
         CALL CPP_LAPACK_cheevr('V','I','U',nsize,largea,
     &     nsize,lb,ub,1,neigd,toler,ne,eig,z,
     &     size(z,1),isuppz,cwork,lwork,work,lrwork,iwork,liwork,info)
      ENDIF
#else
      eig = 0.0
#endif
      IF (info/=0) STOP 'geneigprobl: error in cheevr'
      DEALLOCATE ( isuppz )
      deallocate ( work   )
      deallocate ( iwork  )
      deallocate ( cwork  )

      CALL CPP_LAPACK_ctrtrs('U','N','N',nsize,ne,largeb,
     &                         nsize,z,size(z,1),info)
      IF (info/=0) STOP 'geneigprobl: error in ctrtrs'
#endif
      DEALLOCATE ( largea,largeb )


      END SUBROUTINE geneigprobl
      END MODULE m_geneigprobl

