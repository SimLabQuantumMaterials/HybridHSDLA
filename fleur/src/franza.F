      SUBROUTINE franza(
     >                  nbasfcn,matsize,neigd,
     >                  nsize,lb,ub,
     >                  l_zref,l_J,matind,nred,
     >                  a,b,
     <                  z,eig,ne)

c***********************************************************************
c
c     solves the secular equation a*z=eig*b*z
c
c     nbase,matsize,nval  array dimensions as declared by calling routin
c                         unchanged on exit
c     nsize   actual dimension of a, b
c     a       on entry:   hamiltonian, lower triangle stored row-wise
c                              = upper triangle stored column-wise 
c             on exit :   destroyed
c     b       on entry:   overlap matrix, lower tr. stored row-wise
c             on exit :   destroyed(overwritten by cholesky factor)
c     eig     on entry:   ----
c             on exit :   eigenvalues in ascending rder
c     z       on exit :   corresponding eigenvectors
c     ne      on exit :   number of eigenvalues/vectors actually found
c     nblw    on exit :   number of eigenvalues below e(1)
c     work,iwork          temporary storage of dimensions
c                            nrtmp .ge. 8*nv
c                            nitmp .ge. 5*nv
c
c  This subroutine was written by Gustav Bihlmayer 
c                                 IFF, July 1996
c
c  and it is a LAPACK translation of the EISPACK version  
c  originally written by:
c         m. weinert     1982
c
c  arrays : iwork, work, ifail  are erroneously dimensioned
c           corrected by S.Bl"ugel             BNL 4.Aug.96
c
c***********************************************************************
c lapack routines inserted for the cray
c please note, that if too many ev's are in your interval, you'll get 
c a memory fault. Then change dspevx('V','V' ... to dspevx('V','A'
c to see what's going on .... 
c***********************************************************************
*

#include"cpp_double.h"
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN)    :: nbasfcn,matsize,neigd
      INTEGER, INTENT (INOUT) :: nsize
      INTEGER, INTENT (OUT)   :: ne
      INTEGER, INTENT (IN)    :: nred
      REAL,    INTENT (IN)    :: lb,ub
      LOGICAL, INTENT (IN)    :: l_zref,l_J
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: matind(nbasfcn,2)
      REAL,    INTENT (OUT) :: eig(neigd)
#ifdef CPP_INVERSION
      REAL,    INTENT (INOUT):: a(matsize),b(matsize)
      REAL,    INTENT (OUT) :: z(nbasfcn,neigd)
#else
      COMPLEX, INTENT (INOUT):: a(matsize),b(matsize)
      COMPLEX, INTENT (OUT) :: z(nbasfcn,neigd)
#endif
C     ..
C     .. Local Scalars ..
      REAL toler,sq2i
      INTEGER j,i1,i2,j1,j2,k,nsym,jsym,matsz
      INTEGER info,i,iu,ne_a,nsize_a,n
C     ..
C     .. Local Arrays
#ifdef CPP_INVERSION
      REAL,    ALLOCATABLE :: aa(:), bb(:)
#else
      COMPLEX, ALLOCATABLE :: aa(:), bb(:), cwork(:)
#endif
      REAL,    ALLOCATABLE ::  work(:), etemp(:)
      INTEGER, ALLOCATABLE ::  iwork(:),ifail(:)
      LOGICAL sort(nbasfcn)
C     ..
C     ..
C     .. External Subroutines ..
#ifdef CPP_INVERSION
      EXTERNAL CPP_LAPACK_spptrf,CPP_LAPACK_sspgst,CPP_LAPACK_sspevx,
     +         CPP_LAPACK_stptrs,CPP_BLAS_scopy
#else
      EXTERNAL CPP_LAPACK_cpptrf,CPP_LAPACK_chpgst,CPP_LAPACK_chpevx,
     +         CPP_LAPACK_ctptrs,CPP_BLAS_ccopy
#endif
C     ..
C     .. External Functions (causes problems with ifc)
!      REAL CPP_LAPACK_slamch 
C     ..
c+gu 
      nsym=1
#ifdef DISPLAY_INFO
      print*, 'FRANZA: matsize = ', matsize, 
     +        ', actual (nsize) size = ', nsize
#endif
      IF (l_zref) THEN
c
        ALLOCATE (aa(nred*(nred+1)/2),bb(nred*(nred+1)/2))
c
c separate H and S matrix in symmetric (aa,bb) and antisymmetric (a,b) part
c
        nsym=2
        i2=0
        j2=0
        DO i=1,nred
          DO j=1,i 
            i1=(matind(i,1)-1)*matind(i,1)/2+matind(j,1)
            j1=(matind(i,1)-1)*matind(i,1)/2+matind(j,2)
            i2=i2+1
            aa(i2)=a(i1)+a(j1)
            bb(i2)=b(i1)+b(j1)
            IF ((matind(i,1).NE.matind(i,2)).AND.
     +          (matind(j,1).NE.matind(j,2))) THEN
               j2=j2+1
               a(j2)=a(i1)-a(j1)
               b(j2)=b(i1)-b(j1)
            ENDIF
          ENDDO
        ENDDO
        nsize = nsize - nred
      ENDIF
      DO jsym=1,nsym

        IF (jsym.EQ.2) THEN
c
c second time calculate symmetric part and store the antisym. EV's and EW's
c
          matsz=(nred+1)*nred/2
#ifdef CPP_INVERSION
          CALL CPP_BLAS_scopy(matsz,aa,1,a,1)
          CALL CPP_BLAS_scopy(matsz,bb,1,b,1)
#else
          CALL CPP_BLAS_ccopy(matsz,aa,1,a,1)
          CALL CPP_BLAS_ccopy(matsz,bb,1,b,1)
#endif

          ne_a=ne
          k=1
          nsize_a=nsize
          DO i=1,ne
#ifdef CPP_INVERSION
            aa(i)=eig(i)
#else
            aa(i)=cmplx(eig(i),0.0)
#endif
            DO j=1,nsize
              bb(k)=z(j,i)
              k=k+1
            ENDDO
          ENDDO
          nsize=nred

        ENDIF
c-gu
C --> start with Cholesky factorization of b ( so that b = l * l^t)
C --> b is overwritten by l
c
#ifdef CPP_INVERSION
      CALL CPP_LAPACK_spptrf('U',nsize,b,info)
#else
      CALL CPP_LAPACK_cpptrf('U',nsize,b,info)
#endif
c
      IF (info.ne.0) THEN
#ifdef CPP_INVERSION
        WRITE (*,*) 'Error in spptrf: info =',info
#else
        WRITE (*,*) 'Error in cpptrf: info =',info
#endif
        STOP 'franza1'
      ENDIF
c
C --> now reduce a * z = eig * b * z to the standard form a' * z' = eig * z' 
C --> where a' = (l)^-1 * a * (l^t)^-1 and z' = l^t * z
c
#ifdef CPP_INVERSION
      CALL CPP_LAPACK_sspgst(1,'U',nsize,a,b,info)
#else
      CALL CPP_LAPACK_chpgst(1,'U',nsize,a,b,info)
#endif
c
      IF (info.ne.0) THEN
#ifdef CPP_INVERSION
        WRITE (6,*) 'Error in sspgst: info =',info
#else
        WRITE (6,*) 'Error in chpgst: info =',info
#endif
        STOP 'franza2'
      ENDIF
c
C --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub
c
      iu = neigd
      IF (l_zref) iu = (neigd-1)/2
!      toler=2.0*CPP_LAPACK_slamch('S')
      toler=2.0*TINY(sq2i)
      ALLOCATE ( work(8*nbasfcn),iwork(5*nbasfcn),ifail(nbasfcn) )
      ALLOCATE ( etemp(nbasfcn) )
#ifdef CPP_INVERSION
      IF(l_J)THEN
       CALL CPP_LAPACK_sspevx('N','I','U',nsize,a,lb,ub,1,iu,toler,ne,
     +             etemp,z,nbasfcn,work,iwork,ifail,info)
      ELSE
       CALL CPP_LAPACK_sspevx('V','I','U',nsize,a,lb,ub,1,iu,toler,ne,
     +             etemp,z,nbasfcn,work,iwork,ifail,info)
      ENDIF
#else
       ALLOCATE ( cwork(2*nbasfcn) )
      IF(l_J)THEN
       CALL CPP_LAPACK_chpevx('N','I','U',nsize,a,lb,ub,1,iu,toler,ne,
     +             etemp,z,nbasfcn,cwork,work,iwork,ifail,info)
      ELSE
       CALL CPP_LAPACK_chpevx('V','I','U',nsize,a,lb,ub,1,iu,toler,ne,
     +             etemp,z,nbasfcn,cwork,work,iwork,ifail,info)
      ENDIF
       DEALLOCATE ( cwork )
#endif
c
      IF (ne.gt.neigd) THEN
        WRITE(6,*) 'ne=',ne,' > neigd'
        STOP 'franza : ne.GT.neigd'
      ENDIF
      eig(:) = etemp(1:neigd)
c
      IF(l_J) THEN
      DEALLOCATE ( work,iwork,ifail,etemp )
      ELSE
      IF (info.ne.0) THEN
#ifdef CPP_INVERSION
        WRITE (6,*) 'Error in sspevx: info =',info
#else
        WRITE (6,*) 'Error in chpevx: info =',info
#endif
        WRITE (6,*) 'The following eigenvectors did not converge:'
        WRITE (6,'(30i5)') (ifail(i),i=1,ne)
        STOP 'franza3'
      ENDIF               
      DEALLOCATE ( work,iwork,ifail,etemp )
c
C --> recover the generalized eigenvectors z by solving z' = l^t * z
c
#ifdef CPP_INVERSION
        CALL CPP_LAPACK_stptrs('U','N','N',nsize,ne,b,z,nbasfcn,info)
#else
        CALL CPP_LAPACK_ctptrs('U','N','N',nsize,ne,b,z,nbasfcn,info)
#endif
c
        IF (info.ne.0) THEN
#ifdef CPP_INVERSION
          WRITE (6,*) 'Error in stptrs: info =',info
#else
          WRITE (6,*) 'Error in ctptrs: info =',info
#endif
          STOP 'franza4'
        ENDIF
      ENDIF !l_J
c+gu
      ENDDO
c
c now collect symmetric and antisym. EW's and EV's and sort
c     
      IF (l_zref) THEN

          k=1
          DO i=1,ne
            DO j=1,nsize
              b(k)=z(j,i)
              k=k+1
            ENDDO
          ENDDO
c
c prepare sort-array: even=.true., odd=.false.
c
        i=1
        j=1
        eig(ne+1)    = 99.9e9
#ifdef CPP_INVERSION
        aa(ne_a+1) = 99.9e9
#else
        aa(ne_a+1) = cmplx(99.9e9,0.0)
#endif
        DO k=1,ne+ne_a
#ifdef CPP_INVERSION
          IF (eig(i).lt.aa(j)) THEN
#else
          IF (eig(i).lt.real(aa(j))) THEN
#endif
            sort(k)=.true.
            i=i+1
          ELSE
            sort(k)=.false.
            j=j+1
          ENDIF
        ENDDO
c
c sort EW's and EV's
c
        i=ne
        j=ne_a
        nsize = nsize + nsize_a
        ne = ne + ne_a
        sq2i=1.0/sqrt(2.0)
        DO k=ne,1,-1
          DO n=1,nsize
#ifdef CPP_INVERSION
            z(n,k)=0.0
#else
            z(n,k)=cmplx(0.0,0.0)
#endif
          ENDDO
          IF (sort(k)) THEN
            eig(k)=eig(i)
            i1=nred * (i-1)
            DO n=1,nred
              i1=i1+1
#ifdef CPP_INVERSION
              z(matind(n,1),k) = z(matind(n,1),k)+b(i1)*sq2i
              z(matind(n,2),k) = z(matind(n,2),k)+b(i1)*sq2i
#else
              z(matind(n,1),k) = z(matind(n,1),k)+b(i1)*cmplx(sq2i,0.0)
              z(matind(n,2),k) = z(matind(n,2),k)+b(i1)*cmplx(sq2i,0.0)
#endif
            ENDDO
            i=i-1
          ELSE
#ifdef CPP_INVERSION
            eig(k)=aa(j)
#else
            eig(k)=real(aa(j))
#endif
            j1=nsize_a * (j-1)
            DO n=1,nred
              IF (matind(n,1).NE.matind(n,2)) THEN
                j1=j1+1
#ifdef CPP_INVERSION
                z(matind(n,1),k) = z(matind(n,1),k)+bb(j1)*sq2i
                z(matind(n,2),k) = z(matind(n,2),k)-bb(j1)*sq2i
#else
                z(matind(n,1),k) = z(matind(n,1),k)+bb(j1)
     *                             *cmplx(sq2i,0.0)
                z(matind(n,2),k) = z(matind(n,2),k)-bb(j1)
     *                             *cmplx(sq2i,0.0)
#endif
              ELSE
#ifdef CPP_INVERSION
                z(matind(n,1),k) = 0.0
#else
                z(matind(n,1),k) = cmplx(0.0,0.0)
#endif
              ENDIF
            ENDDO
            j=j-1
          ENDIF
        ENDDO
        DEALLOCATE ( aa,bb )

      ENDIF
c-gu
      END
