      MODULE m_orthoglo
c*********************************************************************
c Each G-vector corresponds to a vector of C-coeff. These vectors must
c be linearly independent. This is checked by this soubroutine for an
c atom that doesn't have an inversion partner.
c Philipp Kurz 99/04
c*********************************************************************
      CONTAINS
      SUBROUTINE orthoglo(
     >                    llod,nlod,
     >                    nkvec,lo,l,linindq,l_lo2,
     X                    cwork,
     <                    linind)
c
c*************** ABBREVIATIONS ***************************************
c cwork   : contains the vectors of C-coeff.
c l_lo2   : changes this routine to old 'orthoglo2': same as orthoglo, 
c           but for a pair of atoms that can be mapped onto eachother 
c           by inversion.
c CF Replaced (unstable) Gram-Schmidt by diagonalization.
c*********************************************************************
c
#include"cpp_double.h"
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: llod,nlod
      INTEGER, INTENT (IN) :: l,lo,nkvec
      REAL,    INTENT (IN) :: linindq
      LOGICAL, INTENT (IN) :: l_lo2
      LOGICAL, INTENT (OUT) :: linind
C     ..
C     .. Array Arguments ..
      COMPLEX,INTENT (INOUT):: cwork(-2*llod:2*llod+1,2*(2*llod+1),nlod)
C     ..
C     .. Local Scalars ..
      INTEGER dim,low,i,j
C     ..
C     .. Local Arrays ..
      REAL eig(nkvec),rwork(3*nkvec)
#ifdef CPP_INVERSION
      REAL olap(nkvec,nkvec)
      EXTERNAL CPP_LAPACK_ssyev
#else
      COMPLEX olap(nkvec,nkvec),work(2*nkvec)
      EXTERNAL CPP_LAPACK_cheev
#endif
      
      IF (l_lo2) THEN
        dim = 2* (2*l+1)
        low = -2*l
      ELSE
        dim = 2*l+1
        low = -l
      ENDIF

      DO i = 1,nkvec
        DO j = 1,nkvec
          olap(i,j) = dot_product(
     &      cwork(low:low+dim-1,i,lo),
     &      cwork(low:low+dim-1,j,lo))
        ENDDO
      ENDDO
#ifdef CPP_INVERSION
      CALL CPP_LAPACK_ssyev('N','U',nkvec,olap,nkvec,eig,
     &                                        rwork,3*nkvec,i)
      IF(i.NE.0) STOP 'orthoglo: (S,D)SYEV failed.'
#else
      CALL CPP_LAPACK_cheev('N','U',nkvec,olap,nkvec,eig,
     +                                   work,2*nkvec,rwork,i)
      IF(i.NE.0) STOP 'orthoglo: (C,Z)HEEV failed.'
#endif
      IF(eig(1).LT.linindq) THEN
        linind = .false.
      ELSE
        linind = .true.
      ENDIF
      RETURN
      
      END SUBROUTINE orthoglo
      END MODULE m_orthoglo
