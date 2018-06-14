      MODULE m_radfun
      CONTAINS
      SUBROUTINE radfun(
     >                  l,e,vr,jri,r0,dx,jmtd,
     <                  f,g,us,dus,uds,duds,ddn,nodeu,noded,wronk)
c*********************************************************************
c     generates the scalar relativistic wavefunctions (function: f;
c     energy derivative: g) at an energy e for angular momentum l.
c     the values on the sphere boundaries are also returned.
c             m. weinert   jan. 1987
c     the solutions r*u(r) are on a log. mesh.
c
c      us ... u(R) if R is the muffin tin Radius
c     dus ... u'(R)   (radial derivative)
c             .               .
c     uds ... u(R)   duds ... u'(R)  (energy derivative)
c              . .
c     ddn ... <u|u>  norm of u-dot
c
c*********************************************************************

      USE m_constants, ONLY : c_light
      USE m_radsra
      USE m_radsrd
      USE binmat

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jmtd,jri,l
      INTEGER, INTENT (OUT):: noded,nodeu
      REAL,    INTENT (IN) :: dx,e,r0
      REAL,    INTENT (OUT):: ddn,duds,dus,uds,us,wronk
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN) :: vr(jmtd)
      REAL,    INTENT (OUT):: f(jmtd,2),g(jmtd,2)
C     ..
C     Local scalars
      REAL c
C     ..
#ifdef DUMP_DATA_OLD
      character(len=50) :: vr_fname
      ! dump vr  and some parameters
!     print*, 'In RADFUN: r0:', r0, ', dx:', dx, ', e:', e,
!    &                   'jri:', jri, ', l:', l
      write(vr_fname, '(A, I1, ".bin")') 'radpot_l-', l
      call write_mat(vr, vr_fname)
#endif
      IF (jri.GT.jmtd) STOP 'radfun'
c
      c = c_light(1.0)
c
c--->    calculate normalized function at e
      CALL radsra(
     >            e,l,vr,r0,dx,jri,jmtd,c,
     <            us,dus,nodeu,f(1,1),f(1,2))
c
c--->    calculate orthogonal energy derivative at e
      CALL radsrd(
     >            e,l,vr,r0,dx,jri,jmtd,c,
     <            uds,duds,ddn,noded,g(1,1),g(1,2),
     >            f(1,1),f(1,2),dus)
c
c--->    calculate wronskian
      wronk = uds*dus - duds*us

      END SUBROUTINE radfun
      END MODULE m_radfun
