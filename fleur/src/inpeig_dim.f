      MODULE m_inpeigdim
c*********************************************************************
c     inputs the necessary quantities for the eigenvalue part (energy
c     parameters, k-points, wavefunction cutoffs, etc.).
c                  m. weinert   jan. 1987
c*********************************************************************
      CONTAINS
      SUBROUTINE inpeig_dim(
     >         rkmax,nwdd,latnam,amat,bmat,film,l_ss,
     X         qss,
     >         odd,l_J,
     <         nkptd,nvd,nv2d,kw1d,kw2d,kw3d)

      USE m_constants, ONLY : pimach
      USE m_dotset
      USE m_od_types, ONLY : od_dim

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: nwdd
      REAL,    INTENT (IN) :: rkmax
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3)
      REAL,    INTENT (INOUT) :: qss(3)
      LOGICAL, INTENT (IN) :: film,l_ss,l_J
      CHARACTER*3, INTENT (IN) :: latnam
      INTEGER, INTENT (OUT) :: nkptd,nvd,nv2d,kw1d,kw2d,kw3d
c-odim
      TYPE (od_dim), INTENT (IN) :: odd
c-odim
      INTEGER nw,nwd,nk,nq,i,nkpt,nqpt,nv,nv2,j,kw1,kw2,kw3
      REAL pi,tpi,scale,s1,s2
      LOGICAL xyu
      REAL aamat(3,3),bbmat(3,3),bk(3)

      EXTERNAL apws_dim
C     ..
      pi = pimach()
      tpi = 2.0 * pi
      nkptd = 0 ; nvd = 0 ; nv2d = 0
      kw1d = 0 ; kw2d = 0 ; kw3d = 0
      nwd = nwdd
      CALL dotset(
     >            amat,bmat,
     <            aamat,bbmat)
c
        nqpt=1
      IF (l_J) THEN
        OPEN (113,file='qpts',form='formatted',status='old')
        READ (113,*) nqpt
      ENDIF
      OPEN (41,file='kpts',form='formatted',status='old')
      DO nw = 1,nwd

c--->    k-mesh: given in units of the reciprocal lattice basis vectors
c--->    scale is a factor to make input easier (default=1.0). k-pt
c--->    weights can be relative weights since they are renormalized.
c--->    input: for bulk - k1,k2,k3,wtkpt
c--->           for film - k1,k2,wtkpt
c--->    read k-points from file 41='kpts'
         IF (film) THEN
            READ (41,fmt=8050) nkpt,scale,xyu
         ELSE
            READ (41,fmt=8040) nkpt,scale
         END IF
 8030    FORMAT (4f10.5)
 8040    FORMAT (i5,f20.10)
 8050    FORMAT (i5,f20.10,3x,l1)

         nkptd = max(nkptd,nkpt)
 8060    FORMAT (i5,f20.10)
         IF (scale.EQ.0.0) scale = 1.0
         DO nq=1,nqpt
           IF(l_J) THEN
             READ (113,fmt=8070) qss(1),qss(2),qss(3)
 8070        FORMAT(2(f14.10,1x),f14.10)
           ENDIF

           DO nk = 1,nkpt
             IF(film) THEN
                READ (41,fmt=8080) (bk(i),i=1,2)
 8080        FORMAT (3f10.5)
             ELSE
                READ (41,fmt=8030) (bk(i),i=1,3)
             ENDIF
             IF (odd%d1) THEN
               bk(1) = 0.
               bk(2) = 0.
             ELSEIF (film .AND. .NOT.odd%d1) THEN
               bk(3) = 0.0
               IF (xyu) THEN
c            transform to cartesian coordinates
                  IF (latnam.EQ.'hex') THEN
                     bk(1) = bk(1)*tpi/amat(2,2)
                     bk(2) = bk(2)*pi/amat(1,1)
                  ELSE
                     bk(1) = bk(1)*pi/amat(1,1)
                     bk(2) = bk(2)*pi/amat(2,2)
                  END IF
c            transform to internal coordinates
                  s1 = 0.0
                  s2 = 0.0
                  DO j = 1,2
                     s1 = s1 + amat(j,1)*bk(j)/tpi
                     s2 = s2 + amat(j,2)*bk(j)/tpi
                  ENDDO
                  bk(1) = s1
                  bk(2) = s2
               END IF
             END IF
             DO i = 1,3
               bk(i) = bk(i)/scale
             ENDDO
             CALL apws_dim(
     >                     bk,bmat,bbmat,rkmax,l_ss,qss,odd,
     <                     nv,nv2,kw1,kw2,kw3)
             kw1d = max(kw1,kw1d)
             kw2d = max(kw2,kw2d)
             kw3d = max(kw3,kw3d)
             nvd = max(nvd,nv)
             nv2d = max(nv2d,nv2)

           ENDDO ! k=pts
         IF (nw == nwd) REWIND(41)
         ENDDO   ! q-pts
      ENDDO      ! windows
      
      IF (l_J) THEN
       CLOSE(113)
      ENDIF
      CLOSE (41)

      END SUBROUTINE inpeig_dim
      END MODULE m_inpeigdim
