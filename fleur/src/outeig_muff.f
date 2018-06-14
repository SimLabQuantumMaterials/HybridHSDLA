      SUBROUTINE outeig_muff(
     >                       neigd,ellow,lpr,film,eonly,
     >                       bkpt,nv,eig,ne,nblw,nv2)
c*********************************************************************
c     outputs the eigenvalues and eigenvectors to unit66.
c     if lpr.gt.0, then also list eigenvectors on output file;
c     otherwise only eigenvalues. form66=.true. gives a formatted
c     eigenvector file.
c             m. weinert
c*********************************************************************

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT (IN) :: neigd
      INTEGER,INTENT (IN) :: lpr,ne,nblw,nv2,nv
      REAL,   INTENT (IN) :: ellow
      LOGICAL eonly,film
C     ..
C     .. Array Arguments ..
      REAL,   INTENT (IN) :: eig(neigd),bkpt(3)
C     ..
C     .. Local Scalars ..
      INTEGER i,m,n4,nbg,nend
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC min
C     ..
      IF (film) THEN
         WRITE (6,FMT=8000) (bkpt(i),i=1,3),nv,nv2,ne
         WRITE (16,FMT=8000) (bkpt(i),i=1,3),nv,nv2,ne
      ELSE
         WRITE (6,FMT=8010) (bkpt(i),i=1,3),nv,ne
         WRITE (16,FMT=8010) (bkpt(i),i=1,3),nv,ne
      END IF
 8000 FORMAT (1x,/,/,' k=(',3f12.6,'):',i10,' 3-d basis functions and',
     +       i6,' 2-d basis functions',/,' the',i4,' eigenvalues are:')
 8010 FORMAT (1x,/,/,' k=(',3f12.6,'):',i10,' 3-d basis functions',/,
     +       ' the',i4,' eigenvalues are:')
      WRITE (6,FMT=8020) (eig(i),i=1,ne)
      WRITE (16,FMT=8020) (eig(i),i=1,ne)
 8020 FORMAT (5x,5f12.6)
c
      IF (nblw.GT.-9999) THEN
        WRITE (6,FMT=8030) nblw,ellow
        WRITE (16,FMT=8030) nblw,ellow
      ENDIF
 8030 FORMAT (5x,60 ('*'),/,10x,i3,
     +       ' eigenvalue(s) below the energy FOR MUFFIN_TIN H',f12.6,/,
     +       5x,60 ('*'),/)
      IF ((lpr.GT.0) .AND. (.NOT.eonly)) THEN
         n4 = ne/4
         IF (4*n4.LT.ne) n4 = n4 + 1
         WRITE (6,FMT=8040)
         nend = 0
         DO 10 m = 1,n4
            nbg = nend + 1
            nend = min(nend+4,ne)
            WRITE (6,FMT=8050) (eig(i),i=nbg,nend)
   10    CONTINUE
      END IF
 8040 FORMAT (1x,/,' eigenvectors:')
 8050 FORMAT (1x,/,' eigenvalues FOR MUFFIN_TIN HAMILTONIAN',
     +       4 (8x,f12.6,5x))
 8060 FORMAT (3i5,5x,4 (3x,2f11.6))
c
      RETURN
      END

