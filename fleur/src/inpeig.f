      MODULE m_inpeig
      CONTAINS
      SUBROUTINE inpeig(
     >                  lmaxd,ntypd,jspd,nkptd,nlod,
     >                  ntype,latnam,nwd,jspins,film,nvac,
     >                  lmax,amat,bmat,sc,neq,nlo,odd,
     <                  bk,wtkpt,nkpt,ello0,llochg,skiplo,
     <                  el0,evac0,lchange,lchg_v,enmix)
c*********************************************************************
c     inputs the necessary quantities for the eigenvalue part (energy
c     parameters, k-points, wavefunction cutoffs, etc.).
c                  m. weinert   jan. 1987
c     modification dec. 1990:
c     dummyline before reading l-dependent energies to make reading of
c     input easier (e.g. insert name of atom etc.)
c     modification dec. 93:
c     for step-forward diagonalization a la wu in case of more
c     than 1 window we read now
c     number of occupied states for EACH window
c*********************************************************************
      USE m_constants, ONLY : pimach
      USE m_enpara,    ONLY : r_enpara
      USE m_od_types,  ONLY : od_dim
      USE m_fleurenv

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,ntypd,jspd,nkptd,nlod
      INTEGER, INTENT (IN) :: ntype,nwd,jspins,nvac
      REAL,    INTENT (IN) :: sc
      LOGICAL, INTENT (IN) :: film
      CHARACTER*3,INTENT (IN) :: latnam
c-odim
      TYPE (od_dim), INTENT (IN) :: odd
c+odim
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: lmax(ntypd),neq(ntypd),nlo(ntypd)
      INTEGER, INTENT (OUT) :: nkpt(nwd),skiplo(ntypd,jspd)
      REAL,    INTENT (IN)  :: amat(3,3),bmat(3,3)
      REAL,    INTENT (OUT) :: el0(0:lmaxd,ntypd,jspd,nwd)
      REAL,    INTENT (OUT) :: evac0(2,jspd,nwd),ello0(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: bk(3,nkptd,nwd),wtkpt(nkptd,nwd)
      LOGICAL, INTENT (OUT) :: lchange(0:lmaxd,ntypd,jspd,nwd)
      LOGICAL, INTENT (OUT) :: lchg_v(2,jspd,nwd)
      LOGICAL, INTENT (OUT) :: llochg(nlod,ntypd,jspd)
      REAL,    INTENT (OUT) :: enmix(jspd,nwd)
C     ..
C     .. Local Scalars ..
      REAL s1,s2,scale,wt,pi,tpi
      INTEGER i,j,nk,nw,jsp
      LOGICAL xyu,l_enpara
C     ..
C     .. External Subroutines ..
      EXTERNAL gkptwgt
c
      pi = pimach()
      tpi = 2 * pi
c
c---> input energy parameters for each atom.
c---> the energy parameters for l.ge.3 have the same value
c---> read from file 40='enpara'  shz Jan.96
c
! YM : it is not nice that fleur crashes every time without enpara
!   let us be nice to the users    
      l_enpara = .FALSE.
      INQUIRE (file='enpara',exist=l_enpara)
      IF (l_enpara) THEN
      OPEN (40,file='enpara',form='formatted',status='old')
      ELSE
      STOP 'enpara is not provided'
      END IF
      DO nw = 1,nwd
         DO jsp = 1,jspins
            CALL r_enpara(     
     >                    lmaxd,nlod,ntype,film,jsp,
     >                    nw,nlo,lmax,neq,
     <                    skiplo(1,jsp),ello0(1,1,jsp),
     <                    el0(0,1,jsp,nw),evac0(1,jsp,nw),
     <                    lchange(0,1,jsp,nw),llochg(1,1,jsp),
     <                    lchg_v(1,jsp,nw),enmix(jsp,nw) )
         ENDDO !jspd
      ENDDO !nw
      CLOSE (40)
!
!---> read k-points from file 41='kpts'
!
      OPEN (41,file='kpts',form='formatted',status='old')
!
c---> k-mesh: given in units of the reciprocal lattice basis vectors
c---> scale is a factor to make input easier (default=1.0). k-pt
c---> weights can be relative weights since they are renormalized.
c---> input: for bulk - k1,k2,k3,wtkpt
c--->        for film - k1,k2,wtkpt
c--->           weights are calculated for films, if wtkpt=0
c     for film calculation k1,k2 may also be read in xy - units : xyu=T
c     1 = boundery of BZ on kx/ky axis
c                                                  shz Feb.96
      DO nw=1,nwd
         READ (41,FMT=8110,ERR=911,END=911) nkpt(nw),scale,xyu
         GOTO 912
  911    CONTINUE
         xyu = .false.
  912    CONTINUE
         
         IF (nkpt(nw).GT.nkptd)  THEN
           CALL fleur_err('nkptd too small')
         ENDIF
 8100    FORMAT (i5,f20.10)
 8110    FORMAT (i5,f20.10,3x,l1)
         IF (scale.EQ.0.0) scale = 1.0
         DO nk = 1,nkpt(nw)
            READ (41,FMT=8040) (bk(i,nk,nw),i=1,3),wtkpt(nk,nw)
 8040       FORMAT (4f10.5)
            IF (film .AND. .NOT.odd%d1) THEN
               wtkpt(nk,nw) = bk(3,nk,nw)
               bk(3,nk,nw) = 0.0
               IF (xyu) THEN
c           transform to cartesian coordinates
                  IF (latnam.EQ.'hex') THEN
                     bk(1,nk,nw) = bk(1,nk,nw)*tpi/amat(2,2)
                     bk(2,nk,nw) = bk(2,nk,nw)*pi/amat(1,1)
                  ELSE
                     bk(1,nk,nw) = bk(1,nk,nw)*pi/amat(1,1)
                     bk(2,nk,nw) = bk(2,nk,nw)*pi/amat(2,2)
                  END IF
c           transform to internal coordinates
                  s1 = 0.0
                  s2 = 0.0
                  DO  j = 1,2
                     s1 = s1 + amat(j,1)*bk(j,nk,nw)/tpi
                     s2 = s2 + amat(j,2)*bk(j,nk,nw)/tpi
                  ENDDO
                  bk(1,nk,nw) = s1
                  bk(2,nk,nw) = s2
               END IF
            ELSEIF (.NOT.film) THEN
              IF (xyu) THEN
                bk(:,nk,nw) = bk(:,nk,nw)*2.0/sc
                bk(:,nk,nw) = matmul( amat, bk(:,nk,nw) )
              ENDIF
            END IF
            DO  i = 1,3
               bk(i,nk,nw) = bk(i,nk,nw)/scale
            ENDDO
c-odim
            IF (odd%d1) THEN
c--> trapezoidal
               bk(1,nk,nw) = 0.0
               bk(2,nk,nw) = 0.0
               IF (bk(3,nk,nw).EQ.0. .OR. bk(3,nk,nw).EQ.0.5) THEN
                  wtkpt(nk,nw) = 1.
               ELSE
                  wtkpt(nk,nw) = 2.
               END IF
c--> 4th Simpson integration
c              IF (bk(3,nk,nw).EQ.0. .OR. bk(3,nk,nw).EQ.0.5) THEN
c                 wtkpt(nk,nw) = 3./8.
c                 WRITE (*,*) '3/8:nkpt=',nk
c              ELSEIF (MOD(nk,4).EQ.0) THEN
c                 wtkpt(nk,nw) = 6./8.
c                 WRITE (*,*) '6/8:nkpt=',nk
c              ELSE
c                 wtkpt(nk,nw) = 9./8.
c                 WRITE (*,*) '9/8:nkpt=',nk
c              END IF
c--> 2nd Simpson integration
c              IF (bk(3,nk,nw).EQ.0. .OR. bk(3,nk,nw).EQ.0.5) THEN
c                 wtkpt(nk,nw) = 1./3.
c                 WRITE (*,*) '1/3:nkpt=',nk
c              ELSEIF (MOD(nk,2).EQ.0) THEN
c                 wtkpt(nk,nw) = 4./3.
c                 WRITE (*,*) '4/3:nkpt=',nk
c              ELSE
c                 wtkpt(nk,nw) = 2./3.
c                 WRITE (*,*) '2/3:nkpt=',nk
c              END IF
            END IF
c-odim
         ENDDO
         wt = 0.0
         DO nk = 1,nkpt(nw)
            wt = wt + wtkpt(nk,nw)
         ENDDO
         IF (wt.EQ.0.0) THEN
            IF (film) THEN
c
c---> generate k-point weights for 2d BZ: squ, rec, cen, hex
c     determine new wt
c
               CALL gkptwgt(
     >                      nkptd,amat,bmat,nkpt(nw),latnam,
     X                      bk(1,1,nw),
     <                      wtkpt(1,nw),wt)
            ELSE
               STOP 'wtkpts'
            END IF
         END IF
c-odim
         DO nk = 1,nkpt(nw)
            IF (odd%d1) THEN
               bk(1,nk,nw) = 0.0
               bk(2,nk,nw) = 0.0
            END IF
         END DO
c+odim
         DO  nk = 1,nkpt(nw)
            wtkpt(nk,nw) = wtkpt(nk,nw)/wt
         ENDDO
         WRITE (6,FMT=8120) nkpt(nw)
         WRITE (16,FMT=8120) nkpt(nw)
         DO  nk = 1,nkpt(nw)
            WRITE (6,FMT=8040) (bk(i,nk,nw),i=1,3),wtkpt(nk,nw)
            WRITE (16,FMT=8040) (bk(i,nk,nw),i=1,3),wtkpt(nk,nw)
         ENDDO
 8120    FORMAT (1x,/,' number of k-points for this window =',i5,/,t12,
     +          'coordinates',t34,'weights')
      ENDDO
      CLOSE (41)

      END SUBROUTINE inpeig
      END MODULE m_inpeig
