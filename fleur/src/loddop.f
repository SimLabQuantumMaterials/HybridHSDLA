      MODULE m_loddop
      CONTAINS
      SUBROUTINE loddop(
     >                  jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >                  jspins,nq3,nq2,nvac,ntype,invs,invs2,film,
     >                  nlh,jri,ntypsd,ntypsy,nu,natd,neq,
     <                  iop,dop,it,fr,fpw,fz,fzxy,name)
c     ***********************************************************
c     reload formatted density or potential   c.l.fu
c     ***********************************************************

      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nu,ntypsd,natd
      INTEGER, INTENT (IN) :: jspd,n3d,n2d,nmzxyd,nmzd,nq3,nq2
      INTEGER, INTENT (IN) :: jmtd,nlhd,ntypd,jspins,nvac,ntype
      INTEGER, INTENT (OUT):: it
      LOGICAL, INTENT (IN) :: invs,invs2,film
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),neq(ntypd)
      COMPLEX, INTENT (OUT):: fpw(n3d,jspd),fzxy(nmzxyd,n2d-1,2,jspd)
      REAL,    INTENT (OUT):: fr(jmtd,0:nlhd,ntypd,jspd),fz(nmzd,2,jspd)
      CHARACTER*8,INTENT (OUT):: dop,iop,name(10)
C     ..
C     .. Local Scalars ..
      REAL delzn,dxn,rmtn,z1n,dummy
      INTEGER i,ivac,ivdummy,j,jrin,jsp,jspdum,k,lh,lhdummy,n,ndum,nlhn,
     +        nmzn,nmzxyn,nn,nq2n,nq3n,ntydum,n_diff,na
      CHARACTER*2 namaux
C     ..
C     .. Local Arrays ..
      REAL, ALLOCATABLE :: fpwr(:,:),fzxyr(:,:,:,:)

      CHARACTER*8 space(10)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx
C     ..
C     .. Data statements ..
      DATA space/10*'        '/
C     ..

      IF (invs) ALLOCATE ( fpwr(n3d,jspd) )
      IF (invs2) ALLOCATE ( fzxyr(nmzxyd,n2d-1,2,jspd) )

      DO 10 i = 1,10
         name(i) = space(i)
   10 CONTINUE
      READ (nu,END=200,ERR=200) name
c      WRITE (*,FMT=8000) name
c 8000 FORMAT (' loddop title:',10a8)
      READ (nu,END=200,ERR=200) iop,dop,it
      DO 130 jsp = 1,jspins
         READ (nu,END=200,ERR=200) jspdum
         READ (nu,END=200,ERR=200) nn
         IF (nn.NE.ntype) STOP 'loddop'
         na = 1
         DO 50 n = 1,nn
            READ (nu,END=200,ERR=200) namaux,ndum,jrin,rmtn,dxn
            READ (nu,END=200,ERR=200) ntydum,nlhn
c+gu
            IF ( nlhn.GT.nlh(ntypsy(na)) ) THEN
              WRITE (*,*) 'nlh (',nlhn,') set to (',nlh(ntypsy(na)),')'
              n_diff = nlhn - nlh(ntypsy(na))
              nlhn = nlh(ntypsy(na))
            ELSE
              n_diff = 0 
            ENDIF
c-gu
            DO 20 lh = 0,nlhn
               READ (nu,END=200,ERR=200) lhdummy
               READ (nu,END=200,ERR=200) (fr(i,lh,n,jsp),i=1,jrin)
   20       CONTINUE
            IF (nlhn.LT.nlh(ntypsy(na))) THEN
               DO lh = nlhn + 1,nlh(ntypsy(na))
                  DO i = 1,jri(n)
                     fr(i,lh,n,jsp) = 0.
                  ENDDO
               ENDDO
            ELSE
               DO lh = 1, n_diff
                 READ (nu,END=200,ERR=200) lhdummy
                 READ (nu,END=200,ERR=200) dummy
               ENDDO 
            ENDIF
            
            na = na + neq(n)
   50    CONTINUE
         READ (nu,END=200,ERR=200) nq3n
c+gu
         IF (nq3n.GT.nq3) THEN
           WRITE (*,*) 'nq3n (',nq3n,') reduced to nq3 (',nq3,')'
           nq3n = nq3
         ENDIF
c-gu
         IF (invs) THEN
            READ (nu,END=200,ERR=200) (fpwr(k,jsp),k=1,nq3n)
            DO 60 k = 1,nq3n
               fpw(k,jsp) = cmplx(fpwr(k,jsp),0.)
   60       CONTINUE
         ELSE
            READ (nu,END=200,ERR=200) (fpw(k,jsp),k=1,nq3n)
         END IF
         IF (nq3n.LT.nq3) THEN
            DO 70 k = nq3n + 1,nq3
               fpw(k,jsp) = (0.,0.)
   70       CONTINUE
         END IF
         IF (film) THEN
            DO 120 ivac = 1,nvac
               READ (nu,END=200,ERR=200) ivdummy
               READ (nu,END=200,ERR=200) nmzn,z1n,delzn
               READ (nu,END=200,ERR=200) (fz(i,ivac,jsp),i=1,nmzn)
               IF (nvac.EQ.1) THEN
                 DO i=1,nmzn
                   fz(i,2,jsp)=fz(i,1,jsp)
                 ENDDO
               ENDIF
               READ (nu,END=200,ERR=200) nq2n,nmzxyn
c+gu
               IF (nq2n.GT.nq2) THEN
                 WRITE (*,*) 'nq2n (',nq2n,') reduced to nq2 (',nq2,')'
                 n_diff = nq2n - nq2
                 nq2n = nq2
               ELSE
                 n_diff = 0
               ENDIF
c-gu
               DO 90 k = 2,nq2n
                  IF (invs2) THEN
                     READ (nu,END=200,ERR=200) 
     +                              (fzxyr(j,k-1,ivac,jsp),j=1,nmzxyn)
                     DO 80 j = 1,nmzxyn
                        fzxy(j,k-1,ivac,jsp) = cmplx(fzxyr(j,k-1,ivac,
     +                                         jsp),0.)
   80                CONTINUE
                  ELSE
                     READ (nu,END=200,ERR=200)  
     +                               (fzxy(j,k-1,ivac,jsp),j=1,nmzxyn)
                  END IF
                  IF (nvac.EQ.1) THEN
                    IF (invs) THEN
                      DO j = 1,nmzxyn
                        fzxy(j,k-1,2,jsp) = conjg(fzxy(j,k-1,1,jsp))
                      ENDDO
                    ELSE
                      DO j = 1,nmzxyn
                        fzxy(j,k-1,2,jsp) = fzxy(j,k-1,1,jsp)
                      ENDDO
                    ENDIF
                  ENDIF
   90          CONTINUE
c+gu
               DO k = 1,n_diff
                  READ (nu,END=200,ERR=200) dummy
               ENDDO
c-gu
               IF (nq2n.LT.nq2) THEN
                  DO 110 k = nq2n + 1,nq2
                     DO 100 j = 1,nmzxyn
                        fzxy(j,k-1,ivac,jsp) = (0.,0.)
  100                CONTINUE
  110             CONTINUE
               END IF
  120       CONTINUE
         END IF
  130 CONTINUE
c
      IF (invs) DEALLOCATE (fpwr)
      IF (invs2) DEALLOCATE ( fzxyr )
      RETURN

 200  WRITE (6,*) 'error reading dop nr.',nu
      IF (nu.NE.98) STOP 'loddop: error reading d/p-file!'

      END SUBROUTINE loddop
      END MODULE m_loddop
