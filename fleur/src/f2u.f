      MODULE m_f2u
c     ****************************************************
c     write unformatted density or potential onto unit '99'
c     e. wimmer   march 1985
c     ****************************************************
      CONTAINS
      SUBROUTINE f2u(
     >               nq3,nq2,jspins,l_noco,jri,nlh,ntype,nmz,nmzxy,
     >               ntypsd,ntypsy,z1,delz,invs,invs2,film,
     >               natd,neq,nvac,zatom,namat,dx,rmt,nwd)
c
      USE m_wrtdop

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nq3,nq2,jspins,natd,nwd
      INTEGER, INTENT (IN) :: ntype,nmz,nmzxy,ntypsd,nvac
      LOGICAL, INTENT (IN) :: invs,invs2,film,l_noco
      REAL,    INTENT (IN) :: z1,delz
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntype),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),neq(ntype)
      REAL,    INTENT (IN) :: rmt(ntype),dx(ntype),zatom(ntype)
      CHARACTER(len=2),INTENT (IN):: namat(0:103)
c
C     .. Local Scalars ..
      INTEGER it,nu,i,jsp,jspdum,nn,n,ndum,jrin,j,ntypsyn,na
      INTEGER nlhn,lh,lhdum,nq3n,k,ivac,ivdum,nmzn,nq2n,nmzxyn
      INTEGER jmtd,nlhd
      REAL    z1n,delzn,rmtn,dxn
      LOGICAL n_exist
      CHARACTER(len=2) namaux
      CHARACTER(len=8) dop,iop
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: n_mmp(:,:,:)
      COMPLEX, ALLOCATABLE :: fpw(:,:),fzxy(:,:,:,:)
      COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)
      REAL,    ALLOCATABLE :: fz(:,:,:),fr(:,:,:,:)
      REAL,    ALLOCATABLE :: fpwr(:,:),fzxyr(:,:,:,:)
      CHARACTER(len=8) name(10),space(10)
C     ..
C     .. Data statements ..
      DATA space/10*'        '/
C     ..
c
      jmtd = maxval(jri(:))
      nlhd = maxval(nlh(:))

      IF (l_noco) THEN
         OPEN (98,file='rhomat_inp',form='unformatted')
         OPEN (99,file='rhomat_form',form='formatted')
      ELSE
         OPEN (98,file='cdn_unf',form='unformatted')
         OPEN (99,file='cdn_form',form='formatted')
      ENDIF
      nu=99
c
      ALLOCATE( fpw(nq3,jspins),fzxy(nmzxy,nq2-1,2,jspins) )
      ALLOCATE( cdom(nq3),cdomvz(nmz,2),cdomvxy(nmzxy,nq2-1,2) )
      ALLOCATE( fr(jmtd,0:nlhd,ntype,jspins),fz(nmz,2,jspins) )
      ALLOCATE( fpwr(nq3,jspins),fzxyr(nmzxy,nq2-1,2,jspins) )

      DO 10 i = 1,10
         name(i) = space(i)
   10 CONTINUE
      READ (nu,FMT=8030) name
      WRITE (6,FMT=8000) name
 8000 FORMAT (' loddop title:',10a8)
      READ (nu,FMT=8050) iop,dop,it
      DO 60 jsp = 1,jspins
         READ (nu,FMT=8060) jspdum
         READ (nu,FMT=8070) nn
         IF (nn.NE.ntype) STOP 'loddop'
         na = 1
         DO 30 n = 1,nn
            READ (nu,FMT=8010) namaux,ndum,jrin,rmtn,dxn
 8010       FORMAT (a2,4x,i3,6x,i3,6x,f10.6,5x,f10.6)
            READ (nu,FMT=8020) ntypsyn,nlhn
            IF (ntypsyn.NE.ntypsy(na)) THEN
              WRITE (6,*) 'ntypsyn = ',ntypsyn, ' =/='
              WRITE (6,*) 'ntypsy(na) = ',ntypsy(na)
              STOP 'f2u: ntypsyn.NE.ntypsy(na)'
            ENDIF
 8020       FORMAT (7x,i2,2x,4x,i3)
            DO 20 lh = 0,nlhn
               READ (nu,FMT=8060) lhdum
               READ (nu,FMT=8100) (fr(i,lh,n,jsp),i=1,jrin)
   20       CONTINUE
            na = na + neq(n)
   30    CONTINUE
         READ (nu,FMT=8110) nq3n
         IF (invs) THEN
            READ (nu,FMT=8100) (fpwr(k,jsp),k=1,nq3n)
            DO k=1,nq3n
              fpw(k,jsp) = cmplx(fpwr(k,jsp),0.0)
            ENDDO
         ELSE
            READ (nu,FMT=8100) (fpw(k,jsp),k=1,nq3n)
         END IF
         IF (film) THEN
            DO 50 ivac = 1,nvac
               READ (nu,FMT=8060) ivdum
               READ (nu,FMT=8120) nmzn,z1n,delzn
               READ (nu,FMT=8040) (fz(i,ivac,jsp),i=1,nmzn)
               READ (nu,FMT=8130) nq2n,nmzxyn
               DO 40 k = 2,nq2n
                  IF (invs2) THEN
                     READ (nu,FMT=8100) (fzxyr(j,k-1,ivac,jsp),j=1,
     +                 nmzxyn)
                     DO j=1,nmzxyn
                       fzxy(j,k-1,ivac,jsp) = 
     +                            cmplx(fzxyr(j,k-1,ivac,jsp),0.0)
                     ENDDO
                  ELSE
                     READ (nu,FMT=8100) (fzxy(j,k-1,ivac,jsp),j=1,
     +                 nmzxyn)
                  END IF
   40          CONTINUE
   50       CONTINUE
         END IF
   60 CONTINUE

      IF (l_noco) THEN
         READ (99,8100) (cdom(k),k=1,nq3)
         IF (film) THEN
            READ (99,8100) ((cdomvz(j,ivac),j=1,nmz),ivac=1,nvac)
            READ (99,8100) (((cdomvxy(j,k-1,ivac),j=1,nmzxy),k=2,nq2),
     +                                              ivac=1,nvac)
         ENDIF
      ENDIF

      CLOSE (99)

 8030 FORMAT (10a8)
 8040 FORMAT (4e20.13)
 8050 FORMAT (a8,1x,a8,11x,i3)
 8060 FORMAT (a9)
 8070 FORMAT (6x,i3)
 8080 FORMAT (15x,i3,6x,f10.6,5x,f10.6)
 8090 FORMAT (15x,i2)
 8100 FORMAT (4e20.13)
 8110 FORMAT (4x,i6)
 8120 FORMAT (4x,i3,5x,f20.13,7x,f8.4)
 8130 FORMAT (4x,i5,8x,i5)
c
      CALL wrtdop(
     >            jspins,nq3,nq2,nmzxy,nmz,jmtd,nlhd,ntype,natd,
     >            jspins,nq3,nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,98,
     >            iop,dop,it,fr,fpw,fz,fzxy,name)

      IF (l_noco) THEN
         WRITE (98) (cdom(k),k=1,nq3)
         IF (film) THEN
            WRITE (98) ((cdomvz(j,ivac),j=1,nmz),ivac=1,nvac)
            WRITE (98) (((cdomvxy(j,k-1,ivac),j=1,nmzxy),k=2,nq2),
     +                                              ivac=1,nvac)
         ENDIF
      ENDIF
      CLOSE (98)
      DEALLOCATE( fpw,fzxy,cdom,cdomvz,cdomvxy,fr,fz,fpwr,fzxyr )
         
      END SUBROUTINE f2u
      END MODULE m_f2u
