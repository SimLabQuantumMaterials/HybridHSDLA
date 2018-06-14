      MODULE m_u2f
c     ****************************************************
c     write formatted density or potential onto unit '99'
c     e. wimmer   march 1985
c     ****************************************************
      CONTAINS
      SUBROUTINE u2f(
     >               nq3,nq2,jspins,l_noco,jri,nlh,ntype,nmz,nmzxy,
     >               ntypsd,ntypsy,z1,delz,invs,invs2,film,
     >               natd,neq,nvac,zatom,namat,dx,rmt,nwd)
c
      USE m_loddop
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nq3,nq2,jspins,ntype,natd
      INTEGER, INTENT (IN) :: nmz,nmzxy,ntypsd,nvac,nwd
      LOGICAL, INTENT (IN) :: invs,invs2,film,l_noco
      REAL,    INTENT (IN) :: z1,delz
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntype),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),neq(ntype)
      REAL,    INTENT (IN) :: rmt(ntype),dx(ntype),zatom(ntype)
      CHARACTER(len=2), INTENT (IN):: namat(0:103)
C     ..
C     .. Local Scalars ..
      INTEGER i,ivac,j,jsp,k,lh,n,izn,na,urec,jmtd,nlhd
      LOGICAL n_exist
C     ..
C     .. Local Arrays ..
      INTEGER it
      CHARACTER(len=8) dop,iop,name(10)
      COMPLEX, ALLOCATABLE :: fpw(:,:),fzxy(:,:,:,:)
      COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)
      REAL,    ALLOCATABLE :: fz(:,:,:),fr(:,:,:,:)
C     ..
C     ..
      jmtd = maxval(jri(:))
      nlhd = maxval(nlh(:))

      IF (l_noco) THEN
         OPEN (98,file='rhomat_inp',form='unformatted')
         OPEN (99,file='rhomat_form',form='formatted')
      ELSE
         OPEN (98,file='f_unf',form='unformatted')
         OPEN (99,file='f_form',form='formatted')
      ENDIF
c
      ALLOCATE( fpw(nq3,jspins),fzxy(nmzxy,nq2-1,2,jspins) )
      ALLOCATE( cdom(nq3),cdomvz(nmz,2),cdomvxy(nmzxy,nq2-1,2) )
      ALLOCATE( fr(jmtd,0:nlhd,ntype,jspins),fz(nmz,2,jspins) )

      CALL loddop(
     >            jspins,nq3,nq2,nmzxy,nmz,jmtd,nlhd,ntype,
     >            jspins,nq3,nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,98,natd,neq,
     <            iop,dop,it,fr,fpw,fz,fzxy,name)
      write(*,*) 'after loddop:', l_noco
      IF (l_noco) THEN
         cdom(:) = (0.0,0.0)
         READ (98,END=101,ERR=101) (cdom(k),k=1,nq3)
         IF (film) THEN
            READ (98) ((cdomvz(j,ivac),j=1,nmz),ivac=1,nvac)
            READ (98) (((cdomvxy(j,k-1,ivac),j=1,nmzxy),k=2,nq2),
     +                                              ivac=1,nvac)
         ENDIF
      ENDIF
         
  101 CLOSE (98)

      WRITE (99,FMT=8000) name
 8000 FORMAT (10a8)
      WRITE (6,FMT=8010) name
 8010 FORMAT (' loddop title:',10a8)
      WRITE (99,FMT=8020) iop,dop,it
 8020 FORMAT (a8,1x,a8,' iteration=',i3)
      DO 60 jsp = 1,jspins
         WRITE (99,FMT=8030) jsp
 8030    FORMAT ('spin=',i1)
         WRITE (99,FMT=8040) ntype
 8040    FORMAT ('ntype=',i3)
         na = 1
         DO 30 n = 1,ntype
            izn = zatom(n) + 0.01
            WRITE (99,FMT=8050) namat(izn),n,jri(n),rmt(n),dx(n)
 8050       FORMAT (a2,2x,'n=',i3,2x,'jri=',i3,2x,'rmt=',f10.6,2x,'dx=',
     +             f10.6)
            WRITE (99,FMT=8060) ntypsy(na),nlh(ntypsy(na))
 8060       FORMAT ('ntypsy=',i2,2x,'nlh=',i3)
            DO 20 lh = 0,nlh(ntypsy(na))
               WRITE (99,FMT=8070) lh
 8070          FORMAT ('lh=',i3)
               WRITE (99,FMT=8080) (fr(i,lh,n,jsp),i=1,jri(n))
 8080          FORMAT (4e20.13)
   20       CONTINUE
            na = na + neq(n)
   30    CONTINUE
         WRITE (99,FMT=8090) nq3
 8090    FORMAT ('nq3=',i6)
         IF (invs) THEN
            WRITE (99,FMT=8080) (real(fpw(k,jsp)),k=1,nq3)
         ELSE
            WRITE (99,FMT=8080) (fpw(k,jsp),k=1,nq3)
         END IF
         IF (film) THEN
            DO 50 ivac = 1,nvac
               WRITE (99,FMT=8100) ivac
 8100          FORMAT ('ivac=',i1)
               WRITE (99,FMT=8110) nmz,z1,delz
 8110          FORMAT ('nmz=',i3,2x,'z1=',f20.13,2x,'delz=',f8.4)
               WRITE (99,FMT=8080) (fz(i,ivac,jsp),i=1,nmz)
               WRITE (99,FMT=8120) nq2,nmzxy
 8120          FORMAT ('nq2=',i5,2x,'nmzxy=',i5)
               DO 40 k = 2,nq2
                  IF (invs2) THEN
                     WRITE (99,FMT=8080) (real(fzxy(j,k-1,ivac,jsp)),
     +                                                    j=1,nmzxy)
                  ELSE
                     WRITE (99,FMT=8080) (fzxy(j,k-1,ivac,jsp),j=1,
     +                 nmzxy)
                  END IF
   40          CONTINUE
   50       CONTINUE
         END IF
   60 CONTINUE

      IF (l_noco) THEN
         WRITE (99,8080) (cdom(k),k=1,nq3)
         IF (film) THEN
            WRITE (99,8080) ((cdomvz(j,ivac),j=1,nmz),ivac=1,nvac)
            WRITE (99,8080) (((cdomvxy(j,k-1,ivac),j=1,nmzxy),k=2,nq2),
     +                                              ivac=1,nvac)
         ENDIF
      ENDIF
      CLOSE (99)
      DEALLOCATE( fpw,fzxy,cdom,cdomvz,cdomvxy,fr,fz)

      END SUBROUTINE u2f
      END MODULE m_u2f
