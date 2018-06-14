      SUBROUTINE unor2or(
     >                   index2,index3,cdnmark,
     >                   n3d,n2d,jspd,jmtd,nlhd,ntypd,nmzd,nmzxyd,natd,
     >                   nq3,nq2,jspins,jri,nlh,ntype,nmz,nmzxy,neq,
     >                   ntypsd,ntypsy,z1,delz,invs,invs2,film,
     >                   nvac,zatom,namat,dx,rmt)
c
c     ****************************************************
c     brings order to the stars
c     ****************************************************
c
      USE m_loddop
      USE m_wrtdop
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      CHARACTER*8 cdnmark
      INTEGER, INTENT (IN) :: n3d,n2d,jspd,jmtd,nlhd
      INTEGER, INTENT (IN) :: ntypd,nmzd,nmzxyd,natd
      INTEGER, INTENT (IN) :: nq3,nq2,jspins,ntype,nmz,nmzxy,ntypsd,nvac
      LOGICAL, INTENT (IN) :: invs,invs2,film
      REAL,    INTENT (IN) :: z1,delz
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: index2(n2d),index3(n3d),neq(ntypd)
      INTEGER, INTENT (IN) :: jri(ntypd),nlh(ntypsd),ntypsy(natd)
      REAL,    INTENT (IN) :: rmt(ntypd),dx(ntypd),zatom(ntypd)
      CHARACTER*2,INTENT (IN):: namat(0:103)
C     ..
C     .. Local Scalars ..
      INTEGER ivac,j,jsp,k
C     ..
C     .. Local Arrays ..
      INTEGER it
      CHARACTER*8 dop,iop,name(10)
      COMPLEX, ALLOCATABLE :: fpw(:,:),fzxy(:,:,:,:)
      COMPLEX, ALLOCATABLE :: fpw_or(:,:),fzxy_or(:,:,:,:)
      REAL,    ALLOCATABLE :: fz(:,:,:),fr(:,:,:,:)
C     ..
c...... cdn1 file
c
      OPEN (98,file='cdn1',form='unformatted')
      CLOSE (99)
      IF (cdnmark.eq.'ordered*') THEN
        OPEN (99,file='cdn_ordered',form='unformatted')
      ELSE
        OPEN (99,file='cdn_unordered',form='unformatted')
      ENDIF

      ALLOCATE( fpw(n3d,jspd),fzxy(nmzxyd,n2d-1,2,jspd) )
      ALLOCATE( fpw_or(n3d,jspd),fzxy_or(nmzxyd,n2d-1,2,jspd) )
      ALLOCATE( fr(jmtd,0:nlhd,ntypd,jspd),fz(nmzd,2,jspd) )

      CALL loddop(
     >            jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,98,natd,neq,
     <            iop,dop,it,fr,fpw,fz,fzxy,name)
c
c order the coefficients
c
      name(10)=cdnmark
      DO jsp = 1,jspins
        DO k = 1,nq3
          fpw_or(k,jsp) = fpw(index3(k),jsp)
        ENDDO
        IF (film) THEN
          DO ivac = 1,nvac
            DO k = 2,nq2
              DO j = 1,nmzxy
                fzxy_or(j,k-1,ivac,jsp) = fzxy(j,index2(k)-1,ivac,jsp)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
c
c write to cdn_(un)ordered
c
      CALL wrtdop(
     >            jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            jspins,nq3,nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,99,
     >            iop,dop,it,fr,fpw_or,fz,fzxy_or,name)
c
      DEALLOCATE( fpw,fzxy,fpw_or,fzxy_or,fr,fz ) 
      RETURN
      END

