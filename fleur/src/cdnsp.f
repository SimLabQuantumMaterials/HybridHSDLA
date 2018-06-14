      MODULE m_cdnsp
c     *******************************************************
c     sets up the starting density for the spin-polarized
c     calculation from a paramagnetic density
c     changed to suit both ferromagnetic and antiferro-
c     magnetic case. changes only in mt-part - r.pentcheva Jan'96
c     *******************************************************
      CONTAINS
      SUBROUTINE cdnsp(
     >                 ntype,jspins,nmz,nmzxy,ntypsd,ntypsy,
     >                 nvac,nq2,nq3,jmtd,natd,neq,
     >                 jri,nlh,zatom,rmt,dx,film,invs,invs2,
     >                 z1,delz,rmsh,msh,bmu,namat,n_u)

      USE m_intgr, ONLY : intgr3
      USE m_constants, ONLY : pimach
      USE m_loddop
      USE m_wrtdop

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntype,jspins,nmzxy,nmz,nvac
      INTEGER, INTENT (IN) :: jmtd,msh,ntypsd,natd,nq2,nq3,n_u
      REAL,    INTENT (IN) :: z1,delz
      LOGICAL, INTENT (IN) :: film,invs,invs2
C     ..
C     .. Array Arguments ..
      INTEGER,     INTENT (IN) :: jri(ntype),neq(ntype)
      INTEGER,     INTENT (IN) :: nlh(ntypsd),ntypsy(natd)
      REAL,        INTENT (IN) :: zatom(ntype),rmt(ntype),dx(ntype)
      REAL,        INTENT (IN) :: rmsh(jmtd,ntype)
      REAL,        INTENT (IN) :: bmu(ntype)
      CHARACTER*2, INTENT (IN) :: namat(0:103)
C     ..
C     .. Local Scalars ..
      REAL dummy,p,pp,qtot1,qtot2,spmtot,qval,sfp
      INTEGER i,iter,ivac,j,k,lh,n,na,nt,jsp_new
      INTEGER ios,nlhd
      CHARACTER(len=8) dop,iop
      LOGICAL n_exist
C     ..
C     .. Local Arrays ..
      REAL rhoc(msh,ntype)
      COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
      REAL   , ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:)
      CHARACTER(len=140), ALLOCATABLE :: clines(:)
      CHARACTER(len=140)              :: lineread
      CHARACTER(len=8) name(10)
C     ..
      sfp = 2 * sqrt(pimach())
      nlhd = maxval(nlh(:))

      IF (jspins.NE.2) STOP 'cdnsp: set jspins=2 and remove fl7para!'
      ALLOCATE ( rho(jmtd,0:nlhd,ntype,jspins),qpw(nq3,jspins) )
      ALLOCATE ( rhtxy(nmzxy,nq2-1,2,jspins),rht(nmz,2,jspins) )

c
      OPEN (17,file='cdnc',form='unformatted',status='old')
      DO n = 1,ntype
         READ (17) (rhoc(j,n),j=1,jri(n))
         READ (17) dummy
      ENDDO
      CLOSE (17)

      nt = 71                   !gs, see sub mix
      OPEN (nt,file='cdn1',form='unformatted',status='old')
c
c     ---> set jspins=1 to read the paramagnetic density
c 
      CALL loddop(
     >            jspins,nq3,nq2,nmzxy,nmz,jmtd,nlhd,ntype,
     >            1,nq3,nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nt,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
c
      qval = 0.
      na = 1
c
c     ---> set jspins=2
      jsp_new = 2
c
      DO n = 1,ntype
         DO j = 1,jri(n)
            rho(j,0,n,1) = rho(j,0,n,1) - rhoc(j,n)/sfp
         ENDDO
!         WRITE (16,FMT='(8f10.4)') (rho(i,0,n,1),i=1,16)
         CALL intgr3(rho(1,0,n,1),rmsh(1,n),dx(n),jri(n),qval)
         p = (bmu(n)+sfp*qval)/ (2.*sfp*qval)
         pp = 1. - p
         DO j = 1,jri(n)
            rho(j,0,n,jsp_new) = pp*rho(j,0,n,1) + rhoc(j,n)/ (2.*sfp)
            rho(j,0,n,1)       =  p*rho(j,0,n,1) + rhoc(j,n)/ (2.*sfp)
         ENDDO
         DO lh = 1,nlh(ntypsy(na))
            DO j = 1,jri(n)
               rho(j,lh,n,jsp_new) = pp*rho(j,lh,n,1)
               rho(j,lh,n,1)       =  p*rho(j,lh,n,1)
            ENDDO
         ENDDO
         na = na + neq(n)
      ENDDO
      DO k = 1,nq3
         qpw(k,jsp_new) = 0.5 * qpw(k,1)
         qpw(k,1)       = qpw(k,jsp_new)
      ENDDO
      IF (film) THEN
         DO ivac = 1,nvac
            DO j = 1, nmz
               rht(j,ivac,jsp_new) = 0.5 * rht(j,ivac,1)
               rht(j,ivac,1)       = rht(j,ivac,jsp_new)
            ENDDO
            DO k = 2, nq2
               DO j = 1,nmzxy
                  rhtxy(j,k-1,ivac,jsp_new) = 0.5 * rhtxy(j,k-1,ivac,1)
                  rhtxy(j,k-1,ivac,1)       = rhtxy(j,k-1,ivac,jsp_new)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
c     ----> write the spin-polarized density on unit nt
      REWIND nt
      CALL wrtdop(
     >            jspins,nq3,nq2,nmzxy,nmz,jmtd,nlhd,ntype,natd,
     >            jsp_new,nq3,nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,nt,
     >            'input   ','density ',iter,rho,qpw,rht,rhtxy,name)
      CLOSE (nt)
c
c     -----> This part is only used for testing th e magnetic moment in 
c     ----->   each sphere
c
      DO n = 1,ntype
          qtot1=0.00
          qtot2=0.00
          CALL intgr3(rho(1,0,n,1),rmsh(1,n),dx(n),jri(n),qtot1)
          CALL intgr3(rho(1,0,n,jsp_new),rmsh(1,n),dx(n),jri(n),qtot2)
          spmtot=sfp*(qtot1-qtot2)
          WRITE (6,'('' moment in sphere '',2x,'':'',f8.4)') spmtot
      ENDDO

c--->   read enpara and then double it
      OPEN(40,file='enpara',status='old',form='formatted')
      REWIND 40
      n=0
      DO
         READ (40,'(a)',iostat=ios) lineread
         IF (ios.ne.0) EXIT
         n = n+1
      ENDDO

      ALLOCATE (clines(n))

      REWIND 40
      DO i=1,n
         READ (40,'(a)') clines(i)
      ENDDO

      REWIND 40
      DO i=1,n
         WRITE (40,'(a)') trim(clines(i))
      ENDDO
      DO i=1,n
         WRITE (40,'(a)') trim(clines(i))
      ENDDO

      DEALLOCATE (clines,rho,qpw,rhtxy,rht)
      CLOSE(40)
!
! for lda+U: flip n-matrix
!
      IF (n_u.GT.0) THEN
        INQUIRE (file='n_mmp_mat',exist=n_exist)
        IF (n_exist) THEN
          OPEN (69,file='n_mmp_mat',status='old',form='formatted')
          REWIND 69

          n=0
          DO
             READ (69,'(a)',iostat=ios) lineread
             IF (ios.ne.0) EXIT
             n = n+1
          ENDDO
          ALLOCATE (clines(n))
          REWIND 69
          DO i=1,n
             WRITE (69,'(a)') trim(clines(i))
          ENDDO
          DO i=1,n
             WRITE (69,'(a)') trim(clines(i))
          ENDDO
          DEALLOCATE (clines)
    
          CLOSE(69)
        ENDIF
      ENDIF


      END SUBROUTINE cdnsp
      END MODULE m_cdnsp
