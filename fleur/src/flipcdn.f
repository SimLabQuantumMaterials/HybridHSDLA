      MODULE m_flipcdn
c     *******************************************************
c     this subroutine reads the charge density and flips the 
c     magnetic moment within the m.t.sphere for each atom 
c     according to the variable nflip. This variable is read in
c     the main program
c             nflip = -1 : flip spin in sphere
c             nflip = -2 : scale spin by bmu(n)
c             nflip = any: no spin flip
c                            r.pentcheva,kfa,Feb'96
c     *******************************************************
      CONTAINS
      SUBROUTINE flipcdn(
     >                   ntype,jspins,nmz,nmzxy,ntypsd,ntypsy,
     >                   nvac,nq2,nq3,film,invs,invs2,l_noco,z1,delz,
     >                   jri,nlh,zatom,rmt,dx,natd,neq,nwd,
     >                   nflip,lda_u,n_u,bmu,namat,nlo)
      USE m_loddop
      USE m_wrtdop

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments from calling program ..
      INTEGER, INTENT (IN) :: ntype,nmzxy,nq2,nq3
      INTEGER, INTENT (IN) :: nmz,ntypsd,natd
      INTEGER, INTENT (IN) :: jspins,nvac,nwd,n_u
      REAL,    INTENT (IN) :: z1,delz
      LOGICAL, INTENT (IN) :: film,invs,invs2,l_noco 
C     ..
C     .. Array Arguments from calling program..
      INTEGER,     INTENT (IN) :: nflip(ntype),jri(ntype),neq(ntype)
      INTEGER,     INTENT (IN) :: nlh(ntypsd),ntypsy(natd),nlo(ntype)
      REAL,        INTENT (IN) :: zatom(ntype),rmt(ntype),dx(ntype)
      CHARACTER*2, INTENT (IN) :: namat(0:103)
      INTEGER,     INTENT (IN) :: lda_u(ntype)
      REAL,        INTENT (IN) :: bmu(ntype)
C     ..
C     .. Local Scalars ..
      REAL    rhodummy,rhodumms
      INTEGER i,iter,n,nt,j,lh,na,m,mp,ispin,n_ldau,nw,urec,itype
      INTEGER jmtd,nlhd
      CHARACTER(len=8) iop,dop
      LOGICAL n_exist
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: n_mmp(:,:,:,:),qpw(:,:),rhtxy(:,:,:,:)
      REAL   , ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:)
      COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)
      CHARACTER(len=80), ALLOCATABLE :: clines(:)
      CHARACTER(len=8) name(10)
C     ..
      jmtd = maxval(jri(:))
      nlhd = maxval(nlh(:))
      ALLOCATE ( rho(jmtd,0:nlhd,ntype,jspins),qpw(nq3,jspins) )
      ALLOCATE ( rhtxy(nmzxy,nq2-1,2,jspins),rht(nmz,2,jspins) )
      IF (l_noco) THEN
        ALLOCATE( cdom(nq3) )
        ALLOCATE( cdomvz(nmz,2),cdomvxy(nmzxy,nq2-1,2) )
      ENDIF 

      nt = 71                   !gs, see sub mix
      IF (l_noco) THEN
        OPEN (nt,file='rhomat_inp',form='unformatted',status='old')
      ELSE
        OPEN (nt,file='cdn1',form='unformatted',status='old')
      ENDIF 
c     ---> read the charge density 
      CALL loddop(
     >            jspins,nq3,nq2,nmzxy,nmz,jmtd,nlhd,ntype,
     >            jspins,nq3,nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nt,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
      IF (l_noco) THEN
         READ (nt) (cdom(j),j=1,nq3)
         IF (film) THEN
           READ (nt) ((cdomvz(n,i),n=1,nmz),i=1,nvac)
           READ (nt) (((cdomvxy(n,j-1,i),n=1,nmzxy),j=2,nq2),i=1,nvac) 
         ENDIF
      ENDIF
c     ---> flip cdn for each atom with nflip=-1
c
      na = 1
      DO n = 1, ntype
         IF (nflip(n).EQ.-1) THEN
c     ---> spherical and non-spherical m.t. charge density
            DO lh = 0,nlh(ntypsy(na))
               DO j = 1,jri(n)
                  rhodummy = rho(j,lh,n,1)
                  rho(j,lh,n,1) = rho(j,lh,n,jspins)
                  rho(j,lh,n,jspins) = rhodummy
               ENDDO
            ENDDO
         ELSEIF (nflip(n).EQ.-2) THEN
            DO lh = 0,nlh(ntypsy(na))
               DO j = 1,jri(n)
                  rhodummy = rho(j,lh,n,1) + rho(j,lh,n,jspins)
                  rhodumms = rho(j,lh,n,1) - rho(j,lh,n,jspins)
                  rho(j,lh,n,1) = 0.5 * ( rhodummy + bmu(n) * rhodumms )
                  rho(j,lh,n,jspins) = 0.5*(rhodummy - bmu(n)*rhodumms )
               ENDDO
            ENDDO 
         END IF
         na = na + neq(n)
      ENDDO
c     ----> write the spin-polarized density on unit nt
      REWIND nt
      CALL wrtdop(
     >            jspins,nq3,nq2,nmzxy,nmz,jmtd,nlhd,ntype,natd,
     >            jspins,nq3,nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,nt,
     >            'input   ','density ',iter,rho,qpw,rht,rhtxy,name)
      IF (l_noco) THEN
         WRITE (nt) (cdom(j),j=1,nq3)
         IF (film) THEN
           WRITE (nt) ((cdomvz(n,i),n=1,nmz),i=1,nvac)
           WRITE (nt) (((cdomvxy(n,j-1,i),n=1,nmzxy),j=2,nq2),i=1,nvac)
         ENDIF
      ENDIF
      CLOSE (nt)
!
! for lda+U: flip n-matrix 
!
      IF (n_u.GT.0) THEN
        INQUIRE (file='n_mmp_mat',exist=n_exist)
        IF (n_exist) THEN
          OPEN (69,file='n_mmp_mat',status='old',form='formatted')
          ALLOCATE (  n_mmp(-3:3,-3:3,n_u,2) )
          DO nw = 1,nwd

            READ (69,9000) n_mmp
c   flip    ...
            n_ldau = 0
            DO n = 1,ntype
              IF (lda_u(n).GE.0) THEN
                n_ldau = n_ldau + 1
                IF (nflip(n).eq.-1) THEN
                  DO m = -3,3
                    DO mp = -3,3
                      rhodummy = n_mmp(m,mp,n_ldau,1)
                      n_mmp(m,mp,n_ldau,1) = n_mmp(m,mp,n_ldau,jspins)
                      n_mmp(m,mp,n_ldau,jspins) = rhodummy
                    ENDDO
                  ENDDO
                ELSEIF (nflip(n).eq.-2) THEN
                  DO m = -3,3
                    DO mp = -3,3
                      rhodummy = n_mmp(m,mp,n_ldau,1) + 
     +                           n_mmp(m,mp,n_ldau,jspins)
                      rhodumms = n_mmp(m,mp,n_ldau,1) - 
     -                           n_mmp(m,mp,n_ldau,jspins)
                      n_mmp(m,mp,n_ldau,1) = 0.5 * ( rhodummy + 
     +                                      bmu(n) * rhodumms )
                      n_mmp(m,mp,n_ldau,jspins) = 0.5*( rhodummy - 
     -                                         bmu(n) * rhodumms )
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
c   flip    ...
            REWIND (69)
            WRITE (69,9000) n_mmp
 9000       FORMAT(7f20.13)
c
          ENDDO
          DEALLOCATE ( n_mmp )
        ENDIF
      ENDIF
c-lda+U
c
c--->   read enpara and  flip lines
c
      OPEN(40,file='enpara',status='old',form='formatted')

      n = 2
      DO itype = 1 , ntype
        n = n + 1
        IF (nlo(itype).GT.0) n = n + 2
      ENDDO 
      IF (film) n = n + 1
      ALLOCATE (clines(2*n))
      DO i=1,2*n
         READ (40,'(a)') clines(i)
      ENDDO

      REWIND 40
      i = 0 
      DO ispin = 1,jspins
        i = i + 2
        WRITE (40,'(a)') trim(clines(i-1))
        WRITE (40,'(a)') trim(clines(i))
        DO itype = 1 , ntype
          i = i + 1
          m = i
          IF (nflip(itype).eq.-1) m = mod(i+n,2*n)
          IF (m.EQ.0) m = 2*n
          WRITE (40,'(a)') trim(clines(m))
          IF (nlo(itype).GT.0) THEN
            WRITE (40,'(a)') trim(clines(m+1))
            WRITE (40,'(a)') trim(clines(m+2))
            i = i + 2
          ENDIF
        ENDDO
        IF (film) THEN
          i = i + 1
          WRITE (40,'(a)') trim(clines(i))
        ENDIF
      ENDDO

      DEALLOCATE (clines,rho,qpw,rhtxy,rht)
      IF (l_noco) THEN
        DEALLOCATE (cdom,cdomvz,cdomvxy)
      ENDIF 
      CLOSE(40)

      END SUBROUTINE flipcdn
      END MODULE m_flipcdn
