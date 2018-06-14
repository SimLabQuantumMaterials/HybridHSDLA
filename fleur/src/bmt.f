      SUBROUTINE bmt(
     > nq3,nq2,jspins,l_noco,jri,nlh,ntype,nmz,nmzxy,
     > ntypsd,ntypsy,z1,delz,invs,invs2,film,
     > natd,neq,nvac,zatom,namat,dx,rmt)
c
      USE m_loddop
      USE m_wrtdop
      IMPLICIT NONE
C     ..
      INTEGER, INTENT (IN) :: nq3,nq2,jspins,ntype,natd
      INTEGER, INTENT (IN) :: nmz,nmzxy,ntypsd,nvac
      LOGICAL, INTENT (IN) :: invs,invs2,film,l_noco
      REAL,    INTENT (IN) :: z1,delz
      INTEGER, INTENT (IN) :: jri(ntype),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ntypsy(natd),neq(ntype)
      REAL,    INTENT (IN) :: rmt(ntype),dx(ntype),zatom(ntype)
      CHARACTER(len=2), INTENT (IN):: namat(0:103)
C     ..
      INTEGER k,i,ivac,jmtd,nlhd,it 
      INTEGER type,typmag 
      CHARACTER(len=8) dop,iop,name(10),filename 
      COMPLEX, ALLOCATABLE :: fpw(:,:),fzxy(:,:,:,:)
      REAL,    ALLOCATABLE :: fz(:,:,:),fr(:,:,:,:)
C     ..
C     ..


      typmag= ntype 
      ! only muffin-tins with type <= typmag remain magnetic  


      IF (jspins/=2) THEN
        STOP 'Stop in bmt:  jspins/=2'
      ENDIF 

      jmtd = maxval(jri(:))
      nlhd = maxval(nlh(:))

      ALLOCATE( fpw(nq3,jspins),fzxy(nmzxy,nq2-1,2,jspins) )
      ALLOCATE( fr(jmtd,0:nlhd,ntype,jspins),fz(nmz,2,jspins) )

      IF (l_noco) THEN 
        OPEN (98,file='rhomat_inp',form='unformatted',action='read')
      ELSE
        OPEN (98,file='cdn1',form='unformatted',action='read')
      ENDIF 
      CALL loddop(
     > jspins,nq3,nq2,nmzxy,nmz,jmtd,nlhd,ntype,
     > jspins,nq3,nq2,nvac,ntype,invs,invs2,film,
     > nlh,jri,ntypsd,ntypsy,98,natd,neq,
     < iop,dop,it,fr,fpw,fz,fzxy,name)
      CLOSE (98)

      IF ( typmag < ntype ) THEN 
        DO type= typmag+1,ntype 
          DO k= 0,nlhd
            DO i= 1,jmtd
              fr(i,k,type,1)= ( fr(i,k,type,1) + fr(i,k,type,2) )/2. 
              fr(i,k,type,2)= fr(i,k,type,1) 
            ENDDO
          ENDDO 
        ENDDO
      ENDIF 

      DO k= 1,nq3
        fpw(k,1)= ( fpw(k,1) + fpw(k,2) )/2.
        fpw(k,2)= fpw(k,1)
      ENDDO
      IF (film) THEN
        DO ivac= 1,nvac
          DO i= 1,nmz
            fz(i,ivac,1)= ( fz(i,ivac,1) + fz(i,ivac,2) )/2.  
            fz(i,ivac,2)= fz(i,ivac,1)
          ENDDO
          DO k= 2,nq2
            DO i= 1,nmzxy
              fzxy(i,k-1,ivac,1)= 
     &         ( fzxy(i,k-1,ivac,1) + fzxy(i,k-1,ivac,2) )/2. 
              fzxy(i,k-1,ivac,2)= fzxy(i,k-1,ivac,1)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      filename= 'cdnbmtXX'
      IF ( typmag < ntype ) THEN
        filename(7:7)= ACHAR(IACHAR('0')+ MOD(typmag,100)/10 )
        filename(8:8)= ACHAR(IACHAR('0')+ MOD(typmag,10) )
        OPEN(98,file=filename(1:8),form='unformatted',status='replace')
      ELSE
        OPEN(98,file=filename(1:6),form='unformatted',status='replace')
      ENDIF 
      CALL wrtdop(
     > jspins,nq3,nq2,nmzxy,nmz,jmtd,nlhd,ntype,natd,
     > jspins,nq3,nq2,nmzxy,nmz,nvac,ntype,neq,
     > invs,invs2,film,delz,z1,dx,rmt,zatom,
     > nlh,jri,ntypsd,ntypsy,namat,98,
     > iop,dop,it,fr,fpw,fz,fzxy,name)
      CLOSE(98) 

      DEALLOCATE( fpw,fzxy,fr,fz)

      END SUBROUTINE bmt 
