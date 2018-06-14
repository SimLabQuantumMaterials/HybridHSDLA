      MODULE m_Ekwritesl
      CONTAINS
      SUBROUTINE Ek_write_sl(
     >                       neigd,nkptd,ntypd,jspd,layerd,natd,nsld,
     >                       nwdd,jspins,nkpt,ntype,nvac,jspin,nw,
     >                       invs,zrfs,bmat,neq,
     >                       nsl,nslat)
c-----------------------------------------------------------------
c-- now write E(k) for all kpts if on T3E
c-- now read data from tmp_dos and write of E(k) in  ek_orbcomp
c-----------------------------------------------------------------
      IMPLICIT NONE
C	..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nkptd,ntypd,jspd,layerd,natd,nsld
      INTEGER, INTENT (IN) :: jspins,nvac,nkpt,nsl,ntype,jspin,nwdd,nw
      LOGICAL, INTENT (IN) :: invs,zrfs
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),nslat(natd,nsld)
      REAL,    INTENT (IN) :: bmat(3,3)
c     ..
c     .. Local Scalars
      INTEGER :: nbands,ikpt,kspin,j,i,n,it,m,na,iband,mt,l
      INTEGER :: ivac
      REAL    :: wk
c     ..
c     .. Local Arrays
      INTEGER  norb(23),iqsl(nsld),iqvacpc(2)
      REAL     bkpt(3),qvact(2)
      REAL, ALLOCATABLE :: eig(:),qvac(:,:,:,:),orbcomp(:,:,:,:,:)
      REAL, ALLOCATABLE :: qintsl(:,:,:,:),qmtsl(:,:,:,:),qmtp(:,:,:,:)
      CHARACTER (len=2) :: chntype
      CHARACTER (len=99) :: chform
c     ..
c     .. Intrinsic Functions
      INTRINSIC  nint 
c
      IF (nsl.GT.nsld)  THEN
          write(6,*)' STOP Ek:   nsl.GT.nsld  '
          STOP 
      ENDIF     
      ALLOCATE(eig(neigd),orbcomp(neigd,23,natd,nkptd,jspd))
      ALLOCATE(qvac(neigd,2,nkptd,jspd),qintsl(nsld,neigd,nkptd,jspd))
      ALLOCATE(qmtsl(nsld,neigd,nkptd,jspd),qmtp(neigd,natd,nkptd,jspd))
c
c  --->     open files for a bandstucture with an orbital composition
c  --->     in the case of the film geometry
c
      IF (nw.EQ.1)  THEN 
          IF (jspin.EQ.1)  OPEN (130,file='ek_orco_11') 
          IF (jspin.EQ.2)  OPEN (130,file='ek_orco_12')
      ENDIF
      IF (nw.EQ.2)  THEN 
	  IF (jspin.EQ.1)  OPEN (130,file='ek_orco_21') 
          IF (jspin.EQ.2)  OPEN (130,file='ek_orco_22')
      ENDIF
      IF (jspin.GT.2)  WRITE (130,FMT=900)
c
c ----->       write bandstructure to ek_orbcomp - file
c 
      WRITE (chntype,'(i2)') nsl
      chform = "('E',i3,'= ',f10.4,4x,'vac ( layers ) vac = ',i3,' ('
     +        ,"//chntype//"(i3,2x),')',i3))"
      WRITE (130,FMT=901) 
      WRITE (130,FMT=902) 
      WRITE (130,FMT=901) 
      WRITE (130,FMT=800) nwdd,nw
      WRITE (130,FMT=903) nsl,nvac,nkpt
      WRITE (130,FMT=904) ntype,(neq(n),n=1,ntype) 
      WRITE (130,FMT=805)  
      DO j=1,nsl
	 WRITE (130,FMT=806) j,(nslat(i,j),i=1,natd)
      ENDDO        
      DO kspin = 1,jspins
	   WRITE (130,FMT=907)  kspin,jspins
c============================================================== 
 800     FORMAT (5X,' nwdd =',i2,'   nw =',i2,/)
 900     FORMAT (5X,'----- !!!!!!!!  nw.GT.2  !!!!!!! -----')
 901     FORMAT (5X,'--------------------------------------')
 902     FORMAT (5X,'-------- E(k) for a film  ------------')
 903     FORMAT (5X,' nsl =',i3,'   nvac =',i2,'   nkpt =',i4,/)
 904     FORMAT (5X,' ntype = ',i3,' neq(n) = ',50i3)
 907     FORMAT (/,5X,' kspin = ',i4,' jspins = ',i4)
 805     FORMAT (5X,'  nsl   nslat(1:nate,nsli)  ')
 806     FORMAT (5X,51i4)
c==============================================================
         DO ikpt=1,nkpt
c                
            READ (129,rec=nkpt*(kspin-1)+ikpt)  bkpt,wk,nbands,eig,
     +            qvac(:,:,ikpt,kspin),qintsl(:,:,ikpt,kspin),
     +            qmtsl(:,:,ikpt,kspin),
     +            orbcomp(:,:,:,ikpt,kspin),qmtp(:,:,ikpt,kspin)
!            write(*,*) kspin,nkpt,qmtp(1,:,ikpt,kspin)
c
            WRITE (130,FMT=8000) (bkpt(i),i=1,3)
 8000       FORMAT (/,3x,'  k =',3f10.5,/)
c
            DO iband = 1,nbands
               qvact = 0.0
               DO ivac = 1,nvac
                  qvact(ivac) = qvac(iband,ivac,ikpt,kspin)
               ENDDO
               IF (invs .OR. zrfs)    qvact(2) = qvact(1)
               iqvacpc(:) = nint(qvact(:)*100.0)
               DO j = 1,nsl
                  iqsl(j) = nint( ( qintsl(j,iband,ikpt,kspin) + 
     +                               qmtsl(j,iband,ikpt,kspin) )*100.0 ) 
               ENDDO         
               WRITE (130,FMT=chform) iband,eig(iband),iqvacpc(2),
     +                               (iqsl(l),l=1,nsl),iqvacpc(1)
       	       WRITE(130,FMT=9) 
	       WRITE(130,FMT=8)
	       WRITE(130,FMT=9) 
 	       DO n = 1,nsl
	          mt=0 
                  DO 1 it=1,ntype
	               DO  m=1,neq(it)
		           mt=mt+1	
	                   na = nslat(mt,n) 
	                   IF (na.EQ.1) THEN
		               DO  j=1,23
                                   norb(j) = 
     +                nint ( orbcomp(iband,j,mt,ikpt,kspin) )
	                       ENDDO	
                               WRITE (130,FMT=5) n,it,m,
     +	   		                  (norb(l),l=1,23),
     +                                    qmtp(iband,mt,ikpt,kspin)
                            ENDIF
	               ENDDO
  1	          CONTINUE  
	       ENDDO              ! over ( n = 1,nsl ) 
	       WRITE(130,FMT=9) 
            ENDDO           ! over ( iband = 1,nbands ) 
         ENDDO        ! over ( ikpt=1,nkpt )
      ENDDO	  ! over ( kspin = 1,jspins )  
        CLOSE (130)
c
! 8040 FORMAT ('E',i3,'= ',f10.4,4x,'vac | layers | vac = ',i3,' | ',50(i3,2x),' | ',i3)
 8    FORMAT('|lyr,tp,at| S | Px  Py  Pz | Dxy  Dyz  Dzx  Dx-y Dz2 |',
     +  ' Fx3  Fy3  Fz3  Fx2y Fy2z Fz2x Fxyz| Fz2x Fz2y Fz3  Fxyz Fx2z',
     +  ' Fx3  Fy3 |  mt  |') 
 5    FORMAT('|',i3,',',i2,',',i2,'|',i3,'|',3(i3,1x),'|',
     +        5(1x,i3,1x),'|',
     +        7(1x,i3,1x),'|',7(1x,i3,1x),'|',f6.1,'|')
 9    FORMAT(133('-'))
c
      DEALLOCATE ( eig,qvac,orbcomp,qintsl,qmtsl,qmtp )

      END SUBROUTINE Ek_write_sl
      END MODULE m_Ekwritesl
