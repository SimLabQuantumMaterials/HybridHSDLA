      MODULE m_cdnread
!------------------------------------------------------------------------
!     Two subroutines to read in the information stored on the 'eig'-file
!     cdn_read0: reads energy-parameters and number of bands stored/pe
!     cdn_read : reads finally the eigenvalues and vectors
!                                                          gb`04
!------------------------------------------------------------------------
      CONTAINS
      SUBROUTINE cdn_read0(
     >                     lmaxd,ntypd,nlod,neigd,jspd,
     >                     irank,isize,jspin,jsp_start,jsp_end,
     >                     l_noco,nrec,ntapwf,
     <                     ello,evac,epar,bkpt,wk,n_bands,nrec1,n_size)
#ifdef CPP_HDF
      USE m_eig66_hdf, ONLY : read_eig,read_enpar
#endif
      IMPLICIT NONE
!
! Arguments ...
!
      INTEGER, INTENT (IN) :: irank,isize,neigd,jsp_start
      INTEGER, INTENT (IN) :: jspd,jspin,ntapwf,jsp_end
      INTEGER, INTENT (IN) :: lmaxd,ntypd,nlod,nrec
      LOGICAL, INTENT (IN) :: l_noco
      INTEGER, INTENT (OUT) :: nrec1,n_size
      REAL,    INTENT (OUT) :: wk

      INTEGER, INTENT (OUT) :: n_bands(0:neigd)
      REAL,    INTENT (OUT) :: bkpt(3)
      REAL,    INTENT (OUT) :: ello(nlod,ntypd,jspd),evac(2,jspd)
      REAL,    INTENT (OUT) :: epar(0:lmaxd,ntypd,jspd)

      INTEGER isp

      nrec1 = nrec
#if ( defined(CPP_MPI) && !defined(CPP_T90) )
      IF (irank.NE.0) THEN 
        IF (nrec==0) THEN
        nrec1 = 0
        ELSE
        nrec1 = INT((nrec-(irank+1))/isize) + 1
        ENDIF
      ENDIF
#endif
      IF (l_noco) THEN
#ifdef CPP_HDF
         CALL read_eig(1,1,
     <                 n_bands(1))
!         DO isp = jsp_start,jsp_end
         DO isp = 1, jspd
           CALL read_enpar(1,isp,
     <                     epar(:,:,isp),ello(:,:,isp),evac(:,isp))
         ENDDO
#else
         READ (ntapwf,rec=nrec1+1) epar,evac,ello,bkpt,wk,n_bands(1)
#endif
      ELSE
#ifdef CPP_HDF
         CALL read_eig(1,jspin,
     <                 n_bands(1))
         CALL read_enpar(1,jspin,
     <                   epar(:,:,jspin),ello(:,:,jspin),evac(:,jspin))
#else
         READ (ntapwf,rec=nrec1+1) epar(:,:,jspin),evac(:,jspin),
     +                             ello(:,:,jspin),bkpt,wk,n_bands(1)
#endif
      ENDIF
c
c n_size is the number of records per k-point,
c n_bands(i) the number of ev's processed on n_rank=0...i-1
c
      n_size = 1
#if ( defined(CPP_MPI) && !defined(CPP_HDF) )
      IF (abs(evac(1,jspin)-999.9).LT.1.e-9) THEN
 97      n_size = n_size + 1
         IF (l_noco) THEN
            READ (ntapwf,rec=nrec1+n_size) epar,evac,ello,
     +                                     bkpt,wk,n_bands(0)
         ELSE
            READ (ntapwf,rec=nrec1+n_size) epar(:,:,jspin),
     +                  evac(:,jspin),ello(:,:,jspin),bkpt,wk,n_bands(0)
         ENDIF
         n_bands(n_size) = n_bands(n_size-1) + n_bands(0)
         IF (abs(evac(1,jspin)-999.9).LT.1.e-9) GOTO 97
      ENDIF
#ifndef CPP_T90
      IF (irank.NE.0) n_size = 1 ! isize
#endif
#endif
      n_bands(0) = 0
      END SUBROUTINE cdn_read0
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
      SUBROUTINE cdn_read(
     >                    lmaxd,ntypd,nlod,neigd,nvd,jspd,irank,isize,
     >                    ikpt,jspin,nbasfcn,nlotot,l_ss,l_noco,
     >                    nrec,i_rec,ntapwf,noccbd,n_start,n_end,
     <                    nmat,nv,ello,evdu,epar,kveclo,
     <                    k1,k2,k3,bkpt,wk,nbands,eig,z,cp_time)

      USE m_cputime
#ifdef CPP_HDF
      USE m_eig66_hdf, ONLY : read_eig,read_kpt,read_basis,
#if ( defined(CPP_SOC) || !defined(CPP_INVERSION) )
     +                        read_vecc
#else
     +                        read_vec
#endif
#endif

      IMPLICIT NONE
!
! Arguments ...
!
      INTEGER, INTENT (IN) :: irank,isize,ikpt,nbasfcn,neigd,nlotot
      INTEGER, INTENT (IN) :: nrec,i_rec,nvd,jspd,jspin,ntapwf
      INTEGER, INTENT (IN) :: lmaxd,ntypd,nlod,noccbd,n_start,n_end
      LOGICAL, INTENT (IN) :: l_ss,l_noco
      INTEGER, INTENT (OUT) :: nbands,nmat
      REAL,    INTENT (OUT) :: wk

      INTEGER, INTENT (OUT) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (OUT) :: nv(jspd),kveclo(nlotot)
      REAL,    INTENT (OUT) :: bkpt(3),eig(neigd)
      REAL,    INTENT (INOUT) :: cp_time(9)
      REAL,    INTENT (OUT) :: ello(nlod,ntypd,jspd),evdu(2,jspd)
      REAL,    INTENT (OUT) :: epar(0:lmaxd,ntypd,jspd)

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, INTENT (OUT) :: z(nbasfcn,noccbd)
#else
      REAL,    INTENT (OUT) :: z(nbasfcn,noccbd)
#endif
!
! Local variables ...
!
      INTEGER :: nrec1,iv,j,isp
      REAL    :: time1,time2
!
! find out the number of the record 
!
#if ( defined(CPP_MPI) && !defined(CPP_T90) && !defined(CPP_HDF) )
      IF (irank.EQ.0) THEN
        nrec1 = nrec + i_rec
      ELSE
        IF (nrec==0) THEN
          nrec1 = 0
        ELSE
          nrec1 = INT((nrec-(irank+1))/isize) + 1
        ENDIF
        nrec1 = nrec1 + (i_rec-irank-1)/isize + 1
      ENDIF
#else
      nrec1 = nrec + i_rec
#endif
!
! For Spin-Spirals
!
      CALL cpu_time(time1)
      IF (l_ss) THEN
#if (defined(CPP_HDF) &&  !defined(CPP_INVERSION) )
        CALL read_kpt(ikpt,1,
     <                bkpt,wk)
        CALL read_eig(ikpt,1,
     <                nbands,eig)
        DO isp = jspd,1,-1
          CALL read_basis(ikpt,isp,nlotot,
     <                    nmat,nv(isp),k1(:,isp),
     <                    k2(:,isp),k3(:,isp),kveclo)
!         write(*,*) kveclo
        ENDDO
        CALL read_vecc(
     >                 ikpt,1,nbands,n_start,n_end,
     <                 nmat,z)
#else
        READ (ntapwf,rec=nrec1) epar,evdu,ello,bkpt,wk,
     +             nbands,nv,nmat,eig,k1,k2,k3,kveclo,
     +             (z(:,j),j=1,noccbd)
#endif
!
! For Non-Collinear, but no Spin-Spirals
!
      ELSEIF (l_noco) THEN
#if (defined(CPP_HDF) &&  !defined(CPP_INVERSION) )
        CALL read_kpt(
     >                ikpt,1,
     <                bkpt,wk)
        CALL read_eig(ikpt,1,
     <                nbands,eig)
        CALL read_basis(ikpt,1,nlotot,
     <                  nmat,nv(1),k1(:,1),
     <                  k2(:,1),k3(:,1),kveclo)
        CALL read_vecc(
     >                 ikpt,1,nbands,n_start,n_end,
     <                 nmat,z)
#else
         READ (ntapwf,rec=nrec1) epar,evdu,ello,bkpt,wk,
     +             nbands,nv(1),nmat,eig,k1(:,1),k2(:,1),
     +             k3(:,1),kveclo,(z(:,j),j=1,noccbd)
#endif
         nv(jspd) = nv(1)
         DO iv = 1,nv(1)
            k1(iv,jspd) = k1(iv,1)
            k2(iv,jspd) = k2(iv,1)
            k3(iv,jspd) = k3(iv,1)
         ENDDO
!
! For Collinear
!
      ELSE
#ifdef CPP_HDF
         CALL read_kpt(ikpt,jspin,
     <                 bkpt,wk)
         CALL read_eig(ikpt,jspin,
     <                 nbands,eig)
         CALL read_basis(ikpt,jspin,nlotot,
     <                   nmat,nv(jspin),k1(:,jspin),
     <                   k2(:,jspin),k3(:,jspin),kveclo)
#if ( defined(CPP_SOC) || !defined(CPP_INVERSION) )
         CALL read_vecc(
     >                  ikpt,jspin,nbands,n_start,n_end,
     <                  nmat,z(:nmat,:noccbd))
#else
         CALL read_vec(
     >                  ikpt,jspin,nbands,n_start,n_end,
     <                  nmat,z(:nmat,:noccbd))
#endif
#else
         READ (ntapwf,rec=nrec1) epar(:,:,jspin),evdu(:,jspin),
     +             ello(:,:,jspin),bkpt,wk,nbands,nv(jspin),nmat,eig,
     +           k1(:,jspin),k2(:,jspin),k3(:,jspin),kveclo,
     +           (z(:,j),j=1,noccbd)
#endif
      ENDIF
      CALL cpu_time(time2)
      cp_time(3) = cp_time(3) + time2 - time1
      ! write(*,*) 'pe',irank,' read ',nrec1,nrec

      IF (nbands.GT.neigd)  STOP 'cdn_read: nbands.GT.neigd'
      IF (nv(jspin).GT.nvd) STOP 'cdn_read: nv.GT.nvd'
      IF (nmat.GT.nbasfcn)  STOP 'cdn_read: nmat.GT.nbasfcn'

      END SUBROUTINE cdn_read
      END MODULE m_cdnread
