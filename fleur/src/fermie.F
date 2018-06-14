      MODULE m_fermie
!-----------------------------------------------------------------------
!     determines the fermi energy by
!            gaussian-integration method                          c.l.fu
!            triangular method (or tetrahedrons)
!            or fermi-function                                    p.kurz
!-----------------------------------------------------------------------
      CONTAINS
      SUBROUTINE fermie(
     > neigd,nkptd,jspd,nwdd,ntypd,lmaxd,nlod,
     > jspins,nwd,ntype,nkpt,ntapwf,irecl,form76,socfile,
     > gauss,delgau,zelec,tkb,tria,film,l_noco,l_ss,
     < ef,seigscv,ts,w,
     X qss,l_J,l_disp,bmat,nkpt_l,eig_l,irank,isize)

!-----------------------------------------------------------------------
!
!     a fist (T=0) approximation to the fermi-energy is determined
!     by:
!           zelec(nwd) = sum { spindg * we }
!                       e=<ef
!
!     TREE STRUCTURE: fermie----sort
!                           ----fergwt
!                           ----fertri---triang
!                                     |--dosint
!                                     |--dosef -- trisrt
!                                     +--doswt
!                           ----ferhis---ef_newton
!
!-----------------------------------------------------------------------

#ifdef CPP_HDF
      USE m_eig66_hdf, ONLY : read_eig, read_kptenpar
#endif
#ifdef CPP_MPI
      USE m_mpi_col_eigJ
#endif
      USE m_fergwt
      USE m_fertri
      USE m_ferhis
      USE m_fergwt
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nkptd,jspd,nwdd,ntypd,lmaxd,nlod
      INTEGER, INTENT (IN) :: jspins,nwd,ntype,irecl,ntapwf,nkpt_l
      LOGICAL, INTENT (IN) :: form76,socfile,gauss,tria,film
      LOGICAL, INTENT (IN) :: l_noco,l_ss,l_disp,l_J
      REAL,    INTENT (IN) :: delgau,tkb
      REAL,    INTENT (OUT):: ef,seigscv,ts
      INTEGER, INTENT (IN) :: irank,isize
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nkpt(nwdd)
      REAL,    INTENT (IN) :: zelec(nwdd)
      REAL,    INTENT (OUT):: w(neigd,nkptd,jspd)
      REAL,    INTENT (IN) :: eig_l(neigd+5,nkpt_l)
C     ..
C     .. Local Scalars ..
      REAL del,seigsc,seigv,spindg,ssc,weight,ws,zc
      INTEGER i,idummy,j,jsp,k,l,n,nbands,nstef,nv,nw,nrec,nmat,nspins
      INTEGER n_help
C     ..
C     .. Local Arrays ..
C
      INTEGER, ALLOCATABLE :: idxeig(:),idxjsp(:),idxkpt(:),index(:)
      REAL,    ALLOCATABLE :: e(:),eig(:,:,:),we(:)
      INTEGER ne(nkptd,jspd)
      REAL bk(3,nkptd),el(0:lmaxd,ntypd,jspd),evac(2,jspd)
      REAL wtkpt(nkptd),ello(nlod,ntypd,jspd)

c--- J constants
      REAL,    INTENT (INOUT) :: qss(3)
      REAL,    INTENT (INOUT) :: bmat(3,3)
c--- J constants

C     ..
C     .. External Subroutines ..
      EXTERNAL sort
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real
C     ..
c***********************************************************************
c                  ABBREVIATIONS
c
c     eig        : array of eigenvalues within all energy-windows
c     wtkpt      : list of the weights of each k-point (from inp-file)
c     e          : linear list of the eigenvalues within the highest
c                  energy-window
c     we         : list of weights of the eigenvalues in e
c     zelec      : number of electrons in a window
c     spindg     : spindegeneracy (2 in nonmagnetic calculations)
c     seigv      : weighted sum of the occupied valence eigenvalues
c     seigsc     : weighted sum of the semi-core eigenvalues
c     seigscv    : sum of seigv and seigsc
C     ts         : entropy contribution to the free energy
c     
c***********************************************************************
C     .. Data statements ..
      DATA del/1.0e-6/
C     ..
      ALLOCATE (idxeig(neigd*nkptd*jspd),idxjsp(neigd*nkptd*jspd),
     +          idxkpt(neigd*nkptd*jspd),index(neigd*nkptd*jspd),
     +          e(nkptd*neigd*jspd),eig(neigd,nkptd,jspd),
     +          we(nkptd*neigd*jspd) )
c
      IF (l_J) THEN
#ifdef CPP_MPI
        CALL mpi_col_eigJ(irank,isize,nkptd,neigd,nkpt(1),nkpt_l,eig_l,
     <                    bk,wtkpt,ne(1,1),eig)
        IF (irank.ne.0) RETURN
#endif
      ENDIF

      IF (socfile) THEN
        IF (nwd/=1) THEN
          STOP 'fermie: eig.soc and multiple windows'
        ENDIF
        OPEN(67,file='eig.soc',form='unformatted',action='read')
      ENDIF
      IF (form76) THEN
        OPEN (76,file='eig.form',form='formatted',status='unknown')
      ENDIF

      WRITE (6,FMT=8000)
 8000 FORMAT (/,/,1x,'fermi energy and band-weighting factors:')
c
c---> READ IN EIGENVALUES
c
      spindg = 2.0/real(jspins)
      n = 0
      nrec = 0
      seigsc = 0.0
      ssc = 0.0
      n_help = 0
c
c---> pk non-collinear
      IF (l_noco) THEN
         nspins = 1
      ELSE
         nspins = jspins
      ENDIF
c---> pk non-collinear
c
      DO 50 nw = 1,nwd
        DO 40 jsp = 1,nspins
          DO 30 k = 1,nkpt(nw)

  100        CONTINUE

             IF (socfile) THEN

               IF (jsp==1) THEN
                 READ(67) i
                 IF (i/=k) THEN
                   STOP 'fermie: error in eig.soc'
                 ENDIF
                 READ(67) bk(1,k), bk(2,k), bk(3,k)
                 READ(67) wtkpt(k)
                 READ(67) ne(k,jsp)
                 DO i= 1,ne(k,jsp)
                   READ(67) eig(i,k,jsp)
                 ENDDO
               ELSE
                 ne(k,jsp)= ne(k,1)
                 DO i= 1,ne(k,1)
                   eig(i,k,jsp)= eig(i,k,1)
                 ENDDO
               ENDIF

               nv= -1

             ELSE

               nrec = nrec + 1
#ifdef CPP_HDF
               IF (nw>1) STOP 'HDF and multiple windows..'
               CALL read_eig(k,jsp,ne(k,jsp),eig(:,k,jsp))
               CALL read_kptenpar(k,jsp,bk(1,k),wtkpt(k),el(:,:,jsp),
     +                            ello(:,:,jsp),evac(1,jsp))
               WRITE (6,'(a2,3f10.5,f12.6)') 'at',bk(:,k),wtkpt(k)
               WRITE (6,'(i5,a14)') ne(k,jsp),' eigenvalues :' 
               WRITE (6,'(8f12.6)') (eig(i,k,jsp),i=1,ne(k,jsp))
               nv= -1 
#else
               IF (l_J) THEN
!                  READ (ntapwf,rec=nrec) (bk(i,k),i=1,3),wtkpt(k),
!     +                          ne(k,jsp),(eig(i,k,jsp),i=1,neigd)
               ELSE
                 IF (l_ss) THEN
                  READ (ntapwf,rec=nrec) el,evac,ello,
     +                               (bk(i,k),i=1,3),wtkpt(k),
     +                                ne(k,jsp),nv,idummy,nmat,
     +                               (eig(i,k,jsp),i=1,neigd)
                 ELSEIF (l_noco) THEN
                  READ (ntapwf,rec=nrec) el,evac,ello,
     +                               (bk(i,k),i=1,3),wtkpt(k),
     +                                ne(k,jsp),nv,nmat,
     +                               (eig(i,k,jsp),i=1,neigd)
                 ELSE
                  READ (ntapwf,rec=nrec) el(:,:,jsp),evac(:,jsp),
     +                               ello(:,:,jsp),
     +                               (bk(i,k),i=1,3),wtkpt(k),
     +                                ne(k,jsp),nv,nmat,
     +                               (eig(i,k,jsp),i=1,neigd)
                 ENDIF
               ENDIF
#endif
             ENDIF ! (socfile) ELSE

!--> form76 start
             IF (form76) THEN
               IF (nv < 0) THEN
                 WRITE (76,FMT='(3f10.6,f12.8,i5)')
     &             bk(:,k),wtkpt(k),ne(k,jsp)
               ELSE
                 WRITE (76,FMT='(3f10.6,f12.8,2i5)')
     &             bk(:,k),wtkpt(k),ne(k,jsp),nv
               ENDIF
               WRITE (76,FMT='(6f15.6)') (eig(i,k,jsp),i=1,ne(k,jsp))
             ENDIF
!--> form76 end

             IF (nw.EQ.nwd) THEN
c
c--->          STORE EIGENVALUES AND WEIGHTS WITHIN THE HIGHEST
C--->          ENERGY WINDOW IN A LINEAR LIST. AND MEMORIZE  THE CONECTION
C--->          CONECTION TO THE ORIGINAL ARRAYS
c
               DO 10 j = 1,ne(k,jsp)
                 e(n+j) = eig(j,k,jsp)
                 we(n+j) = wtkpt(k)
                 idxeig(n+j) = j+n_help
                 idxkpt(n+j) = k
                 idxjsp(n+j) = jsp
   10          CONTINUE
C--->          COUNT THE NUMBER OF EIGENVALUES IN THE HIGHEST WINDOW
               n = n + ne(k,jsp)
             ELSE
               nbands = zelec(nw)/2.e0 + del
               IF (l_noco) nbands = nint(zelec(nw)) 
#ifdef CPP_MPI
                 DO 20 j = 1,min(nbands,ne(k,jsp))
#else
                 DO 20 j = 1,nbands
                   IF (j.GT.ne(k,jsp)) THEN
                     WRITE(6,8030) nw,k,ne(k,jsp)
                     STOP 'not enough eigenvalues in semicore'
                   ENDIF
#endif
                   ssc = ssc + spindg*wtkpt(k)
                   seigsc = seigsc + eig(j,k,jsp)*spindg*wtkpt(k)
   20            CONTINUE
               END IF
c
#ifdef CPP_MPI
               IF (.NOT.socfile) THEN 
                 n_help = n_help + ne(k,jsp)
                 IF (abs(evac(1,jsp)-999.9).LT.1.e-9) THEN
c
c                  obviously this record was only one part of a k-point
c
                   IF (gauss.OR.tria) STOP 'EV-parallel & gauss || tria' 
                   GOTO 100
                 ENDIF
                 ne(k,jsp) = n_help
                 n_help = 0
c                write(*,'(8f10.5)') (e(n-ne(k,jsp)+i),i=1,ne(k,jsp))
               ENDIF
#endif
c
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
 8030 FORMAT('There are not enough eigenvalues in a semicore window.',
     +     /,'In window=',i3,' k-point ',i5,' there are only ',i4,
     +       ' states.')

      IF (socfile) THEN
        CLOSE(67)
      ENDIF
      IF (form76) THEN
        CLOSE(76)
        STOP 'form66 complete'
      ENDIF

      IF (n.GT.nkptd*neigd*jspd) THEN
        WRITE (6,*) 'n = ',n,' > nkptd*neigd*jspd = ',nkptd*neigd*jspd
        STOP
      ENDIF
      CALL sort(n,e,index)
c
c---> DETERMINE EF BY SUMMING WEIGHTS
c
      weight = zelec(nwd)/spindg
      seigv = 0.0e0
      ws = 0.0e0
      l = 0
c      WRITE (6,*) zelec(nwd),spindg
c      WRITE (6,*) weight
c      WRITE (6,*)
c      DO 125 i = 1,n
c         WRITE (6,FMT='(i10,f12.8,f12.8)') i,e(index(i)),we(index(i))
c 125  CONTINUE
      DO WHILE ((ws+del).LT.weight)
         l = l + 1
         IF (l.GT.n) THEN
            WRITE (16,FMT=8010) n,ws,weight
            WRITE (6,FMT=8010) n,ws,weight
            STOP 'fermi'
 8010       FORMAT (/,10x,'error: not enough wavefunctions.',i10,
     +             2d20.10)
         END IF
         ws = ws + we(index(l))
         seigv = seigv + e(index(l))*we(index(l))*spindg
c         WRITE (6,FMT='(2f10.7)') e(index(l)),we(index(l))
      END DO
      ef = e(index(l))
      nstef = l
      zc = zelec(nwd)
      WRITE (6,FMT=8020) ef,nstef,seigv,ws,seigsc,ssc
 8020 FORMAT (/,'FERMIE:',/,
     +       10x,'first approx. to ef    (T=0)  :',f10.6,' htr',
     +       '   (energy of the highest occ. eigenvalue)',/,
     +       10x,'number of occ. states  (T=0)  :',i10,/,
     +       10x,'first approx. to seigv (T=0)  :',f10.6,' htr',/,
     +       10x,'sum of weights of occ. states :',f10.6,/,
     +       10x,'sum of semicore eigenvalues   :',f10.6,' htr',/,
     +       10x,'sum of semicore charge        :',f10.6,' e',/)
c+po
      ts = 0.0
c-po
      IF (gauss) THEN
         CALL fergwt(
     >               neigd,nkptd,jspd,nkpt(nwd),jspins,ne,
     >               delgau,zc,eig,wtkpt,
     X               ef,
     <               seigv,w)
      ELSE IF (tria) THEN
         CALL fertri(
     >               neigd,nkptd,jspd,
     >               ne,nkpt(nwd),nspins,zc,eig,bk,spindg,
     X               ef,
     <               seigv,w)
      ELSE
      nspins = jspins
      IF (l_noco) nspins = 1
         CALL ferhis(nspins,nkpt(nwd),tkb,index,idxeig,idxkpt,idxjsp,n,
     +               nstef,ws,weight,spindg,ef,seigv,ts,
     >               neigd,nkptd,jspd,e,ne,wtkpt,
     X               we,
     <               w,
     >               qss,l_J,l_disp,bmat)

      END IF
c     7.12.95 r.pentcheva seigscv must be calculated outside if (gauss)
      seigscv = seigsc + seigv
c
      DEALLOCATE ( idxeig,idxjsp,idxkpt,index,e,eig,we )
c
      RETURN
      END SUBROUTINE fermie
      END MODULE m_fermie
