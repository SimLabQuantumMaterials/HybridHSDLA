      MODULE m_evaldos
      CONTAINS
      SUBROUTINE evaldos(
     >                   film,vacdos,nstars,nvac,nkpt,nkptd,ntype,layers
     >                   ,jspins,neigd,ntypd,layerd,neq,eminarg,
     >                emaxarg,sigmaarg,efermiarg,nop2,invs,l_noco,odi,
     >                   l_mcd,ncored,ncore,e_mcd,ndir,natd,nsld,bmat)
c----------------------------------------------------------------------
c
c     vk: k-vectors
c     wk weight of k-point (not used)
c     nevk no of eigenvalues
c     ev eigenvalue
c     qal partial charges
c             partial charge interstitial: qal(lmax*ntype+1...
c             partial charge vacuum      : qal(lmax*ntype+2...
c     qlay,qstar read in vacuum charges
c     qval partial charges in the vacuum
c             qval(m*nstars,neigd,nkptd):charge in m'th layer
c             qval(m*nstars+nstars,... ):star-resolved charge of that layer
c             qval(layers*nstars+....  ):same for second vacuum
c     ntb=max(nevk)
c
c----------------------------------------------------------------------
      USE m_triang
      USE m_maketetra
      USE m_tetrados
      USE m_dosbin
      USE m_ptdos
      USE m_smooth
      USE m_od_types, ONLY : od_inp
      USE m_cotra, ONLY : cotra3

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: neigd,nkptd,layerd,ntypd,ncored,nop2
      INTEGER, INTENT(IN) :: ndir,natd,nsld
      INTEGER, INTENT(IN) :: nstars,nvac,nkpt,ntype,layers,jspins
      REAL,    INTENT(IN) :: eminarg,emaxarg,sigmaarg,efermiarg
      LOGICAL, INTENT(IN) :: vacdos,film,l_noco,l_mcd,invs

      INTEGER, INTENT(IN) :: neq(ntype),ncore(ntype)
      REAL,    INTENT(IN) :: e_mcd(ntypd,jspins,ncored),bmat(3,3)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim 
c    locals
      INTEGER, PARAMETER :: lmax = 4, ned = 1301
      REAL,    PARAMETER :: factor = 27.2
      INTEGER  i,s,v,index,jspin,k,l,l1,l2,ln,n,nl,ntb,ntria,ntetra
      INTEGER  icore,qdim,n_orb
      REAL     as,de,efermi,emax,emin,qmt,sigma,totdos
      REAL     e_up,e_lo,e_test1,e_test2,fac,sumwei,dk
      LOGICAL  l_tria,l_orbcomp

      INTEGER  itria(3,2*nkpt),nevk(nkpt),itetra(4,6*nkpt)
      INTEGER, ALLOCATABLE :: ksym(:),jsym(:)
      REAL     vk(3,nkptd),wt(nkptd),voltet(6*nkpt),kx(nkpt),vkr(3,nkpt)
      REAL     ev(neigd,nkptd),e(ned),gpart(ned,ntype),atr(2*nkpt)
      REAL     e_grid(ned+1),spect(ned,3*ntypd),ferwe(neigd,nkptd)
      REAL,    ALLOCATABLE :: qal(:,:,:),qval(:,:,:),qlay(:,:,:),g(:,:)
      REAL,    ALLOCATABLE :: mcd(:,:,:),orbcomp(:,:,:),qmtp(:,:)
      REAL,    ALLOCATABLE :: qintsl(:,:),qmtsl(:,:),qvac(:,:)
      COMPLEX, ALLOCATABLE :: qstars(:,:,:,:)
      CHARACTER(len=2) :: spin12(2),ch_mcd(3)
      CHARACTER(len=8) :: chntype*2,chform*19
      DATA spin12/'.1' , '.2'/
      DATA ch_mcd/'.+' , '.-' , '.0'/
c
      IF  (ndir.NE.-3) THEN
        qdim = lmax*ntype+3
      ELSE
        qdim = 2*nsld 
        n_orb = 0
        l_orbcomp = .false.
        INQUIRE(file='orbcomp',exist=l_orbcomp)
        IF (l_orbcomp) THEN
          OPEN (1,file='orbcomp',form='formatted')
          READ (1,*) n_orb
          WRITE (*,*) 'DOS: orbcomp',n_orb
          CLOSE (1)
          qdim = 23
        ENDIF
      ENDIF
      ALLOCATE( qal(qdim,neigd,nkptd),
     +          qval(nstars*layers*nvac,neigd,nkptd),
     +          qlay(neigd,layerd,2),qstars(nstars,neigd,layerd,2))
      IF (l_mcd) ALLOCATE( mcd(3*ntypd*ncored,neigd,nkptd) )
c
c scale energies
      sigma = sigmaarg*factor
      emin = eminarg*factor
      emax = emaxarg*factor
      efermi = efermiarg*factor
 
      WRITE (6,'(a)') 'DOS-Output is generated!'

      IF ( NINT((emax - emin)/sigma) > ned ) THEN
        WRITE(6,*) 'sig_dos too small for DOS smoothing:'   
        WRITE(6,*) 'Reduce energy window or enlarge sig_dos!'
        WRITE(6,*) 'For now: setting sigma to zero !'
        sigma = 0.0
      ENDIF

      WRITE (6,*) 'sigma=   ' , sigma
      WRITE (6,*) 'emax=   ' , emax
      WRITE (6,*) 'emin=   ' , emin
      WRITE (6,*) 'ef_inp=   ' , efermi
c
c     create energy grid
      emax = emax - efermi
      emin = emin - efermi
      de = (emax-emin)/(ned-1)
      DO i=1,ned
         e(i) = emin + (i-1)*de
      ENDDO
 
      IF ( l_mcd ) THEN ! create an energy grid for mcd-spectra
        e_lo =  9.9d+9 
        e_up = -9.9d+9     
        DO jspin = 1,jspins
          DO n = 1,ntype
            DO icore = 1 , ncore(n)
              e_lo = min(e_mcd(n,jspin,icore),e_lo)
              e_up = max(e_mcd(n,jspin,icore),e_up)
            ENDDO
          ENDDO
        ENDDO
        e_lo = e_lo*factor - efermi - emax 
        e_up = e_up*factor - efermi
        de = (e_up-e_lo)/(ned-1)
        DO i=1,ned
          e_grid(i) = e_lo + (i-1)*de
          spect(i,:) = 0.0
        ENDDO
        e_grid(ned+1) = e_lo + ned*de
      ENDIF

      DO jspin = 1,jspins
         ntb = 0
         DO k = 1,nkpt
c
c     initialize arrays
c
            DO n = 1,neigd
               DO i = 1,qdim
                  qal(i,n,k) = 0.
               ENDDO
               DO i = 1,nstars*layers*nvac
                  qval(i,n,k) = 0.
               ENDDO
            ENDDO
c
c     read data from dos_tmp file!
c
            IF (( .not.l_mcd ).AND.(ndir.NE.-3)) THEN
            READ(84,rec=nkpt*(jspin-1)+k) vk(:,k),wt(k),nevk(k),ev(:,k),
     +                                    qal(1:lmax*ntype,:,k),          ! atoms 1...ntype
     +                                    qal(lmax*ntype+2,:,k),          ! vacuum 1
     +                                    qal(lmax*ntype+3,:,k),          ! vacuum 2
     +                                    qal(lmax*ntype+1,:,k),          ! interstitial
     +                                    qlay,qstars
            ELSEIF (ndir.NE.-3) THEN
            ALLOCATE( ksym(neigd),jsym(neigd) )
            READ(84,rec=nkpt*(jspin-1)+k) vk(:,k),wt(k),nevk(k),ev(:,k),
     +                                    qal(1:lmax*ntype,:,k),          ! atoms 1...ntype
     +                                    qal(lmax*ntype+2,:,k),          ! vacuum 1
     +                                    qal(lmax*ntype+3,:,k),          ! vacuum 2
     +                                    qal(lmax*ntype+1,:,k),          ! interstitial
     +                                    qlay,qstars,ksym,jsym,
     +                                    mcd(:,:,k)
            DEALLOCATE( ksym,jsym )
            ELSE  
             ALLOCATE( orbcomp(neigd,23,natd),qintsl(nsld,neigd))
             ALLOCATE( qmtsl(nsld,neigd),qmtp(neigd,natd),qvac(neigd,2))
             READ (129,rec=nkpt*(jspin-1)+k) vk(:,k),wt(k),nevk(k),
     +                      ev(:,k),qvac,qintsl,qmtsl,orbcomp,qmtp
             IF (n_orb == 0) THEN
               qal(1:nsld,:,k)        = qintsl(:,:)
               qal(nsld+1:2*nsld,:,k) = qmtsl(:,:)
             ELSE 
               DO i = 1, 23
                 DO l = 1, nevk(k)
                   qal(i,l,k) = orbcomp(l,i,n_orb)*qmtp(l,n_orb)/10000.
                 ENDDO
                 DO l = nevk(k)+1, neigd
                  qal(i,l,k) = 0.0
                 ENDDO 
               ENDDO
             ENDIF
             DEALLOCATE( orbcomp,qintsl,qmtsl,qmtp,qvac)
            ENDIF 
            ntb = max(ntb,nevk(k))
c
c     set vacuum partial charge zero, if bulk calculation
c     otherwise, write vacuum charge in correct arrays
c
            IF ((.NOT.film).AND.(ndir.NE.-3)) THEN
               DO n = 1,neigd
                  qal(lmax*ntype+2,n,k) = 0.0
                  qal(lmax*ntype+3,n,k) = 0.0
               ENDDO
            ELSEIF ( vacdos .and. film ) THEN
               DO i = 1,nevk(k)
                  DO v = 1,nvac
                     DO l = 1,layers
                        index = (l-1)*nstars + (v-1)*(nstars*layers) + 1
                        qval(index,i,k) = qlay(i,l,v)
                        DO s = 1,nstars - 1
                           qval(index+s,i,k) = real(qstars(s,i,l,v))
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
c
c     calculate interstitial dos if not noco
c     in the noco case, qis has been calculated in pwden and is read in from tmp_dos
c
            IF ((.NOT.l_noco).AND.(ndir.NE.-3)) THEN
               DO i = 1 , neigd
                  qal(lmax*ntype+1,i,k) = 1.
                  DO nl = 1 , ntype
                     l1 = lmax*(nl-1) + 1
                     l2 = lmax*nl
                     qmt=0.0
                     DO l = l1 , l2
                        qmt = qmt + qal(l,i,k)*neq(nl)
                     ENDDO
                     qal(lmax*ntype+1,i,k) = qal(lmax*ntype+1,i,k) - qmt
                  ENDDO
                  qal(lmax*ntype+1,i,k) = qal(lmax*ntype+1,i,k)
     &                          -qal(lmax*ntype+2,i,k)*(3-nvac)
     &                          -qal(lmax*ntype+3,i,k)*(nvac-1) 
               ENDDO
            ENDIF
c
c---- >     convert eigenvalues to ev and shift them by efermi
c
            DO i = 1 , nevk(k)
               ev(i,k) = ev(i,k)*factor - efermi
            ENDDO
            DO i = nevk(k) + 1, neigd
               ev(i,k) = 9.9e+99
            ENDDO
c
c
         ENDDO                                                 ! end of k-point loop
c
c     calculate the triangles!
c
         IF ( jspin.EQ.1 ) THEN
           l_tria=.true.
           IF (film .AND. .NOT.odi%d1) THEN
             CALL triang(vk,nkpt,itria,ntria,atr,as,l_tria)
             IF (invs) THEN
               IF (abs(nop2*as-0.5).GT.0.000001) l_tria=.false.
             ELSE
               IF (abs(nop2*as-1.0).GT.0.000001) l_tria=.false.
             ENDIF
             write(*,*) as,nop2,l_tria
             l_tria=.true.
           ELSE
             OPEN (41,file='kpts',FORM='formatted',STATUS='old')
             DO i = 1, nkpt+1
                READ (41,*)
             ENDDO
             READ (41,'(i5)',END=66,ERR=66) ntetra
             READ (41,'(4(4i6,4x))') ((itetra(i,k),i=1,4),k=1,ntetra)
             READ (41,'(4f20.13)') (voltet(k),k=1,ntetra)
             CLOSE(41)
             voltet(1:ntetra) = voltet(1:ntetra) / ntetra
             l_tria=.true.
             GOTO 67
 66          CONTINUE                       ! no tetrahedron-information of file
             CALL triang(vk,nkpt,itria,ntria,atr,as,l_tria)
             l_tria=.true.
c YM: tetrahedrons is not the way in 1D
             IF (odi%d1) as = 0.0         
             IF (invs) THEN
               IF (abs(nop2*as-1.0).GT.0.000001) l_tria=.false.
             ELSE
               IF (abs(nop2*as-0.5).GT.0.000001) l_tria=.false.
             ENDIF

             IF (l_tria) THEN
               CALL make_tetra(
     >                         nkptd,nkpt,vk,ntria,itria,atr,
     <                         ntetra,itetra,voltet)
             ELSE
               WRITE (6,*) 'no tetrahedron method with these k-points!'
               WRITE (6,*) nop2,as
             ENDIF
 67          CONTINUE                       ! tetrahedron-information read or created
           ENDIF
         ENDIF
c
        IF ( .not.l_mcd ) THEN
         ALLOCATE (g(ned,qdim))
        ELSE
         ALLOCATE (g(ned,3*ntypd*ncored))
        ENDIF
c
         IF ( l_tria.and.(.not.l_mcd).and.(ndir.NE.-3) ) THEN
c
c     DOS calculation: use triangular method!!
c
            IF ( film ) THEN
!             CALL ptdos(
!    >                  emin,emax,jspins,ned,qdim,neigd,
!    >                  ntria,as,atr,2*nkpt,itria,nkpt,ev,qal,e,
!    <                  g)
              CALL ptdos(
     >                   emin,emax,jspins,ned,qdim,ntb,ntria,as,
     >                   atr,2*nkpt,itria,nkpt,ev(1:ntb,1:nkpt),
     >                   qal(:,1:ntb,1:nkpt),e,
     <                   g)
            ELSE
              write(*,*) efermi
              CALL tetra_dos(
     >                       lmax,ntype,neigd,ned,ntetra,nkpt,
     >                       itetra,efermi,voltet,e,nevk,
     >                       ev,qal,
     <                       g)
              IF (jspins.EQ.1) g(:,:) = 2 * g(:,:)
            ENDIF
         ELSE
c
c     DOS calculation: use histogram method
c
            IF ( .not.l_mcd ) THEN
            CALL dos_bin(
     >                   jspins,qdim,ned,emin,emax,neigd,nkpt,
     >                   nevk,wt,ev,qal,
     <                   g)
            ELSE
            CALL dos_bin(
     >                   jspins,3*ntypd*ncored,ned,emin,emax,ntb,nkpt,
     >                   nevk(1:nkpt),wt(1:nkpt),ev(1:ntb,1:nkpt),
     >                   mcd(1:3*ntypd*ncored,1:ntb,1:nkpt),
     <                   g)
            ENDIF
         ENDIF
c
c---- >     smoothening
c
         IF ( .not.l_mcd ) THEN
            IF ( sigma.GT.0.0 ) THEN
              DO ln = 1 , qdim
                CALL smooth(e,g(1,ln),sigma,ned)
              ENDDO
            ENDIF
 
c*** sum up for all atoms
 
         IF (ndir.NE.-3) THEN
            DO l = 1 , ntype
               l1 = lmax*(l-1) + 1
               l2 = lmax*l
               DO i = 1 , ned
                  gpart(i,l) = 0.0
                  DO nl = l1 , l2
                     gpart(i,l) = gpart(i,l) + g(i,nl)
                  ENDDO
               ENDDO
            ENDDO
         ELSEIF (n_orb == 0) THEN
            DO l = 1 , nsld
               nl = nsld+l
               DO i = 1 , ned
                  gpart(i,l) = g(i,l) + g(i,nl)
               ENDDO
            ENDDO
         ENDIF
    
c**** write out DOS
         OPEN (18,FILE='DOS'//spin12(jspin))

         DO i = 1 , ned
           totdos = 0.0
           IF (ndir.NE.-3) THEN
             DO nl = 1 , ntype
                totdos = totdos + gpart(i,nl)*neq(nl)
             ENDDO
             totdos = totdos + g(i,lmax*ntype+1) + g(i,lmax*ntype+2) *
     +                        (3 - nvac) + g(i,lmax*ntype+3)*(nvac - 1)
             IF (ntype < 20) THEN
             WRITE (18,99001)  e(i),totdos,g(i,lmax*ntype+1), 
     +                        g(i,lmax*ntype+2),g(i,lmax*ntype+3),
     +                        (gpart(i,l),l=1,ntype), 
     +                        (g(i,l),l=1,ntype*lmax)
             ELSE
             WRITE (18,99001)  e(i),totdos,g(i,lmax*ntype+1), 
     +                        g(i,lmax*ntype+2),g(i,lmax*ntype+3),
     +                        (gpart(i,l),l=1,ntype)
             ENDIF
           ELSEIF (n_orb == 0) THEN
             DO nl = 1 , nsld
                totdos = totdos + gpart(i,nl)
             ENDDO
             WRITE (18,99001)  e(i),totdos,(gpart(i,nl),nl=1,nsld),
     +                                     (g(i,l),l=1,2*nsld)
           ELSE
             DO nl = 1 , 23
                totdos = totdos + g(i,nl)
             ENDDO
             WRITE (18,99001)  e(i),totdos,(g(i,l),l=1,23)
           ENDIF
         ENDDO
         CLOSE (18)

         ELSE
           write(*,'(4f15.8)') ((e_mcd(n,jspin,i),n=1,ntype),i=1,ncored)
           write(*,*)
           write(*,'(4f15.8)') (g(800,n),n=1,3*ntype*ncored)
           write(*,*)
           write(*,'(4f15.8)') (mcd(n,10,8),n=1,3*ntype*ncored)
           DO n = 1,ntype
             DO l = 1 , ned
               DO icore = 1 , ncore(n)
                 DO i = 1 , ned-1
                   IF (e(i).GT.0) THEN     ! take unoccupied part only
                   e_test1 = -e(i) - efermi +e_mcd(n,jspin,icore)*factor
                   e_test2 = -e(i+1)-efermi +e_mcd(n,jspin,icore)*factor
                   IF ((e_test2.LE.e_grid(l)).AND.
     +                 (e_test1.GT.e_grid(l))) THEN
                     fac = (e_grid(l)-e_test1)/(e_test2-e_test1)
                     DO k = 3*(n-1)+1,3*(n-1)+3
                       spect(l,k) = spect(l,k)+ g(i,3*ntypd*(icore-1)+k)
     +                      *(1.-fac) + fac * g(i+1,3*ntypd*(icore-1)+k)
                     ENDDO
                   ENDIF
                   ENDIF
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
           CLOSE (18) 
         ENDIF
         DEALLOCATE (g)
c         
c------------------------------------------------------------------------------
c     now calculate the VACOS
c------------------------------------------------------------------------------
            
         IF ( vacdos .and. film ) THEN
            ALLOCATE(g(ned,nstars*layers*nvac))
!            CALL ptdos(
!     >                 emin,emax,jspins,ned,nstars*nvac*layers,neigd,
!     >                 ntria,as,atr,2*nkpt,itria,nkptd,ev,qval,e,
!     <                 g)
            CALL ptdos(emin,emax,jspins,ned,nstars*nvac*layers,ntb,ntria
     &           ,as,atr,2*nkpt,itria,nkptd,ev(1:ntb,1:nkptd),
     &           qval(:,1:ntb,1:nkptd),e,g)
            
c---- >     smoothening
            IF ( sigma.GT.0.0 ) THEN
               DO ln = 1 , nstars*nvac*layers
                  CALL smooth(e,g(1,ln),sigma,ned)
               ENDDO
            ENDIF
            
c     write VACDOS
            
            OPEN (18,FILE='VACDOS'//spin12(jspin))
c            WRITE (18,'(i2,25(2x,i3))') Layers , (Zlay(l),l=1,Layers)
            DO i = 1 , ned
             WRITE (18,99001) e(i) , (g(i,l),l=1,Layers*Nstars*Nvac)
            ENDDO
            CLOSE (18)
            DEALLOCATE(g)
         ENDIF
c
c------------------------------------------------------------------------------
c     for bandstructures
c------------------------------------------------------------------------------

         IF (ndir == -4) THEN
            OPEN (18,FILE='bands'//spin12(jspin))
            ntb = minval(nevk(:))    
            kx(1) = 0.0
            CALL cotra3(vk(1,1),vkr(1,1),bmat)
            DO k = 2, nkpt
              CALL cotra3(vk(1,k),vkr(1,k),bmat)
              dk = ( vkr(1,k) - vkr(1,k-1) )**2 +
     +             ( vkr(2,k) - vkr(2,k-1) )**2 +
     +             ( vkr(3,k) - vkr(3,k-1) )**2
              kx(k) = kx(k-1) + sqrt(dk)
            ENDDO
            DO i = 1, ntb
              DO k = 1, nkpt
                write(18,'(2f12.6)') kx(k),ev(i,k)
              ENDDO
            ENDDO
            CLOSE (18)
         ENDIF

      ENDDO
c         
c------------------------------------------------------------------------------
c     for MCD calculations ...
c------------------------------------------------------------------------------

      IF (l_mcd) THEN
        WRITE (chntype,'(i2)') ntype+1
        chform = '('//chntype//'f15.8)'
        IF ( sigma.GT.0.0 ) THEN
           IF ( l_mcd ) THEN
             DO ln = 1 , 3*ntypd
               CALL smooth(e_grid,spect(1,ln),sigma,ned)
             ENDDO
           ENDIF
        ENDIF
        DO l = 1,3
          OPEN (18,FILE='MCD_SPEC'//ch_mcd(l))
          DO i = 1 , ned
          WRITE (18,FMT=chform) e_grid(i),(spect(i,3*(n-1)+l),n=1,ntype)
          ENDDO
          CLOSE (18)
        ENDDO
      ENDIF

      DEALLOCATE(qal,qval,qlay,qstars)
      IF (l_mcd) DEALLOCATE( mcd )
99001 FORMAT (f10.5,110(1x,d10.3))

      END SUBROUTINE evaldos
      END MODULE m_evaldos
