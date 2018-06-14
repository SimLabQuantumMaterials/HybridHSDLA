      MODULE m_setinp
!---------------------------------------------------------------------
!  Check muffin tin radii and determine a reasonable choice for MTRs.
!  Derive also other parameters for the input file, to provide some
!  help in the out-file.                                        gb`02
!---------------------------------------------------------------------
      CONTAINS
      SUBROUTINE set_inp(
     >                   ntype,nop,nop2,nat,omtil,title,
     >                   film,invs,invs2,zrfs,amat,dvac,l_soc,
     >                   neq,zatom,taual,l_ss,a1,a2,a3,theta,phi,
     <                   rmt1)

      USE m_chkmt
      USE m_rwinp
      USE m_types, ONLY : t_utype
      USE m_od_types, ONLY : od_dim
      IMPLICIT NONE

      INTEGER, INTENT (INOUT) :: ntype,nop,nat,nop2
      LOGICAL, INTENT (IN) :: film,invs,invs2,zrfs,l_ss,l_soc
      INTEGER, INTENT (INOUT) :: neq(ntype)
      REAL,    INTENT (IN) :: dvac,theta,phi,omtil 
      REAL,    INTENT (IN) :: amat(3,3),zatom(ntype)
      REAL,    INTENT (INOUT) :: a1(3),a2(3),a3(3),taual(3,nat)
      REAL,    INTENT (OUT) :: rmt1(ntype)
      CHARACTER(len=80), INTENT (IN) :: title
 
      INTEGER, PARAMETER :: nwdd = 1
      INTEGER nel,i,j
      REAL    kmax,dtild,dvac1,n1,n2,gam
      LOGICAL l_test,l_gga,strho
      REAL    pos(3,nat),rmt(ntype),dx(ntype),zelec(nwdd)
      INTEGER nz(ntype),nlo(ntype),llo(2,ntype),ncst(ntype),lmax(ntype)
      INTEGER jri(ntype)
      CHARACTER(len=1) :: ch_rw
      CHARACTER(len=3) :: namat(0:103),noel(ntype),latnam
      CHARACTER(len=4) :: namgrp,namex
      CHARACTER(len=8) :: name(10)
      CHARACTER(len=12) :: relcor
      LOGICAL stden,dos,secvar,vchk,cdinf,pot8,form66,ctail
      LOGICAL l_noco,l_J,l_u2f,l_f2u,l_bmt
      LOGICAL eonly,gauss,tria,frcor,slice,disp,starcoeff
      LOGICAL swsp,lflip,vacdos,integ,iplot,score,plpot,pallst,lwb,l_f
      INTEGER layerd,nlod,nu,iofile,nwd,jspins,lpr,lepr,kcrel,nkpt
      INTEGER itmax,imix,maxiter,kk,nnne,nstars,nstm,isec1,ndir
      INTEGER igrd,ndvgrd,idsprs,iggachk,idsprs0,idsprsl,idsprsi,idsprsv
      INTEGER ntypd,natd,layers,n,gw,gw_neigd,iostat
      REAL    alpha,spinf,xa,thetad,epsdisp,epsforce,gmax,gmaxxc
      REAL    e1s,e2s,e1_dos,e2_dos,sig_dos,tworkf,chng,area
      REAL    scale,scpos,delgau,zc,tkb,rkm(nwdd)

      INTEGER lnonsph(ntype),nflip(ntype),izlay(1,2),relax(3,ntype)
      REAL    locx(2),locy(2),bmu(ntype),ellow(nwdd),elup(nwdd)
      LOGICAL l_geo(ntype),eig66(2),soc_opt(ntype+2)
      TYPE (t_utype) :: lda_u(ntype)
c-odim
      TYPE (od_dim) :: odd
c+odim
      REAL, PARAMETER :: eps=0.00000001
C     ..
C     .. DATA statements ..
      DATA namat/'va',' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     +     'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',
     +     ' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',
     +     'Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
     +     'Ag','Cd','In','Sn','Sb','Te',' J','Xe','Cs','Ba','La','Ce',
     +     'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     +     'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb',
     +     'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu',
     +     'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lw'/

      l_test = .false.
      l_gga  = .true.
      nz(:) = NINT(zatom(:))
      noel(:) = namat(nz(:))
      rmt(:) = 999.9
      pos(:,:) = matmul( amat , taual(:,:) )
      relcor = 'non-relativi' 
      ch_rw = 'w'
      namgrp= 'any ' ; namex = 'pbe '
      stden= .true.  ; dos   = .false. ; secvar = .false.
      vchk = .false. ; cdinf = .false. ; ctail  = .true.
      pot8 = .false. ; eig66(1) = .false. ; eig66(2) = .true.  
      form66 = .false. ; l_u2f= .false. ; l_f2u = .false. 
      l_bmt= .false. ; eonly  = .false.
      gauss= .false. ; tria  = .false. ; frcor  = .false.
      slice= .false. ; disp  = .false. ; swsp  = .false.
      lflip= .false. ; vacdos= .false. ; integ = .false.
      iplot= .false. ; score = .false. ; plpot = .false.
      pallst = .false. ; lwb = .false. ; starcoeff = .false.
      strho  = .true.  ; l_f = .false. ; l_geo(:) = .false.
      l_noco = l_ss ; l_J = .false. ; soc_opt(:) = .false. ; jspins = 1  
      lpr = 0 ; itmax = 9 ; maxiter = 99 ; imix = 7 ; alpha = 0.05 
      spinf = 2.0 ; igrd = 1 ; ndvgrd = 6 ; idsprs = 0 ; lepr = 0 
      kcrel = 0 ; kk = 0 ; nnne = 0 ; nwd = 1 ; nstars = 0 ; nstm = 0 
      isec1 = 99 ; nu = 5 ; layerd = 1 ; nlod = 2  ; iofile = 6  
      ndir = 0 ; layers = 0 ; nflip(:) = 1 ; izlay(:,:) = 0 
      lda_u%l = -1 ; relax(1:2,:) = 0 ; relax(3,:) = 1
      epsdisp = 0.00001 ; epsforce = 0.00001 ; xa = 2.0 ; thetad = 330.0
      e1s = 0.0 ; e2s = 0.0 ; e1_dos = 0.5 ; e2_dos = -0.5 ; tkb = 0.001
      sig_dos = 0.015 ; tworkf = 0.0 ; scale = 1.0 ; scpos = 1.0 
      zc = 0.0 ; chng = -1.0e-12 ; locx(:) = 0.0 ;  locy(:) = 0.0 

!+odim
      odd%mb = 0 ; odd%M = 0 ; odd%m_cyl = 0 ; odd%chi = 0 ; odd%rot = 0
      odd%k3 = 0 ; odd%n2d= 0 ; odd%nq2 = 0 ; odd%nn2d = 0 
      odd%nop = 0 ; odd%kimax2 = 0 ; odd%nat = 0
      odd%invs = .false. ; odd%zrfs = .false. ; odd%d1 = .false.
!-odim
! check for magnetism
      bmu(:) = 0.0
      DO n = 1, ntype
        IF (nz(n).EQ.24) bmu(n) = 1.0  ! Cr - Ni
        IF (nz(n).EQ.25) bmu(n) = 3.5
        IF (nz(n).EQ.26) bmu(n) = 2.2
        IF (nz(n).EQ.27) bmu(n) = 1.6
        IF (nz(n).EQ.28) bmu(n) = 1.1
        IF (nz(n).EQ.59) bmu(n) = 2.1  ! Pr - Tm
        IF (nz(n).EQ.60) bmu(n) = 3.1
        IF (nz(n).EQ.61) bmu(n) = 4.1
        IF (nz(n).EQ.62) bmu(n) = 5.1
        IF (nz(n).EQ.63) bmu(n) = 7.1
        IF (nz(n).EQ.64) bmu(n) = 7.1 
        IF (nz(n).EQ.65) bmu(n) = 6.1
        IF (nz(n).EQ.66) bmu(n) = 5.1
        IF (nz(n).EQ.67) bmu(n) = 4.1
        IF (nz(n).EQ.68) bmu(n) = 3.1
        IF (nz(n).EQ.69) bmu(n) = 2.1
      ENDDO
      IF ( ANY(bmu(:) > 0.0) ) jspins=2 

      delgau = tkb ; ntypd = ntype ; natd = nat
      DO i = 1, 10
        j = (i-1) * 8 + 1
        name(i) = title(j:j+7)
      ENDDO 
      IF (l_noco) jspins = 2
       
      a1(:) = amat(:,1) ; a2(:) = amat(:,2) ; a3(:) = amat(:,3) 

      CALL chkmt(
     >           nat,ntype,neq,film,pos,dvac,rmt,amat,nz,
     >           l_gga,noel,l_test,odd,jspins,
     <           nel,kmax,dtild,dvac1,
     <           nlo,llo,ncst,lmax,jri,rmt1,dx)

      IF ( ANY(nlo(:).NE.0) ) THEN
        ellow(1) = -1.8
      ELSE
        ellow(1) = -0.8  
      ENDIF
      IF (film) THEN
         elup(1) = 0.5
      ELSE
         elup(1) = 1.0
      ENDIF
      gmax = 3.0 * kmax ; gmaxxc = 2.5 * kmax ; rkm(1) = kmax
      lnonsph(:) = min( max( (lmax(:)-2),3 ), 8 ) 
      zelec(1) = nel
      IF (film) taual(3,:) = taual(3,:) * a3(3) / dtild

      IF (.not.film) THEN
         dvac1 = a3(3) ; dtild = dvac1
      ENDIF
      IF ( (abs(a1(3)).GT.eps).OR.(abs(a2(3)).GT.eps).OR.
     +     (abs(a3(1)).GT.eps).OR.(abs(a3(2)).GT.eps) ) THEN          
        latnam = 'any'
      ELSE
        IF ( (abs(a1(2)).LT.eps).AND.(abs(a2(1)).LT.eps) ) THEN
          IF (abs(a1(1)-a2(2)).LT.eps) THEN
            latnam = 'squ'
          ELSE
            latnam = 'p-r'
          ENDIF
        ELSE
          n1 = sqrt(a1(1)**2 + a1(2)**2); n2 = sqrt(a2(1)**2 + a2(2)**2)
          IF (abs(n1-n2).LT.eps) THEN
            gam = ( a1(1)*a2(1) + a1(2)*a2(2) ) / (n1 * n2)
            gam = 57.295779512*acos(gam)
            IF (abs(gam-60.).LT.eps) THEN
               latnam = 'hex'
               a1(2) = n1 * 0.5
               a1(1) = a1(2) * sqrt(3.0)
            ELSEIF (abs(gam-120.).LT.eps) THEN
               latnam = 'hx3'
               a1(1) = n1 * 0.5
               a1(2) = a1(1) * sqrt(3.0)
            ELSE
               latnam = 'c-r'
               gam = 0.5 * gam / 57.295779512
               a1(1) = n1 * cos(gam)
               a1(2) = n1 * sin(gam)
            ENDIF
            a2(1) = a1(1)
            a2(2) = - a1(2)
          ELSE
            latnam = 'obl'
          ENDIF
        ENDIF
      ENDIF

! rounding
      rmt1(:) = real(INT(rmt1(:)*100)/100.)
      dx(:) = real(INT(dx(:)*1000)/1000.)
      gmax = real(INT(gmax*10)/10.)
      rkm(:) = real(INT(rkm(:)*10)/10.)
      gmaxxc = real(INT(gmaxxc*10)/10.)
!
      CLOSE (6)
      iofile = 6
      OPEN (iofile,file='inp',form='formatted',status='new',
     &      iostat=iostat)
      IF (iostat /= 0) THEN
        STOP 'Cannot create new file "inp". Maybe it already exists?'
      ENDIF
      nu = 8 
      gw = 0 ; gw_neigd = 0
      CALL rw_inp(
     >            nu,ntypd,natd,nwdd,layerd,nlod,iofile,ch_rw,
     >  dvac1,dtild,scale,scpos,gmax,rkm,delgau,zc,tkb,alpha,spinf,
     X  e1s,e2s,isec1,ndir,jspins,lpr,nwd,lepr,kcrel,itmax,imix,layers,
     X  l_u2f,l_f2u,l_bmt,kk,nnne,maxiter,latnam,noel,namex,namgrp,
     X  relcor,
     X  strho,film,dos,secvar,invs,zrfs,invs2,vchk,cdinf,pot8,form66,
     X  l_f,eonly,eig66,gauss,tria,frcor,slice,ctail,disp,swsp,lflip,
     X  vacdos,integ,iplot,score,plpot,pallst,a1,a2,a3,rmt1,dx,ellow,
     X  elup,zelec,bmu,ncst,jri,lmax,lnonsph,nflip,izlay,
     X  name,igrd,ndvgrd,idsprs,lwb,chng,gmaxxc,l_soc,soc_opt,theta,phi,
     X  taual,
     X  ntype,neq,nz,xa,relax,thetad,epsdisp,epsforce,nlo,llo,tworkf,
     X  nstars,nstm,starcoeff,locx,locy,l_noco,l_J,l_geo,e1_dos,e2_dos,
     X  sig_dos,lda_u,gw,gw_neigd,odd)

      IF (film) THEN
        area = omtil / dvac1
        nkpt = nint((3600/area)/nop2)
      ELSE
        nkpt = nint((216000/omtil)/nop) 
      ENDIF
      WRITE (iofile,'(a5,i5)') 'nkpt=',nkpt
      CLOSE (iofile)

      END SUBROUTINE set_inp
      END MODULE m_setinp
