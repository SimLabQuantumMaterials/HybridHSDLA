      MODULE m_inped
c     *******************************************************
c     read in input parameters
c     modified to include also empty spheres (z=1.e-10)
c     r.p. aug. 1990
c     now includes readin of k-points          * shz Apr.96 *
c     modified to include all exchange correlation potentials
c     and relativistic correction to vxc
c     r.pentcheva dec. 1995
c
cta+
c     igrd=0: no gradient correction.
c     igrd=1: pw-91. icorr=6.

c     ndvgrd: number of points used to calculate the numerical
cc           derivatives. 6 is biggest and presumably the best.
c     ntimec: if ntime ge ntimec, iteration stops with code=2.
c     distc : distance of convergence in charge criterion.
c     tendfc: read-in in mhtree.
cc            if tendf (total energy diff. in mili-hartree from former
cc            tenrg) becomes less than tendfc, ntime=ntime+1.
c     chng  : charge-negative.
cc             if(ro.lt.chng) ineg=1 and stop.
c     idsprs: if eq 1, gga contributions are forced to vanish.
cta-
c     *******************************************************
c
      CONTAINS
      SUBROUTINE inped( 
     > ntypd,natd,nwdd,layerd,nmzxyd,nmzd,jmtd,lmaxd,
     > nlod,llod,ivers,namat,noinpfile,
     < strho,dos,isec1,ndir,secvar,icorr,krla,namgrp,
     < cdinf,ncst,pot8,vchk,latnam,name,
     < igrd,ndvgrd,idsprs,isprsv,idsprs0,idsprsl,
     < idsprsi,idsprsv,chng,sprsv,lwb,gmaxxc,
     < gauss,tria,delgau,tkb,ellow,elup,nwd,zelec,
     < lpr,form66,l_f,eonly,eig66,lnonsph,rkmax,
     < frcor,slice,ctail,disp,kcrel,lepr,
     < itmax,maxiter,imix,alpha,spinf,l_u2f,l_f2u,l_bmt,
     < film,invs,invs2,zrfs,omtil,vol,volint,volmts,
     < area,z1,delz,lmax,l_soc,soc_opt,theta,phi,l_geo,
     < jspins,nvac,nmz,nmzxy,dvac,scale,ntype,neq,jri,gmax,
     < swsp,lflip,bmu,nflip,vacdos,layers,integ,
     < izlay,iplot,score,plpot,kk,e1s,e2s,nnne,pallst,
     < amat,bmat,dx,pos,taual,rmt,rmsh,zatom,lapw_l,ulo_der,
     < nlo,llo,nlol,lo1l,l_dulo,tworkf,nstars,nstm,starcoeff,locx,
     < locy,l_noco,l_ss,l_mperp,l_constr,mix_b,alph,beta,
     < b_con,l_relax,qss,e1_dos,e2_dos,sig_dos,lda_u,n_u,
     < gw,gw_neigd,odd,l_J,thetaJ,l_magn,nmagn,M,mtypes,
     < magtype,nmagtype,nsh,l_disp,sso_opt)
      USE m_rwinp
      USE m_chkmt
      USE m_inpnoco
      USE m_cotra,     ONLY : cotra0
      USE m_constants, ONLY : pimach
      USE m_types,     ONLY : t_utype
      USE m_od_types,  ONLY : od_dim
      USE m_inv3
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER,    INTENT (IN)  :: ntypd,natd,nwdd,layerd,nmzd,nmzxyd
      INTEGER,    INTENT (IN)  :: jmtd,lmaxd,nlod,llod,noinpfile
      INTEGER,    INTENT (OUT) :: gw,gw_neigd
      LOGICAL,    INTENT (OUT) :: cdinf,pot8,vchk,strho,dos,secvar
      LOGICAL,    INTENT (OUT) :: gauss,tria,frcor,slice,ctail,disp
      LOGICAL,    INTENT (OUT) :: swsp,lflip,vacdos,integ
      LOGICAL,    INTENT (OUT) :: film,invs,invs2,zrfs,l_soc
      LOGICAL,    INTENT (OUT) :: l_noco,l_ss,l_mperp,l_constr
      LOGICAL,    INTENT (OUT) :: iplot,score,plpot,pallst,l_f
      LOGICAL,    INTENT (OUT) :: eonly,form66,l_u2f,l_f2u,l_bmt
      INTEGER,    INTENT (OUT) :: icorr,isec1,ndir,krla,kcrel,lpr,lepr
      INTEGER,    INTENT (OUT) :: itmax,maxiter,imix,layers,kk,nnne
      INTEGER,    INTENT (OUT) :: nwd,jspins,nvac,nmz,nmzxy,ntype
      REAL,       INTENT (OUT) :: delgau,tkb,alpha,spinf,e1s,e2s,gmax
      REAL,       INTENT (OUT) :: z1,delz,omtil,vol,volint,area,dvac
      REAL,       INTENT (OUT) :: scale
      REAL,       INTENT (OUT) :: theta,phi,mix_b,e1_dos,e2_dos,sig_dos
      CHARACTER(len=3), INTENT (OUT) :: latnam
      CHARACTER(len=4), INTENT (OUT) :: namgrp
      CHARACTER(len=9), INTENT (IN)  :: ivers
c--- J<
      INTEGER,    INTENT(OUT)   ::  nmagn,nsh,mtypes
      LOGICAL,    INTENT(OUT)   ::  l_J,l_disp
      REAL,       INTENT(OUT)   ::  thetaJ
c--- J>

c+dw new arguments for STM
      REAL,       INTENT (OUT) :: tworkf
      INTEGER,    INTENT (OUT) :: nstars,nstm
      LOGICAL,    INTENT (OUT) :: starcoeff
      REAL,       INTENT (OUT) :: locx(2),locy(2)
c-dw
c---> pk non-collinear
c---> new arrays for non-collinear calculations
      REAL,       INTENT (OUT) :: alph(ntypd),b_con(2,ntypd)
      REAL,       INTENT (OUT) :: beta(ntypd)
      REAL,       INTENT (OUT) :: qss(3)
c---> pk non-collinear
c--- J<
      LOGICAL,    INTENT(INOUT)   ::  l_magn(ntypd)
      INTEGER,    INTENT(OUT)     ::  magtype(ntypd),nmagtype(ntypd)
      REAL,       INTENT (INOUT)  :: M(ntypd)
c--- J>

c+lda+u
      TYPE (t_utype), INTENT (OUT) :: lda_u(ntypd)
      INTEGER,        INTENT (OUT) :: n_u
c-lda+u
C     ..
c-odim
      TYPE (od_dim), INTENT (INOUT) :: odd
c+odim
C     .. Array Arguments ..
      INTEGER,    INTENT (OUT) :: ncst(ntypd),nflip(ntypd)
      INTEGER,    INTENT (OUT) :: izlay(layerd,2),lnonsph(ntypd)
      INTEGER,    INTENT (OUT) :: jri(ntypd),neq(ntypd),lmax(ntypd)
      INTEGER,    INTENT (OUT) :: nlo(ntypd),llo(nlod,ntypd)
      INTEGER,    INTENT (OUT) :: lapw_l(ntypd),ulo_der(nlod,ntypd)
      INTEGER,    INTENT (OUT) :: lo1l(0:llod,ntypd),nlol(0:llod,ntypd)
      REAL,       INTENT (OUT) :: amat(3,3),bmat(3,3),dx(ntypd)
      REAL,       INTENT (OUT) :: rmsh(jmtd,ntypd),rmt(ntypd)
      REAL,       INTENT (OUT) :: taual(3,natd),pos(3,natd),zatom(ntypd)
      REAL,       INTENT (OUT) :: bmu(ntypd),zelec(nwdd),volmts(ntypd)
      REAL,       INTENT (OUT) :: rkmax(nwdd),ellow(nwdd),elup(nwdd)
      CHARACTER*2,INTENT (IN)  :: namat(0:103)
      CHARACTER*8,INTENT (OUT) :: name(10)
      LOGICAL,    INTENT (OUT) :: l_geo(ntypd),l_relax(ntypd)
      LOGICAL,    INTENT (OUT) :: l_dulo(nlod,ntypd)
      LOGICAL,    INTENT (OUT) :: soc_opt(ntypd+2),eig66(2),sso_opt(2)  
C     ..
C     .. Arguments (for GGA use )
      INTEGER igrd,ndvgrd,idsprs,isprsv,
     +        iggachk,idsprs0,idsprsl,idsprsi,idsprsv
      REAL    chng,sprsv,gmaxxc
      LOGICAL lwb
C     ..
C     .. Local Scalars ..
      REAL dr,dtild,r,s,scpos,sss,sum,kmax1,dtild1,dvac1
      REAL epsdisp,epsforce,thetad,zc,xa,rscale,tpi,fpi,zp
      INTEGER i,iz,j,n,n1,na,nw,nat,nel1,ntst,nn,ios
      LOGICAL l_gga,l_test,l_vca
      CHARACTER(len=2)  :: str_up,str_do
      CHARACTER(len=4)  :: namex
      CHARACTER(len=12) :: relcor
C     ..
C     .. Local Arrays ..
      CHARACTER*3 noel(ntypd)
      CHARACTER*8 llr(0:1)
      CHARACTER*11 pmod(0:1)
      INTEGER relax(3,ntypd),nz(ntypd),jri1(ntypd)
      INTEGER nlo1(ntypd),llo1(nlod,ntypd),ncst1(ntypd),lmax1(ntypd)
      REAL    rmt1(ntypd),dx1(ntypd)
      REAL    rkm(nwdd),a1(3),a2(3),a3(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL setlomap
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp
C     ..
C     .. Data statements ..
      DATA llr(0)/'absolute'/,llr(1)/'floating'/
      DATA pmod(0)/'not printed'/,pmod(1)/'printed    '/
c
      tpi = 2 * pimach()
      fpi = 2 * tpi
C     ..
      DO i = 1,3
         a1(i) = 0
         a2(i) = 0
         a3(i) = 0
      ENDDO

      na = 0
      rscale = 1.0

      WRITE (6,FMT=8000) ivers
      WRITE (16,FMT=8000) ivers
 8000 FORMAT ('1   f l a p w  version ',a9,/,/)
      CALL rw_inp(
     >            16,ntypd,natd,nwdd,layerd,nlod,5,'R',
     < dvac,dtild,scale,scpos,gmax,rkm,delgau,zc,tkb,alpha,spinf,
     < e1s,e2s,isec1,ndir,jspins,lpr,nwd,lepr,kcrel,itmax,imix,layers,
     < l_u2f,l_f2u,l_bmt,kk,nnne,maxiter,latnam,noel,namex,namgrp,
     < relcor,
     < strho,film,dos,secvar,invs,zrfs,invs2,vchk,cdinf,pot8,form66,
     < l_f,eonly,eig66,gauss,tria,frcor,slice,ctail,disp,swsp,lflip,
     < vacdos,integ,iplot,score,plpot,pallst,a1,a2,a3,rmt,dx,ellow,
     < elup,zelec,bmu,ncst,jri,lmax,lnonsph,nflip,izlay,
     < name,igrd,ndvgrd,idsprs,lwb,chng,gmaxxc,l_soc,soc_opt,theta,phi,
     < taual,
     < ntype,neq,nz,xa,relax,thetad,epsdisp,epsforce,nlo,llo,tworkf,
     < nstars,nstm,starcoeff,locx,locy,l_noco,l_J,l_geo,e1_dos,e2_dos,
     < sig_dos,lda_u,gw,gw_neigd,odd)

c---> pk non-collinear
c---> read the angle information from nocoinf
      qss(1) = 0.0
      qss(2) = 0.0
      qss(3) = 0.0
      IF (l_noco) THEN 
         CALL inpnoco(
     >    natd,noinpfile,itmax,isec1,secvar,gauss,tria,
     >    ctail,score,ntype,neq,taual,nstm,starcoeff,l_f,l_J,l_soc,
     <    l_ss,l_mperp,l_constr,mix_b,alph,beta,b_con,
     <    l_relax,qss,thetaJ,l_magn,nmagn,M,mtypes,
     <    magtype,nmagtype,nsh,l_disp,sso_opt)
      ELSE
         l_ss = .false.
         l_mperp = .false.
         l_constr = .false.
         mix_b = 0.0
         thetaJ = 0.0
         nmagn=1
         DO n = 1,ntype
           l_relax(n) = .false.
           l_magn(n) = .false.
           alph(n) = 0.0
           beta(n) = 0.0
           b_con(1,n) = 0.0
           b_con(2,n) = 0.0
         ENDDO
      ENDIF
c---> pk non-collinear

      WRITE (6,FMT=8010) name
      WRITE (16,FMT=8010) name
 8010 FORMAT (/,/,4x,10a8,/,/)
 8020 FORMAT (a3,1x,a4,6x,l1,6x,l1,7x,l1,8x,i1)
c--->    the menu for namgrp can be found in subroutine spgset
      WRITE (6,FMT=8030) latnam,namgrp,invs,zrfs,invs2,jspins
      WRITE (16,FMT=8030) latnam,namgrp,invs,zrfs,invs2,jspins
 8030 FORMAT (' lattice=',a3,/,' name of space group=',a4,/,
     +       ' inversion symmetry=   ',l1,/,' z-reflection symmetry=',
     +       l1,/,' vacuum-inversion symm=',l1,/,' jspins=',i1)

      IF (film.AND.(invs.OR.zrfs)) THEN
        IF ( (invs.AND.zrfs).NEQV.invs2 ) THEN
          WRITE (6,*) 'Settings of inversion and z-reflection symmetry='
          WRITE (6,*) 'are inconsistent with vacuum-inversion symmetry!'
          STOP 'inped: invs, zrfs and invs2 do not match!'
        ENDIF
      ENDIF

      IF ((a1(1).EQ.0.).AND.(a1(2).EQ.0.).AND.(a1(3).EQ.0.)) THEN
        WRITE (6,'(a4,3f10.5,a8,a4)') 'a1 =',a1(:),' latnam=',latnam
        STOP 'inped: latnam'
      ENDIF
      a3(3) = dtild
      IF (scale.EQ.0.) scale = 1.
      dvac = scale*dvac
      dtild = scale*dtild
c+odim
      IF (.NOT.odd%d1) THEN 
         IF ((dtild-dvac.LT.0.0).AND.film) THEN
             write(6,'(2(a7,f10.5))') 'dtild:',dtild,' dvac:',dvac
             STOP 'inped: dtild < dvac'
         ENDIF
      ELSE
         IF (dvac.GE.sqrt(a1(1)**2 + a1(2)**2).OR.
     &       dvac.GE.sqrt(a2(1)**2 + a2(2)**2)) THEN
            STOP 'one-dim:inped: dvac >= amat(1,1) or amat(2,2)'
         END IF
      ENDIF
c-odim
      nvac = 2
      IF (zrfs .OR. invs) nvac = 1
      IF (odd%d1) nvac = 1
      z1 = dvac/2
      nmz = nmzd
      delz = 25.0/nmz
      IF (odd%d1) delz = 20.0/nmz
      IF (nmz.GT.nmzd) STOP 'nmzd'
      nmzxy = nmzxyd
      IF (nmzxy.GT.nmzxyd) STOP 'nmzxyd'
      DO 20 j = 1,3
         a1(j) = scale*a1(j)
         a2(j) = scale*a2(j)
         a3(j) = scale*a3(j)
   20 CONTINUE
      WRITE (6,FMT=8050) scale
      WRITE (16,FMT=8050) scale
 8050 FORMAT (' unit cell scaled by    ',f10.6)
      WRITE (6,FMT=8060) z1
      WRITE (16,FMT=8060) z1
 8060 FORMAT (' the vacuum begins at z=',f10.6)
      WRITE (6,FMT=8070) dtild/2.
      WRITE (16,FMT=8070) dtild/2.
 8070 FORMAT (' dtilda/2=              ',f10.6)
c     set up bravais matrices of real and reciprocal lattices
      DO 30 i = 1,3
         amat(i,1) = a1(i)
         amat(i,2) = a2(i)
         amat(i,3) = a3(i)
   30 CONTINUE
      CALL inv3(amat,bmat,omtil)

      omtil = abs(omtil)

      IF (film .AND. .NOT.odd%d1) THEN
         vol = omtil/amat(3,3)*dvac
         area = omtil/amat(3,3)
c-odim
      ELSEIF (odd%d1) THEN
         area = tpi*amat(3,3)
         vol = pimach()*(dvac**2)*amat(3,3)/4.
c+odim
      ELSE
         vol = omtil
         area = amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)
         IF (area.lt.1.0e-7) THEN
           IF (latnam.EQ.'any') THEN
            area = 1.
           ELSE
            STOP 'area = 0'
           ENDIF
         ENDIF
      ENDIF

      DO 50 j = 1,3
         DO 40 i = 1,3
            bmat(i,j) = tpi*bmat(i,j)
   40    CONTINUE
   50 CONTINUE
      WRITE (6,FMT=8080)
 8080 FORMAT (/,/,1x,'bravais matrices of real and reciprocal lattices',
     +       /)
      DO 60 i = 1,3
         WRITE (6,FMT=8090) (amat(i,j),j=1,3), (bmat(i,j),j=1,3)
   60 CONTINUE
 8090 FORMAT (3x,3f10.6,3x,3f10.6)
      WRITE (6,FMT=8100) omtil,vol,area
 8100 FORMAT (/,4x,'the volume of the unit cell omega-tilda=',f12.6,/,
     +       10x,'the volume of the unit cell omega=',f12.6,/,2x,
     +       'the area of the two-dimensional unit cell=',f12.6)
      WRITE (6,FMT=8120) namex,relcor
 8120 FORMAT (1x,'exchange-correlation: ',a4,2x,a12,1x,'correction')
      icorr = -99

c     l91: lsd(igrd=0) with dsprs=1.d-19 in pw91.
      IF (namex.EQ.'l91')  icorr =-1
      IF (namex.EQ.'x-a')  icorr = 0
      IF (namex.EQ.'wign') icorr = 1
      IF (namex.EQ.'mjw')  icorr = 2
      IF (namex.EQ.'hl')   icorr = 3
      IF (namex.EQ.'bh')   icorr = 3
      IF (namex.EQ.'vwn')  icorr = 4
      IF (namex.EQ.'pz')   icorr = 5
      IF (namex.EQ.'pw91') icorr = 6
c     pbe: easy_pbe [Phys.Rev.Lett. 77, 3865 (1996)]
c     rpbe: rev_pbe [Phys.Rev.Lett. 80, 890 (1998)]
c     Rpbe: Rev_pbe [Phys.Rev.B 59, 7413 (1999)]
      IF (namex.eq.'pbe')  icorr = 7
      IF (namex.eq.'rpbe') icorr = 8
      IF (namex.eq.'Rpbe') icorr = 9
      IF (namex.eq.'wc')   icorr =10
!     wc: Wu & Cohen, [Phys.Rev.B 73, 235116 (2006)]
      IF (namex.eq.'PBEs') icorr =11
!     PBEs: PBE for solids ( arXiv:0711.0156v2 )

      IF (icorr == -99) THEN
        WRITE(6,*) 'Name of XC-potential not recognized. Use one of:'
        WRITE(6,*) 'x-a,wign,mjw,hl,bh,vwn,pz,l91,pw91,pbe,rpbe,Rpbe'
        STOP 'Wrong name of XC-potential!'
      ENDIF
      igrd = 0
      IF (icorr.GE.6) igrd = 1
      krla = 0
      IF (relcor.EQ.'relativistic') THEN 
         krla = 1    
         IF (igrd.EQ.1) THEN
           WRITE(6,'(18a,a4)') 'Use XC-potential: ',namex
           WRITE(6,*) 'only without relativistic corrections !'
           STOP 'relativistic corrections + GGA not implemented'
         ENDIF
      ENDIF

      IF (icorr.eq.0) WRITE(6,*) 'WARNING: using X-alpha for XC!'
      IF (icorr.eq.1) WRITE(6,*) 'INFO   : using Wigner  for XC!'
      IF ((icorr.eq.2).and.(namex.NE.'mjw')) 
     +                WRITE(6,*) 'CAUTION: using MJW(BH) for XC!'

c+guta
      IF ((icorr.EQ.-1).OR.(icorr.GE.6)) THEN

         isprsv=0
         sprsv=0.0
         ndvgrd = max(ndvgrd,3)
         IF ((igrd.NE.0).AND.(igrd.NE.1)) THEN 
           WRITE (6,*) 'selecting l91 or pw91 as XC-Potental you should'
           WRITE (6,*) ' have 2 lines like this in your inp-file:'
           WRITE (6,*) 'igrd=1,lwb=F,ndvgrd=4,idsprs=0,chng= 1.000E-16'
           WRITE (6,*)
     +          'iggachk=1,idsprs0=1,idsprsl=1,idsprsi=1,idsprsv=1'
           STOP 'inped: igrd =/= 0 or 1'
         ENDIF

c        iggachk: removed; triggered via idsprs (see below)
c                 idsprs-0(mt,l=0),-l(nmt),-i(interstitial),-v(vacuum)
c                 enable to make gga partially enactive if corresponding
c                 idsprs set to be zero.

         IF (idsprs.EQ.0) THEN
           idsprs0=0 ; idsprsl=0 ; idsprsi=0 ; idsprsv=0
         ELSEIF (idsprs.EQ.1) THEN
           idsprs0=1 ; idsprsl=0 ; idsprsi=0 ; idsprsv=0
         ELSEIF (idsprs.EQ.2) THEN
           idsprs0=1 ; idsprsl=1 ; idsprsi=0 ; idsprsv=0
         ELSEIF (idsprs.EQ.3) THEN
           idsprs0=0 ; idsprsl=0 ; idsprsi=1 ; idsprsv=0
         ELSEIF (idsprs.EQ.4) THEN
           idsprs0=0 ; idsprsl=0 ; idsprsi=0 ; idsprsv=1
         ELSEIF (idsprs.EQ.5) THEN
           idsprs0=0 ; idsprsl=0 ; idsprsi=1 ; idsprsv=1
         ELSEIF (idsprs.EQ.6) THEN
           idsprs0=1 ; idsprsl=1 ; idsprsi=1 ; idsprsv=1
         ELSE
           idsprs0=0 ; idsprsl=0 ; idsprsi=0 ; idsprsv=0
         ENDIF

         WRITE (16,FMT=8122) igrd,lwb,ndvgrd,idsprs,chng
         IF (idsprs.EQ.1) WRITE (16,8125) 'MT (l=0) '
         IF (idsprs.EQ.2) WRITE (16,8125) 'MT (l>0) '
         IF (idsprs.EQ.3) WRITE (16,8125) 'interst. '
         IF (idsprs.EQ.4) WRITE (16,8125) 'vacuum   '
         IF (idsprs.EQ.5) WRITE (16,8125) 'int.+vac.'
         IF (idsprs.EQ.6) WRITE (16,8125) 'all parts'
         WRITE (16,'(/)')
         WRITE (6,FMT=8122) igrd,lwb,ndvgrd,idsprs,chng
         IF (idsprs.EQ.1) WRITE (6,8125) 'MT (l=0) '
         IF (idsprs.EQ.2) WRITE (6,8125) 'MT (l>0) '
         IF (idsprs.EQ.3) WRITE (6,8125) 'interst. '
         IF (idsprs.EQ.4) WRITE (6,8125) 'vacuum   '
         IF (idsprs.EQ.5) WRITE (6,8125) 'int.+vac.'
         IF (idsprs.EQ.6) WRITE (6,8125) 'all parts'
         WRITE (6,'(/)')
 8122    FORMAT ('igrd=',i1,',lwb=',l1,',ndvgrd=',i1,',idsprs=',i1,
     +           ',chng=',d10.3)
 8125    FORMAT ('GGA switched off in the ',a9) 

      ENDIF
c-guta
c     specification of atoms
      IF (ntype.GT.ntypd) THEN
         WRITE (6,FMT='(a)') 'ntype > ntypd !!!'
         WRITE (16,FMT='(a)') 'ntype > ntypd !!!'
         STOP 'ntypd'
      END IF
      volint = vol

      DO  n = 1,ntype
       IF (TRIM(ADJUSTL(noel(n))).NE.TRIM(ADJUSTL(namat(nz(n))))) THEN
         CALL trans(namat(n),str_up,str_do)
         IF ( (TRIM(ADJUSTL(noel(n))).NE.TRIM(ADJUSTL(str_up))) .OR.
     &        (TRIM(ADJUSTL(noel(n))).NE.TRIM(ADJUSTL(str_do))) ) THEN
           WRITE(16,*) 'Element ',noel(n),' does not match Z = ',nz(n)
           WRITE( 6,*) 'Element ',noel(n),' does not match Z = ',nz(n)
           STOP 'Element name and nuclear number do not match!'
         ENDIF
       ENDIF
       WRITE (6,8140) noel(n),nz(n),ncst(n),lmax(n),jri(n),rmt(n),dx(n)
       WRITE (16,8140) noel(n),nz(n),ncst(n),lmax(n),jri(n),rmt(n),dx(n)
 8140  FORMAT (a3,i3,3i5,2f10.6)
         IF (jri(n).GT.jmtd) STOP 'jmtd'
         zatom(n) = nz(n)
         IF (nz(n).EQ.0) zatom(n) = 1.e-10
!
! check for virtual crystal approximation
!
         l_vca = .false.
         INQUIRE (file="vca.in", exist=l_vca)
         IF (l_vca) THEN
           OPEN (17,file='vca.in',form='formatted')
           DO nn = 1, n
            READ (17,*,IOSTAT=ios) ntst,zp
            IF (ios /= 0) EXIT
            IF (ntst == n) THEN
              zatom(n) = zatom(n) + zp
            ENDIF
           ENDDO 
           CLOSE (17)
         ENDIF
!
         r = rmt(n)*exp(dx(n)* (1-jri(n)))
         dr = exp(dx(n))
         DO i = 1,jri(n)
            rmsh(i,n) = r
            r = r*dr
         ENDDO  
         volmts(n) = fpi/3.*rmt(n)**3
         volint = volint - volmts(n)*neq(n)

         DO n1 = 1,neq(n)
            na = na + 1
            IF (na.GT.natd) STOP 'natd too small'
!
!--->    the in-plane coordinates refer to the lattice vectors a1 and a2,
!--->    i.e. they are given in internal units scaled by 'scpos'
!
            WRITE (6,FMT=8170) (taual(i,na),i=1,3),scpos
            WRITE (16,FMT=8170) (taual(i,na),i=1,3),scpos
 8170       FORMAT (4f10.6)
!
!--->   for films, the z-coordinates are given in absolute values:
!
            IF (film) taual(3,na) = scale*taual(3,na)/a3(3)
!
! Transform intern coordinates to cartesian:
!
            CALL cotra0(taual(1,na),pos(1,na),amat)

         ENDDO  ! l.o. equivalent atoms (n1)
      ENDDO     ! loop over atom-types (n)
c
c  check muffin tin radii
c
      l_gga = .false.
      IF (icorr.GE.6) l_gga = .true.
      l_test = .true.                  ! only checking, dont use new parameters
      CALL chkmt(
     >           natd,ntype,neq,film,pos,dvac,rmt,amat,nz,
     >           l_gga,noel,l_test,odd,jspins,
     <           nel1,kmax1,dtild1,dvac1,
     <           nlo1,llo1,ncst1,lmax1,jri1,rmt1,dx1)

      WRITE (6,FMT=8180) volint
 8180 FORMAT (13x,' volume of interstitial region=',f12.6)
      nat = na
c--->    evaluate cartesian coordinates of positions
      WRITE (6,FMT=8190) ntype,nat
 8190 FORMAT (/,/,' number of atom types=',i3,/,
     +       ' total number of atoms=',i2,/,/,t3,'no.',t10,'type',
     +       t21,'int.-coord.',t49,'cart.coord.',t76,'rmt',t84,
     +       'jri',t92,'dx',t98,'lmax',/)
      na = 0
      DO 140 n = 1,ntype
         DO 130 n1 = 1,neq(n)
            na = na + 1
            iz = zatom(n) + 0.001
            WRITE (6,FMT=8200) na,namat(iz),n,
     +        (taual(i,na),i=1,3), (pos(i,na),i=1,3),rmt(n),jri(n),
     +        dx(n),lmax(n)
 8200       FORMAT (1x,i3,4x,a2,t12,i3,2x,3f6.2,3x,3f10.6,3x,
     +             f10.6,i6,3x,f6.4,3x,i2)
  130    CONTINUE
  140 CONTINUE
c
c--->    input various parameters for eigenvalue parts: see intro. to
c--->    eigen for the various values:
c--->    lpr=0,form66=f,l_f=f,eonly=f   is an example.
      lpr = max(lpr,0)
      lpr = min(lpr,1)
c
c--->    lnonsph(n): max. l for H -setup in each atom type;
c
      DO n = 1,ntype
c+APW
         lapw_l(n) = (lnonsph(n) - mod(lnonsph(n),10) )/10
         lnonsph(n) = mod(lnonsph(n),10)
c-APW
         IF (lnonsph(n).EQ.0) lnonsph(n) = lmax(n)
         lnonsph(n) = min(lnonsph(n),lmax(n))
      ENDDO
c--->    nwd = number of energy windows; lepr = 0 (1) for energy
c--->    parameters given on absolute (floating) scale
      WRITE (16,FMT=*) 'nwd=',nwd,'lepr=',lepr
      IF (nwd.GT.nwdd) STOP 'nwdd'
      lepr = max(lepr,0)
      lepr = min(lepr,1)
      WRITE (6,FMT=8320) pmod(lpr),form66,l_f,eonly,nwd,llr(lepr)
      WRITE (16,FMT=8320) pmod(lpr),form66,l_f,eonly,nwd,llr(lepr)
      WRITE (6,FMT=8330) ntype, (lnonsph(n),n=1,ntype)
      WRITE (16,FMT=8330) ntype, (lnonsph(n),n=1,ntype)
 8300 FORMAT (4x,i1,8x,l1,5x,l1,7x,l1,7x,l1)
 8320 FORMAT (1x,/,/,/,' input of parameters for eigenvalues:',/,t5,
     +       'eigenvectors are ',a11,/,t5,
     +       'formatted eigenvector file = ',l1,/,t5,
     +       'calculate Pulay-forces = ',l1,/,t5,'eigenvalues ',
     +       'only = ',l1,/,t5,'number of energy windows =',i2,/,t5,
     +       'energy parameter mode: ',a8,/,/)
 8330 FORMAT (t5,'max. l value in wavefunctions for atom type(s) 1 to',
     +       i3,':',16i3,/, (t59,16i3,/))
c
c--->    input information  for each window
c
      IF (lepr.eq.1) THEN
         WRITE ( 6,'(//,''  Floating energy parameters: relative'',
     +                  '' window(s):'')')
         WRITE (16,'(//,''  Floating energy parameters: relative'',
     +                  '' window(s):'')')
      ENDIF
      DO nw = 1,nwd
c--->    energy window

c--->    for floating energy parameters, the window will be given relative
c--->    to the highest/lowest energy parameters. a sanity check is made here
         IF (lepr.eq.1) THEN
            IF (nw.lt.nwd) THEN  ! semicore should have narrow range
               ellow(nw) = min( ellow(nw) , -0.2 )
               elup(nw)  = min( max( elup(nw), 0.04 ) , 0.15 )
            ELSE
               ellow(nw) = min( ellow(nw) , -0.2 )
               elup(nw)  = max( elup(nw)  ,  0.15 )
            ENDIF
         ENDIF
c
         WRITE (6,FMT=8350) nw,nwd,ellow(nw),elup(nw),zelec(nw)
         WRITE (16,FMT=8350) nw,nwd,ellow(nw),elup(nw),zelec(nw)
 8350    FORMAT (1x,/,/,' energy window',i3,' of',i3,' from',f8.3,' to',
     +          f8.3,' hartrees; nr. of electrons=',f6.1)
c--->    input of wavefunction cutoffs: input is a scaled quantity
c--->    related to the absolute value by rscale (e.g. a muffin-tin
c--->    radius)
         rkmax(nw) = rkm(nw)/rscale
         WRITE (6,FMT=8290) rkm(nw),rkmax(nw)
         WRITE (16,FMT=8290) rkm(nw),rkmax(nw)
 8290    FORMAT (1x,/,' wavefunction cutoff (relative) =',f10.5,/,
     +          '                     (absolute) =',f10.5)
      ENDDO
c
      IF ((tria) .AND. (gauss)) THEN
         WRITE (6,FMT='(a)') 'choose: either gaussian or triangular!'
         WRITE (16,FMT='(a)') 'choose: either gaussian or triangular!'
         STOP 'integration method'
      END IF
      WRITE (6,FMT=8230) gauss,delgau
      WRITE (6,FMT=8240) zc,tkb
 8230 FORMAT (/,10x,'gauss-integration is used  =',3x,l1,/,10x,
     +       'gaussian half width        =',f10.5)
 8240 FORMAT (/,10x,'number of valence electrons=',f10.5,/,10x,
     +       'temperature broadening     =',f10.5)
      WRITE (6,FMT=*) 'itmax=',itmax,' broy_sv=',maxiter,' imix=',imix
      WRITE (6,FMT=*) 'alpha=',alpha,' spinf=',spinf
      WRITE (16,FMT=*) 'itmax=',itmax,' broy_sv=',maxiter,' imix=',imix
      WRITE (16,FMT=*) 'alpha=',alpha,' spinf=',spinf

!      IF ((.NOT.invs).AND.(isec1.LE.itmax)) THEN
!         WRITE(6,*)'The Wu-diagonalization is not implemented in the'
!         WRITE(6,*)'complex version of the program.'
!         STOP 'Wu not implemented in complex version'
!      ENDIF
      IF ((.NOT.invs).AND.secvar) THEN
         WRITE(6,*)'The second variation is not implemented in the'
         WRITE(6,*)'complex version of the program.'
         STOP 'second variation not implemented in complex version'
      ENDIF

#ifdef CPP_INVERSION
      IF (.NOT.invs .AND. .NOT.odd%d1)
     &             STOP 'recompile without -D CPP_INVERSION'
#endif
#ifndef CPP_INVERSION
      IF (secvar) STOP 'Second variation only with  -D CPP_INVERSION'
!      IF (isec1.LE.itmax)
!     +            STOP 'Wu-diagonalization only with  -D CPP_INVERSION'
#endif
#ifndef CPP_SOC
      IF (l_soc .and. (.not. l_noco)) STOP 'recompile with -D CPP_SOC'
#else
#ifdef CPP_INVERSION
      IF (.NOT.l_soc) STOP 'recompile without -D CPP_SOC'
#endif
#endif
      IF ( (jspins.EQ.1).AND.(kcrel.EQ.1) )  THEN
        WRITE (6,*) 'WARNING : in a non-spinpolarized calculation the'
        WRITE (6,*) 'coupled-channel relativistic coreprogram (kcrel=1)'
        WRITE (6,*) 'makes no sense; **** setting kcrel = 0 ****'
        kcrel = 0
      ENDIF

      WRITE (6,'(a7,l1)') 'swsp = ',swsp
      WRITE (6,'(15f6.2)') (bmu(i),i=1,ntype)
C     ******** read new variables vacdos,layers
C     ******* and izlay(i=1,layers)   shz Jan.96
      IF (layers.gt.layerd) STOP 'too many layers'
      IF (slice) THEN
         cdinf = .false.
         WRITE (6,FMT=8390) kk,e1s,e2s
         WRITE (16,FMT=8390) kk,e1s,e2s
      END IF
 8390 FORMAT (' slice: k=',i3,' e1s=',f10.6,' e2s=',f10.6)
c
c Check the LO stuff:
c
      DO n=1,ntype
        IF (nlo(n).GE.1) THEN
#ifdef CPP_INVERSION
          IF (l_soc.AND.(neq(n).GT.1)) THEN
!           STOP 'for LO + SOC use complex version in this case!'
          ENDIF
#endif
          IF (secvar)         STOP 'LO + sevcar not implemented'
          IF (isec1.LT.itmax) STOP 'LO + Wu not implemented'
          IF (nlo(n).GT.nlod) THEN
            WRITE (6,*) 'nlo(n) =',nlo(n),' > nlod =',nlod 
            STOP 'inped: nlo(n).GT.nlod'
          ENDIF
          DO j=1,nlo(n)
#ifndef CPP_APW
            IF (llo(j,n).LT.0) THEN !STOP 'llo<0 ; compile with DCPP_APW!'
              WRITE(6,'(A)') 'Info: Not compiled with CPP_APW.'
              WRITE(6,'(A,I2,A,I2,A)') '      LO #',j,' at atom type',n,
     &          ' is an e-derivative.'
            ENDIF
#endif
            IF ( (llo(j,n).GT.llod).OR.(mod(-llod,10)-1).GT.llod ) THEN
              WRITE (6,*) 'llo(j,n) =',llo(j,n),' > llod =',llod
              STOP 'inped: llo(j,n).GT.llod'
            ENDIF
          ENDDO
          CALL setlomap(
     >                  ntypd,nlod,llod,
     >                  n,nlo,llo,
     <                  lo1l,nlol,l_dulo,ulo_der)
#ifdef CPP_APW
        WRITE (6,*) 'lapw_l(n) = ',lapw_l(n)
#endif
        ENDIF
      ENDDO
c 
c Check for LDA+U:
c
      n_u = 0
      DO  n = 1,ntype
         IF (lda_u(n)%l.GE.0)  THEN
            n_u = n_u + 1
            IF (nlo(n).GE.1) THEN
               DO j = 1, nlo(n)
                 IF ((abs(llo(j,n)).EQ.lda_u(n)%l) .AND.
     +               (.not.l_dulo(j,n)) ) WRITE (*,*)
     +               'LO and LDA+U for same l not implemented'
               ENDDO
            ENDIF
         ENDIF 
      ENDDO
      IF (n_u.GT.0) THEN
          IF (secvar)         STOP 'LDA+U and sevcar not implemented'
          IF (isec1.LT.itmax) STOP 'LDA+U and Wu not implemented'
          IF (l_mperp)        STOP 'LDA+U and l_mperp not implemented'
      ENDIF
c
C     check all the dos-related switches!
c
      IF (ndir.lt.0.and..not.dos) THEN
         WRITE(6,*) 'STOP DOS: the inbuild dos-program ndir<0
     +               can only be used if dos=true'
         STOP 'DOS'
      ENDIF

      IF (ndir.lt.0.and.dos) THEN
         IF (e1_dos-e2_dos.lt.1e-3) THEN
            WRITE(6,*) 'STOP DOS: no valid energy window for
     +                    internal dos-program'
            STOP 'DOS'
         ENDIF
         IF (sig_dos.lt.0) THEN
            WRITE(6,*) 'STOP DOS: no valid broadening (sig_dos) for
     +                      internal dos-program'
            STOP 'DOS'
         ENDIF
      ENDIF

      IF (vacdos) THEN
         IF (.not. dos) THEN
            WRITE(6,*) 'STOP DOS: only set vacdos=.true. if dos=.true.'
            STOP 'DOS'
         ENDIF
         IF (.not.starcoeff.and.(nstars.ne.1))THEN
            WRITE(6,*) 'STOP DOS: if stars=f set nstars=1'
            STOP 'DOS'
         ENDIF
         IF (layers.lt.1) THEN
             WRITE(6,*) 'STOP DOS: specify layers if vacdos=true'
            STOP 'DOS'
         ENDIF
         DO i=1,layers
            IF (izlay(i,1).lt.1) THEN
               WRITE(6,*) 'STOP DOS: all layers must be at z>0'
               STOP 'DOS'
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE inped
!--------------------------------------------------------------
      SUBROUTINE trans(
     >                 string,
     <                 str_up,str_do)

      IMPLICIT NONE
      CHARACTER(len=2), INTENT(IN)  :: string
      CHARACTER(len=2), INTENT(OUT) :: str_up,str_do

      INTEGER offs,i,n,j
      CHARACTER(len=2) :: str_in
      CHARACTER(len=1) :: st(2)

      str_up='  ' ; str_do='  ' ; st(:)=' '
      offs = iachar('A') - iachar('a')
      str_in = trim(adjustl(string))
      n = len_trim(str_in)
      st = (/(str_in(i:i),i=1,n)/)
      DO i=1,n
        IF (iachar(st(i)) > iachar('Z')) THEN ! lowercase
          str_up(i:i) = char( iachar(st(i)) + offs)
        ELSE
          str_up(i:i) = st(i)
        ENDIF
      ENDDO
      DO i=1,n
        str_do(i:i) = char( iachar(str_up(i:i)) - offs)
      ENDDO
      END SUBROUTINE trans

      END MODULE
