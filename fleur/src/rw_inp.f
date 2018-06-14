      MODULE m_rwinp
      CONTAINS
      SUBROUTINE rw_inp(
     >                  nu,ntypd,natd,nwdd,layerd,nlod,iofile,ch_rw,
     X  dvac,dtild,scale,scpos,gmax,rkm,delgau,zc,tkb,alpha,spinf,
     X  e1s,e2s,isec1,ndir,jspins,lpr,nwd,lepr,kcrel,itmax,imix,layers,
     X  l_u2f,l_f2u,l_bmt,kk,nnne,maxiter,latnam,noel,namex,namgrp,
     X  relcor,
     X  strho,film,dos,secvar,invs,zrfs,invs2,vchk,cdinf,pot8,form66,
     X  l_f,eonly,eig66,gauss,tria,frcor,slice,ctail,disp,swsp,lflip,
     X  vacdos,integ,iplot,score,plpot,pallst,a1,a2,a3,rmt,dx,ellow,
     X  elup,zelec,bmu,ncst,jri,lmax,lnonsph,nflip,izlay,
     X  name,igrd,ndvgrd,idsprs,lwb,chng,gmaxxc,l_soc,soc_opt,theta,phi,
     X  taual,
     X  ntype,neq,nz,xa,relax,thetad,epsdisp,epsforce,nlo,llo,tworkf,
     X  nstars,nstm,starcoeff,locx,locy,l_noco,l_J,l_geo,e1_dos,e2_dos,
     X  sig_dos,lda_u,gw,gw_neigd,odd)

*********************************************************************
* This subroutine reads or writes an inp - file on unit iofile      *
* for ch_rw = 'R' read, ch_rw = 'W' write.                          *
*                                                           Gustav  *
*********************************************************************
      USE m_calculator
      USE m_types, ONLY : t_utype
      USE m_od_types, ONLY : od_dim
      IMPLICIT NONE
C ..
C ..  Scalar Arguments ..
      CHARACTER*1,INTENT (IN) :: ch_rw
      INTEGER, INTENT (IN)    :: ntypd,natd,nwdd,layerd,nlod
      INTEGER, INTENT (IN)    :: iofile,nu
      INTEGER, INTENT (INOUT) :: ntype
      REAL,    INTENT (INOUT) :: xa,thetad,epsdisp,epsforce
      INTEGER gw,gw_neigd
      REAL dvac,dtild,scale,scpos,gmax,delgau,zc,tkb,alpha,spinf
      REAL e1s,e2s,chng,gmaxxc,theta,phi,tworkf,e1_dos,e2_dos,sig_dos
      INTEGER isec1,ndir,jspins,lpr,nwd,lepr,kcrel,itmax
      INTEGER imix,layers,kk,nnne,maxiter,nw,nstars,nstm
      INTEGER igrd,ndvgrd,idsprs
      CHARACTER*3 latnam
      CHARACTER*4 namex,namgrp
      CHARACTER*12 relcor
      LOGICAL strho,film,dos,secvar,invs,zrfs,invs2,vchk,cdinf,pot8
      LOGICAL form66,l_f,eonly,gauss,tria,frcor,slice,ctail,disp
      LOGICAL swsp,lflip,vacdos,integ,iplot,score,plpot,pallst,lwb
      LOGICAL l_u2f,l_f2u,l_bmt,l_soc,starcoeff,l_noco,l_J
C ..
C ..  Array Arguments ..
      INTEGER, INTENT (INOUT) :: nz(ntypd),neq(ntypd),relax(3,ntypd)
      INTEGER, INTENT (INOUT) :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (INOUT) :: taual(3,natd),rkm(nwdd)
      REAL a1(3),a2(3),a3(3),rmt(ntypd),dx(ntypd),locx(2),locy(2)
      REAL ellow(nwdd),elup(nwdd),zelec(nwdd),bmu(ntypd)
      INTEGER ncst(ntypd),jri(ntypd),lmax(ntypd)
      INTEGER lnonsph(ntypd),nflip(ntypd),izlay(layerd,2)
      CHARACTER*3 noel(ntypd)
      CHARACTER*8 name(10)
      LOGICAL l_geo(ntypd),eig66(2),soc_opt(ntypd+2)
c+lda+u
      TYPE (t_utype), INTENT (INOUT) :: lda_u(ntypd)
      REAL    u,j
      INTEGER l
      LOGICAL l_amf
      CHARACTER*3 ch_test
      NAMELIST /ldaU/ l,u,j,l_amf
c-lda+u
c+odim
      TYPE (od_dim), INTENT (INOUT) :: odd
      INTEGER MM,vM,m_cyl
      LOGICAL invs1,zrfs1
      INTEGER chi,rot
      LOGICAL d1,band
      NAMELIST /odim/ d1,MM,vM,m_cyl,chi,rot,invs1,zrfs1
c-odim
C ..
C ..  Local Variables
      INTEGER ieq,i,k,na,n,ilo
      REAL s3,ah,a,hs2,rest
      LOGICAL l_sym,ldum
      INTEGER :: ierr
C ..
C...  Local Arrays
      CHARACTER :: helpchar(ntypd)
      CHARACTER(len=  3) :: chntype
      CHARACTER(len= 40) :: chform
      CHARACTER(len=100) :: line
c
c---------------------------------------------------------------------
      IF (ch_rw.eq.'R') THEN
c---------------------------------------------------------------------

      OPEN (iofile,file='inp',form='formatted',status='old')

      WRITE (nu,*) '-------- dump of inp-file ------------'
c
      !<-- Added possibility to define variables here

      DO
         READ (UNIT = iofile,FMT = 7182,END=77,ERR=77) ch_test
         BACKSPACE(iofile)
         IF (ch_test   /="def") EXIT
         READ(unit = iofile,FMT="(4x,a)") line
         n = INDEX(line,"=")
         IF (n == 0.OR.n>len_TRIM(line)-1) STOP
     $        "Error in variable definitions"
         CALL ASSIGN_var(line(:n-1),evaluate(line(n+1:)))
      ENDDO

      !>
      READ (UNIT=iofile,FMT=8000,END=99,ERR=99) 
     +                strho,film,dos,isec1,ndir,secvar
      WRITE (nu,9000) strho,film,dos,isec1,ndir,secvar
 8000 FORMAT (6x,l1,6x,l1,5x,l1,7x,i2,6x,i2,8x,l1)
c
      READ (UNIT=iofile,FMT=7000,END=99,ERR=99) name
      WRITE (nu,9010) name
 7000 FORMAT (10a8)
c
      READ (UNIT=iofile,FMT=7020,END=99,ERR=99)
     +     latnam,namgrp,invs,zrfs,invs2,jspins,l_noco,l_J
      WRITE (nu,9020)
     +     latnam,namgrp,invs,zrfs,invs2,jspins,l_noco,l_J
 7020 FORMAT (a3,1x,a4,6x,l1,6x,l1,7x,l1,8x,i1,8x,l1,5x,l1)
c
      IF ((latnam.EQ.'squ').OR.(latnam.EQ.'hex').OR.
     +    (latnam.EQ.'c-b').OR.(latnam.EQ.'hx3').OR.
     +    (latnam.EQ.'fcc').OR.(latnam.EQ.'bcc')) THEN
          READ (UNIT=iofile,FMT=*,END=99,ERR=99) a1(1)
          WRITE (nu,9030) a1(1)
      ELSEIF ((latnam.EQ.'c-r').OR.(latnam.EQ.'p-r')) THEN
          READ (UNIT=iofile,FMT=*,END=99,ERR=99) a1(1),a2(2)
          WRITE (nu,9030) a1(1),a2(2)
      ELSEIF (latnam.EQ.'obl') THEN
          READ (UNIT=iofile,FMT=*,END=99,ERR=99) a1(1),a1(2)
          READ (UNIT=iofile,FMT=*,END=99,ERR=99) a2(1),a2(2)
          WRITE (nu,9030) a1(1),a1(2)
          WRITE (nu,9030) a2(1),a2(2)
      ELSEIF (latnam.EQ.'any') THEN
          READ (UNIT=iofile,FMT=*,END=99,ERR=99) a1
          READ (UNIT=iofile,FMT=*,END=99,ERR=99) a2
          WRITE (nu,9030) a1(1),a1(2),a1(3)
          WRITE (nu,9030) a2(1),a2(2),a2(3)
      ELSE
          WRITE (6,*) 'rw_inp: latnam ',latnam,' unknown'
          STOP 'rw_inp: latnam'
      ENDIF
c
c
      IF (latnam.EQ.'squ') THEN
         a2(2) = a1(1)
      END IF
c
c     Centered rectangular, special case for bcc(110)
c
      IF (latnam.EQ.'c-b') THEN
         a = a1(1)
         hs2 = sqrt(2.)*0.5e0
         a1(1) = a*hs2
         a1(2) = -a*0.5e0
         a2(1) = a*hs2
         a2(2) = a*0.5e0
      END IF
c
c     Centered rectangular, general case
c     on input: a ---> half of long diagonal
c               b ---> half of short diagonal
c
      IF (latnam.EQ.'c-r') THEN
         a1(2) = -a2(2)
         a2(1) =  a1(1)
      END IF
      IF (latnam.EQ.'hex') THEN
         s3 = sqrt(3.)
         ah = a1(1)/2.
         a1(1) = ah*s3
         a1(2) = -ah
         a2(1) = a1(1)
         a2(2) = ah
      END IF
      IF (latnam.EQ.'hx3') THEN
         s3 = sqrt(3.)
         ah = a1(1)/2.
         a1(1) = ah
         a1(2) = -ah*s3
         a2(1) = a1(1)
         a2(2) = -a1(2)
      END IF

      IF (namgrp.EQ.'any ') THEN
        INQUIRE (file='sym.out',exist=l_sym)
        IF (.not.l_sym)
     +     STOP 'for namgrp="any " please provide a sym-file!'
      ENDIF
      IF (latnam.EQ.'any') THEN
!       STOP 'please specify lattice type (squ,p-r,c-r,hex,hx3,obl)'
        READ (UNIT=iofile,FMT=*,END=99,ERR=99) a3(1),a3(2),a3(3),
     +                                            dvac,scale
        WRITE (nu,9031) a3(1),a3(2),a3(3),dvac,scale
        dtild = a3(3)
      ELSE
        READ (UNIT=iofile,FMT=*,END=99,ERR=99) dvac,dtild,scale
        WRITE (nu,9030) dvac,dtild,scale
        a3(3) = dtild
      ENDIF
 7040 FORMAT (3f10.6)
 7041 FORMAT (5f10.5)
c
      READ (UNIT=iofile,FMT=7110,END=99,ERR=99) namex,relcor
      WRITE (nu,9040) namex,relcor
      namex = TRIM(ADJUSTL(namex))
 7110 FORMAT (a4,3x,a12)
      IF ((namex.EQ.'pw91').OR.(namex.EQ.'l91').OR.
     +    (namex.eq.'pbe').OR.(namex.eq.'rpbe').OR.
     +    (namex.EQ.'Rpbe').OR.(namex.eq.'wc')) THEN                    ! some defaults
        igrd=1 ; lwb=.false. ; ndvgrd=6; idsprs=0 ; chng=-0.1e-11
      ENDIF
!
! look what comes in the next two lines
!
      READ (UNIT=iofile,FMT=7182,END=77,ERR=77) ch_test
      IF (ch_test.EQ.'igr') THEN                          ! GGA input
         BACKSPACE (iofile)
         READ (UNIT=iofile,FMT=7121,END=99,ERR=99) 
     +                   igrd,lwb,ndvgrd,idsprs,chng
         WRITE (nu,9121) igrd,lwb,ndvgrd,idsprs,chng
 7121    FORMAT (5x,i1,5x,l1,8x,i1,8x,i1,6x,d10.3)

         READ (UNIT=iofile,FMT=7182,END=77,ERR=77) ch_test
         IF (ch_test.EQ.'igg') THEN                      ! GGA 2nd line
           GOTO 76
         ELSEIF (ch_test.EQ.'&od') THEN                  ! continue with 1D
           GOTO 78
         ELSE
           GOTO 77
         ENDIF
      ELSEIF (ch_test.EQ.'&od') THEN                  ! continue with 1D
         GOTO 78
      ELSE
         GOTO 77
      ENDIF
c-odim
      READ (UNIT=iofile,FMT=7182,END=99,ERR=99) ch_test
 7182 FORMAT (a3)
   78 IF (ch_test.EQ.'&od') THEN
        BACKSPACE (iofile)
        READ (iofile,odim)
        odd%d1 = d1 ; odd%mb = vM ; odd%M = MM ; odd%m_cyl = m_cyl
        odd%chi = chi ; odd%rot = rot
        odd%invs = invs1 ; odd%zrfs = zrfs1
        WRITE (nu,8182) d1,MM,vM,m_cyl,chi,rot,invs1,zrfs1
      END IF
c+odim
      GOTO 76
   77 BACKSPACE (iofile)                                ! continue with atoms
   76 IF (ch_test /= '&od') THEN
        WRITE (nu,*) '   '
        odd%d1 = .false.
        odd%M = 1 ; odd%mb = 1 ; odd%m_cyl = 1
        odd%chi = 1 ; odd%rot = 1
        odd%invs = .FALSE. ; odd%zrfs = .FALSE.
      END IF
      READ (UNIT=iofile,FMT=*,END=99,ERR=99) ntype
      WRITE (nu,9050) ntype
 7130 FORMAT (i3)
c
      na = 0
      READ (UNIT=iofile,FMT=7110,END=99,ERR=99)
      WRITE (nu,9060)
      DO n=1,ntype
c
         READ (UNIT=iofile,FMT=7140,END=99,ERR=99) noel(n),nz(n),
     +                      ncst(n),lmax(n),jri(n),rmt(n),dx(n)
         WRITE (nu,9070) noel(n),nz(n),ncst(n),lmax(n),jri(n),
     +                      rmt(n),dx(n)
 7140    FORMAT (a3,i3,3i5,2f10.6)
c
c+lda+u
         READ (UNIT=iofile,FMT=7180,END=199,ERR=199) ch_test
 7180    FORMAT (a3)
         IF (ch_test.EQ.'&ld') THEN
           l=0 ; u=0.0 ; j=0.0 ; l_amf = .false.
           BACKSPACE (iofile)
           READ (iofile,ldaU)
           lda_u(n)%l = l ; lda_u(n)%u = u ; lda_u(n)%j = j 
           lda_u(n)%l_amf= l_amf
           WRITE (nu,8180) l,u,j,l_amf
         ELSE
            WRITE (nu,*) '   '
            lda_u(n)%l = -1
         ENDIF
 199     CONTINUE
c-lda+u
c
c---> read extra info for local orbitals, and l_geo. if l_geo=T
c---> calculate force on this atom.
c---> p.kurz 97-06-05
c
         READ (UNIT=iofile,FMT=7160,END=99,ERR=99) neq(n),
     +                    l_geo(n),nlo(n),(llo(ilo,n),ilo=1,nlo(n))
 7160    FORMAT (i2,8x,l1,5x,i2,5x,30i3)
         WRITE (nu,9090) neq(n),l_geo(n),nlo(n),
     +                                    (llo(ilo,n),ilo=1,nlo(n))
c
         DO ieq=1,neq(n)
            na = na + 1
            READ (UNIT = iofile,FMT = 7170,iostat = ierr) 
     +                      (taual(i,na),i=1,3),scpos
            IF (ierr.NE.0) THEN
               BACKSPACE(iofile)
               !<-- read positions with new format
               READ (UNIT = iofile,FMT ="(a)",END = 99,ERR = 99) line
               taual(1,na) = evaluatefirst(line)
               taual(2,na) = evaluatefirst(line)
               taual(3,na) = evaluatefirst(line)
               scpos = evaluatefirst(line)
               IF (scpos == 0.0)  scpos          = 1.0
               !>
            ENDIF
            WRITE (nu,9100) (taual(i,na),i=1,3),scpos
 7170       FORMAT (4f10.6)
            IF (scpos.EQ.0.) scpos = 1.
            DO i = 1,2
               taual(i,na) = taual(i,na)/scpos
            ENDDO
            IF (.not.film) taual(3,na) = taual(3,na)/scpos
c+odim
c in 1D case all the coordinates are given cartesian'ly
            IF (odd%d1) THEN
               taual(1,na) = taual(1,na)/a1(1)
               taual(2,na) = taual(2,na)/a2(2)
            END IF
c-odim
         ENDDO
         READ (iofile,*)
         WRITE (nu,9060)
      ENDDO
c
      READ (UNIT=iofile,FMT=7210,END=99,ERR=99) gmax,gmaxxc
      WRITE (nu,9110) gmax,gmaxxc
 7210 FORMAT (2f10.6)
c
      INQUIRE(file='fl7para',exist=ldum)  ! fl7para must not exist for gw=2
      IF (gw.eq.-1.and.ldum) THEN         ! in the first run of rw_inp
        ldum = .true.                     ! (then, gw=-1 at this point).
      ELSE                                !
        ldum = .false.                    !
      ENDIF                               !
      gw = 0
      READ (UNIT=iofile,FMT=7220,END=99,ERR=7215)
     &                                       vchk,cdinf,pot8,gw,gw_neigd
 7215 IF(strho) gw=0
      IF(gw.eq.2) THEN
        IF(ldum) STOP 'rw_inp: Remove fl7para before run with gw=2!'
        IF(gw_neigd.eq.0) STOP 'rw_inp: No gw_neigd-value given.'
        IF(.not.pot8)     STOP 'rw_inp: pot8 must be set if gw=2!'
      ELSE
        INQUIRE(file='QGpsi',exist=ldum)
        IF(ldum)          STOP 'rw_inp: QGpsi exists but gw/=2 in inp.'
      ENDIF
      BACKSPACE(iofile)                                         ! Make sure that vchk,cdinf,pot8 are all given.
      READ (UNIT=iofile,FMT=7220,END=99,ERR=99) vchk,cdinf,pot8 !
      WRITE (nu,9120) vchk,cdinf,pot8,gw,gw_neigd
 7220 FORMAT (5x,l1,1x,6x,l1,1x,5x,l1,1x,3x,i1,1x,9x,i3)
c
      DO i=1,100 ; line(i:i)=' ' ; ENDDO
      READ (UNIT=iofile,fmt='(A)',END=99,ERR=99) line 
      eig66(2)= ( line(38:44)=='soc66=T' ).or.( line(38:44)=='soc66=t' )
      BACKSPACE (UNIT=iofile)
      READ (UNIT=iofile,FMT=6000,END=99,ERR=99) 
     +                lpr,form66,l_f,eonly,eig66(1)
      WRITE (nu,9130) lpr,form66,l_f,eonly,eig66(1),eig66(2)
 6000 FORMAT (4x,i1,8x,l1,5x,l1,7x,l1,7x,l1)
c
c+roa
      WRITE (chntype,'(i3)') ntype
      chform = '('//chntype//'i3 )'
      READ (UNIT=iofile,FMT=chform,END=99,ERR=99) 
     +                (lnonsph(n),n=1,ntype)
      WRITE (nu,FMT=chform) (lnonsph(n),n=1,ntype)
 6010 FORMAT (25i3)
c
      READ (UNIT=iofile,FMT=6010,END=99,ERR=99) nwd,lepr
      WRITE (nu,9140) nwd,lepr
c
      zc=0.0
      DO nw=1,nwd
         READ (UNIT=iofile,FMT=*,END=99,ERR=99)
         WRITE (nu,'(a8,i2)') 'Window #',nw
c
         READ (UNIT=iofile,FMT=6040,END=99,ERR=99) 
     +                   ellow(nw),elup(nw),zelec(nw)
         WRITE (nu,9150) ellow(nw),elup(nw),zelec(nw) 
 6040    FORMAT (4f10.5)
         zc = zc + zelec(nw)
c
         READ (UNIT=iofile,FMT='(f10.5)',END=99,ERR=99) rkm(nw)
         WRITE (nu,FMT='(f10.5,1x,A)') rkm(nw), '=kmax' 
      ENDDO
c
      READ (UNIT=iofile,FMT=8010,END=99,ERR=99) gauss,delgau,tria
      WRITE (nu,9160) gauss,delgau,tria
 8010 FORMAT (6x,l1,f10.5,5x,l1)
c
      tkb=delgau
c
      READ(iofile,fmt='(27x,l1)',END=99,ERR=99) l_soc  
      DO i=1,100 ; line(i:i)=' ' ; ENDDO
      BACKSPACE(iofile)
      READ(iofile,fmt='(A)',END=99,ERR=99) line
      BACKSPACE(iofile)
      IF (line(9:10)=='pi') THEN
        READ(iofile,fmt='(f8.4)') theta
        theta= theta*4.*ATAN(1.)
      ELSE
        READ(iofile,fmt='(f10.6)',END=99,ERR=99) theta
      ENDIF
      BACKSPACE(iofile)
      IF (line(19:20)=='pi') THEN
        READ(iofile,fmt='(10x,f8.4)',END=99,ERR=99) phi 
        phi= phi*4.*ATAN(1.)
      ELSE
        READ(iofile,fmt='(10x,f10.6)',END=99,ERR=99) phi
      ENDIF
      IF ( line(30:34)=='spav=' ) THEN
        BACKSPACE(iofile)
        READ(iofile,fmt='(34x,l1)',END=99,ERR=99) soc_opt(ntype+2) 
      ELSE
        soc_opt(ntype+2)= .false. 
      ENDIF 
      IF ( line(37:40)=='off=' ) THEN 
        BACKSPACE(iofile) 
        chform= '(40x,l1,1x,'//chntype//'a1)'
        READ(iofile,fmt=chform,END=99,ERR=99) 
     &   soc_opt(ntype+1),(helpchar(i),i=1,ntype) 
        DO i= 1,ntype
          soc_opt(i)= (helpchar(i)=='1') 
        ENDDO
      ELSE
        DO i= 1,ntype+1 
          soc_opt(i)= .false.
        ENDDO 
      ENDIF 
      IF (soc_opt(ntype+1)) THEN
        DO i= 1,ntype
          IF (soc_opt(i)) THEN
            helpchar(i)= '1'
          ELSE
            helpchar(i)= '0'
          ENDIF
        ENDDO  
        chform= '(2f10.6,a7,l1,a6,l1,a5,l1,a1,'//chntype//'a1)'
        WRITE(nu,fmt=chform) 
     &   theta,phi,',l_soc=',l_soc,',spav=',soc_opt(ntype+2),
     &   ',off=',soc_opt(ntype+1),',',(helpchar(i),i=1,ntype) 
      ELSE
        WRITE(nu,fmt='(2f10.6,a7,l1,a6,l1,a5,l1)') 
     &   theta,phi,',l_soc=',l_soc,',spav=',soc_opt(ntype+2),
     &   ',off=',soc_opt(ntype+1)
      ENDIF 
c
      l_u2f=.false.
      l_f2u=.false.
      READ (UNIT=iofile,FMT=8050,END=99,ERR=99) 
     +                 frcor,slice,ctail,disp,kcrel,l_u2f,l_f2u
      BACKSPACE(iofile)
      READ (UNIT=iofile,fmt='(A)') line
      l_bmt= ( line(52:56)=='bmt=T' ).or.( line(52:56)=='bmt=t' )
      WRITE (nu,9170)  frcor,slice,ctail,disp,kcrel,l_u2f,l_f2u,l_bmt
 8050 FORMAT (6x,l1,7x,l1,7x,l1,6x,l1,7x,i1,5x,l1,5x,l1)
      READ (UNIT=iofile,FMT=8060,END=99,ERR=99) 
     +                 itmax,maxiter,imix,alpha,spinf
      WRITE (nu,9180)  itmax,maxiter,imix,alpha,spinf
 8060 FORMAT (6x,i2,9x,i3,6x,i2,7x,f6.2,7x,f6.2)
      chform = '(5x,l1,'//chntype//'f6.2)'
c      chform = '(5x,l1,23f6.2)'
      READ (UNIT=iofile,FMT=chform,END=99,ERR=99) 
     +                                   swsp, (bmu(i),i=1,ntype)
      chform = '(6x,l1,'//chntype//'i3 )'
c      chform = '(6x,l1,23i3 )'
      READ (UNIT=iofile,FMT=chform,END=99,ERR=99) 
     +                                   lflip, (nflip(i),i=1,ntype)
c-
      chform = '("swsp=",l1,'//chntype//'f6.2)'
c      chform = '("swsp=",l1,23f6.2)'
      WRITE (nu,FMT=chform) swsp, (bmu(i),i=1,ntype)
      chform = '("lflip=",l1,'//chntype//'i3 )'
c      chform = '("lflip=",l1,23i3 )'
      WRITE (nu,FMT=chform) lflip, (nflip(i),i=1,ntype)
c-roa
c+stm
      READ (UNIT=iofile,FMT=8075,END=99,ERR=99) 
     +      vacdos,layers,integ,starcoeff,nstars,
     +      locx(1),locy(1),locx(2),locy(2),nstm,tworkf
      WRITE (nu,9210) vacdos,layers,integ,starcoeff,nstars,
     +      locx(1),locy(1),locx(2),locy(2),nstm,tworkf
 8075 FORMAT (7x,l1,8x,i2,7x,l1,6x,l1,8x,i2,4(4x,f5.2),6x,i1,8x,f10.6)
c-stm
      IF (vacdos) THEN
        IF (integ) THEN
          READ (UNIT=iofile,FMT=8076,END=99,ERR=99) 
     +                    ((izlay(i,k),k=1,2),i=1,layers)
          WRITE (nu,9220) ((izlay(i,k),k=1,2),i=1,layers)
 8076     FORMAT (10(2(i3,1x),1x))
        ELSE
          READ (UNIT=iofile,FMT=8077,END=99,ERR=99) 
     +                    (izlay(i,1),i=1,layers)
          WRITE (nu,9230) (izlay(i,1),i=1,layers)
 8077     FORMAT (20(i3,1x))
        END IF
      ELSE
        READ (UNIT=iofile,FMT=*,END=99,ERR=99)
        WRITE (nu,fmt='(x)') 
      END IF
c
      band = .false.
      READ (UNIT=iofile,FMT=8050,END=992,ERR=992) iplot,score,plpot,band
      WRITE (nu,9240) iplot,score,plpot,band
      IF (band) THEN
        dos=.true. ; ndir = -4
      ENDIF
      GOTO 993
 992  BACKSPACE(iofile)
      READ (UNIT=iofile,FMT=8050,END=99,ERR=99) iplot,score,plpot
      WRITE (nu,9240) iplot,score,plpot,band
c
 993  READ (UNIT=iofile,FMT='(i3,2f10.6,6x,i3,8x,l1)',END=99,ERR=99) 
     +                kk,e1s,e2s,nnne,pallst
      WRITE (nu,9250) kk,e1s,e2s,nnne,pallst
c
      READ (UNIT=iofile,FMT=8090,END=99,ERR=99) 
     +                xa,thetad,epsdisp,epsforce
      WRITE (nu,9260) xa,thetad,epsdisp,epsforce
 8090 FORMAT (3x,f10.5,8x,f10.5,9x,f10.5,10x,f10.5)
c

c+/-odim YM : changed to '70' in the format, sometimes caused probl.
      chform = '(6x,'//chntype//'(3i1,1x))'
      READ (UNIT=iofile,FMT=chform,END=99,ERR=99) 
     +      ((relax(i,k),i=1,3),k=1,ntype)
      chform = '("relax ",'//chntype//'(3i1,1x))'
      WRITE (nu,FMT=chform) ((relax(i,k),i=1,3),k=1,ntype)

c read dos_params! These will be set automatically if not present!
      e1_dos=0.0
      e2_dos=-1.0
      sig_dos=1e-4
      READ (UNIT=iofile,FMT='(9x,f10.5,10x,f10.5,9x,f10.5)',
     +     END=98,ERR=98) e2_dos,e1_dos,sig_dos
 98   Continue
      WRITE (nu,'(a,f10.5,a,f10.5,a,f10.5)')
     +     'emin_dos=',e2_dos,',emax_dos=',e1_dos,',sig_dos=',sig_dos      
      CLOSE (iofile)

c---------------------------------------------------------------------
      ELSEIF ((ch_rw.eq.'W').OR.(ch_rw.eq.'w'))  THEN
c---------------------------------------------------------------------

      IF (ch_rw.eq.'W') THEN
      OPEN (iofile,file='inp_new',form='formatted',status='unknown')
      REWIND (iofile)
      ENDIF

      WRITE (iofile,9000) strho,film,dos,isec1,ndir,secvar
 9000 FORMAT ('strho=',l1,',film=',l1,',dos=',l1,',isec1=',i2,
     +        ',ndir=',i2,',secvar=',l1)
      WRITE (iofile,9010) name
 9010 FORMAT (10a8)
      WRITE(iofile,9020) latnam,namgrp,invs,zrfs,invs2,jspins,l_noco,l_J
 9020 FORMAT (a3,1x,a4,',invs=',l1,',zrfs=',l1,',invs2=',l1,
     +       ',jspins=',i1,',l_noco=',l1,',l_J=',l1)
c
      IF (latnam.EQ.'c-b') THEN
         a1(1) = sqrt(2.)* a1(1)
      END IF
      IF (latnam.EQ.'hex') THEN
         s3 = sqrt(3.)
         a1(1) = 2*a1(1)/sqrt(3.)
      END IF
      IF (latnam.EQ.'hx3') THEN
         a1(1) = 2*a1(1) 
      END IF
c
      IF ((latnam.EQ.'squ').OR.(latnam.EQ.'hex').OR.
     +    (latnam.EQ.'c-b').OR.(latnam.EQ.'hx3')) THEN
          WRITE (iofile,9030) a1(1)
      ELSEIF ((latnam.EQ.'c-r').OR.(latnam.EQ.'p-r')) THEN
          WRITE (iofile,9030) a1(1),a2(2)
      ELSEIF (latnam.EQ.'obl') THEN
          WRITE (iofile,9030) a1(1),a1(2)
          WRITE (iofile,9030) a2(1),a2(2)
      ELSEIF (latnam.EQ.'any') THEN
          WRITE (iofile,9030) a1(1),a1(2),a1(3)
          WRITE (iofile,9030) a2(1),a2(2),a2(3)
      ELSE
          WRITE (6,*) 'rw_inp: latnam ',latnam,' unknown'
          STOP 'rw_inp: latnam'
      ENDIF
c
      IF (latnam.EQ.'any') THEN
        WRITE (iofile,9031)  a3(1),a3(2),a3(3),dvac,scale
        dtild = a3(3)
      ELSE
        WRITE (iofile,9030) dvac,dtild,scale
        a3(3) = scale * dtild
      ENDIF
 9030 FORMAT (3f15.8)
 9031 FORMAT (5f15.8)
      WRITE (iofile,9040) namex,relcor
 9040 FORMAT (a4,3x,a12)
      IF ((namex.EQ.'pw91').OR.(namex.EQ.'l91').OR.
     +    (namex.eq.'pbe').OR.(namex.eq.'rpbe').OR.
     +    (namex.EQ.'Rpbe').OR.(namex.eq.'wc') ) THEN
         WRITE (iofile,FMT=9121) igrd,lwb,ndvgrd,idsprs,chng
 9121    FORMAT ('igrd=',i1,',lwb=',l1,',ndvgrd=',i1,',idsprs=',i1,
     +           ',chng=',d10.3)
      ENDIF
c-odim
      IF (odd%d1) THEN
        WRITE (iofile,8182) odd%d1,odd%M,odd%mb,
     &                      odd%m_cyl,odd%chi,odd%rot,odd%invs,odd%zrfs
 8182   FORMAT ('&odim d1=',l1,',MM=',i3,',vM=',i3,
     &          ',m_cyl=',i3,',chi=',i3,',rot=',i3,
     &          ',invs1=',l1,',zrfs1=',l1,'/')
      ELSE
        WRITE (iofile,*) '   '
      END IF
c+odim
      WRITE (iofile,9050) ntype
 9050 FORMAT (i3)
      na = 0
      WRITE (iofile,9060)
 9060 FORMAT ('**********************************')
      DO n=1,ntype
         WRITE (iofile,9070) noel(n),nz(n),ncst(n),lmax(n),jri(n),
     +                       rmt(n),dx(n)
 9070    FORMAT (a3,i3,3i5,2f10.6)
c+lda_u
         IF (lda_u(n)%l.GE.0) THEN
           WRITE (iofile,8180) lda_u(n)%l,lda_u(n)%u,lda_u(n)%j,
     +                                               lda_u(n)%l_amf
 8180      FORMAT ('&ldaU l=',i1,',u=',f4.2,',j=',f4.2,',l_amf=',l1,'/')
         ELSE
            WRITE (iofile,*) '   '
         ENDIF
c-lda_u
         WRITE (iofile,9090) neq(n),l_geo(n),nlo(n),
     +                                     (llo(ilo,n),ilo=1,nlo(n))
 9090    FORMAT (i2,',force =',l1,',nlo=',i2,',llo=',30i3)
         DO ieq=1,neq(n)
            na = na + 1
            DO i = 2,9
              rest = ABS(i*taual(1,na) - NINT(i*taual(1,na)) )
     +             + ABS(i*taual(2,na) - NINT(i*taual(2,na)) )
              IF (rest.LT.(i*0.000001)) EXIT
            ENDDO
            IF (i.LT.10) scpos = real(i)  ! common factor found (x,y)
            IF (.not.film) THEN           ! now check z-coordinate
              DO i = 2,9
                rest = ABS(i*scpos*taual(3,na) -
     +                NINT(i*scpos*taual(3,na)) )
                IF (rest.LT.(i*scpos*0.000001)) THEN
                  scpos = i*scpos
                  EXIT
                ENDIF
              ENDDO
            ENDIF
            DO i = 1,2
               taual(i,na) = taual(i,na)*scpos
            ENDDO
            IF (.not.film) taual(3,na) = taual(3,na)*scpos
            IF (film) taual(3,na) = a3(3)*taual(3,na)/scale
c+odim in 1D case all the coordinates are given in cartesian YM
            IF (odd%d1) THEN
               taual(1,na) = taual(1,na)*a1(1)
               taual(2,na) = taual(2,na)*a2(2)
            END IF
c-odim
            WRITE (iofile,9100) (taual(i,na),i=1,3),scpos
 9100       FORMAT (4f10.6)
         ENDDO
         WRITE (iofile,9060)
      ENDDO
      IF ((gmaxxc.LE.0).OR.(gmaxxc.GT.gmax)) gmaxxc=gmax
      WRITE (iofile,9110) gmax,gmaxxc
 9110 FORMAT (2f10.6)
      WRITE (iofile,9120) vchk,cdinf,pot8,gw,gw_neigd
 9120 FORMAT ('vchk=',l1,',cdinf=',l1,',pot8=',l1,',gw=',i1,
     &        ',gw_neigd=',i3)
      WRITE (iofile,9130) lpr,form66,l_f,eonly,eig66(1),eig66(2)
 9130 FORMAT ('lpr=',i1,',form66=',l1,',l_f=',l1,',eonly=',l1,
     +        ',eig66=',l1,',soc66=',l1)
      WRITE (chntype,'(i3)') ntype
      chform = '('//chntype//'i3 )'
      WRITE (iofile,FMT=chform) (lnonsph(n),n=1,ntype)
 9140 FORMAT (25i3)
      WRITE (iofile,9140) nwd,lepr
      DO nw=1,nwd
         WRITE (iofile,'(a8,i2)') 'Window #',nw
         WRITE (iofile,9150) ellow(nw),elup(nw),zelec(nw)
 9150    FORMAT (4f10.5)
         WRITE (iofile,fmt='(f10.5,1x,A)') rkm(nw), '=kmax'
      ENDDO
      WRITE (iofile,9160) gauss,delgau,tria
 9160 FORMAT ('gauss=',l1,f10.5,'tria=',l1)
      IF (soc_opt(ntype+1)) THEN
        DO i= 1,ntype
          IF (soc_opt(i)) THEN
            helpchar(i)= '1'
          ELSE
            helpchar(i)= '0'
          ENDIF
        ENDDO 
        chform= '(2f10.6,a7,l1,a6,l1,a5,l1,a1,'//chntype//'a1)'
        WRITE(iofile,fmt=chform) 
     &   theta,phi,',l_soc=',l_soc,',spav=',soc_opt(ntype+2),
     &   ',off=',soc_opt(ntype+1),',',(helpchar(i),i=1,ntype)     
      ELSE
        WRITE(iofile,fmt='(2f10.6,a7,l1,a6,l1,a5,l1)') 
     &   theta,phi,',l_soc=',l_soc,',spav=',soc_opt(ntype+2),
     &   ',off=',soc_opt(ntype+1)
      ENDIF
      WRITE (iofile,9170) frcor,slice,ctail,disp,kcrel,l_u2f,l_f2u,l_bmt
 9170 FORMAT ('frcor=',l1,',slice=',l1,',ctail=',l1,',disp=',
     +        l1,',kcrel=',i1,',u2f=',l1,',f2u=',l1,',bmt=',l1)
      WRITE (iofile,9180) itmax,maxiter,imix,alpha,spinf
 9180 FORMAT ('itmax=',i2,',maxiter=',i3,',imix=',i2,',alpha=',
     +        f6.2,',spinf=',f6.2)
c+roa
      WRITE (chntype,'(i3)') ntype
      chform = '("swsp=",l1,'//chntype//'f6.2)'
      WRITE (iofile,FMT=chform) swsp, (bmu(i),i=1,ntype)
      chform = '("lflip=",l1,'//chntype//'i3 )'
      WRITE (iofile,FMT=chform) lflip, (nflip(i),i=1,ntype)
c-roa
c+stm
      WRITE (iofile,9210) vacdos,layers,integ,starcoeff,nstars,
     +      locx(1),locy(1),locx(2),locy(2),nstm,tworkf
 9210 FORMAT ('vacdos=',l1,',layers=',i2,',integ=',l1,',star=',l1,
     + ',nstars=',i2,4(4x,f5.2),',nstm=',i1,',tworkf=',f10.6)
c-stm
      IF (vacdos) THEN
        IF (integ) THEN
          WRITE (iofile,9220) ((izlay(i,k),k=1,2),i=1,layers)
 9220     FORMAT (10(2(i3,1x),1x))
        ELSE
          WRITE (iofile,9230) (izlay(i,1),i=1,layers)
 9230     FORMAT (20(i3,1x))
        END IF
      ELSE
        WRITE (iofile,*)
      END IF
      band = .false.
      WRITE (iofile,9240) iplot,score,plpot,band
 9240 FORMAT ('iplot=',l1,',score=',l1,',plpot=',l1,',band=',l1)
      WRITE (iofile,9250) kk,e1s,e2s,nnne,pallst
 9250 FORMAT (i3,2f10.6,',nnne=',i3,',pallst=',l1)
      WRITE (iofile,9260) xa,thetad,epsdisp,epsforce
 9260 FORMAT ('xa=',f10.5,',thetad=',f10.5,',epsdisp=',f10.5,
     +        ',epsforce=',f10.5)
c+/-gb
      chform = '("relax ",'//chntype//'(3i1,1x))'
      WRITE (iofile,FMT=chform) ((relax(i,k),i=1,3),k=1,ntype)
      WRITE (iofile,'(a,f10.5,a,f10.5,a,f10.5)') 
     +     'emin_dos=',e2_dos,',emax_dos=',e1_dos,',sig_dos=',sig_dos

      IF (ch_rw.eq.'W') CLOSE (iofile)
      ELSE 
        WRITE (6,*) 'specify either W to write or R to read!'
      ENDIF

      RETURN
  99  CLOSE (nu)
      WRITE (6,*) 'Error reading inp-file'
      STOP 'error reading inp in rw_inp'

      END SUBROUTINE rw_inp
      END MODULE m_rwinp
