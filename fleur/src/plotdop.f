      MODULE m_plotdop
!     +++++++++++++++++++++++++++++++++++++++++++++++++
!     plot the charge density for fleur  code output
!     
!     if twodim = .false. a 3-D plot with nplot layers in z-direction
!     is constructed; the 3x3 matrix gives the 3 vectors of the cell ..
!     .gustav
!
!    Changed the input/output for easier use. 
!    This subroutine uses the file plot_inp for input. 
!    The old plotin-file is still supported by the old subroutine at the end of the module
!                      Juelich, 21.1.06 DW
!
!     +++++++++++++++++++++++++++++++++++++++++++++++++
      PRIVATE
      PUBLIC plotdop
      CONTAINS
      SUBROUTINE plotdop(odi,ods,
     >     n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,ntype,
     >     lmaxd,jmtd,ntypd,natd,nmzd,jspins,neq,
     >     nq3,nvac,nmz,nmzxy,nq2,nop,nop2,volint,film,
     >     plpot,score,slice,symor,invs,invs2,z1,delz,
     >     clnu,llh,nmem,mlh,nlh,ngopr,ntypsy,jri,pos,zatom,
     >     nstr,nstr2,lmax,kv2,kv3,mrot,tau,rmsh,invtab,
     >     amat,bmat,l_noco,cdnfname)
!    *****************************************************
      USE m_constants, ONLY : pimach
      USE m_od_types, ONLY : od_inp, od_sym
      USE m_outcdn
      USE m_loddop
      USE m_xsf_io
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd
      INTEGER, INTENT (IN) :: lmaxd,jmtd,ntypd,natd,nmzd,jspins
      INTEGER, INTENT (IN) :: nq3,nvac,nmz,nmzxy,nq2,nop,nop2,ntype
      LOGICAL, INTENT (IN) :: plpot,symor,invs,score,slice,invs2,film
      LOGICAL, INTENT (IN) :: l_noco
      REAL,    INTENT (IN) :: z1,delz,volint
      CHARACTER*10, INTENT (IN) :: cdnfname
!     ..
!     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ngopr(natd),ntypsy(natd),lmax(ntypd)
      INTEGER, INTENT (IN) :: nstr(n3d),nstr2(n2d),jri(ntypd),neq(ntypd)
      INTEGER, INTENT (IN) :: kv2(2,n2d),kv3(3,n3d),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: invtab(nop)
      REAL,    INTENT (IN) :: zatom(:)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3),pos(3,natd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),tau(3,nop)
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
!     ..
!     .. Local Scalars ..
      REAL          :: s,tec,qint,sfp,xdnout
      INTEGER       :: i,i1,i2,i3,ii3,ix,iy,iz,jsp,na,nplo
      INTEGER       :: nplot,nq,nt,jm,jspin,iter,ii1,ii2
      CHARACTER*8   :: dop,iop
      LOGICAL       :: twodim,oldform,newform
!     ..
!     .. Local Arrays ..
      COMPLEX :: qpw(n3d,jspd),rhtxy(nmzxyd,n2d-1,2,jspd)
      REAL    :: rho(jmtd,0:nlhd,ntypd,jspd),rht(nmzd,2,jspd)
      REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3)
      INTEGER :: grid(3)
      LOGICAL :: cartesian,xsf
      REAL    :: rhocc(jmtd)
      REAL    :: point(3)
      CHARACTER (len=30) :: filename
      CHARACTER (len=8)  :: name(10)
      CHARACTER (len=7)  :: append
      NAMELIST /plot/twodim,cartesian,vec1,vec2,vec3,grid,zero,filename


      oldform = .false.
      INQUIRE(file ="plotin",exist = oldform) 
      IF ( oldform ) THEN 
         CALL priv_old_plot(odi,ods,
     >                n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,ntype,
     >                lmaxd,jmtd,ntypd,natd,nmzd,jspins,neq,
     >                nq3,nvac,nmz,nmzxy,nq2,nop,nop2,volint,film,
     >                plpot,score,slice,symor,invs,invs2,z1,delz,
     >                clnu,llh,nmem,mlh,nlh,ngopr,ntypsy,jri,pos,
     >                nstr,nstr2,lmax,kv2,kv3,mrot,tau,rmsh,invtab,
     >                amat,bmat,l_noco,cdnfname)
         RETURN
      ENDIF
      INQUIRE(file ="plot_inp",exist= newform)
      IF (.NOT.newform) THEN !no input file exists, create a template and
                            !exit
         OPEN(20,file ="plot_inp")
         WRITE(20,'(i2,a5,l1)') 2,",xsf=",.true.
         WRITE(20,*) "&PLOT twodim=t,cartesian=t"
         WRITE(20,*) "  vec1(1)=10.0 vec2(2)=10.0"
         WRITE(20,*) "  filename='plot1' /"
         WRITE(20,*) "&PLOT twodim=f,cartesian=f"
         WRITE(20,*) "  vec1(1)=1.0 vec1(2)=0.0 vec1(3)=0.0 "
         WRITE(20,*) "  vec2(1)=0.0 vec2(2)=1.0 vec2(3)=0.0 "
         WRITE(20,*) "  vec3(1)=0.0 vec3(2)=0.0 vec3(3)=1.0 "
         WRITE(20,*) "  grid(1)=30  grid(2)=30  grid(3)=30  "
         WRITE(20,*) "  zero(1)=0.0 zero(2)=0.0 zero(3)=0.5 "
         WRITE(20,*) "  filename ='plot2' /"
         CLOSE(20)
         WRITE(*,*) "No plot_inp file found. Created a template"
         STOP "Missing input for plot; modify plot_inp"
      ENDIF

      ! new input
      !<-- Open the charge/potential file
      OPEN (20,file = cdnfname,form='unformatted',status='old')
      CALL loddop(
     >            jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,20,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
      !<--Perhaps only the core charge should be plotted
      IF (score) THEN
         OPEN (17,file='cdnc',form='unformatted',status='old')
         REWIND 17
         DO jspin = 1,jspins
            DO nt = 1,ntype
               jm = jri(nt)
               READ (17) (rhocc(i),i=1,jm)
               DO i = 1,jri(nt)
                  rho(i,0,nt,jspin) = rho(i,0,nt,jspin) - rhocc(i)/2.0
     $                 /SQRT( pimach() )
               ENDDO
               READ (17) tec
            ENDDO
            READ (17) qint
            qpw(1,jspin) = qpw(1,jspin) - qint/volint
         ENDDO
         CLOSE (17)
      END IF
      !>
      !>
      !<-- Open the plot_inp file for input
      OPEN (18,file='plot_inp')
      READ(18,'(i2,5x,l1)')    nplot,xsf
      ! If xsf is specified we create an input file for xcrysden
      IF (xsf) THEN
         IF (l_noco) THEN
             append = '_pl.xsf'
             OPEN (55,file = trim(cdnfname)//append,form='formatted')
         ELSE
             OPEN(55,file="plot.xsf")
         ENDIF
         CALL xsf_WRITE_atoms(
     >                        55,film,odi%d1,amat,neq(:ntype),
     >                        zatom(:ntype),pos)
      ENDIF
      !<-- Loop over all plots
      DO nplo=1,nplot
         ! the defaults
         twodim = .TRUE.;cartesian=.TRUE.;grid=(/100,100,100/)
         vec1 = (/0.,0.,0./);vec2=(/0.,0.,0./);vec3=(/0.,0.,0./)
         zero = (/0.,0.,0./);filename="default"
         READ(18,plot)
         IF (twodim.AND.ANY(grid(1:2)<1)) 
     +        STOP "Illegal grid size in plot"
         IF (.NOT.twodim.AND.ANY(grid<1)) 
     +        STOP "Illegal grid size in plot"
         IF (twodim) grid(3) = 1
         !calculate cartesian coordinates if needed
         IF (.NOT.cartesian) THEN
            vec1=matmul(amat,vec1)
            vec2=matmul(amat,vec2)
            vec3=matmul(amat,vec3)
            zero=matmul(amat,zero)
         ENDIF
         !Open the file
         IF (filename =="default") WRITE(filename,'(a,i2)') "plot",nplo
         IF (xsf) THEN
            CALL xsf_WRITE_header(55,twodim,filename,vec1,vec2,vec3,zero
     $           ,grid)
         ELSE
            IF (l_noco) THEN
               OPEN (55,file = filename//cdnfname,form='formatted')
            ELSE
               OPEN (55,file = filename,form='formatted')
            ENDIF
         ENDIF
         !loop over spins
         DO jsp = 1,jspins
            !loop over all points
            DO iz = 0,grid(3)-1
               DO iy = 0,grid(2)-1
                  xloop:DO ix = 0,grid(1)-1
                    point = zero + vec1*REAL(ix)/(grid(1)-1) +
     +                             vec2*REAL(iy)/(grid(2)-1)
                    IF (.NOT.twodim) point = point +
     +                             vec3*REAL(iz)/(grid(3)-1)
                    !Check if the point is in MT-sphere
                    ii1 = 3
                    ii2 = 3
                    ii3 = 3
                    IF (film .AND. .NOT.odi%d1) ii3 = 0
                    IF (odi%d1) THEN
                       ii1 = 0 ; ii2 = 0
                    END IF
                    DO  i1 = -ii1,ii1
                       DO  i2 = -ii2,ii2
                          DO  i3 = -ii3,ii3
                             pt = point+MATMUL(amat,(/i1,i2,i3/))
                             na = 0
                             DO nt = 1,ntype
                                DO nq = 1,neq(nt)
                                   na   = na + 1
                                   s  = SQRT(dot_PRODUCT(pos(:,na)
     $                                  -pt,pos(:,na)-pt))
                                   IF (s<rmsh(jri(nt),nt)) THEN
                                      CALL outcdn(
     >                                     pt,nt,na,0,1,jsp,plpot,kv2
     $                                     ,kv3,nstr,nstr2,n3d,jspd
     $                                     ,nmzxyd,n2d,memd,nlhd,ntypsd
     $                                     ,lmaxd,jmtd,natd,ntypd,nmzd
     $                                     ,nop,nop2,mrot,tau,symor
     $                                     ,invtab,nq3,nvac,invs,z1,delz
     $                                     ,nmz,nmzxy,nq2,mlh,llh,clnu
     $                                     ,nmem,nlh,lmax,rmsh,jri,pos
     $                                     ,ngopr,ntypsy,amat,bmat,qpw,
     $                                     rhtxy,rho,rht,odi,ods,xdnout)
                                      IF (xsf) THEN
                                         write(55,*) xdnout
                                      ELSE
                                         IF (twodim) THEN
                                            WRITE(55,'(3e15.7)')
     $                                           point(1:2),xdnout
                                         ELSE
                                            WRITE(55,'(4e15.7)') point
     $                                           ,xdnout
                                         ENDIF
                                      ENDIF
                                      CYCLE xloop
                                   ENDIF
                                ENDDO
                             ENDDO !nt
                          ENDDO
                       ENDDO
                    ENDDO !i1
                    !Check for point in vacuum

                    IF (film.AND..NOT.odi%d1.AND.ABS(point(3))>=z1) THEN
                       CALL outcdn(
     >                      point,0,0,1,0,jsp,plpot,kv2,kv3,nstr,nstr2,
     >                      n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,lmaxd
     $                      ,jmtd,natd,ntypd,nmzd,nop,nop2,mrot,tau
     $                      ,symor,invtab,nq3,nvac,invs,z1,delz,nmz
     $                      ,nmzxy,nq2,mlh,llh,clnu,nmem,nlh,lmax,rmsh
     $                      ,jri,pos,ngopr,ntypsy,amat,bmat,qpw,rhtxy
     $                      ,rho,rht,odi,ods,xdnout)
                       IF (xsf) THEN
                          write(55,*) xdnout
                       ELSE
                          IF (twodim) THEN
                             WRITE(55,'(3e15.7)') point(1:2),xdnout
                          ELSE
                             WRITE(55,'(4e15.7)') point(:),xdnout
                          ENDIF
                       ENDIF
                       CYCLE xloop
                    END IF
                    IF (odi%d1) THEN
                       IF (SQRT((pt(1))**2 + (pt(2))**2)>=z1) THEN
                          CALL outcdn(
     >                         pt,0,0,1,0,jsp,plpot,kv2,kv3,nstr,nstr2
     $                         ,n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd
     $                         ,lmaxd,jmtd,natd,ntypd,nmzd,nop,nop2,mrot
     $                         ,tau,symor,invtab,nq3,nvac,invs,z1,delz
     $                         ,nmz,nmzxy,nq2,mlh,llh,clnu,nmem,nlh,lmax
     $                         ,rmsh,jri,pos,ngopr,ntypsy,amat,bmat,qpw
     $                         ,rhtxy,rho,rht,odi,ods,xdnout)
                          IF (xsf) THEN
                             WRITE(55,*) xdnout
                          ELSE
                             IF (twodim) THEN
                                WRITE (55,'(3e15.7)') point(1:2),xdnout
                             ELSE
                                WRITE (55,'(4e15.7)') point(:),xdnout
                             ENDIF
                          ENDIF
                          CYCLE xloop
                       END IF
                    END IF
                    CALL outcdn(
     >                   point,0,0,0,2,jsp,plpot,kv2,kv3,nstr,nstr2,
     >                   n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,lmaxd,jmtd
     $                   ,natd,ntypd,nmzd,nop,nop2,mrot,tau,symor,invtab
     $                   ,nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,mlh,llh
     $                   ,clnu,nmem,nlh,lmax,rmsh,jri,pos,ngopr,ntypsy
     $                   ,amat,bmat,qpw,rhtxy,rho,rht,odi,ods,xdnout)
                    IF (xsf) THEN
                       WRITE(55,*) xdnout
                    ELSE
                       IF (twodim) THEN
                          WRITE(55,'(3e15.7)') point(1:2),xdnout
                       ELSE
                          WRITE(55,'(4e15.7)') point(:),xdnout
                       ENDIF
                    ENDIF
                 ENDDO xloop
              ENDDO
           ENDDO !z-loop
           IF (xsf.AND.jsp /= jspins) CALL xsf_WRITE_newblock(55,twodim
     $          ,vec1,vec2,vec3,zero,grid)
        ENDDO !Spin-loop
        IF (xsf) THEN
           CALL xsf_WRITE_endblock(55,twodim)
        ELSE
           CLOSE(55)
        ENDIF
      ENDDO   !nplot      
      CLOSE(18)
      IF (xsf) CLOSE(55)
      RETURN
      END SUBROUTINE plotdop
!------------------------------------------
!     The old subroutine from Fleur is here
!------------------------------------------
      SUBROUTINE priv_old_plot(odi,ods,
     >                n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,ntype,
     >                lmaxd,jmtd,ntypd,natd,nmzd,jspins,neq,
     >                nq3,nvac,nmz,nmzxy,nq2,nop,nop2,volint,film,
     >                plpot,score,slice,symor,invs,invs2,z1,delz,
     >                clnu,llh,nmem,mlh,nlh,ngopr,ntypsy,jri,pos,
     >                nstr,nstr2,lmax,kv2,kv3,mrot,tau,rmsh,invtab,
     >                amat,bmat,l_noco,cdnfname)
c
      USE m_constants, ONLY : pimach
      USE m_od_types, ONLY : od_inp, od_sym
      USE m_outcdn
      USE m_loddop
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd
      INTEGER, INTENT (IN) :: lmaxd,jmtd,ntypd,natd,nmzd,jspins
      INTEGER, INTENT (IN) :: nq3,nvac,nmz,nmzxy,nq2,nop,nop2,ntype
      LOGICAL, INTENT (IN) :: plpot,symor,invs,score,slice,invs2,film
      LOGICAL, INTENT (IN) :: l_noco
      REAL,    INTENT (IN) :: z1,delz,volint
      CHARACTER*10, INTENT (IN) :: cdnfname
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: ngopr(natd),ntypsy(natd),lmax(ntypd)
      INTEGER, INTENT (IN) :: nstr(n3d),nstr2(n2d),jri(ntypd),neq(ntypd)
      INTEGER, INTENT (IN) :: kv2(2,n2d),kv3(3,n3d),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: invtab(nop)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3),pos(3,natd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),tau(3,nop)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      REAL rx,ry,s,sl,sm,su,x,xm,y,ym,tec,qint,sfp,xdnout
      INTEGER i,i1,i2,i3,ii3,imshx,imshy,ix,iy,j,jsp,na,nfile,nplo,
     +        nplot,nq,nt,nplott,jm,jspin,iter,ii1,ii2
      CHARACTER*8 dop,iop
      LOGICAL twodim
C     ..
C     .. Local Arrays ..
      COMPLEX qpw(n3d,jspd),rhtxy(nmzxyd,n2d-1,2,jspd)
      REAL rho(jmtd,0:nlhd,ntypd,jspd),rht(nmzd,2,jspd)
      REAL ptp(3),rngl(3),rngm(3),rngu(3),tl(3),tm(3)
      REAL tu(3),vx1(3),vx2(3),tl_r(3),tm_r(3),tu_r(3),rhocc(jmtd)
      REAL pt(3),a3(3)
      REAL, ALLOCATABLE :: cdn(:,:)
      CHARACTER*10, ALLOCATABLE :: plotname(:)
      CHARACTER*8  name(10)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,sqrt
C     ..
      sfp = 2.0 * sqrt( pimach() )
      a3(3) = amat(3,3)
C     ..
      OPEN (20,file=cdnfname,form='unformatted',status='old')
      CALL loddop(
     >            jspd,n3d,n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,20,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
c
      IF (score) THEN
         OPEN (17,file='cdnc',form='unformatted',status='old')
         REWIND 17
c
         DO jspin = 1,jspins
            DO nt = 1,ntype
               jm = jri(nt)
               READ (17) (rhocc(i),i=1,jm)
               DO i = 1,jri(nt)
                  rho(i,0,nt,jspin) = rho(i,0,nt,jspin) - rhocc(i)/sfp
               ENDDO
               READ (17) tec
             ENDDO
            READ (17) qint
            qpw(1,jspin) = qpw(1,jspin) - qint/volint
         ENDDO
c
         CLOSE (17)
      END IF
      OPEN (18,file='plotin')
      READ (18,FMT='(7x,l1)') twodim

      READ (18,FMT=8000) nplot
      ALLOCATE ( plotname(nplot) )
 8000 FORMAT (6x,i2)
      nplott=nplot
      if (nplot.eq.1) nplott=2
      DO 140 nplo = 1,nplot
         IF (twodim.OR.(nplo.eq.1)) THEN
           nfile = 55 + nplo
           READ (18,FMT='(a,3x,a)') plotname(nplo)
           IF (l_noco) THEN
              OPEN (nfile,file=plotname(nplo)//cdnfname,
     +                    form='formatted')
           ELSE
              OPEN (nfile,file=plotname(nplo),form='formatted')
           ENDIF
           READ (18,FMT=8010) (tu_r(i),i=1,3)
           READ (18,FMT=8010) (tm_r(i),i=1,3)
           READ (18,FMT=8010) (tl_r(i),i=1,3)
 8010      FORMAT (4f10.6)
           READ (18,FMT=8020) imshx,imshy
 8020      FORMAT (6x,i5,7x,i5)
           ALLOCATE (cdn(imshx,imshy))
         ENDIF
         IF (twodim) THEN
           DO i=1,3
             tu(i) = tu_r(i)
             tm(i) = tm_r(i)
             tl(i) = tl_r(i)
           ENDDO 
         ELSE
           DO i=1,2
             tu(i) = tu_r(i)
             tm(i) = tm_r(i)
             tl(i) = tl_r(i)
           ENDDO 
           tu(3) = tu_r(3) + a3(3)*(nplo-1)/(2*nplott-2)
           tm(3) = tm_r(3) + a3(3)*(nplo-1)/(2*nplott-2)
           tl(3) = tl_r(3) + a3(3)*(nplo-1)/(2*nplott-2)
         ENDIF

         tu(3) = tu(3)/a3(3)
         tm(3) = tm(3)/a3(3)
         tl(3) = tl(3)/a3(3)
c--->    evaluate cartesian coordinates of positions
         DO 20 i = 1,3
            su = 0.
            sm = 0.
            sl = 0.
            DO 10 j = 1,3
               su = su + amat(i,j)*tu(j)
               sm = sm + amat(i,j)*tm(j)
               sl = sl + amat(i,j)*tl(j)
   10       CONTINUE
            rngu(i) = su
            rngm(i) = sm
            rngl(i) = sl
   20    CONTINUE
         DO 30 i = 1,3
            vx1(i) = rngu(i) - rngm(i)
            vx2(i) = rngl(i) - rngm(i)
   30    CONTINUE
         rx = sqrt(vx1(1)*vx1(1)+vx1(2)*vx1(2)+vx1(3)*vx1(3))
         ry = sqrt(vx2(1)*vx2(1)+vx2(2)*vx2(2)+vx2(3)*vx2(3))
         DO 130 jsp = 1,jspins
            WRITE (16,FMT=8030) rx,ry
 8030       FORMAT (2f10.6)
            WRITE (nfile,FMT=8050) imshy,imshx,ry,rx
            WRITE (16,FMT=8050) imshx,imshy
            xm = imshx - 1
            ym = imshy - 1
            DO 120 ix = 1,imshx
               DO 110 iy = 1,imshy
                  x = ix - 1
                  y = iy - 1
                  pt(1) = rngm(1) + vx1(1)*x/xm + vx2(1)*y/ym
                  pt(2) = rngm(2) + vx1(2)*x/xm + vx2(2)*y/ym
                  pt(3) = rngm(3) + vx1(3)*x/xm + vx2(3)*y/ym
                  ii1 = 3
                  ii2 = 3
                  ii3 = 3
                  IF (film .AND. .NOT.odi%d1) ii3 = 0
                  IF (odi%d1) THEN
                     ii1 = 0 ; ii2 = 0
                  END IF
                  DO 100 i1 = -ii1,ii1
                     DO 90 i2 = -ii2,ii2
                        DO 80 i3 = -ii3,ii3
                           DO 40 i = 1,3
                              ptp(i) = pt(i) + i1*amat(i,1) +
     +                                 i2*amat(i,2) + i3*amat(i,3)
   40                      CONTINUE
                           na = 0
                           DO 70 nt = 1,ntype
                              DO 60 nq = 1,neq(nt)
                                 na = na + 1
                                 s = 0.
                                 DO 50 i = 1,3
                                    s = s + (ptp(i)-pos(i,na))**2
   50                            CONTINUE
                                 s = sqrt(s)
                                 IF (s.LT.rmsh(jri(nt),nt)) THEN
                                    CALL outcdn(
     >                  ptp,nt,na,0,1,jsp,plpot,kv2,kv3,nstr,nstr2,
     >                  n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,lmaxd,jmtd,
     >                  natd,ntypd,nmzd,nop,nop2,mrot,tau,symor,invtab,
     >                  nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,mlh,llh,
     >                  clnu,nmem,nlh,lmax,rmsh,jri,pos,ngopr,ntypsy,
     >                  amat,bmat,qpw,rhtxy,rho,rht,odi,ods,
     <                  xdnout)
                                    cdn(ix,iy) = xdnout
                                    GO TO 110
                                 END IF
   60                         CONTINUE
   70                      CONTINUE
   80                   CONTINUE
   90                CONTINUE
  100             CONTINUE
                  IF (film .AND. .NOT.odi%d1) THEN
                    IF (abs(pt(3)).GE.z1) THEN
                      CALL outcdn(
     >                  pt,0,0,1,0,jsp,plpot,kv2,kv3,nstr,nstr2,
     >                  n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,lmaxd,jmtd,
     >                  natd,ntypd,nmzd,nop,nop2,mrot,tau,symor,invtab,
     >                  nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,mlh,llh,
     >                  clnu,nmem,nlh,lmax,rmsh,jri,pos,ngopr,ntypsy,
     >                  amat,bmat,qpw,rhtxy,rho,rht,odi,ods,
     <                  xdnout)
                        cdn(ix,iy) = xdnout
                      GO TO 110
                    END IF
                  END IF
c-odim
                  IF (odi%d1) THEN
                     IF (sqrt((pt(1))**2 + (pt(2))**2).GE.z1) THEN
                      CALL outcdn(
     >                  pt,0,0,1,0,jsp,plpot,kv2,kv3,nstr,nstr2,
     >                  n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,lmaxd,jmtd,
     >                  natd,ntypd,nmzd,nop,nop2,mrot,tau,symor,invtab,
     >                  nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,mlh,llh,
     >                  clnu,nmem,nlh,lmax,rmsh,jri,pos,ngopr,ntypsy,
     >                  amat,bmat,qpw,rhtxy,rho,rht,odi,ods,
     <                  xdnout)
                        cdn(ix,iy) = xdnout
                      GO TO 110
                     END IF
                  END IF
c+odim
                  CALL outcdn(
     >                  pt,0,0,0,2,jsp,plpot,kv2,kv3,nstr,nstr2,
     >                  n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,lmaxd,jmtd,
     >                  natd,ntypd,nmzd,nop,nop2,mrot,tau,symor,invtab,
     >                  nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,mlh,llh,
     >                  clnu,nmem,nlh,lmax,rmsh,jri,pos,ngopr,ntypsy,
     >                  amat,bmat,qpw,rhtxy,rho,rht,odi,ods,
     <                  xdnout)
                  cdn(ix,iy) = xdnout
  110          CONTINUE
  120       CONTINUE
            WRITE (nfile,FMT=8040) ((cdn(ix,iy),ix=1,imshx),iy=1,imshy)
 8040       FORMAT (5e14.6)
 8050       FORMAT (2i5,2f12.6)
  130    CONTINUE
         IF (twodim) THEN
            DEALLOCATE (cdn)
            CLOSE (nfile)
         ENDIF
  140 CONTINUE
      DEALLOCATE ( plotname )
      IF (.not.twodim) CLOSE (nfile)
      CLOSE (18)
      CLOSE (20)
      END SUBROUTINE priv_old_plot

      END MODULE m_plotdop
