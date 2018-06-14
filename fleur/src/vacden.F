      MODULE m_vacden
c     *************************************************************
c     determines the 2-d star function expansion coefficients of
c     vacuum charge density. speed up by r. wu 1992
c     *************************************************************
      CONTAINS
      SUBROUTINE vacden(
     >                  nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,nvd,jspd,
     >                  nkptd,neigd,layerd,nvac,jspins,nmz,nmzxy,
     >                  mx1,mx2,z1,area,delz,omtil,bmat,recno,wkpts,
     >                  nstr2,ig2,rgphs,ig,gvac1,gvac2,nlotot,
     >                  we,ikpt,jspin,l_noco,l_ss,npotmatfile,vz,vz0,
     >                  ne,z,bbmat,bkpt,k1,k2,k3,nv,nbasfcn,nmat,
     >                  layers,izlay,vacdos,integ,evac,eig,
     >                  tworkf,nstars,nstm,starcoeff,locx,locy,qss,
     >                  odi,ods,sk2,phi2,
     X                  rhtxy,rht,qvac,qvlay,emax,n2max,
     <                  stcoeff,cdomvz,cdomvxy)

c***********************************************************************
C     ****** change vacden(....,q) for vacuum density of states shz Jan.96
C     ****** change vacden(......,stcoeff) for starcoefficients, shz. Jan.99
c     ****** changed for fleur dw
c     In non-collinear calculations the density becomes a hermitian 2x2
c     matrix. This subroutine generates this density matrix in the 
c     vacuum region. The diagonal elements of this matrix (n_11 & n_22)
c     are store in rht and rhtxy, while the real and imaginary part
c     of the off-diagonal element are store in cdomvz and cdomvxy. 
c
c     Philipp Kurz 99/07
c***********************************************************************

c******** ABBREVIATIONS ************************************************
c     qvac     : vacuum charge of each eigenstate, needed in in cdnval
c                to determine the vacuum energy parameters
c     vz       : non-warping part of the vacuum potential (matrix)
c                collinear    : 2. index = ivac (# of vaccum)
c                non-collinear: 2. index = ipot (comp. of pot. matr.)
c     rht      : non-warping part of the vacuum density matrix,
c                diagonal elements n_11 and n_22
c     rhtxy    : warping part of the vacuum density matrix,
c                diagonal elements n_11 and n_22
c     cdomvz   : non-warping part of the vacuum density matrix,
c                off-diagonal elements n_21
c     cdomvxy  : warping part of the vacuum density matrix,
c                off-diagonal elements n_21
c***********************************************************************
c
      USE m_constants, ONLY : pimach
      USE m_grdchlh
      USE m_qsf
      USE m_cylbes
      USE m_dcylbs
      USE m_dotir, ONLY : dotirp
      USE m_od_abvac
      USE m_od_types, ONLY : od_inp,od_sym

      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d
      INTEGER, INTENT (IN) :: jspin,nvac,jspins,nmz,nmzxy,mx1,mx2
      INTEGER, INTENT (IN) :: layers,ne,layerd,nvd,jspd,neigd
      INTEGER, INTENT (IN) :: nbasfcn,nmat,npotmatfile,nkptd,ikpt,n2max
      INTEGER, INTENT (IN) :: nlotot
      REAL,    INTENT (IN) :: z1,area,delz,omtil,emax
      LOGICAL, INTENT (IN) :: vacdos,integ,l_noco,l_ss
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: izlay(layerd,2),nstr2(n2d)
      INTEGER, INTENT (INOUT):: nv(jspd)
      INTEGER, INTENT (INOUT):: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: bkpt(3),bbmat(3,3),bmat(3,3)
      REAL,    INTENT (IN) :: evac(2,jspd)
      REAL,    INTENT (IN) :: sk2(n2d),phi2(n2d)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, INTENT (IN):: z(nbasfcn,neigd)
#else
      REAL,    INTENT (IN):: z(nbasfcn,neigd)
#endif
      COMPLEX, INTENT (INOUT):: rhtxy(nmzxyd,odi%n2d-1,2,jspd)
      REAL,    INTENT (INOUT):: rht(nmzd,2,jspd)
      REAL,    INTENT (OUT)  :: qvlay(neigd,layerd,2,nkptd,jspd)
      REAL qvac(neigd,2,nkptd,jspd),we(neigd),vz(nmzd,2),vz0(2)
      COMPLEX, INTENT (INOUT):: cdomvz(nmzd,2)
      COMPLEX, INTENT (INOUT):: cdomvxy(nmzxyd,odi%n2d-1,2)
C
C     STM-Arguments
      REAL,    INTENT (IN) :: tworkf,wkpts
      REAL,    INTENT (IN) :: eig(neigd),locx(2),locy(2)
      REAL,    INTENT (IN) :: qss(3)
      INTEGER, INTENT (IN) :: nstars,nstm,recno
      INTEGER, INTENT (IN) :: gvac1(nv2d),gvac2(nv2d)
      LOGICAL, INTENT (IN) :: starcoeff
      COMPLEX, INTENT (OUT):: stcoeff(nstars,neigd,layerd,2)
c
c     local STM variables
      REAL    tpi
      INTEGER nv2(jspd)
      INTEGER kvac1(nv2d,jspd),kvac2(nv2d,jspd),map2(nvd,jspd)
      INTEGER kvac3(nv2d,jspd),map1(nvd,jspd)
      INTEGER mapg2k(nv2d)
C     .. Local Scalars ..
      COMPLEX aa,ab,av,ba,bb,bv,t1,aae,bbe,abe,bae,aaee,bbee,abee,baee,
     +     factorx,factory,c_1,aa_1,ab_1,ba_1,bb_1,ic,av_1,bv_1,d 

      REAL arg,const,ddui,dduj,dduei,dduej,eps,ev,evacp,phs,phsp,qout,
     +     scale,sign,uei,uej,ui,uj,wronk,zks,result(1),ui2,uei2,
     +     k_diff,k_d1,k_d2,ui_1,uei_1,uj_1,uej_1,wronk_1

      INTEGER i,ii,i1,i2,i3,ig3,ig3p,ik,ind2,ind2p,istar,m,m1,m3,
     +        ivac,j,jj,jz,k,l,ll,l1,n,n2,ispin,kspin,jsp_start,jsp_end,
     +        ipot,ie,imz,isp_start,isp_end,
     +        ind1,ind1p,irec2,irec3
C
C     .. Local Arrays ..
      REAL qssbti(3,2),qssbtii
      REAL bess(-odi%mb:odi%mb),dbss(-odi%mb:odi%mb)
      COMPLEX, ALLOCATABLE :: ac(:,:,:),bc(:,:,:)
      REAL,    ALLOCATABLE :: dt(:),dte(:),du(:),ddu(:,:),due(:)
      REAL,    ALLOCATABLE :: ddue(:,:),t(:),te(:),tei(:,:),dummy(:)
      REAL,    ALLOCATABLE :: u(:,:,:),ue(:,:,:),v(:),yy(:)
c-odim
      COMPLEX, ALLOCATABLE :: ac_1(:,:,:,:),bc_1(:,:,:,:)
      REAL,    ALLOCATABLE :: dt_1(:,:),dte_1(:,:)
      REAL,    ALLOCATABLE :: du_1(:,:),ddu_1(:,:,:),due_1(:,:)
      REAL,    ALLOCATABLE :: ddue_1(:,:,:),t_1(:,:)
      REAL,    ALLOCATABLE :: tei_1(:,:,:),te_1(:,:)
      REAL,    ALLOCATABLE :: u_1(:,:,:,:),ue_1(:,:,:,:)
c+odim
C     ..
C     .. External Subroutines ..
      EXTERNAL vacudz,vacuz
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,conjg,exp,iabs,real,sqrt
C     ..

c     *******************************************************************************
c

c
c    layers: no. of layers to be calculated (in vertical direction with z-values as given by izlay)
c    izlay : defines vertical position of layers in delz (=0.1 a.u.) units from begining of vacuum region
c    vacdos: =T: calculate vacuum dos in layers as given by the above
c    integ : =T: integrate in vertical position between izlay(layer,1)..izlay(layer,2)
c    nstm  : 0: s-Tip, 1: p_z-Tip, 2: d_z^2-Tip (following Chen's derivative rule) ->rhzgrd.f is used 
c                 to calculate derivatives 
c    tworkf: Workfunction of Tip (in hartree units) is needed for d_z^2-Orbital)
c    starcoeff: =T: star coefficients are calculated at values of izlay for 0th (=q) to nstars-1. star 
c                (stcoeff(1..nstars-1))
c    nstars: number of star functions to be used (0th star is given by value of q=charge integrated in 2D) 
c
c    further possibility: (readin of locx, locy has to be implemented in flapw7.f or they have to be set explicitly)
c
c     locx and locy can be used to calculate local DOS at a certain vertical position z (or integrated in z)
c     within a restricted area of the 2D unit cell, the corners of this area is given by locx and locy
c     they are defined in internal coordinates, i.e. \vec{r}_1=locx(1)*\vec{a}_1+locy(1)*\vec{a}_2
c                                                    \vec{r}_2=locx(2)*\vec{a}_1+locy(2)*\vec{a}_2
c                 \vec{a}_1,2 are the 2D lattice vectors           
c  
c     **************************************************************************************************
c 
      ALLOCATE ( ac(nv2d,neigd,jspd),bc(nv2d,neigd,jspd),dt(nv2d),
     +           dte(nv2d),du(nmzd),ddu(nmzd,nv2d),due(nmzd),
     +           ddue(nmzd,nv2d),t(nv2d),te(nv2d),
     +           tei(nv2d,jspd),u(nmzd,nv2d,jspd),ue(nmzd,nv2d,jspd),
     +           v(3),yy(nmzd),
     +                  ac_1(nv2d,-odi%mb:odi%mb,neigd,jspd),
     +                  bc_1(nv2d,-odi%mb:odi%mb,neigd,jspd),
     +                  dt_1(nv2d,-odi%mb:odi%mb),
     +                 dte_1(nv2d,-odi%mb:odi%mb),
     +                  du_1(nmzd,-odi%mb:odi%mb),
     +            ddu_1(nmzd,nv2d,-odi%mb:odi%mb),
     +                 due_1(nmzd,-odi%mb:odi%mb),
     +           ddue_1(nmzd,nv2d,-odi%mb:odi%mb),
     +                   t_1(nv2d,-odi%mb:odi%mb),
     +                  te_1(nv2d,-odi%mb:odi%mb),
     +                 tei_1(nv2d,-odi%mb:odi%mb,jspd),
     +              u_1(nmzd,nv2d,-odi%mb:odi%mb,jspd),
     +             ue_1(nmzd,nv2d,-odi%mb:odi%mb,jspd) )
c

      eps=0.01
      tpi = 2 * pimach()
      ic = cmplx(0.,1.)
c    ------------------
c     WRITE (16,'(a,i2)') 'nstars=',nstars
      DO n = 1,neigd
        DO ll = 1,layerd  
          DO i = 1,2
            IF (starcoeff .and. vacdos) THEN
              DO istar=1,nstars-1
                stcoeff(istar,n,ll,i) = cmplx(0.0,0.0)
              ENDDO
            END IF
          ENDDO
        ENDDO
      ENDDO

c     -----> set up mapping arrays
      IF (l_ss) THEN
         jsp_start = 1
         jsp_end   = 2
      ELSE
         jsp_start = jspin
         jsp_end   = jspin
      ENDIF
      DO ispin = jsp_start,jsp_end
         IF (odi%d1) THEN
            n2 = 0
            DO 35 k = 1,nv(ispin)
               DO 45 j = 1,n2
                  IF (k3(k,ispin).EQ.kvac3(j,ispin)) THEN
                     map1(k,ispin) = j
                     GOTO 35
                  END IF
 45            CONTINUE
               n2 = n2 + 1
               IF (n2.GT.nv2d) STOP 'vacden0'
               kvac3(n2,ispin) =  k3(k,ispin)
               map1(k,ispin) = n2
 35         CONTINUE
            nv2(ispin) = n2
         ELSE
            n2 = 0
            DO 40 k = 1,nv(ispin)
               DO 30 j = 1,n2
                  IF ( k1(k,ispin).EQ.kvac1(j,ispin) .AND.
     +                 k2(k,ispin).EQ.kvac2(j,ispin) ) THEN
                     map2(k,ispin) = j
                     GO TO 40
                  END IF
 30            CONTINUE
               n2 = n2 + 1
               IF (n2.GT.nv2d) STOP 'vacden0'
               kvac1(n2,ispin) = k1(k,ispin)
               kvac2(n2,ispin) = k2(k,ispin)
               map2(k,ispin) = n2
 40         CONTINUE
            nv2(ispin) = n2
         END IF
      ENDDO
      IF ( l_noco .AND. (.NOT. l_ss) ) THEN
         nv(2)  = nv(1)
         nv2(2) = nv2(1)
         DO k = 1,nv2(1)
            kvac1(k,2) = kvac1(k,1)
            kvac2(k,2) = kvac2(k,1)
            if(odi%d1) kvac3(k,2) =  kvac3(k,1)
         ENDDO
         DO k = 1,nv(1)
            k3(k,2) = k3(k,1)
            map2(k,2) = map2(k,1)
            if(odi%d1)then
               k1(k,2) = k1(k,1)
               k2(k,2) = k2(k,1)
               map1(k,2) = map1(k,1)
            endif
         ENDDO
      ENDIF

c+dw
c    if tunneling current should be calculated we need to write out 
c     info on electronic structure: --> mapping from kvac to gvac by mapg2k
c                                             shz, Jan.99
      IF (nstm.EQ.3) THEN
         DO j=1, n2max
            mapg2k(j)=j
            DO i=1, nv2(jspin)
               IF ( kvac1(i,jspin).EQ.gvac1(j) .AND.
     +              kvac2(i,jspin).EQ.gvac2(j) ) 
     +              mapg2k(j)=i
            END DO
         END DO
      END IF
c
c-dw
      IF (l_noco) THEN
         OPEN (npotmatfile,FILE='potmat',FORM='unformatted',
     +         STATUS='old')
c--->    skip the four components of the interstitial potential matrix
         DO ipot = 1,3
            READ (npotmatfile)
         ENDDO
      ENDIF

      wronk = 2.0
      const = 1.0 / ( sqrt(omtil)*wronk )
c-odim
      ac_1(:,:,:,:) = cmplx(0.0,0.0)
      bc_1(:,:,:,:) = cmplx(0.0,0.0)
c+odim
      DO 140 ivac = 1,nvac
         DO ispin = 1,jspins
            DO 60 n = 1,neigd
               DO 50 i = 1,nv2d
                  ac(i,n,ispin) = cmplx(0.0,0.0)
                  bc(i,n,ispin) = cmplx(0.0,0.0)
 50            CONTINUE
 60         CONTINUE
         ENDDO
         sign = 3. - 2.*ivac
         IF (l_noco) THEN
c--->       read the non-warping potential matrix, it is needed to 
c--->       calculate the vacuum basis functions
            READ (npotmatfile)((vz(imz,ipot),imz=1,nmzd),ipot=1,2)
c--->       skip the warping potential matrix
            DO ipot = 1,3
               READ (npotmatfile)
            ENDDO
         ENDIF
         IF (l_noco) THEN
c--->    In a non-collinear calculation vacden is only called once.
c--->    Thus, the vaccum wavefunctions and the A- and B-coeff. (ac bc)
c--->    have to be calculated for both spins on that call.
c--->       setup the spin-spiral q-vector
            qssbti(1,1) = - qss(1)/2
            qssbti(2,1) = - qss(2)/2
            qssbti(1,2) = + qss(1)/2
            qssbti(2,2) = + qss(2)/2
            qssbti(3,1) = - qss(3)/2
            qssbti(3,2) = + qss(3)/2
            DO ispin = 1,jspins
c     -----> set up vacuum wave functions
            IF (odi%d1) THEN
              evacp = evac(ivac,ispin)
              CALL od_abvac(
     >           z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,
     >           ig,odi%ig,tpi,qssbti(3,ispin),
     >           nmzxy,nmz,delz,ig2,odi%n2d,
     >           bbmat,wronk,evacp,bkpt,odi%M,odi%mb,
     >           vz(1,ispin),kvac3(1,ispin),nv2(ispin),
     <           t_1(1,-odi%mb),dt_1(1,-odi%mb),u_1(1,1,-odi%mb,ispin),
     <           te_1(1,-odi%mb),dte_1(1,-odi%mb),
     <           tei_1(1,-odi%mb,ispin),
     <           ue_1(1,1,-odi%mb,ispin))
              DO k = 1,nv(ispin)
                 kspin = (nv(1)+nlotot)*(ispin-1) + k
                 l = map1(k,ispin) 
                 irec3 = ig(k1(k,ispin),k2(k,ispin),k3(k,ispin))
                 IF (irec3.ne.0) THEN
                    irec2 = ig2(irec3)
                    zks = sk2(irec2)*z1
                    arg = phi2(irec2)
                    CALL cylbes(odi%mb,zks,bess)
                    CALL dcylbs(odi%mb,zks,bess,dbss)
                    DO m = -odi%mb,odi%mb
                       wronk_1 = t_1(l,m)*dte_1(l,m) -
     -                      te_1(l,m)*dt_1(l,m)
                       av_1 = exp(-cmplx(0.0,m*arg))*(ic**m)*
     *                      cmplx(dte_1(l,m)*bess(m) -
     +                      te_1(l,m)*sk2(irec2)*dbss(m),0.0)/
     /                      ((wronk_1)*sqrt(omtil))
                       bv_1 = exp(-cmplx(0.0,m*arg))*(ic**m)*
     *                      cmplx(-dt_1(l,m)*bess(m) +
     -                      t_1(l,m)*sk2(irec2)*dbss(m),0.0)/
     /                      ((wronk_1)*sqrt(omtil))
                       DO n = 1,ne
                           ac_1(l,m,n,ispin) = ac_1(l,m,n,ispin) +
     +                        z(kspin,n)*av_1
c     +                         conjg(z(k,n))*av_1
                           bc_1(l,m,n,ispin) = bc_1(l,m,n,ispin) +
     +                        z(kspin,n)*bv_1
c     +                         conjg(z(k,n))*bv_1
                       END DO
                    END DO      ! -mb:mb
                 END IF
              END DO
            ELSE ! 1-dimensional
            vz0(ispin) = vz(nmz,ispin)
            evacp = evac(ivac,ispin)
            DO ik = 1,nv2(ispin)
               v(1) = bkpt(1) + kvac1(ik,ispin) + qssbti(1,ispin)
               v(2) = bkpt(2) + kvac2(ik,ispin) + qssbti(2,ispin)
               v(3) = 0.
               ev = evacp - 0.5*dotirp(v,v,bbmat)
               CALL vacuz(ev,vz(1,ispin),vz0(ispin),nmz,delz,t(ik),
     +              dt(ik),u(1,ik,ispin))
               CALL vacudz(ev,vz(1,ispin),vz0(ispin),nmz,delz,te(ik),
     +              dte(ik),tei(ik,ispin),ue(1,ik,ispin),dt(ik),
     +              u(1,ik,ispin))
               scale = wronk/ (te(ik)*dt(ik)-dte(ik)*t(ik))
               te(ik) = scale*te(ik)
               dte(ik) = scale*dte(ik)
               tei(ik,ispin) = scale*tei(ik,ispin)
               DO j = 1,nmz
                  ue(j,ik,ispin) = scale*ue(j,ik,ispin)
               ENDDO
            ENDDO
c     -----> construct a and b coefficients
            DO k = 1,nv(ispin)
c--->          the coefficients of the spin-down basis functions are
c--->          stored in the second half of the eigenvector
               kspin = (nv(1)+nlotot)*(ispin-1) + k
               l = map2(k,ispin)
               zks = k3(k,ispin)*bmat(3,3)*sign
               arg = zks*z1
               c_1 = cmplx(cos(arg),sin(arg)) * const
               av = -c_1 * cmplx( dte(l),zks*te(l) ) 
               bv =  c_1 * cmplx(  dt(l),zks* t(l) ) 
c     -----> loop over basis functions
               DO n = 1,ne
                  ac(l,n,ispin) = ac(l,n,ispin) + z(kspin,n)*av
                  bc(l,n,ispin) = bc(l,n,ispin) + z(kspin,n)*bv
               ENDDO
            ENDDO
c--->       end of spin loop
            ENDIF
c---Y       end of geometry 1d/film
            ENDDO
c--->       output for testing
c            DO k = 1,10
c               DO n = 1,5
c                  DO ispin = 1,jspins
c                     write(*,9000)k,n,ispin,ac(k,n,ispin),bc(k,n,ispin)
c                  ENDDO
c               ENDDO
c            ENDDO
c 9000       FORMAT('k=',i3,' ie=',i3,' isp=',i3,
c     +             ' ac= (',e12.6,',',e12.6,')',
c     +             ' bc= (',e12.6,',',e12.6,')')
         ELSE
c     -----> set up vacuum wave functions
           IF (odi%d1) THEN
              qssbtii = 0.
              evacp = evac(ivac,jspin)
              CALL od_abvac(
     >           z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,
     >           ig,odi%ig,tpi,qssbtii,
     >           nmzxy,nmz,delz,ig2,odi%n2d,
     >           bbmat,wronk,evacp,bkpt,odi%M,odi%mb,
     >           vz(1,ivac),kvac3(1,jspin),nv2(jspin),
     <           t_1(1,-odi%mb),dt_1(1,-odi%mb),u_1(1,1,-odi%mb,jspin),
     <           te_1(1,-odi%mb),dte_1(1,-odi%mb),
     <           tei_1(1,-odi%mb,jspin),
     <           ue_1(1,1,-odi%mb,jspin))
              DO k = 1,nv(jspin)
                 l = map1(k,jspin)
                 irec3 = ig(k1(k,jspin),k2(k,jspin),k3(k,jspin))
                 IF (irec3.ne.0) THEN
                    irec2 = ig2(irec3)
                    zks = sk2(irec2)*z1
                    arg = phi2(irec2)
                    CALL cylbes(odi%mb,zks,bess)
                    CALL dcylbs(odi%mb,zks,bess,dbss)
                    DO m = -odi%mb,odi%mb
                       wronk_1 = t_1(l,m)*dte_1(l,m) -
     -                      te_1(l,m)*dt_1(l,m)
                       av_1 = exp(-cmplx(0.0,m*arg))*(ic**m)*
     *                      cmplx(dte_1(l,m)*bess(m) -
     +                      te_1(l,m)*sk2(irec2)*dbss(m),0.0)/
     /                      ((wronk_1)*sqrt(omtil))
                       bv_1 = exp(-cmplx(0.0,m*arg))*(ic**m)*
     *                      cmplx(-dt_1(l,m)*bess(m) +
     -                      t_1(l,m)*sk2(irec2)*dbss(m),0.0)/
     /                      ((wronk_1)*sqrt(omtil))
                       DO n = 1,ne
                           ac_1(l,m,n,jspin) = ac_1(l,m,n,jspin) +
     +                        z(k,n)*av_1
c     +                         conjg(z(k,n))*av_1
                           bc_1(l,m,n,jspin) = bc_1(l,m,n,jspin) +
     +                        z(k,n)*bv_1
c     +                         conjg(z(k,n))*bv_1
                       END DO
                    END DO      ! -mb:mb
                 END IF
              END DO         ! k = 1,nv
           ELSE     !odi%d1
            evacp = evac(ivac,jspin)
            DO ik = 1,nv2(jspin)
               v(1) = bkpt(1) + kvac1(ik,jspin)
               v(2) = bkpt(2) + kvac2(ik,jspin)
               v(3) = 0.
               ev = evacp - 0.5*dotirp(v,v,bbmat)
               CALL vacuz(ev,vz(1,ivac),vz0(ivac),nmz,delz,t(ik),dt(ik),
     +              u(1,ik,jspin))
               CALL vacudz(ev,vz(1,ivac),vz0(ivac),nmz,delz,te(ik),
     +              dte(ik),tei(ik,jspin),ue(1,ik,jspin),dt(ik),
     +              u(1,ik,jspin))
               scale = wronk/ (te(ik)*dt(ik)-dte(ik)*t(ik))
               te(ik) = scale*te(ik)
               dte(ik) = scale*dte(ik)
               tei(ik,jspin) = scale*tei(ik,jspin)
               DO j = 1,nmz
                  ue(j,ik,jspin) = scale*ue(j,ik,jspin)
               ENDDO
            ENDDO
c     -----> construct a and b coefficients
            DO k = 1,nv(jspin)
               l = map2(k,jspin)
               zks = k3(k,jspin)*bmat(3,3)*sign
               arg = zks*z1
               c_1 = cmplx(cos(arg),sin(arg)) * const
               av = -c_1 * cmplx( dte(l),zks*te(l) ) 
               bv =  c_1 * cmplx(  dt(l),zks* t(l) ) 
c     -----> loop over basis functions
               DO n = 1,ne
                  ac(l,n,jspin) = ac(l,n,jspin) + z(k,n)*av
                  bc(l,n,jspin) = bc(l,n,jspin) + z(k,n)*bv
               ENDDO
            ENDDO
           END IF ! D1
         ENDIF
c
c   ----> calculate first and second derivative of u,ue 
c        in order to simulate p_z or d_z^2 Tip in Chen's model , shz. 97
c 
         IF (nstm.gt.0) THEN
            DO 89 ik = 1,nv2(jspin)
!               CALL rhzgrd(nmz,delz,u(1,ik,jspin),4,du,ddu(1,ik))
!               CALL rhzgrd(nmz,delz,ue(1,ik,jspin),4,due,ddue(1,ik))   

               ALLOCATE ( dummy(nmz) )
               CALL grdchlh(
     >                      0,1,nmz,delz,dummy,u(1,ik,jspin),4,
     <                      du,ddu(1,ik))
               CALL grdchlh(
     >                      0,1,nmz,delz,dummy,ue(1,ik,jspin),4,
     <                      due,ddue(1,ik))
               DEALLOCATE ( dummy )

               IF (nstm.eq.1) THEN 
                  DO 87 j= 1,nmz  
                     u(j,ik,jspin) = du(j)
                     ue(j,ik,jspin) = due(j)
 87               CONTINUE  
               END IF
 89         CONTINUE
         END IF  

c+dw

c
c       --> to calculate Tunneling Current between two systems
c           within Bardeens Approach one needs ac(l,n), bc(l,n);
c           they are written to the file vacwave 
c                           IF nstm=3
c                              tworkf is then the fermi energy (in hartree)
c
         IF (nstm.EQ.3) THEN
#ifdef CPP_MPI
            WRITE (*,*) 'Rec',recno,'was written'
            WRITE (88,rec=recno) wkpts,ne,bkpt(1),bkpt(2),
     +                   (eig(n),n=1,neigd),
     +                   ((ac(mapg2k(i),n,jspin),i=1,n2max),n=1,neigd), 
     +                   ((bc(mapg2k(i),n,jspin),i=1,n2max),n=1,neigd)
#else
            i=0
            DO n = 1, ne
               IF (ABS(eig(n)-tworkf).LE.emax) i=i+1
            END DO
            WRITE (87,FMT=990) bkpt(1), bkpt(2), i, n2max
            DO n = 1, ne
               IF (ABS(eig(n)-tworkf).LE.emax) THEN
                  WRITE (87,FMT=1000) eig(n)
                  DO j=1,n2max
                     WRITE (87,FMT=1010) ac(mapg2k(j),n,jspin),
     +                                   bc(mapg2k(j),n,jspin)
                  END DO
               END IF
            END DO
#endif



         END IF
 990     FORMAT(2(f8.4,1x),i3,1x,i3)
 1000    FORMAT(e10.4)
 1010    FORMAT(2(2e20.8,1x))
c
c        ------------------------------------------------------------
c-dw
c
c---->   non-warping part of the density (g=0 star)
c
         IF (nstm.EQ.2) THEN
c
c  ----> d_z^2-Tip needs: |d^2(psi)/dz^2 - kappa^2/3 psi|^2
c
            DO l = 1,nv2(jspin)
               aa = 0.0
               bb = 0.0
               ba = 0.0
               ab = 0.0
               DO n = 1,ne
                  aa = aa + we(n)*conjg(ac(l,n,jspin))*ac(l,n,jspin)
                  bb = bb + we(n)*conjg(bc(l,n,jspin))*bc(l,n,jspin)
                  ab = ab + we(n)*conjg(ac(l,n,jspin))*bc(l,n,jspin)
                  ba = ba + we(n)*conjg(bc(l,n,jspin))*ac(l,n,jspin)
                  qout = real(conjg(ac(l,n,jspin))*ac(l,n,jspin)+
     +                 tei(l,jspin)*conjg(bc(l,n,jspin))*bc(l,n,jspin))
                  qvac(n,ivac,ikpt,jspin) = qvac(n,ivac,ikpt,jspin) + 
     +                                                      qout*area
               END DO
               aae=-aa*tworkf*2/3
               bbe=-bb*tworkf*2/3
               abe=-ab*tworkf*2/3
               bae=-ba*tworkf*2/3
               aaee=aa*tworkf*tworkf*4/9
               bbee=bb*tworkf*tworkf*4/9
               abee=ab*tworkf*tworkf*4/9
               baee=ba*tworkf*tworkf*4/9  
               DO 105 jz = 1,nmz
                  ui = u(jz,l,jspin)
                  uei = ue(jz,l,jspin)
                  ddui = ddu(jz,l) 
                  dduei = ddue(jz,l)
                  rht(jz,ivac,jspin) = rht(jz,ivac,jspin) +
     +                 real(aaee*ui*ui+bbee*uei*uei+
     +                 (abee+baee)*ui*uei+aa*ddui*ddui+
     +                 bb*dduei*dduei+(ab+ba)*ddui*dduei+
     +                 2*aae*ui*ddui+2*bbe*uei*dduei+
     +                 (abe+bae)*(ui*dduei+uei*ddui))

 105           CONTINUE   
            END DO
c
c    -----> s-Tip: |psi|^2 and p-Tip: |d(psi)/dz|^2 
c
         ELSE
            IF (l_noco) THEN
c--->          diagonal elements of the density matrix, n_11 and n_22
c--->          the non-warping part of n_21 is calculated together with
c--->          the warping part of n_21
               DO ispin = 1,jspins
               IF (odi%d1) THEN
                  DO 16 l = 1,nv2(ispin)
                   DO 17 m = -odi%mb,odi%mb
                    aa = cmplx(0.0,0.0)
                    bb = cmplx(0.0,0.0)
                    ba = cmplx(0.0,0.0)
                    ab = cmplx(0.0,0.0)
                    DO 18 n = 1,ne
                     aa = aa + we(n)*
     *                   conjg(ac_1(l,m,n,ispin))*ac_1(l,m,n,ispin)
                     bb = bb + we(n)*
     *                   conjg(bc_1(l,m,n,ispin))*bc_1(l,m,n,ispin)
                     ab = ab + we(n)*
     *                   conjg(ac_1(l,m,n,ispin))*bc_1(l,m,n,ispin)
                     ba = ba + we(n)*
     *                   conjg(bc_1(l,m,n,ispin))*ac_1(l,m,n,ispin)
                     qout = real(conjg(ac_1(l,m,n,ispin))*
     *                                     ac_1(l,m,n,ispin) +
     +                      tei_1(l,m,ispin)*conjg(bc_1(l,m,n,ispin))*
     *                                     bc_1(l,m,n,ispin))
                     qvac(n,ivac,ikpt,ispin) = qvac(n,ivac,ikpt,ispin)+
     +                      qout*area
 18                 END DO
                    DO 19 jz = 1,nmz
                     ui = u_1(jz,l,m,ispin)
                     uei = ue_1(jz,l,m,ispin)
                     rht(jz,ivac,ispin) = rht(jz,ivac,ispin) +
     +                    real(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)
 19                 END DO
 17                END DO
 16               END DO
               ELSE 
                  DO l = 1,nv2(ispin)
                     aa = 0.0
                     bb = 0.0
                     ba = 0.0
                     ab = 0.0
                     DO n = 1,ne
                        aa=aa + we(n)*conjg(ac(l,n,ispin))*ac(l,n,ispin)
                        bb=bb + we(n)*conjg(bc(l,n,ispin))*bc(l,n,ispin)
                        ab=ab + we(n)*conjg(ac(l,n,ispin))*bc(l,n,ispin)
                        ba=ba + we(n)*conjg(bc(l,n,ispin))*ac(l,n,ispin)
                        qout = real(conjg(ac(l,n,ispin))*ac(l,n,ispin)+
     +                  tei(l,ispin)*conjg(bc(l,n,ispin))*bc(l,n,ispin))
                        qvac(n,ivac,ikpt,ispin) = 
     +                               qvac(n,ivac,ikpt,ispin) + qout*area
                     END DO
                     DO jz = 1,nmz
                        ui = u(jz,l,ispin)
                        uei = ue(jz,l,ispin)
                        rht(jz,ivac,ispin) = rht(jz,ivac,ispin) +
     +                       real(aa*ui*ui+bb*uei*uei+
     +                       (ab+ba)*ui*uei)
                     ENDDO
                  ENDDO
               END IF ! one-dimensional
               ENDDO
            ELSE
              IF (odi%d1) THEN
               DO 106 l = 1,nv2(jspin)
                DO 107 m = -odi%mb,odi%mb
                 aa = cmplx(0.0,0.0)
                 bb = cmplx(0.0,0.0)
                 ba = cmplx(0.0,0.0)
                 ab = cmplx(0.0,0.0)
                 DO 108 n = 1,ne
                  aa = aa + we(n)*
     *                conjg(ac_1(l,m,n,jspin))*ac_1(l,m,n,jspin)
                  bb = bb + we(n)*
     *                conjg(bc_1(l,m,n,jspin))*bc_1(l,m,n,jspin)
                  ab = ab + we(n)*
     *                conjg(ac_1(l,m,n,jspin))*bc_1(l,m,n,jspin)
                  ba = ba + we(n)*
     *                conjg(bc_1(l,m,n,jspin))*ac_1(l,m,n,jspin)
                  qout = real(conjg(ac_1(l,m,n,jspin))*
     *                                    ac_1(l,m,n,jspin) +
     +                    tei_1(l,m,jspin)*conjg(bc_1(l,m,n,jspin))*
     *                                    bc_1(l,m,n,jspin))
                  qvac(n,ivac,ikpt,jspin) = qvac(n,ivac,ikpt,jspin)+
     +                      qout*area
 108             END DO
                 DO 109 jz = 1,nmz
                  ui = u_1(jz,l,m,jspin)
                  uei = ue_1(jz,l,m,jspin)
                  rht(jz,ivac,jspin) = rht(jz,ivac,jspin) +
     +                 real(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)
 109             END DO
 107            END DO
 106           END DO
              ELSE          ! D1
               DO l = 1,nv2(jspin)
                  aa = cmplx(0.0,0.0)
                  bb = cmplx(0.0,0.0)
                  ba = cmplx(0.0,0.0)
                  ab = cmplx(0.0,0.0)
                  DO n = 1,ne
                     aa = aa + we(n)*conjg(ac(l,n,jspin))*ac(l,n,jspin)
                     bb = bb + we(n)*conjg(bc(l,n,jspin))*bc(l,n,jspin)
                     ab = ab + we(n)*conjg(ac(l,n,jspin))*bc(l,n,jspin)
                     ba = ba + we(n)*conjg(bc(l,n,jspin))*ac(l,n,jspin)
                     qout = real(conjg(ac(l,n,jspin))*ac(l,n,jspin)+
     +                  tei(l,jspin)*conjg(bc(l,n,jspin))*bc(l,n,jspin))
                     qvac(n,ivac,ikpt,jspin) = qvac(n,ivac,ikpt,jspin) + 
     +                                                         qout*area
                  END DO
                  DO 110 jz = 1,nmz
                     ui = u(jz,l,jspin)
                     uei = ue(jz,l,jspin)
                     rht(jz,ivac,jspin) = rht(jz,ivac,jspin) +
     +                    real(aa*ui*ui+bb*uei*uei+
     +                    (ab+ba)*ui*uei)
 110              CONTINUE
               END DO
              END IF ! D1
            ENDIF
         END IF 
c
C     ****************** change for vacuum density of states shz Jan.96 ***
c
          IF (vacdos) THEN                       
c
c  ----> d_z^2-Tip needs: |d^2(psi)/dz^2 - kappa^2/3 psi|^2
c
             IF (nstm.eq.2) THEN
                DO l=1,nv2(jspin)
                   DO n = 1,ne
                      aa = conjg(ac(l,n,jspin))*ac(l,n,jspin)
                      bb = conjg(bc(l,n,jspin))*bc(l,n,jspin)
                      ab = conjg(ac(l,n,jspin))*bc(l,n,jspin)
                      ba = conjg(bc(l,n,jspin))*ac(l,n,jspin)
                      aae = -tworkf*aa*2/3
                      bbe = -tworkf*bb*2/3
                      abe = -tworkf*ab*2/3
                      bae = -tworkf*ba*2/3
                      aaee = aa*tworkf*tworkf*4/9
                      bbee = bb*tworkf*tworkf*4/9
                      abee = ab*tworkf*tworkf*4/9
                      baee = ba*tworkf*tworkf*4/9
                      DO jj = 1,layers
c
c     ----> either integrated LDOS(z1,z2) or LDOS(z1)
c     
                         IF (integ) THEN
                            ll = 1
                            DO ii = izlay(jj,1),izlay(jj,2)
                               ui = u(ii,l,jspin)
                               uei = ue(ii,l,jspin)
                               ddui = ddu(ii,l)
                               dduei = ddue(ii,l)
                               yy(ll) = real(aaee*ui*ui+bbee*uei*uei+
     +                             (abee+baee)*ui*uei+aa*ddui*ddui+bb*d
     +                             duei*dduei+(ab+ba)*ddui*dduei+2*aae*
     +                             ui*ddui+2*bbe*uei*dduei+
     +                             (abe+bae)*(ui*dduei+uei*ddui))*area
                               ll = ll+1
                            END DO
                            CALL qsf(delz,yy,result,ll-1,0)
                            qvlay(n,jj,ivac,ikpt,jspin) =
     +                           qvlay(n,jj,ivac,ikpt,jspin) + result(1)
                         ELSE
                            ui = u(izlay(jj,1),l,jspin)       
                            uei = ue(izlay(jj,1),l,jspin)
                            ddui = ddu(izlay(jj,1),l)
                            dduei = ddue(izlay(jj,1),l)
                            yy(1) = real(aaee*ui*ui+bbee*uei*uei+
     +                           (abee+baee)*ui*uei+aa*ddui*ddui+
     +                           bb*dduei*dduei+(ab+ba)*ddui*dduei+
     +                           2*aae*ui*ddui+2*bbe*uei*dduei+
     +                           (abe+bae)*(ui*dduei+uei*ddui))
                            qvlay(n,jj,ivac,ikpt,jspin) = 
     +                             qvlay(n,jj,ivac,ikpt,jspin) +yy (1)
                         END IF
                      END DO
                   END DO
                END DO
c     
c     ----> s-Tip = calculate LDOS and(!) p_z-Tip (since u->du/dz, ue->due/dz)
c     
             ELSE 
                IF (abs(locx(1)-locx(2)).le.eps) THEN
c     
c     ----> integrated over 2D-unit cell
c
                  IF (l_noco) THEN
                     isp_start = 1
                     isp_end   = jspins
                  ELSE
                     isp_start = jspin
                     isp_end   = jspin
                  ENDIF
                  DO ispin = isp_start, isp_end
                   DO l=1,nv2(ispin)
                      DO n = 1,ne
                         aa = conjg(ac(l,n,ispin))*ac(l,n,ispin)
                         bb = conjg(bc(l,n,ispin))*bc(l,n,ispin)
                         ab = conjg(ac(l,n,ispin))*bc(l,n,ispin)
                         ba = conjg(bc(l,n,ispin))*ac(l,n,ispin) 
                         DO jj = 1,layers
c     
c     ---> either integrated (z1,z2) or slice (z1)
c     
                            IF (integ) THEN
                               ll = 1
                               DO ii = izlay(jj,1),izlay(jj,2)
                                  ui = u(ii,l,ispin)
                                  uei = ue(ii,l,ispin)
                                  yy(ll) = real(aa*ui*ui+bb*uei*uei+
     +                                 (ab+ba)*ui*uei)
                                  ll = ll+1
                               END DO
                               CALL qsf(delz,yy,result,ll-1,0)
                               qvlay(n,jj,ivac,ikpt,ispin) = 
     +                           qvlay(n,jj,ivac,ikpt,ispin) + result(1)
                            ELSE
                               ui = u(izlay(jj,1),l,ispin)       
                               uei = ue(izlay(jj,1),l,ispin)
                               qvlay(n,jj,ivac,ikpt,ispin) = 
     +                              qvlay(n,jj,ivac,ikpt,ispin) + real(
     +                              aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)
                               
                            END IF
                         END DO
                      END DO
                   END DO
                  ENDDO
                ELSE
c     
c     ----> if LDOS should be calculated over restricted area of the 2D-unit cell
c     lower left corner: (locx(1), locy(1))   }  in internal
c     upper right corner: (locx(2), locy(2))  }  coordinates
c     
                   DO l=1, nv2(jspin)
                      DO l1=1, nv2(jspin)
                         IF (kvac1(l,jspin).eq.kvac1(l1,jspin)) THEN
                            factorx = cmplx((locx(2)-locx(1)), 0.)
                         ELSE
                            k_diff=tpi*(kvac1(l,jspin)-kvac1(l1,jspin))
                            k_d1 = k_diff*locx(1)
                            k_d2 = k_diff*locx(2)
                            factorx=( cmplx( cos(k_d2), sin(k_d2)) -
     +                                cmplx( cos(k_d1), sin(k_d1)) ) /
     +                                cmplx( 0.,k_diff )
                         END IF
                         IF (kvac2(l,jspin).eq.kvac2(l1,jspin)) THEN
                            factory = cmplx((locy(2)-locy(1)), 0.)
                         ELSE
                            k_diff=tpi*(kvac2(l,jspin)-kvac2(l1,jspin))
                            k_d1 = k_diff*locy(1)
                            k_d2 = k_diff*locy(2)
                            factory=( cmplx( cos(k_d2), sin(k_d2)) -
     +                                cmplx( cos(k_d1), sin(k_d1)) ) /
     +                                cmplx( 0.,k_diff )
                         END IF
                         DO n=1, ne
                            aa = conjg(ac(l1,n,jspin))*ac(l,n,jspin)
                            bb = conjg(bc(l1,n,jspin))*bc(l,n,jspin)
                            ab = conjg(ac(l1,n,jspin))*bc(l,n,jspin)
                            ba = conjg(bc(l1,n,jspin))*ac(l,n,jspin)   
                            DO jj = 1,layers
c     
c     ---> either integrated (z1,z2) or slice (z1)
c     
                               IF (integ) THEN
                                  ll = 1
                                  DO ii = izlay(jj,1),
     +                                 izlay(jj,2)
                                     ui = u(ii,l,jspin)
                                     uei = ue(ii,l,jspin)
                                     uj = u(ii,l1,jspin)
                                     uej = ue(ii,l1,jspin)
                                     yy(ll) = real((aa*ui*uj+bb*uei*uej+
     +                                    ab*uei*uj+ba*ui*uej)
     +                                    *factorx*factory)
                                     ll = ll+1
                                  END DO
                                  CALL qsf(delz,yy,result,ll-1,0)
                                  qvlay(n,jj,ivac,ikpt,jspin) = 
     +                           qvlay(n,jj,ivac,ikpt,jspin) + result(1)
                               ELSE
                                  ui = u(izlay(jj,1),l,jspin)       
                                  uei = ue(izlay(jj,1),l,jspin)
                                  uj = u(izlay(jj,1),l1,jspin)
                                  uej = ue(izlay(jj,1),l1,jspin)
                                  qvlay(n,jj,ivac,ikpt,jspin) = real((
     +                                 aa*ui*uj + bb*uei*uej+ab*uei*uj
     +                              +ba*ui**uej)*factorx*factory)
                               END IF                                       
                            END DO
                         END DO 
                      END DO
                   END DO
                END IF 
             END IF 
          END IF  
      
c
c     **********************************************************************
c
c--->    warping part of the density (g.ne.0 stars)
c
c   ---> d_z^2-Tip
c
         IF (nstm.eq.2) THEN
            DO l = 1,nv2(jspin)
               DO 115 l1 = 1,l - 1
                  i1 = kvac1(l,jspin) - kvac1(l1,jspin)
                  i2 = kvac2(l,jspin) - kvac2(l1,jspin)
                  i3 = 0
                  IF (iabs(i1).GT.mx1) GO TO 115
                  IF (iabs(i2).GT.mx2) GO TO 115
                  ig3 = ig(i1,i2,i3)
                  IF (ig3.EQ.0)  GO TO 115
                  phs = rgphs(i1,i2,i3)
                  ig3p = ig(-i1,-i2,i3)
                  phsp = rgphs(-i1,-i2,i3)
                  ind2 = ig2(ig3)
                  ind2p = ig2(ig3p)
                  aa = 0.0
                  bb = 0.0
                  ba = 0.0
                  ab = 0.0
                  DO n = 1,ne
                     aa = aa + we(n)*conjg(ac(l1,n,jspin))*ac(l,n,jspin)
                     bb = bb + we(n)*conjg(bc(l1,n,jspin))*bc(l,n,jspin)
                     ab = ab + we(n)*conjg(ac(l1,n,jspin))*bc(l,n,jspin)
                     ba = ba + we(n)*conjg(bc(l1,n,jspin))*ac(l,n,jspin)
                  END DO
                  aae=-aa*2/3*tworkf
                  bbe=-bb*2/3*tworkf
                  abe=-ab*2/3*tworkf
                  bae=-ba*2/3*tworkf
                  aaee=aa*4/9*tworkf*tworkf
                  bbee=bb*4/9*tworkf*tworkf
                  abee=ab*4/9*tworkf*tworkf
                  baee=ba*4/9*tworkf*tworkf
                  DO 113 jz = 1,nmzxy
                     ui = u(jz,l,jspin)
                     uj = u(jz,l1,jspin)
                     ddui = ddu(jz,l)
                     dduj = ddu(jz,l1)
                     uei = ue(jz,l,jspin)
                     uej = ue(jz,l1,jspin)
                     dduei = ddue(jz,l)
                     dduej = ddue(jz,l1)
                     t1=aaee*ui*uj+bbee*uei*uej+baee*ui*uej+abee*uei*uj
     +                 + aae*(ui*dduj+uj*ddui)+bbe*(uei*dduej+uej*dduei)
     +                 + abe*(ui*dduej+uj*dduei)+bae*(ddui*uej+dduj*uei)
     +                 + aa*ddui*dduj+bb*dduei*dduej+ba*ddui*dduej
     +                 + ab*dduei*dduj  
                     rhtxy(jz,ind2-1,ivac,jspin) = rhtxy(jz,ind2-1,ivac,
     +                    jspin) + t1*phs/nstr2(ind2)
                     rhtxy(jz,ind2p-1,ivac,jspin)=rhtxy(jz,ind2p-1,ivac,
     +                    jspin) + conjg(t1)*phsp/nstr2(ind2p)
 113              CONTINUE
 115           CONTINUE
            END DO
c
c ---> s-Tip and p_z-Tip
c
         ELSE
!=============================================================
!           continuation of vacden....
!=============================================================
            IF (l_noco) THEN
c--->       diagonal elements of the density matrix, n_11 and n_22
               DO ispin = 1,jspins
                IF (odi%d1) THEN
                 DO l = 1,nv2(ispin)
                  DO m = -odi%mb,odi%mb
                   lprimee: DO l1 = 1,l-1
                    mprimee: DO m1 = -odi%mb,m-1
                       i3 = kvac3(l,ispin) - kvac3(l1,ispin)
                       m3 = m-m1
                       IF (m3.EQ.0 .and. i3.EQ.0) CYCLE mprimee
                       IF (iabs(m3).GT.odi%M) CYCLE mprimee
                       IF (iabs(i3).GT.k3d) CYCLE lprimee
                       ind1 = odi%ig(i3,m3)
                       ind1p = odi%ig(-i3,-m3)
                       IF (ind1.NE.0 .OR. ind1p.NE.0) THEN
                       aa = cmplx(0.,0.)
                       bb = cmplx(0.,0.)
                       ba = cmplx(0.,0.)
                       ab = cmplx(0.,0.)
                       DO n = 1,ne 
                        aa=aa+we(n)*conjg(ac_1(l1,m1,n,ispin))*
     *                                      ac_1(l,m,n,ispin)
                        bb=bb+we(n)*conjg(bc_1(l1,m1,n,ispin))*
     *                                      bc_1(l,m,n,ispin)
                        ab=ab+we(n)*conjg(ac_1(l1,m1,n,ispin))*
     *                                      bc_1(l,m,n,ispin)
                        ba=ba+we(n)*conjg(bc_1(l1,m1,n,ispin))*
     *                                      ac_1(l,m,n,ispin)
                       END DO
                       xys1: DO jz = 1,nmzxy
                        ui = u_1(jz,l,m,ispin)
                        uj = u_1(jz,l1,m1,ispin)
                        uei = ue_1(jz,l,m,ispin)
                        uej = ue_1(jz,l1,m1,ispin)
                        t1 = aa*ui*uj + bb*uei*uej +
     +                       ba*ui*uej + ab*uei*uj
                        IF (ind1.NE.0) THEN
                         rhtxy(jz,ind1-1,ivac,ispin) =
     =                   rhtxy(jz,ind1-1,ivac,ispin) + t1/
     /                        odi%nst2(ind1)
                        END IF
                        IF (ind1p.NE.0) THEN
                         rhtxy(jz,ind1p-1,ivac,ispin) =
     =                   rhtxy(jz,ind1p-1,ivac,ispin) + conjg(t1)/
     /                        odi%nst2(ind1p)
                        END IF

                       END DO xys1
                      END IF   ! ind1 and ind1p =0
                     END DO mprimee
                   END DO lprimee
                  END DO  ! m
                 END DO   ! l
                ELSE
                 DO l = 1,nv2(ispin)
                  DO 131 l1 = 1,l - 1
                     i1 = kvac1(l,ispin) - kvac1(l1,ispin)
                     i2 = kvac2(l,ispin) - kvac2(l1,ispin)
                     i3 = 0
                     IF (iabs(i1).GT.mx1) GO TO 131
                     IF (iabs(i2).GT.mx2) GO TO 131
                     ig3 = ig(i1,i2,i3)
                     IF (ig3.EQ.0)  GO TO 131
                     phs = rgphs(i1,i2,i3)
                     ig3p = ig(-i1,-i2,i3)
                     phsp = rgphs(-i1,-i2,i3)
                     ind2 = ig2(ig3)
                     ind2p = ig2(ig3p)
                     aa = 0.0
                     bb = 0.0
                     ba = 0.0
                     ab = 0.0
                     DO n = 1,ne
                      aa=aa+we(n)*conjg(ac(l1,n,ispin))*ac(l,n,ispin)
                      bb=bb+we(n)*conjg(bc(l1,n,ispin))*bc(l,n,ispin)
                      ab=ab+we(n)*conjg(ac(l1,n,ispin))*bc(l,n,ispin)
                      ba=ba+we(n)*conjg(bc(l1,n,ispin))*ac(l,n,ispin)
                     END DO
                     DO jz = 1,nmzxy
                        ui = u(jz,l,ispin)
                        uj = u(jz,l1,ispin)
                        uei = ue(jz,l,ispin)
                        uej = ue(jz,l1,ispin)
                        t1 = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                        rhtxy(jz,ind2-1,ivac,ispin) =
     +                       rhtxy(jz,ind2-1,ivac,ispin)
     +                       + t1*phs/nstr2(ind2)
                        rhtxy(jz,ind2p-1,ivac,ispin) =
     +                       rhtxy(jz,ind2p-1,ivac,ispin)
     +                       + conjg(t1)*phsp/nstr2(ind2p)
                     ENDDO
 131              CONTINUE
                 ENDDO
                END IF ! odi%d1
               END DO
c--->          off-diagonal element of the density matrix, n_21
               IF (odi%d1) THEN
                 DO l = 1,nv2(1)
                  DO m = -odi%mb,odi%mb
                   lprimea: DO l1 = 1,nv2(2)
                    mprimea: DO m1 = -odi%mb,odi%mb
                       i3 = kvac3(l,1) - kvac3(l1,2)
                       m3 = m-m1
                       IF (iabs(m3).GT.odi%M) CYCLE mprimea
                       IF (iabs(i3).GT.k3d) CYCLE lprimea
                       ind1 = odi%ig(i3,m3)
                      IF (ind1.NE.0) THEN
                       IF (m3.EQ.0 .and. i3.EQ.0) THEN
                        aa = cmplx(0.,0.)
                        bb = cmplx(0.,0.)
                        ba = cmplx(0.,0.)
                        ab = cmplx(0.,0.)
                        DO n = 1,ne
                         aa=aa+we(n)*conjg(ac_1(l1,m1,n,2))*
     *                                      ac_1(l,m,n,1)
                         bb=bb+we(n)*conjg(bc_1(l1,m1,n,2))*
     *                                      bc_1(l,m,n,1)
                         ab=ab+we(n)*conjg(ac_1(l1,m1,n,2))*
     *                                      bc_1(l,m,n,1)
                         ba=ba+we(n)*conjg(bc_1(l1,m1,n,2))*
     *                                      ac_1(l,m,n,1)
                        END DO
                        xys3: DO jz = 1,nmzxy
                         ui = u_1(jz,l,m,1)
                         uj = u_1(jz,l1,m1,2)
                         uei = ue_1(jz,l,m,1)
                         uej = ue_1(jz,l1,m1,2)
                         cdomvz(jz,ivac) =
     =                   cdomvz(jz,ivac) + 
     +                        aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                        END DO xys3
                       ELSE ! the warped part (ind1 > 1)
                        aa = cmplx(0.,0.)
                        bb = cmplx(0.,0.)
                        ba = cmplx(0.,0.)
                        ab = cmplx(0.,0.)
                        DO n = 1,ne
                         aa=aa+we(n)*conjg(ac_1(l1,m1,n,2))*
     *                                      ac_1(l,m,n,1)
                         bb=bb+we(n)*conjg(bc_1(l1,m1,n,2))*
     *                                      bc_1(l,m,n,1)
                         ab=ab+we(n)*conjg(ac_1(l1,m1,n,2))*
     *                                      bc_1(l,m,n,1)
                         ba=ba+we(n)*conjg(bc_1(l1,m1,n,2))*
     *                                      ac_1(l,m,n,1)
                        END DO
                        xys2: DO jz = 1,nmzxy
                         ui = u_1(jz,l,m,1)
                         uj = u_1(jz,l1,m1,2)
                         uei = ue_1(jz,l,m,1)
                         uej = ue_1(jz,l1,m1,2)
                         t1 = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                         cdomvxy(jz,ind1-1,ivac) =
     =                   cdomvxy(jz,ind1-1,ivac) + t1/odi%nst2(ind1)
                        END DO xys2
                       END IF  ! the non-warped (ind1 = 1)
                      END IF   ! ind1.ne.0
                     END DO mprimea
                   END DO lprimea
                  END DO  ! m
                 END DO   ! l
               ELSE ! odi%d1
                DO l = 1,nv2(1)
                  DO 132 l1 = 1,nv2(2)
                     i1 = kvac1(l,1) - kvac1(l1,2)
                     i2 = kvac2(l,1) - kvac2(l1,2)
                     i3 = 0
c--->                treat only the warping part
                     IF (iabs(i1).GT.mx1) GO TO 132
                     IF (iabs(i2).GT.mx2) GO TO 132
                     ig3 = ig(i1,i2,i3)
                     IF (ig3.EQ.0)  GO TO 132
                     phs = rgphs(i1,i2,i3)
                     ind2 = ig2(ig3)
                     IF ( ind2.EQ.1) THEN
c--->                non-warping part (1st star G=0)
                        aa = 0.0
                        bb = 0.0
                        ba = 0.0
                        ab = 0.0
                        DO ie = 1,ne
                           aa=aa+we(ie)*conjg(ac(l1,ie,2))*ac(l,ie,1)
                           bb=bb+we(ie)*conjg(bc(l1,ie,2))*bc(l,ie,1)
                           ab=ab+we(ie)*conjg(ac(l1,ie,2))*bc(l,ie,1)
                           ba=ba+we(ie)*conjg(bc(l1,ie,2))*ac(l,ie,1)
                        END DO
                        DO jz = 1,nmz
                           ui = u(jz,l,1)
                           ui2 = u(jz,l1,2)
                           uei = ue(jz,l,1)
                           uei2 = ue(jz,l1,2)
                           cdomvz(jz,ivac) = cdomvz(jz,ivac) + aa*ui2*ui
     +                         + bb*uei2*uei + ab*ui2*uei + ba*uei2*ui
                        ENDDO
                     ELSE
c--->                warping part
                        aa = 0.0
                        bb = 0.0
                        ba = 0.0
                        ab = 0.0
                        DO ie = 1,ne
                           aa=aa + we(ie)*conjg(ac(l1,ie,2))*ac(l,ie,1)
                           bb=bb + we(ie)*conjg(bc(l1,ie,2))*bc(l,ie,1)
                           ab=ab + we(ie)*conjg(ac(l1,ie,2))*bc(l,ie,1)
                           ba=ba + we(ie)*conjg(bc(l1,ie,2))*ac(l,ie,1)
                        END DO
                        DO  jz = 1,nmzxy
                           ui = u(jz,l,1)
                           uj = u(jz,l1,2)
                           uei = ue(jz,l,1)
                           uej = ue(jz,l1,2)
                           t1 = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                           cdomvxy(jz,ind2-1,ivac) = cdomvxy(jz,
     +                          ind2-1,ivac) + t1*phs/nstr2(ind2)
                        ENDDO
                     ENDIF
 132              CONTINUE
                END DO
               END IF ! odi%d1
            ELSE                                ! collinear part

              IF (odi%d1) THEN
               DO l = 1,nv2(jspin)
                DO m = -odi%mb,odi%mb
                 lprime: DO l1 = 1,l-1
                  mprime: DO m1 = -odi%mb,m-1
                     i3 = kvac3(l,jspin) - kvac3(l1,jspin)
                     m3 = m-m1
                     IF (m3.EQ.0 .and. i3.EQ.0) CYCLE mprime
                     IF (iabs(m3).GT.odi%M) CYCLE mprime
                     IF (iabs(i3).GT.k3d) CYCLE lprime
                     ind1 = odi%ig(i3,m3)
                     ind1p = odi%ig(-i3,-m3)
                     IF (ind1.NE.0 .OR. ind1p.NE.0) THEN
                     aa = cmplx(0.,0.)
                     bb = cmplx(0.,0.)
                     ba = cmplx(0.,0.)
                     ab = cmplx(0.,0.)
                     DO n = 1,ne
                      aa=aa+we(n)*conjg(ac_1(l1,m1,n,jspin))*
     *                                      ac_1(l,m,n,jspin)
                      bb=bb+we(n)*conjg(bc_1(l1,m1,n,jspin))*
     *                                      bc_1(l,m,n,jspin)
                      ab=ab+we(n)*conjg(ac_1(l1,m1,n,jspin))*
     *                                      bc_1(l,m,n,jspin)
                      ba=ba+we(n)*conjg(bc_1(l1,m1,n,jspin))*
     *                                      ac_1(l,m,n,jspin)
                     END DO
                     xys: DO jz = 1,nmzxy
                      ui = u_1(jz,l,m,jspin)
                      uj = u_1(jz,l1,m1,jspin)
                      uei = ue_1(jz,l,m,jspin)
                      uej = ue_1(jz,l1,m1,jspin)
                      t1 = aa*ui*uj + bb*uei*uej +
     +                     ba*ui*uej + ab*uei*uj
                      IF (ind1.NE.0) THEN
                       rhtxy(jz,ind1-1,ivac,jspin) =
     =                 rhtxy(jz,ind1-1,ivac,jspin) + t1/
     /                        odi%nst2(ind1)
                      END IF
                      IF (ind1p.NE.0) THEN
                       rhtxy(jz,ind1p-1,ivac,jspin) =
     =                 rhtxy(jz,ind1p-1,ivac,jspin) + conjg(t1)/
     /                        odi%nst2(ind1p)
                      END IF

                     END DO xys
                    END IF   ! ind1 and ind1p =0
                   END DO mprime
                 END DO lprime
                END DO  ! m
               END DO   ! l

              ELSE         !D1

               DO l = 1,nv2(jspin)
                  DO 130 l1 = 1,l - 1
                     i1 = kvac1(l,jspin) - kvac1(l1,jspin)
                     i2 = kvac2(l,jspin) - kvac2(l1,jspin)
                     i3 = 0
                     IF (iabs(i1).GT.mx1) GO TO 130
                     IF (iabs(i2).GT.mx2) GO TO 130
                     ig3 = ig(i1,i2,i3)
                     IF (ig3.EQ.0)  GO TO 130
                     phs = rgphs(i1,i2,i3)
                     ig3p = ig(-i1,-i2,i3)
                     phsp = rgphs(-i1,-i2,i3)
                     ind2 = ig2(ig3)
                     ind2p = ig2(ig3p)
                     aa = 0.0
                     bb = 0.0
                     ba = 0.0
                     ab = 0.0
                     DO n = 1,ne
                        aa=aa+we(n)*conjg(ac(l1,n,jspin))*ac(l,n,jspin)
                        bb=bb+we(n)*conjg(bc(l1,n,jspin))*bc(l,n,jspin)
                        ab=ab+we(n)*conjg(ac(l1,n,jspin))*bc(l,n,jspin)
                        ba=ba+we(n)*conjg(bc(l1,n,jspin))*ac(l,n,jspin)
                     END DO
                     DO 120 jz = 1,nmzxy
                        ui = u(jz,l,jspin)
                        uj = u(jz,l1,jspin)
                        uei = ue(jz,l,jspin)
                        uej = ue(jz,l1,jspin)
                        t1 = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                        rhtxy(jz,ind2-1,ivac,jspin) = rhtxy(jz,ind2-1,
     +                       ivac,jspin) + t1*phs/nstr2(ind2)
                        rhtxy(jz,ind2p-1,ivac,jspin) = rhtxy(jz,ind2p-1,
     +                       ivac,jspin) + conjg(t1)*phsp/nstr2(ind2p)
 120                 CONTINUE
 130              CONTINUE
               END DO
              END IF ! D1
            ENDIF
         END IF
!=============================================================
c
c       calculate 1. to nstars. starcoefficient for each k and energy eigenvalue 
c           to stcoeff(ne,layer,ivac) if starcoeff=T (the star coefficient values are written to vacdos)
c
        IF (starcoeff .and. vacdos) THEN
           DO 139 n=1,ne
              DO l = 1,nv2(jspin)
                 DO 135 l1 = 1,l - 1
                    i1 = kvac1(l,jspin) - kvac1(l1,jspin)
                    i2 = kvac2(l,jspin) - kvac2(l1,jspin)
                    i3 = 0
                    IF (iabs(i1).GT.mx1) GO TO 135
                    IF (iabs(i2).GT.mx2) GO TO 135
                    ig3 = ig(i1,i2,i3)
                    IF (ig3.EQ.0)  GO TO 135
                    ind2 = ig2(ig3)
                    ig3p = ig(-i1,-i2,i3)
                    ind2p = ig2(ig3p)
                    IF ((ind2.GE.2.AND.ind2.LE.nstars).OR.
     +                   (ind2p.GE.2.AND.ind2p.LE.nstars)) THEN
                       phs = rgphs(i1,i2,i3)
                       phsp = rgphs(-i1,-i2,i3)
                       aa = conjg(ac(l1,n,jspin))*ac(l,n,jspin)
                       bb = conjg(bc(l1,n,jspin))*bc(l,n,jspin)
                       ab = conjg(ac(l1,n,jspin))*bc(l,n,jspin)
                       ba = conjg(bc(l1,n,jspin))*ac(l,n,jspin)
                       DO jj = 1,layers
                          ui = u(izlay(jj,1),l,jspin)
                          uj = u(izlay(jj,1),l1,jspin)
                          uei = ue(izlay(jj,1),l,jspin)
                          uej = ue(izlay(jj,1),l1,jspin)
                          t1 = aa*ui*uj + bb*uei*uej
     +                         +ba*ui*uej + ab*uei*uj
                          IF (ind2.GE.2.AND.ind2.LE.nstars) 
     +                      stcoeff(ind2-1,n,jj,ivac) =
     +                      stcoeff(ind2-1,n,jj,ivac)+ 
     +                         t1*phs/nstr2(ind2)
                          IF (ind2p.GE.2.AND.ind2p.LE.nstars) 
     +                stcoeff(ind2p-1,n,jj,ivac) = 
     +                         stcoeff(ind2p-1,n,jj,ivac) 
     +                                       +conjg(t1)*phs/nstr2(ind2p)
                       END DO
                    END IF
 135             CONTINUE  
              END DO
 139       CONTINUE
         END IF
  140 CONTINUE
      IF (l_noco) CLOSE (npotmatfile)
      DEALLOCATE (ac,bc,dt,dte,du,ddu,due,ddue,t,te,tei,u,ue,v,yy )
      DEALLOCATE (ac_1,bc_1,dt_1,dte_1,du_1,ddu_1,due_1,ddue_1)
      DEALLOCATE (t_1,te_1,tei_1,u_1,ue_1)

      END SUBROUTINE vacden
      END MODULE m_vacden
