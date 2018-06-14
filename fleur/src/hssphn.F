      MODULE m_hssphn
      CONTAINS
      SUBROUTINE hssphn(
     >                  nvd,lmaxd,ntypd,nlhd,jmtd,natd,nop,jspd,nwdd,
     >                  llod,nlod,loplod,matsize,lmd,lmplmd,nlotot,
     >                  ntype,neq,lmax,jri,rmt,rmsh,dx,vr,taual,
     >                  bmat,mrot,ngopr,omtil,
     >                  nb0,jsp,jspins,nw,secvar,el,ello,lnonsph,
     >                  us,dus,uds,duds,ddn,nlo,llo,l_dulo,ulo_der,
     >                  kveclo,ulos,dulos,uulon,dulon,uloulopn,lapw_l,
     >                  bkpt,k1,k2,k3,rk,nv,invsat,invtab,
     >                  l_noco,alph,beta,l_constr,b_con,l_ss,qss,
     >                  l_soc,soc_opt,n_u,vs_mmp,lmaxb,lda_u,odi,ods,
     >                  tuu,tud,tdu,tdd,ind,tuulo,tdulo,tuloulo,mlotot,
     >                  mlolotot,n_size,n_rank,SUB_COMM,gwc,irank,
     X                  aa,bb,nmat,cp_sph,cp_nsph)
c*********************************************************************
c     updates the hamiltonian and overlap matrices with the
c     contributions from the spheres, both spherical and non-
c     spherical.
c                m. weinert  1986
c     modified for inversion symmetry (real version), rp.p nov. 92
c
c     Since this non-collinear version allows only one definite 
c     direction of magnetic field in each muffin-tin, no matrix 
c     potential has to be used inside the MT's. However, the boundary
c     conditions change, because the local (spin-) frame is different
c     from the global frame. As a result, the usual spin-up and -down
c     Hamitonian- and overlapp-matrix elements have to be
c     premultiplied by a phasefactor, depending on the direction of 
c     the magentic field of each MT, to form the elements of the full
c     complex matrices.
c
c     Philipp Kurz 98/01/27
c*********************************************************************
c------------------------------------------------------------------+
c Note for ev-parallelization:                                     |
c                             unlike in hsint and hsvac, we do not |
c move through the H-matrix block by block (i.e. first the up/up,  |
c then the down/down, and finally the up/down spin blocks) but     |
c go  through one block only and update all blocks simultaniously. |
c Therefore, virtually we move through the whole matrix (up to     |
c nv(1)+nv(2)) and project back to the first block.                |
c                                                   gb-00          |
c------------------------------------------------------------------+
c******** ABBREVIATIONS ***********************************************
c     aa       : hamitonian matrix
c     bb       : overlapp  matrix
c     alph,beta: Euler angles of the local magnetic field direction of
c                each atom (-type). 
c     chi      : Pauli spinors of the local spin-coordinate-frames
c     chinn    : prefactors to be multiplied with the Hamiltonian- and
c                overlappmatrix elements
c**********************************************************************
c
c******** Spin-orbit interaction *************************************
c   When l_soc=true & l_noco=true & l_ss=false, spin-orbit interaction 
c   is added in the first variation
c   The formula by Youn et al. J.Comp.Phys. 172, 387 (2001) is used
c   (note the wrong sign in the paper): 
c   <G|Hso|G'>= -i sigma . (G x G') sum_l (2l+1)/4pi Rso P'_l(G.G')
c   
c   Jussi Enkovaara 2004, Juelich
c********************************************************************

#include"cpp_double.h"

      USE m_constants, ONLY : pimach
      USE m_matmul   , ONLY : matmul3,matmul3r
      USE m_sphbes
      USE m_dsphbs
      USE m_uham
      USE m_ylm
      USE m_cputime
      USE m_soinit
      USE m_abccoflo
      USE m_hlomat
      USE m_slomat
      USE m_vecforlo
      USE m_setabc1lo
      USE m_od_types, ONLY : od_inp, od_sym
#ifdef CPP_MPI
#if ( defined(CPP_T90) || defined(CPP_HDF) )
      USE m_mingeselle
#endif
#endif
#ifdef DUMP_DATA
      USE binmat ! for matrix dumping
#endif

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nvd,lmaxd,ntypd,nlhd,jmtd,natd,nop,jspd
      INTEGER, INTENT (IN) :: nwdd,SUB_COMM,n_size,n_rank,nlotot
      INTEGER, INTENT (IN) :: jsp,jspins,nb0,nw,ntype,matsize,irank
      INTEGER, INTENT (IN) :: lmd,lmplmd,llod,nlod,loplod,n_u,lmaxb
      INTEGER, INTENT (IN) :: mlotot,mlolotot,gwc
      LOGICAL, INTENT (IN) :: secvar,l_noco,l_constr,l_ss,l_soc
      REAL,    INTENT (IN) :: omtil
      REAL,    INTENT (INOUT) :: cp_sph,cp_nsph
      INTEGER, INTENT (INOUT) :: nmat
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: lnonsph(ntypd),lda_u(ntypd),lapw_l(ntypd)
      INTEGER, INTENT (IN) :: invsat(natd),nlo(ntypd),llo(nlod,ntypd)
      INTEGER, INTENT (IN) :: invtab(nop),ulo_der(nlod,ntypd)
      INTEGER, INTENT (IN) :: nv(jspd),mrot(3,3,nop),jri(ntypd)
      INTEGER, INTENT (IN) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (IN) :: lmax(ntypd),neq(ntypd),ngopr(natd)
      REAL,    INTENT (IN) :: bkpt(3),rk(nvd,jspd)
      REAL,    INTENT (IN) :: rmt(ntypd),taual(3,natd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: dx(ntypd),vr(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (IN) :: bmat(3,3)
      REAL,    INTENT (IN) :: el(0:lmaxd,ntypd,jspd,nwdd)
      REAL,    INTENT (IN) :: ello(nlod,ntypd,jspd) 
      REAL,    INTENT (INOUT) :: us(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (INOUT) :: uds(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (INOUT) :: dus(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (INOUT) :: duds(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (INOUT) :: ulos(nlod,ntypd,jspd)
      REAL,    INTENT (INOUT) :: dulos(nlod,ntypd,jspd)
      REAL,    INTENT (INOUT) :: uulon(nlod,ntypd,jspd)
      REAL,    INTENT (INOUT) :: dulon(nlod,ntypd,jspd)
      REAL,    INTENT (IN) :: uloulopn(nlod,nlod,ntypd,jspd)
      REAL,    INTENT (INOUT) :: ddn(0:lmaxd,ntypd,jspd)
      REAL,    INTENT (IN) :: alph(ntypd),beta(ntypd),qss(3)
      REAL,    INTENT (IN) :: b_con(2,ntypd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      LOGICAL, INTENT (IN) :: soc_opt(ntypd+2) 
      INTEGER, INTENT (OUT) :: kveclo(nlotot)
#ifdef CPP_INVERSION
      REAL,    INTENT (INOUT) :: aa(matsize),bb(matsize)
#else
      COMPLEX, INTENT (INOUT) :: aa(matsize),bb(matsize)
#endif
      COMPLEX,INTENT(IN):: vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,n_u,jspins)
      COMPLEX,INTENT(IN):: tdulo(0:lmd,-llod:llod,mlotot,*)
      COMPLEX,INTENT(IN):: tuulo(0:lmd,-llod:llod,mlotot,*)
      COMPLEX,INTENT(IN):: tuloulo(-llod:llod,-llod:llod,mlolotot,*)
      COMPLEX :: tdd(0:lmplmd,ntypd,*),tdu(0:lmplmd,ntypd,*)
      COMPLEX :: tud(0:lmplmd,ntypd,*),tuu(0:lmplmd,ntypd,*)
      INTEGER :: ind(0:lmd,0:lmd,ntypd,*)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim

C     ..
C     .. Local Scalars ..
      REAL con1,ff,gg,gs,termi,termr,th,ws,
     +     t1nn,t2nn,t3nn,time1,time2,invsfct,doti,
     +     elall,fct,fjkiln,gjkiln,ddnln,ski1,ski2,ski3
      REAL tpi,fpi,sfp,apw_lo1,apw_lo2,apw1,w1
      REAL qss1,qss2,qss3
      COMPLEX ci,chi11,chi21,chi22,chihlp,capw1
      INTEGER ii,iii,ij,im,in,k,ki,kj,i,lmplm,
     +        l,ll1,lm,lmaxx,lmp,lp,lwn,lo,jd,
     +        m,mp,n,n0,n1,na,nn,np,mlo,mlolo,
     +        iiloh,iilos,nkvecprevath,nkvecprevats,kjmax,
     +        hlpmsize,isp,nsp,jsp_start,jsp_end,nintsp,iintsp,jintsp
      INTEGER nc,locolh,locols,nv_tot,kii,nkvecprevatu,iilou,locolu
      INTEGER ab_dim,nkvec_sv
      LOGICAL enough,l_lo1
#ifdef DISPLAY_INFO
      INTEGER :: xxx ! find some unused variable name...
#endif
#ifdef DUMP_DATA_OLD
C     .. for matrix dumping ..
      character(len=32) :: misc_name
#endif
#ifdef DUMP_ENABLE_SPHERICAL
      logical :: spherical_only
#endif
C     ..
C     .. Local Arrays ..
      INTEGER kvec(2*(2*llod+1),nlod)
      REAL alo1(nlod),blo1(nlod),clo1(nlod)
      REAL bmrot(3,3),dtd(2),dtu(2),fb(0:lmaxd),
     +     fleg1(0:lmaxd),fleg2(0:lmaxd),fl2p1(0:lmaxd),
     +     fl2p1bt(0:lmaxd),gb(0:lmaxd),gkrot(nvd,3),
     +     utd(2),utu(2),vmult(3),v(3)
      REAL qssbti(3),qssbtj(3)
      REAL, ALLOCATABLE :: ar(:,:,:),ai(:,:,:),br(:,:,:),bi(:,:,:)
      REAL, ALLOCATABLE :: axr(:),axi(:),bxr(:),bxi(:)
      REAL, ALLOCATABLE :: fj(:,:,:,:),gj(:,:,:,:),plegend(:,:)
      REAL, ALLOCATABLE :: rph(:,:),cph(:,:),gk(:,:,:),vk(:,:,:)
      REAL, ALLOCATABLE :: uun21(:,:),udn21(:,:),dun21(:,:),ddn21(:,:)
      COMPLEX ylm( (lmaxd+1)**2 )
      COMPLEX cwork(-2*llod:2*llod+1,2*(2*llod+1),nlod)
      COMPLEX chi(2,2),chj(2,2,2,ntypd),aawa(nvd),bbwa(nvd)
      COMPLEX, ALLOCATABLE :: aahlp(:),bbhlp(:)
      COMPLEX, ALLOCATABLE :: alo(:,:,:,:),blo(:,:,:,:),clo(:,:,:,:)
      INTEGER, ALLOCATABLE :: nkvec(:,:)
      LOGICAL apw(0:lmaxd)

      REAL, PARAMETER :: eps = 1.0e-30
c for Spin-orbit...
      LOGICAL :: l_socfirst
      LOGICAL, SAVE :: first_k = .TRUE.
      REAL, ALLOCATABLE :: dplegend(:,:)
      REAL, ALLOCATABLE :: cross_k(:,:)
      INTEGER :: j1,j2
      COMPLEX :: isigma(2,2,3),isigma_x(2,2),isigma_y(2,2),isigma_z(2,2)
      REAL    :: rsopp  (ntypd,lmaxd,2,2)
      REAL    :: rsoppd (ntypd,lmaxd,2,2)
      REAL    :: rsopdp (ntypd,lmaxd,2,2)
      REAL    :: rsopdpd(ntypd,lmaxd,2,2)
      REAL    :: rsoplop (ntypd,nlod,2,2)
      REAL    :: rsoplopd(ntypd,nlod,2,2)
      REAL    :: rsopdplo(ntypd,nlod,2,2)
      REAL    :: rsopplo (ntypd,nlod,2,2)
      REAL    :: rsoploplop(ntypd,nlod,nlod,2,2)
      COMPLEX :: chi11so(2,2),chi21so(2,2),chi22so(2,2),angso(nvd,2,2)
      COMPLEX :: chis(2,2)
      CHARACTER*3 :: chntype
#ifdef SOC_DEBUG
      COMPLEX :: radso
      INTEGER :: m1,m2,is1,is2,lm1,lm2
      COMPLEX :: soangl(0:lmaxd,-lmaxd:lmaxd,2,-lmaxd:lmaxd,2),angso2
      COMPLEX :: ylm1( (lmaxd+1)**2 ), ylm2( (lmaxd+1)**2 ),aawa2(nvd)
      INTEGER :: ispjsp(2)
      DATA ispjsp/1,-1/
      COMPLEX :: anglso
      EXTERNAL anglso
#endif
C     ..
C     .. External Subroutines ..
      EXTERNAL CPP_BLAS_caxpy
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cos,sin,sqrt
c
      ab_dim = 1
      IF (l_ss) ab_dim = 2
      ALLOCATE ( axr(nvd),axi(nvd),bxr(nvd),bxi(nvd),nkvec(nlod,ab_dim),
     +           ar(nvd,0:lmd,ab_dim),ai(nvd,0:lmd,ab_dim),
     +           br(nvd,0:lmd,ab_dim),bi(nvd,0:lmd,ab_dim),
     +           fj(nvd,0:lmaxd,ntypd,ab_dim),rph(nvd,ab_dim),
     +           gj(nvd,0:lmaxd,ntypd,ab_dim),cph(nvd,ab_dim),
     +           vk(nvd,3,ab_dim),gk(nvd,3,ab_dim),plegend(nvd,0:lmaxd),
     +           alo(-llod:llod,2*(2*llod+1),nlod,ab_dim),
     +           blo(-llod:llod,2*(2*llod+1),nlod,ab_dim),
     +           clo(-llod:llod,2*(2*llod+1),nlod,ab_dim) )
      IF (l_constr) THEN
         DEALLOCATE ( fj,gj )
         ALLOCATE ( fj(nvd,0:lmaxd,ntypd,2),gj(nvd,0:lmaxd,ntypd,2),
     +              uun21(0:lmaxd,ntypd),udn21(0:lmaxd,ntypd),
     +              dun21(0:lmaxd,ntypd),ddn21(0:lmaxd,ntypd) )
      ENDIF

      l_socfirst= l_soc .AND. l_noco .AND. (.NOT. l_ss)
      IF (l_socfirst) THEN
         DEALLOCATE ( fj,gj )
         ALLOCATE ( dplegend(nvd,0:lmaxd),cross_k(nvd,3) )
         ALLOCATE ( fj(nvd,0:lmaxd,ntypd,2),gj(nvd,0:lmaxd,ntypd,2))
c     isigma= -i * sigma, where sigma is Pauli matrix
         isigma=cmplx(0.0,0.0)
         isigma(1,2,1)=cmplx(0.0,-1.0)
         isigma(2,1,1)=cmplx(0.0,-1.0)
         isigma(1,2,2)=cmplx(-1.0,0.0)
         isigma(2,1,2)=cmplx(1.0,0.0)
         isigma(1,1,3)=cmplx(0.0,-1.0)
         isigma(2,2,3)=cmplx(0.0,1.0)

         CALL soinit(
     >            jmtd,jspd,ntypd,nwdd,lmaxd,nlhd,nlod,jri,jspins,
     >            ntype,nw,lmax,dx,rmsh,el,ello,nlo,llo,l_dulo,
     >            ulo_der,vr,soc_opt(ntype+2),
     <            rsopp,rsoppd,rsopdp,rsopdpd,ddn,
     <            rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     <            us,dus,uds,duds,ulos,dulos,uulon,dulon)
         IF (soc_opt(ntype+1)) THEN
           DO n= 1,ntype
             IF (.not. soc_opt(n)) THEN 
               rsopp(n,:,:,:) = 0.0 
               rsopdp(n,:,:,:) = 0.0 
               rsoppd(n,:,:,:) = 0.0 
               rsopdpd(n,:,:,:) = 0.0 
               rsoplop(n,:,:,:) = 0.0 
               rsoplopd(n,:,:,:) = 0.0 
               rsopdplo(n,:,:,:) = 0.0 
               rsopplo(n,:,:,:) = 0.0  
               rsoploplop(n,:,:,:,:) = 0.0
             ENDIF
           ENDDO
         ENDIF 

         IF ((first_k).AND.(irank.EQ.0)) THEN
            DO n = 1,ntype
               WRITE (6,FMT=8000)
               WRITE (6,FMT=9000)
               WRITE (6,FMT=8001) (2*rsopp(n,l,1,1),l=1,3)
               WRITE (6,FMT=8001) (2*rsopp(n,l,2,2),l=1,3)
               WRITE (6,FMT=8001) (2*rsopp(n,l,2,1),l=1,3)
            ENDDO
            IF (soc_opt(ntype+1)) THEN
              WRITE (chntype,'(i3)') ntype
              WRITE (6,fmt='(A,2x,'//chntype//'l1)') 'SOC contributi'
     &         //'on of certain atom types included in Hamiltonian:',
     &         (soc_opt(n),n=1,ntype)
            ELSE
              WRITE(6,fmt='(A,1x,A)') 'SOC contribution of all atom'//
     &         ' types inculded in Hamiltonian.'
            ENDIF
            IF (soc_opt(ntype+2)) THEN
              WRITE(6,fmt='(A)')
     &         'SOC Hamiltonian is constructed by neglecting B_xc.'
            ENDIF
            first_k=.FALSE.
         ENDIF
 8000 FORMAT (' spin - orbit parameter HR  ')
 8001 FORMAT (8f8.4)
 9000 FORMAT (5x,' p ',5x,' d ', 5x, ' f ')

      ENDIF

#ifdef CPP_MPI
c
c determine size of close-packed matrix
c
      nc = 0
      k = 0
      DO  i=1+n_rank, nv(1)+nlotot, n_size
        nc = nc + 1
        k = k + n_size*(nc-1) + n_rank + 1
      ENDDO
      hlpmsize = k
#else
      hlpmsize = (nvd+nlotot)*(nvd+nlotot+1)/2
#endif
      IF ( l_noco .AND. (.NOT. l_ss) )
     +     ALLOCATE ( aahlp(hlpmsize),bbhlp(hlpmsize) )
C     ..
      tpi = 2 * pimach()
      fpi = 2 * tpi
      sfp = sqrt( fpi )
      con1 = fpi/sqrt(omtil)
      ci = cmplx(0.0,1.0)
      DO l = 0,lmaxd
         fleg1(l) = real(l+l+1)/real(l+1)
         fleg2(l) = real(l)/real(l+1)
         fl2p1(l) = real(l+l+1)/fpi
         fl2p1bt(l) = fl2p1(l)*0.5
      END DO

      
c---> pk non-collinear
c---> calculate the overlap of the radial basis functions with different
c---> spin directions.
      IF (l_constr) CALL rad_ovlp(
     >           jspd,ntypd,lmaxd,nlhd,jmtd,
     >           jspins,ntype,lmax,jri,rmsh,dx,vr,el(0,1,1,nw),
     <           uun21,udn21,dun21,ddn21)
c---> pk non-collinear
c---> loop over the local spins and interstitial spins
      nintsp = 1
      IF (l_noco) THEN
         jsp_start = 1
         jsp_end   = 2
c--->    in a spin-spiral calculations the a- and b-coeff and other
c--->    quantities depend on the interstitial spin. therefore, an
c--->    additional loop over the interstitial spin is needed in some
c--->    places.
         IF (l_ss) nintsp = 2
      ELSE
         jsp_start = jsp
         jsp_end   = jsp
      ENDIF

      DO iintsp = 1,nintsp
         IF (iintsp .EQ. 1) THEN
            qss1= - qss(1)/2
            qss2= - qss(2)/2
            qss3= - qss(3)/2
         ELSE
            qss1= + qss(1)/2
            qss2= + qss(2)/2
            qss3= + qss(3)/2
         ENDIF
         DO k = 1,nv(iintsp)
            v(1) = bkpt(1) + k1(k,iintsp) + qss1
            v(2) = bkpt(2) + k2(k,iintsp) + qss2
            v(3) = bkpt(3) + k3(k,iintsp) + qss3
            vk(k,1,iintsp) = v(1)
            vk(k,2,iintsp) = v(2)
            vk(k,3,iintsp) = v(3)
            doti = max (rk(k,iintsp),eps)
            plegend(k,0) = 1.e0
            IF (l_socfirst) THEN 
               dplegend(k,0)=0.e0
               dplegend(k,1)=1.e0
            ENDIF
            gk(k,1,iintsp) = ( bmat(1,1)*v(1) + bmat(2,1)*v(2)
     +                       + bmat(3,1)*v(3) )/doti
            gk(k,2,iintsp) = ( bmat(1,2)*v(1) + bmat(2,2)*v(2)
     +                       + bmat(3,2)*v(3) )/doti
            gk(k,3,iintsp) = ( bmat(1,3)*v(1) + bmat(2,3)*v(2)
     +                       + bmat(3,3)*v(3) )/doti
         END DO
      ENDDO
#ifdef DISPLAY_INFO
      print*, 'HSSPHN: gk = '
      print '(3G15.5)',(gk(xxx, 1, 1), gk(xxx, 2, 1), gk(xxx, 3, 1),
     +                      xxx=1, nv(1))
#endif

      DO isp = jsp_start,jsp_end

c---> loop over each atom type
      lmaxx = 0
      DO iintsp = 1,nintsp
         DO n = 1,ntype
            lwn = lmax(n)
            lmaxx = max(lmaxx,lwn)
            DO l = 0,lwn
              apw(l) = .false.
              DO lo = 1,nlo(n)
                IF (l_dulo(lo,n)) apw(l) = .true.
              ENDDO
#ifdef CPP_APW
              IF (lapw_l(n).GE.l) apw(l) = .false.
#endif
            ENDDO
            DO lo = 1,nlo(n)
              IF (l_dulo(lo,n)) apw(llo(lo,n)) = .true.
            ENDDO
#ifdef DISPLAY_INFO
      print*, 'rk = '
      print '(3G10.3)', (rk(k, 1), k=1,nv(1))
#endif
            DO k = 1,nv(iintsp)
               gs = rk(k,iintsp)*rmt(n)
               CALL sphbes(
     >                     lwn,gs,
     <                     fb)
               CALL dsphbs(
     >                     lwn,gs,fb,
     <                     gb)
               DO l = 0,lwn
c---> set up wronskians for the matching conditions for each ntype
                  ws = con1/(uds(l,n,isp)*dus(l,n,isp) 
     +                    - us(l,n,isp)*duds(l,n,isp))
                  ff = fb(l)
                  gg = rk(k,iintsp)*gb(l)
                  IF ( apw(l) ) THEN
                    fj(k,l,n,iintsp) = 1.0*con1 * ff / us(l,n,isp) 
                    gj(k,l,n,iintsp) = 0.0d0 
                  ELSE
                    IF (l_constr.or.l_socfirst) THEN
c--->                 in a constrained calculation fj and gj are needed for
c--->                 both local spin directions (isp) at the same time
                      fj(k,l,n,isp) = ws * ( uds(l,n,isp)*gg - 
     +                                      duds(l,n,isp)*ff )
                      gj(k,l,n,isp) = ws * ( dus(l,n,isp)*ff -
     +                                        us(l,n,isp)*gg )
                    ELSE
c--->                 in a spin-spiral calculation fj and gj are needed for
c--->                 both interstitial spin directions at the same time.
c--->                 In all other cases iintsp runs from 1 to 1.
                      fj(k,l,n,iintsp) = ws * ( uds(l,n,isp)*gg - 
     +                                         duds(l,n,isp)*ff )
                      gj(k,l,n,iintsp) = ws * ( dus(l,n,isp)*ff -
     +                                           us(l,n,isp)*gg )
#ifdef DISPLAY_INFO
! print fj, gj (which are basically AB coeff without factors ylm etc.)
      print '("fj(", 3I3, ")=", G13.5)', k, l, n, fj(k, l, n, iintsp)
#endif
                    ENDIF
                  ENDIF
               ENDDO
            ENDDO ! k = 1, nv
         ENDDO    ! n = 1,ntype
      ENDDO       ! iintsp = 1,nintsp
c
      CALL cpu_time(time1)

      IF ( l_noco .AND. (.NOT. l_ss) ) THEN
c---> pk non-collinear
c--->    initialize hamiltonian help array
         DO ii = 1,hlpmsize
            aahlp(ii)=cmplx(0.0,0.0)
            bbhlp(ii)=cmplx(0.0,0.0)
         ENDDO
c---> pk non-collinear
      ENDIF

      DO iintsp = 1,nintsp
        IF (iintsp.EQ.1) THEN
          qssbti(1) = - qss(1)/2
          qssbti(2) = - qss(2)/2
          qssbti(3) = - qss(3)/2
        ELSE
          qssbti(1) = + qss(1)/2
          qssbti(2) = + qss(2)/2
          qssbti(3) = + qss(3)/2
        ENDIF
        DO jintsp = 1,iintsp
         IF (jintsp.EQ.1) THEN
            qssbtj(1) = - qss(1)/2
            qssbtj(2) = - qss(2)/2
            qssbtj(3) = - qss(3)/2
         ELSE
            qssbtj(1) = + qss(1)/2
            qssbtj(2) = + qss(2)/2
            qssbtj(3) = + qss(3)/2
         ENDIF

         nc = 0                                    ! maybe IF (iintsp.EQ.1) nc = 0
         IF ( l_noco .AND. (n_size.GT.1) ) THEN
c--->       for EV-parallelization & noco ( see comments at top )
            nv_tot = nv(1) + nv(2)
            IF (l_ss) STOP 'ev-|| & spin-spiral !' ! maybe  nv_tot = nv(iintsp) works
         ELSE
            nv_tot = nv(iintsp)
         ENDIF
         DO 40 kii =  n_rank, nv_tot-1, n_size
            ki = mod(kii,nv(iintsp)) + 1  
            nc = nc + 1
            iii = nc*(nc-1)/2 * n_size - (nc-1)*(n_size - n_rank - 1)
c            ii = nc*(nc+1)/2 * n_size -  nc   *(n_size - n_rank - 1)

            IF (l_ss.OR.l_constr.OR.l_socfirst) THEN
               kjmax = nv(jintsp)
            ELSE
               kjmax = ki
            ENDIF
            ski1 = k1(ki,iintsp) + qssbti(1)
            ski2 = k2(ki,iintsp) + qssbti(2)
            ski3 = k3(ki,iintsp) + qssbti(3)
c--->       legendre polynomials
            DO kj = 1,kjmax
               plegend(kj,1) = gk(kj,1,jintsp)*gk(ki,1,iintsp) +
     +                         gk(kj,2,jintsp)*gk(ki,2,iintsp) +
     +                         gk(kj,3,jintsp)*gk(ki,3,iintsp)
               IF (l_socfirst) THEN
c#ifdef CPP_SOCFIRST
                  cross_k(kj,1)=gk(ki,2,jintsp)*gk(kj,3,iintsp)-
     +                          gk(ki,3,jintsp)*gk(kj,2,iintsp)
                  cross_k(kj,2)=gk(ki,3,jintsp)*gk(kj,1,iintsp)-
     +                          gk(ki,1,jintsp)*gk(kj,3,iintsp)
                  cross_k(kj,3)=gk(ki,1,jintsp)*gk(kj,2,iintsp)-
     +                          gk(ki,2,jintsp)*gk(kj,1,iintsp)
               ENDIF
c#endif
            END DO 
            DO l = 1,lmaxx - 1
               DO kj = 1,kjmax
                  plegend(kj,l+1) = fleg1(l)*plegend(kj,1)*
     +                           plegend(kj,l) -
     +                           fleg2(l)*plegend(kj,l-1)
                  IF (l_socfirst) THEN
c#ifdef CPP_SOCFIRST
                  dplegend(kj,l+1)=REAL(l+1)*plegend(kj,l)+
     +                  plegend(kj,1)*dplegend(kj,l)
c#endif
                  ENDIF
               END DO
            END DO
c--->       loop over equivalent atoms
            n1 = 0
            DO 31 n = 1,ntype

               IF (l_noco) THEN
c--->          pk non-collinear
c--->             set up the spinors of this atom within global
c--->             spin-coordinateframe            
                  chi(1,1) =  exp(-ci*alph(n)/2)*cos(beta(n)/2)
                  chi(1,2) = -exp(-ci*alph(n)/2)*sin(beta(n)/2)
                  chi(2,1) =  exp(ci*alph(n)/2)*sin(beta(n)/2)
                  chi(2,2) =  exp(ci*alph(n)/2)*cos(beta(n)/2)
c--->             and determine the prefactors for the Hamitonian- and
c--->             overlapp-matrix elements
                  chi11 = chi(1,isp)*conjg(chi(1,isp))
                  chi21 = chi(2,isp)*conjg(chi(1,isp))
                  chi22 = chi(2,isp)*conjg(chi(2,isp))
c--->             determine the prefactors of the Hamiltonian elements
c--->             that are off-diag in the local spin (isp).
c--->             (only in a constrained calculation)
                  IF (l_constr) THEN
                     nsp = 3 - isp
                     chj(isp,1,1,n) = chi(1,nsp)*conjg(chi(1,isp))
                     chj(isp,2,1,n) = chi(2,nsp)*conjg(chi(1,isp))
                     chj(isp,2,2,n) = chi(2,nsp)*conjg(chi(2,isp))
                  ENDIF
                  IF ((isp.EQ.2).AND.l_socfirst) THEN
c                     chis(1,1) =  exp(-ci*phi/2)*cos(theta/2)
c                     chis(1,2) = -exp(-ci*phi/2)*sin(theta/2)
c                     chis(2,1) =  exp(ci*phi/2)*sin(theta/2)
c                     chis(2,2) =  exp(ci*phi/2)*cos(theta/2)
                     isigma_x=MATMUL(conjg(TRANSPOSE(chi)),
     +                    MATMUL(isigma(:,:,1),chi))
                     isigma_y=MATMUL(conjg(TRANSPOSE(chi)),
     +                    MATMUL(isigma(:,:,2),chi))
                     isigma_z=MATMUL(conjg(TRANSPOSE(chi)),
     +                    MATMUL(isigma(:,:,3),chi))
                     DO j1=1,jspins
                        DO j2=1,jspins
                           chi11so(j1,j2)=chi(1,j1)*conjg(chi(1,j2))
                           chi21so(j1,j2)=chi(2,j1)*conjg(chi(1,j2))
                           chi22so(j1,j2)=chi(2,j1)*conjg(chi(2,j2))
                           DO kj = 1,kjmax
                              angso(kj,j1,j2)=
     +                            isigma_x(j1,j2)*cross_k(kj,1)+
     +                            isigma_y(j1,j2)*cross_k(kj,2)+
     +                            isigma_z(j1,j2)*cross_k(kj,3)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
c---> pk non-collinear

               n0 = n1 + 1
               n1 = n1 + neq(n)
               lwn = lmax(n)
               DO kj = 1,kjmax
                  rph(kj,1) = 0.0
                  cph(kj,1) = 0.0
               END DO
               DO nn = n0,n1
                  t1nn = tpi*taual(1,nn)
                  t2nn = tpi*taual(2,nn)
                  t3nn = tpi*taual(3,nn)
c--->             set up phase factors
                  DO kj = 1,kjmax
                     rph(kj,1) = rph(kj,1) +
     +                    cos(  
     *                    ( ski1 - (k1(kj,jintsp)+qssbtj(1)) )*t1nn +
     +                    ( ski2 - (k2(kj,jintsp)+qssbtj(2)) )*t2nn +
     +                    ( ski3 - (k3(kj,jintsp)+qssbtj(3)) )*t3nn ) 
#ifndef CPP_INVERSION
c--->                if the system does not posses inversion symmetry
c--->                the complex part of the exponential is needed.
                     cph(kj,1) = cph(kj,1) +
     +                    sin(
     *                    ( k1(kj,jintsp)+qssbtj(1) - ski1 )*t1nn +
     +                    ( k2(kj,jintsp)+qssbtj(2) - ski2 )*t2nn +
     +                    ( k3(kj,jintsp)+qssbtj(3) - ski3 )*t3nn )
#endif
                  END DO
               END DO
#ifdef DISPLAY_INFO
! show phase factors, i.e. (rph, cph)
      print '(A, 2I4)', '(rph,cph) for ki - kj=1..ki. ki=', ki
      print '(2G10.3)', ( (rph(xxx,1), cph(xxx,1)), xxx=1,kjmax )
#endif
#ifdef DUMP_PHASES
! save rph, cph as complex numbers
      write(misc_name, '(A, "_", I3.3, ".bin")') 'phases', ki
      call write_cvec(rph(1:kjmax, 1), cph(1:kjmax, 1), misc_name)
#endif

c--->          update overlap and l-diagonal hamiltonian matrix
               DO 30 l = 0,lwn
                  IF (l_constr.or.l_socfirst) THEN
                     fjkiln = fj(ki,l,n,isp)
                     gjkiln = gj(ki,l,n,isp)
                  ELSE
                     fjkiln = fj(ki,l,n,iintsp)
                     gjkiln = gj(ki,l,n,iintsp)
                  ENDIF
c
                  w1 = 0.5 * ( uds(l,n,isp)*dus(l,n,isp) + 
     +                         us(l,n,isp)*duds(l,n,isp) )
                  apw_lo1 = fl2p1(l) * 0.5 * rmt(n)**2 * ( gjkiln * w1 +
     +                       fjkiln * us(l,n,isp) * dus(l,n,isp) )
                  apw_lo2 = fl2p1(l) * 0.5 * rmt(n)**2 * ( fjkiln * w1 +
     +                     gjkiln * uds(l,n,isp) * duds(l,n,isp) )
c
                  ddnln =  ddn(l,n,isp)
                  elall = el(l,n,isp,nw)

                  IF ( l_noco .AND. (.NOT. l_ss) ) THEN
c--->             pk non-collinear
#ifndef CPP_INVERSION
                     IF (l_constr.or.l_socfirst) THEN
                       DO kj = 1,ki
                        fct  = plegend(kj,l)*fl2p1(l)*
     +                         ( fjkiln*fj(kj,l,n,isp) +
     +                           gjkiln*gj(kj,l,n,isp)*ddnln )

                        bbwa(kj) = cmplx(rph(kj,1),cph(kj,1))*fct
                        aawa(kj) = cmplx(rph(kj,1),cph(kj,1)) * ( 
     +                             fct*elall + plegend(kj,l)*fl2p1bt(l)*
     +                             ( fjkiln*gj(kj,l,n,isp) +
     +                               gjkiln*fj(kj,l,n,isp) ) )
#ifdef CPP_APW
                        capw1 = cmplx(rph(kj,1),cph(kj,1))*plegend(kj,l)
     +                                 * ( apw_lo1 * fj(kj,l,n,isp) +
     +                                     apw_lo2 * gj(kj,l,n,isp) )
                        aawa(kj) = aawa(kj) + capw1
#endif
                       ENDDO
                     ELSE
                       DO kj = 1,ki
                        fct  = plegend(kj,l)*fl2p1(l)*
     +                         ( fjkiln*fj(kj,l,n,jintsp) +
     +                           gjkiln*gj(kj,l,n,jintsp)*ddnln )

                        bbwa(kj) = cmplx(rph(kj,1),cph(kj,1))*fct
                        aawa(kj) = cmplx(rph(kj,1),cph(kj,1)) * (
     +                             fct*elall + plegend(kj,l)*fl2p1bt(l)*
     +                             ( fjkiln*gj(kj,l,n,jintsp) +
     +                               gjkiln*fj(kj,l,n,jintsp) ) )
#ifdef CPP_APW
                        capw1 = cmplx(rph(kj,1),cph(kj,1))*plegend(kj,l)
     +                                 * ( apw_lo1 * fj(kj,l,n,jintsp) +
     +                                     apw_lo2 * gj(kj,l,n,jintsp) )
                        aawa(kj) = aawa(kj) + capw1
#endif
                       ENDDO
                     ENDIF
c+||
                     IF ( kii+1.LE.nv(1) ) THEN
c--->                spin-up spin-up part
                        CALL CPP_BLAS_caxpy(ki,chi11,bbwa,1,bb(iii+1),1)
                        CALL CPP_BLAS_caxpy(ki,chi11,aawa,1,aa(iii+1),1)
c--->                spin-down spin-up part, upper triangle.
c--->                the help array is used to allow vectorization on
c--->                Cray PVP systems. it is mapped onto the hamiltonian
c--->                matrix after the setup is complete.
                        DO kj = 1,ki - 1
                           ii = iii + kj
                           aahlp(ii)=aahlp(ii)+conjg(aawa(kj))*chi21
                           bbhlp(ii)=bbhlp(ii)+conjg(bbwa(kj))*chi21
                        END DO
                     ENDIF
                     IF ( (kii+1.GT.nv(1)).OR.(n_size.EQ.1) ) THEN
                        IF (n_size.EQ.1) THEN
                          ii = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2 
                        ELSE
                          ii = iii
                        ENDIF
c--->                spin-down spin-up part, lower triangle
                        CALL CPP_BLAS_caxpy(ki,chi21,bbwa,1,bb(ii+1),1)
                        CALL CPP_BLAS_caxpy(ki,chi21,aawa,1,aa(ii+1),1)
c--->                spin-down spin-down part
                        ii = ii + nv(1)+nlotot
                        CALL CPP_BLAS_caxpy(ki,chi22,bbwa,1,bb(ii+1),1)
                        CALL CPP_BLAS_caxpy(ki,chi22,aawa,1,aa(ii+1),1)
                     ENDIF
c-||
c--->                when fj and gj are available for both local spins
c--->                (isp), add the contribution from the constraint
c--->                B-field.
                     IF (l_constr .AND. (isp .EQ. 2)) THEN
                      DO nsp = 1,jspins
                       IF (nsp.EQ.1) THEN
                       DO kj = 1,nv(1)
                        aawa(kj) = (-0.5)*cmplx(b_con(1,n),b_con(2,n))*
     *                             cmplx(rph(kj,1),cph(kj,1))*
     *                             plegend(kj,l)*fl2p1(l)*(
     +                           + fj(ki,l,n,2)*fj(kj,l,n,1)*uun21(l,n)
     +                           + fj(ki,l,n,2)*gj(kj,l,n,1)*udn21(l,n)
     +                           + gj(ki,l,n,2)*fj(kj,l,n,1)*dun21(l,n)
     +                           + gj(ki,l,n,2)*gj(kj,l,n,1)*ddn21(l,n)
     +                                                    )
                       ENDDO
                       ELSE
                       DO kj = 1,nv(1)
                        aawa(kj) = (-0.5)*cmplx(b_con(1,n),-b_con(2,n))*
     *                             cmplx(rph(kj,1),cph(kj,1))*
     *                             plegend(kj,l)*fl2p1(l)*(
     +                           + fj(ki,l,n,1)*fj(kj,l,n,2)*uun21(l,n)
     +                           + fj(ki,l,n,1)*gj(kj,l,n,2)*dun21(l,n)
     +                           + gj(ki,l,n,1)*fj(kj,l,n,2)*udn21(l,n)
     +                           + gj(ki,l,n,1)*gj(kj,l,n,2)*ddn21(l,n)
     +                                                    )
                       ENDDO
                       ENDIF
c--->                  spin-up spin-up part
                       ii = (ki-1)*(ki)/2
                       DO kj = 1,ki
                        aa(ii+kj) = aa(ii+kj) + aawa(kj)*chj(nsp,1,1,n)
                       ENDDO
c--->                  spin-down spin-down part
                       ii = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2 + 
     +                       nv(1)+nlotot
                       DO kj = 1,ki
                        aa(ii+kj) = aa(ii+kj) + aawa(kj)*chj(nsp,2,2,n)
                       ENDDO
c--->                  spin-down spin-up part
                       ii = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2
                       DO kj = 1,nv(1)
                        aa(ii+kj) = aa(ii+kj) + aawa(kj)*chj(nsp,2,1,n)
                       ENDDO
                      ENDDO
                     ENDIF

!                 Add spin-orbit coupling
                     IF (isp.EQ.2) THEN
                        IF ((l.GT.0).AND.l_socfirst) THEN

                     DO j1 = 1,jspins
                      DO j2 = 1,jspins
                        DO kj = 1,nv(1)
                         aawa(kj)=cmplx(rph(kj,1),cph(kj,1))*(
     +                 + fj(ki,l,n,j1)*fj(kj,l,n,j2)*rsopp(n,l,j1,j2)
     +                 + fj(ki,l,n,j1)*gj(kj,l,n,j2)*rsopdp(n,l,j1,j2)
     +                 + gj(ki,l,n,j1)*fj(kj,l,n,j2)*rsoppd(n,l,j1,j2)
     +                 + gj(ki,l,n,j1)*gj(kj,l,n,j2)*rsopdpd(n,l,j1,j2))
     +                 *dplegend(kj,l)*fl2p1(l)*angso(kj,j1,j2)
#ifdef SOC_DEBUG
c  For debugging purposes it is possible to perform numerically
c  the summation over m and m' 
c  Warning!!! Works only if there are now equivalent atoms (spherical 
c  harmonics are not rotated). The current implementation is also slow...
                         radso=cmplx(rph(kj,1),cph(kj,1))*(
     +                 + fj(ki,l,n,j1)*fj(kj,l,n,j2)*rsopp(n,l,j1,j2)
     +                 + fj(ki,l,n,j1)*gj(kj,l,n,j2)*rsopdp(n,l,j1,j2)
     +                 + gj(ki,l,n,j1)*fj(kj,l,n,j2)*rsoppd(n,l,j1,j2)
     +                 + gj(ki,l,n,j1)*gj(kj,l,n,j2)*rsopdpd(n,l,j1,j2))      
                         CALL ylm4(l,gk(ki,:,1),ylm1)
                         CALL ylm4(l,gk(kj,:,1),ylm2)
                         angso2=0.0
                         is1=ispjsp(j1)
                         is2=ispjsp(j2)
                         DO m1=-l,l
                            lm1=l*(l+1)+m1+1
                            DO m2=-l,l
                              lm2=l*(l+1)+m2+1
                              angso2=angso2+ylm1(lm1)*conjg(ylm2(lm2))*
c     +                             anglso(theta,phi,l,m1,is1,l,m2,is2)
     +                        anglso(beta(n),alph(n),l,m1,is1,l,m2,is2)
                           ENDDO
                        ENDDO
c                         radso=0.0
                         aawa2(kj)=aawa(kj)
                         aawa(kj)=angso2*radso
c  As gamma=pi is used in the spin rotations in anglso, off-diagonal (in spin)
c  elements of aawa and aawa2 have different sign...
c                         IF (abs(aawa(kj)-aawa2(kj)).GT.0.00001) THEN
c                         WRITE(6,666) ki,kj,j1,j2,l,angso1,angso2,radso
c 666      FORMAT('ang parts',2i4,2i2,i3,2f8.5,1x,2f8.5,1x,2f8.5)
c                         ENDIF
#endif
                        ENDDO
c                         IF (ki.LT.11) THEN
c                         DO kj=1,min(10,nv(1))
c                            WRITE(6,666) ki,kj,j1,j2,chi11so*aawa(kj),
c     +                      aawa2(kj)*chi11so

c                         ENDDO
c                      ENDIF

                      IF (n_size.EQ.1) THEN

c--->                  spin-up spin-up part
                       ii = (ki-1)*(ki)/2
                       DO kj = 1,ki
                        aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi11so(j1,j2)
                       ENDDO
c--->                  spin-down spin-down part
                       ii = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2 +
     +                       nv(1)+nlotot
                       DO kj = 1,ki
                        aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi22so(j1,j2)
                       ENDDO
c--->                  spin-down spin-up part
                       ii = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2
                       DO kj = 1,nv(1)
                        aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi21so(j1,j2)
                       ENDDO

                       ELSE  ! eigenvalue parallelization

                       IF ( kii+1.LE.nv(1) ) THEN
c--->                    spin-up spin-up part
                         CALL CPP_BLAS_caxpy(ki,chi11so(j1,j2),aawa,1,
     <                                                   aa(iii+1),1)

c--->                    spin-down spin-up part, upper triangle.
                         DO kj = 1,ki - 1
                            ii = iii + kj
                            aahlp(ii) = aahlp(ii) +
     +                                  conjg(aawa(kj))*chi21so(j2,j1)
                         END DO
                       ENDIF
                       IF  (kii+1.GT.nv(1)) THEN
                         ii = iii
c--->                    spin-down spin-up part, lower triangle
                         CALL CPP_BLAS_caxpy(ki,chi21so(j1,j2),aawa,1,
     <                                                     aa(ii+1),1)
c--->                    spin-down spin-down part
                         ii = ii + nv(1)+nlotot
                         CALL CPP_BLAS_caxpy(ki,chi22so(j1,j2),aawa,1,
     <                                                     aa(ii+1),1)
                       ENDIF

                      ENDIF ! eigenvalue par.

                    ENDDO ! j2
                 ENDDO    ! j1
              ENDIF       ! ( l > 0 ) & socfirst
           ENDIF          ! (isp = 2)
!               End spin-orbit
c#endif
                  ELSEIF ( l_noco .AND. l_ss ) THEN
                     IF ( iintsp.EQ.2 .AND. jintsp.EQ.1 ) THEN
                        kjmax = nv(1)
                     ELSE
                        kjmax = ki
                     ENDIF
                     DO kj = 1,kjmax
                        fct  = plegend(kj,l)*fl2p1(l)*
     +                         ( fjkiln*fj(kj,l,n,jintsp) +
     +                           gjkiln*gj(kj,l,n,jintsp)*ddnln )

                        bbwa(kj) = cmplx(rph(kj,1),cph(kj,1))*fct
                        aawa(kj) = cmplx(rph(kj,1),cph(kj,1)) * ( 
     +                             fct*elall + plegend(kj,l)*fl2p1bt(l)*
     +                             ( fjkiln*gj(kj,l,n,jintsp) +
     +                               gjkiln*fj(kj,l,n,jintsp) ) )
                     ENDDO
                     IF ( iintsp.EQ.1 .AND. jintsp.EQ.1 ) THEN
c--->                   spin-up spin-up part
                        ii = (ki-1)*(ki)/2
                        DO kj = 1,ki
                           bb(ii+kj) = bb(ii+kj) + bbwa(kj)*chi11
                           aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi11
                        ENDDO
                     ELSEIF ( iintsp.EQ.2 .AND. jintsp.EQ.2 ) THEN
c--->                   spin-down spin-down part
                        ii = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2 +
     +                        nv(1)+nlotot
                        DO kj = 1,ki
                           bb(ii+kj) = bb(ii+kj) + bbwa(kj)*chi22
                           aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi22
                        ENDDO
                     ELSE
c--->                   spin-down spin-up part
                        ii = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2
                        DO kj = 1,nv(1)
                           bb(ii+kj) = bb(ii+kj) + bbwa(kj)*chi21
                           aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi21
                        ENDDO
                     ENDIF
#endif
c--->             pk non-collinear
                  ELSE
                     DO kj = 1,ki
                        fct  = plegend(kj,l)*fl2p1(l)*
     +                     ( fjkiln*fj(kj,l,n,jintsp) +
     +                       gjkiln*gj(kj,l,n,jintsp)*ddnln )

                        ij = iii + kj
#ifdef CPP_INVERSION
                        bb(ij) = bb(ij) + rph(kj,1) * fct
                        aa(ij) = aa(ij) + rph(kj,1) * ( fct * elall +
     +                               plegend(kj,l) * fl2p1bt(l) *
     +                               ( fjkiln*gj(kj,l,n,jintsp) +
     +                                 gjkiln*fj(kj,l,n,jintsp) ) )
c+APW
#ifdef CPP_APW
                        apw1 = rph(kj,1) * plegend(kj,l)  * 
     +                                   ( apw_lo1 * fj(kj,l,n,iintsp) +
     +                                     apw_lo2 * gj(kj,l,n,iintsp) )
                        aa(ij) = aa(ij) + apw1       
#endif
c-APW
#else
                        bb(ij) = bb(ij) + cmplx(rph(kj,1),cph(kj,1))*fct

                        aa(ij) = aa(ij) + cmplx(rph(kj,1),cph(kj,1)) * (
     +                           fct*elall + plegend(kj,l)*fl2p1bt(l) *
     +                                  ( fjkiln*gj(kj,l,n,jintsp) +
     +                                    gjkiln*fj(kj,l,n,jintsp) ) )
#ifdef CPP_APW
                        capw1 = cmplx(rph(kj,1),cph(kj,1))*plegend(kj,l)
     +                                 * ( apw_lo1 * fj(kj,l,n,iintsp) +
     +                                     apw_lo2 * gj(kj,l,n,iintsp) )
                       aa(ij) = aa(ij) + capw1
#endif
#endif
                     END DO
                  ENDIF

c--->          end loop over l
 30            CONTINUE

c--->       end loop over atom types (ntype)
 31         CONTINUE

c--->    end loop over ki
 40      CONTINUE

c---> end loops over interstitial spins
        ENDDO ! jintsp = 1,iintsp
      ENDDO   ! iintsp = 1,nintsp

c---> pk non-collinear
      IF ( l_noco .AND. (.NOT. l_ss) ) THEN
c--->    map the hamiltonian help array onto the hamitonian matrix
         IF (n_size.EQ.1) THEN
            DO ki = 1,nv(1)
               ii = (ki-1)*(ki)/2
               DO kj = 1,ki-1
                  ij = (nv(1)+nlotot+kj-1)*(nv(1)+nlotot+kj)/2 + ki
                  aa(ij) = aa(ij) + aahlp(ii+kj)
                  bb(ij) = bb(ij) + bbhlp(ii+kj)
               ENDDO
            ENDDO
         ELSE
#ifdef CPP_MPI
#if ( defined(CPP_T90) || defined(CPP_HDF) )
            CALL mingeselle(
     >                      SUB_COMM,n_size,n_rank,nv,
     >                      hlpmsize,matsize,aahlp,
     X                      aa)
            CALL mingeselle(
     >                      SUB_COMM,n_size,n_rank,nv,
     >                      hlpmsize,matsize,bbhlp,
     X                      bb)
#endif
#endif
         ENDIF
      ENDIF

      CALL cpu_time(time2)
      cp_sph = cp_sph + time2 - time1
#ifdef DUMP_ENABLE_SPHERICAL
! skip the whole nonspherical loop
      inquire(file=".SPHERICAL_ONLY", exist=spherical_only)
      if (spherical_only) then
        print*, 'DUMP_DATA: skipping nonspherical part'
        goto 77
      endif
#endif
      IF (secvar) GOTO 77
      CALL cpu_time(time1)
      na = 0
      nkvecprevats = 0
      nkvecprevath = 0
      nkvec_sv = 0
      locols = nv(1)
      locolh = nv(1)
#ifdef CPP_MPI
      nc = 0
      k = 0
      DO  i=1+n_rank, nv(1), n_size
        nc = nc + 1
        k = k + n_size*(nc-1) + n_rank + 1
      ENDDO
      iilos = k
      iiloh = k
#else
      iilos = nv(1)* (nv(1)+1)/2 
      iiloh = nv(1)* (nv(1)+1)/2
#endif

      mlo = 1; mlolo = 1
      DO 130 n=1,ntype

         IF (l_noco) THEN
c--->    pk non-collinear
            IF (.NOT. l_ss) THEN
c--->          initialize hamiltonian help array
               DO ii = 1,hlpmsize
                  aahlp(ii)=cmplx(0.0,0.0)
                  bbhlp(ii)=cmplx(0.0,0.0)
               ENDDO
            ENDIF
c--->       set up the spinors of this atom within global
c--->       spin-coordinateframe            
            chi(1,1) =  exp(-ci*alph(n)/2)*cos(beta(n)/2)
            chi(1,2) = -exp(-ci*alph(n)/2)*sin(beta(n)/2)
            chi(2,1) =  exp(ci*alph(n)/2)*sin(beta(n)/2)
            chi(2,2) =  exp(ci*alph(n)/2)*cos(beta(n)/2)
c--->       and determine the prefactors for the Hamitonian- and
c--->       overlapp-matrix elements
            chi11 = chi(1,isp)*conjg(chi(1,isp))
            chi21 = chi(2,isp)*conjg(chi(1,isp))
            chi22 = chi(2,isp)*conjg(chi(2,isp))
c--->    pk non-collinear
         ENDIF

         lwn = lmax(n)

         DO 120 nn = 1,neq(n)
         na = na + 1
         IF (lnonsph(n).LT.0) GOTO 121
         IF ((invsat(na).EQ.0) .OR. (invsat(na).EQ.1)) THEN
            IF (invsat(na).EQ.0) invsfct = 1
            IF (invsat(na).EQ.1) invsfct = 2
            np = invtab(ngopr(na))
            IF (odi%d1) THEN
               np = ods%ngopr(na)
            END IF
            CALL vec_for_lo(
     >                      llod,nlod,ntypd,natd,nvd,jspd,nintsp,nop,na,
     >                      n,np,lmaxd,con1,tpi,eps,l_ss,nlo,llo,invsat,
     >                      lnonsph,nv,mrot,k1,k2,k3,qss,taual,rmt,bmat,
     >                      rk,gk,vk,
     <                      nkvec,kvec)
            DO lo = 1,nlo(n)
               kveclo(nkvec_sv+1:nkvec_sv+nkvec(lo,1)) = 
     +                                            kvec(1:nkvec(lo,1),lo)
               nkvec_sv = nkvec_sv+nkvec(lo,1)
               nkvec(lo,:) = 0
            ENDDO

c--->       loop over interstitial spins
            DO iintsp = 1,nintsp
               IF (iintsp.EQ.1) THEN
                  qssbti(1) = - qss(1)/2
                  qssbti(2) = - qss(2)/2
                  qssbti(3) = - qss(3)/2
               ELSE
                  qssbti(1) = + qss(1)/2
                  qssbti(2) = + qss(2)/2
                  qssbti(3) = + qss(3)/2
               ENDIF

c--->          set up phase factors
               t1nn =  tpi*taual(1,na) 
               t2nn =  tpi*taual(2,na) 
               t3nn =  tpi*taual(3,na) 
               DO k = 1,nv(iintsp)
                  th = ( k1(k,iintsp) + qssbti(1) )*t1nn +
     +                 ( k2(k,iintsp) + qssbti(2) )*t2nn +
     +                 ( k3(k,iintsp) + qssbti(3) )*t3nn
                  rph(k,iintsp) = cos(th)
                  cph(k,iintsp) = -sin(th)
               END DO

c--->          set up the a,b and c  coefficients
c--->          for the local orbitals, if necessary.
               IF (nlo(n).GE.1) THEN
                  enough = .false.
                  CALL setabc1lo(
     >                        ntypd,lmaxd,nlod,
     >                        n,us(0,1,isp),dus(0,1,isp),
     >                        uds(0,1,isp),duds(0,1,isp),ddn(0,1,isp),
     >                        nlo,llo,l_dulo,
     >                        ulos(1,1,isp),dulos(1,1,isp),
     >                        uulon(1,1,isp),dulon(1,1,isp),
     <                        alo1,blo1,clo1)
               ELSE
                  enough = .true.
               END IF

c--->          set up the a and b coefficients

               IF (np.EQ.1) THEN
                  DO k = 1,nv(iintsp)
                     gkrot(k,1) = gk(k,1,iintsp)
                     gkrot(k,2) = gk(k,2,iintsp)
                     gkrot(k,3) = gk(k,3,iintsp)
                  END DO
               ELSE 
                  IF (odi%d1) THEN
                     CALL matmul3r(ods%mrot(:,:,np),bmat,bmrot)
                  ELSE
                    CALL matmul3(mrot(1,1,np),bmat,bmrot)
                  END IF
                  DO k = 1,nv(iintsp)
c-->                 apply the rotation that brings this atom into the
c-->                 representative (this is the definition of ngopr(na))
c-->                 and transform to cartesian coordinates
                     v(1) = vk(k,1,iintsp)
                     v(2) = vk(k,2,iintsp)
                     v(3) = vk(k,3,iintsp)
                     gkrot(k,1) = bmrot(1,1)*v(1) + bmrot(2,1)*v(2) +
     +                            bmrot(3,1)*v(3)
                     gkrot(k,2) = bmrot(1,2)*v(1) + bmrot(2,2)*v(2) +
     +                            bmrot(3,2)*v(3) 
                     gkrot(k,3) = bmrot(1,3)*v(1) + bmrot(2,3)*v(2) +
     +                            bmrot(3,3)*v(3)
                  END DO
               END IF
               DO 20 k = 1,nv(iintsp)
c-->              generate spherical harmonics
                  vmult(1) =  gkrot(k,1)
                  vmult(2) =  gkrot(k,2)
                  vmult(3) =  gkrot(k,3)
                  CALL ylm4(
     >                     lnonsph(n),vmult,
     <                     ylm)
c+lo
                  IF (.NOT.enough) THEN
                     l_lo1=.false. 
                     IF ((rk(k,iintsp).LT.eps).AND.(.not.l_ss)) THEN
                         l_lo1=.true.
                     ELSE
                         l_lo1=.false. 
                     ENDIF
                     CALL abccoflo(
     >                          llod,nlod,lmaxd,ntypd,natd,con1,rmt(n),
     >                          rph(k,iintsp),cph(k,iintsp),ylm,n,
     >                          na,k,nv(iintsp),
     >                          l_lo1,sfp,llo,nlo,invsat,alo1,blo1,clo1,
     X                          nkvec(1,iintsp),cwork,enough,
     <                      alo(-llod,1,1,iintsp),blo(-llod,1,1,iintsp),
     <                      clo(-llod,1,1,iintsp),kvec)
c-lo
                  ENDIF
c-->              synthesize the complex conjugates of a and b
                  DO l = 0,lnonsph(n)
                     ll1 = l* (l+1)
                     DO m = -l,l
                        termr = rph(k,iintsp)* real(ylm(ll1+m+1)) -
     +                          cph(k,iintsp)*aimag(ylm(ll1+m+1))
                        termi = rph(k,iintsp)*aimag(ylm(ll1+m+1)) +
     +                          cph(k,iintsp)* real(ylm(ll1+m+1))
                        IF (l_constr.or.l_socfirst) THEN
                           ar(k,ll1+m,iintsp) = fj(k,l,n,isp)*termr
                           ai(k,ll1+m,iintsp) = fj(k,l,n,isp)*termi
                           br(k,ll1+m,iintsp) = gj(k,l,n,isp)*termr
                           bi(k,ll1+m,iintsp) = gj(k,l,n,isp)*termi
                        ELSE
                           ar(k,ll1+m,iintsp) = fj(k,l,n,iintsp)*termr
                           ai(k,ll1+m,iintsp) = fj(k,l,n,iintsp)*termi
c                           IF (.not.apw(l)) THEN
                             br(k,ll1+m,iintsp) = gj(k,l,n,iintsp)*termr
                             bi(k,ll1+m,iintsp) = gj(k,l,n,iintsp)*termi
c                           ENDIF
                        ENDIF
                     END DO
                  END DO
 20            CONTINUE
c--->       end loop over interstitial spin
            ENDDO

c--->       add the local orbital contribution to the overlap and
c--->       hamiltonian matrix, if they are used for this atom.
            IF (nlo(n).GE.1) THEN
              IF (n_u.GT.0) THEN
                nkvecprevatu = nkvecprevats
                iilou = iilos ; locolu = locols
              ENDIF
              IF ( l_noco .AND. (.not.l_ss) ) THEN
#ifndef CPP_INVERSION
               CALL slomat(
     >              llod,nlod,ntypd,natd,nvd,lmaxd,hlpmsize,
     >              n,na,nv,invsat,fpi,con1,rmt,n_size,n_rank,
     >              ddn(0,1,isp),gk,rph,cph,fj,gj,
     >              kvec,llo,nlo,dulon(1,1,isp),uulon(1,1,isp),
     >              uloulopn(1,1,1,isp),alo1,blo1,clo1,nlotot,l_ss,
     >              ab_dim,1,1,chi11,chi22,chi21,
     X              iilos,locols,nkvecprevats,bbhlp)
               CALL hlomat(
     >              llod,nlod,ntypd,natd,nvd,
     >              hlpmsize,loplod,lmd,lmplmd,n_size,n_rank,
     >              n,na,nv,lnonsph(n),invsat,
     >              ar(1,0,1),br(1,0,1),ai(1,0,1),bi(1,0,1),
     >              tuu(0,n,isp),tud(0,n,isp),tdu(0,n,isp),tdd(0,n,isp),
     >              ind(0,0,n,isp),el(0,n,isp,nw),ddn(0,n,isp),
     >              llo,nlo,alo,blo,clo,
     >              tuulo(0,-llod,mlo,isp),tdulo(0,-llod,mlo,isp),
     >              tuloulo(-llod,-llod,mlolo,isp),lmaxd,
     >              us(0,1,isp),dus(0,1,isp),uds(0,1,isp),duds(0,1,isp),
     >              ulos(1,1,isp),dulos(1,1,isp),rmt(n),nlotot,l_ss,
     >              ab_dim,1,1,chi11,chi22,chi21,
     X              iiloh,locolh,nkvecprevath,aahlp)
#endif
              ELSE
               jd = 1 ; IF (l_noco) jd = isp
               DO iintsp = 1,nintsp
               DO jintsp = 1,nintsp
               CALL slomat(
     >              llod,nlod,ntypd,natd,nvd,lmaxd,matsize,
     >              n,na,nv,invsat,fpi,con1,rmt,n_size,n_rank,
     >              ddn(0,1,isp),gk,rph,cph,fj,gj,
     >              kvec,llo,nlo,dulon(1,1,isp),uulon(1,1,isp),
     >              uloulopn(1,1,1,isp),alo1,blo1,clo1,nlotot,l_ss,
     >              ab_dim,iintsp,jintsp,chi11,chi22,chi21,
     X              iilos,locols,nkvecprevats,bb)
               CALL hlomat(
     >              llod,nlod,ntypd,natd,nvd,
     >              matsize,loplod,lmd,lmplmd,n_size,n_rank,
     >              n,na,nv,lnonsph(n),invsat,ar(1,0,jintsp),
     >              br(1,0,jintsp),ai(1,0,jintsp),bi(1,0,jintsp),
     >              tuu(0,n,jd),tud(0,n,jd),tdu(0,n,jd),tdd(0,n,jd),
     >              ind(0,0,n,jd),el(0,n,isp,nw),ddn(0,n,isp),
     >              llo,nlo,alo,blo,clo,
     >              tuulo(0,-llod,mlo,jd),tdulo(0,-llod,mlo,jd),
     >              tuloulo(-llod,-llod,mlolo,jd),lmaxd,
     >              us(0,1,isp),dus(0,1,isp),uds(0,1,isp),duds(0,1,isp),
     >              ulos(1,1,isp),dulos(1,1,isp),rmt(n),nlotot,l_ss,
     >              ab_dim,iintsp,jintsp,chi11,chi22,chi21,
     X              iiloh,locolh,nkvecprevath,aa)
               ENDDO
               ENDDO
              ENDIF
            END IF
 
            IF (lda_u(n).GE.0.AND.gwc.eq.1) THEN
               IF ( l_noco .AND. (.not.l_ss) ) THEN
#ifndef CPP_INVERSION
                  CALL u_ham(
     >                       jspd,nvd,lmaxd,ntypd,lmd,hlpmsize,llod,
     >                       nlod,jspins,nv,lda_u,isp,n,invsfct,
     >                       ar,ai,br,bi,ddn,vs_mmp,n_u,lmaxb,
     >                       llo,nlo,alo,blo,clo,uulon(1,1,isp),
     >                       dulon(1,1,isp),uloulopn(1,1,1,isp),
     >                       n_size,n_rank,nlotot,l_ss,
     >                       ab_dim,1,1,chi11,chi22,chi21,
     X                       nkvecprevatu,iilou,locolu,aahlp)
#endif
               ELSE
                  DO iintsp = 1,nintsp
                    DO jintsp = 1,iintsp
                      CALL u_ham(
     >                       jspd,nvd,lmaxd,ntypd,lmd,matsize,llod,
     >                       nlod,jspins,nv,lda_u,isp,n,invsfct,
     >                       ar,ai,br,bi,ddn,vs_mmp,n_u,lmaxb,
     >                       llo,nlo,alo,blo,clo,uulon(1,1,isp),
     >                       dulon(1,1,isp),uloulopn(1,1,1,isp),
     >                       n_size,n_rank,nlotot,l_ss,
     >                       ab_dim,iintsp,jintsp,chi11,chi22,chi21,
     X                       nkvecprevatu,iilou,locolu,aa)
                    ENDDO
                  ENDDO
               ENDIF
            ENDIF

c--->       loops over the interstitial spin
            DO iintsp = 1,nintsp
            DO jintsp = 1,iintsp

            jd = 1 ; IF (l_noco) jd = isp
c--->       loop over l',m'
            DO 110 lp = 0,lnonsph(n)
               DO 100 mp = -lp,lp
                  lmp = lp* (lp+1) + mp

c--->             initialize ax and bx
                  DO k = 1,nv(jintsp)
                     axr(k) = 0.0
                     axi(k) = 0.0
                     bxr(k) = 0.0
                     bxi(k) = 0.0
                  END DO
c--->             loop over l,m
                  DO 70 l = 0,lnonsph(n)
                     DO 60 m = -l,l
                        lm = l* (l+1) + m
                        in = ind(lmp,lm,n,jd)
                        IF (in.NE.-9999) THEN
                           IF (in.GE.0) THEN
                              utu(1) =   real(tuu(in,n,jd))*invsfct
                              utu(2) = -aimag(tuu(in,n,jd))*invsfct
                              dtu(1) =   real(tdu(in,n,jd))*invsfct
                              dtu(2) = -aimag(tdu(in,n,jd))*invsfct
                              utd(1) =   real(tud(in,n,jd))*invsfct
                              utd(2) = -aimag(tud(in,n,jd))*invsfct
                              dtd(1) =   real(tdd(in,n,jd))*invsfct
                              dtd(2) = -aimag(tdd(in,n,jd))*invsfct
                           ELSE
                              im = -in
                              utu(1) =  real(tuu(im,n,jd))*invsfct
                              utu(2) = aimag(tuu(im,n,jd))*invsfct
                              dtu(1) =  real(tud(im,n,jd))*invsfct
                              dtu(2) = aimag(tud(im,n,jd))*invsfct
                              utd(1) =  real(tdu(im,n,jd))*invsfct
                              utd(2) = aimag(tdu(im,n,jd))*invsfct
                              dtd(1) =  real(tdd(im,n,jd))*invsfct
                              dtd(2) = aimag(tdd(im,n,jd))*invsfct
                           END IF
c--->    update ax, bx
                           DO 50 k = 1,nv(jintsp)
                              axr(k) = axr(k) +
     +                                  utu(1)*ar(k,lm,jintsp) -
     +                                  utu(2)*ai(k,lm,jintsp) +
     +                                  utd(1)*br(k,lm,jintsp) -
     +                                  utd(2)*bi(k,lm,jintsp)
                              axi(k) = axi(k) -
     +                                  utu(1)*ai(k,lm,jintsp) -
     +                                  utu(2)*ar(k,lm,jintsp) -
     +                                  utd(1)*bi(k,lm,jintsp) -
     +                                  utd(2)*br(k,lm,jintsp)
c                              IF (.not.apw(l)) THEN
                               bxr(k) = bxr(k) +
     +                                  dtu(1)*ar(k,lm,jintsp) -
     +                                  dtu(2)*ai(k,lm,jintsp) +
     +                                  dtd(1)*br(k,lm,jintsp) -
     +                                  dtd(2)*bi(k,lm,jintsp)
                               bxi(k) = bxi(k) -
     +                                  dtu(1)*ai(k,lm,jintsp) -
     +                                  dtu(2)*ar(k,lm,jintsp) -
     +                                  dtd(1)*bi(k,lm,jintsp) -
     +                                  dtd(2)*br(k,lm,jintsp)
c                              ENDIF
   50                      CONTINUE
                        END IF
   60                CONTINUE
   70             CONTINUE
c
c--->             update hamiltonian and overlap matrices
                  nc = 0 
                  IF ( l_noco .AND. (n_size.GT.1) ) THEN
                     nv_tot = nv(1) + nv(2)
                  ELSE
                     nv_tot = nv(iintsp)
                  ENDIF

                  DO 90 kii =  n_rank, nv_tot-1, n_size
                     ki = mod(kii,nv(iintsp)) + 1

c--->                pk non-collinear
                     IF ( l_noco .AND. (.NOT. l_ss) ) THEN
#ifndef CPP_INVERSION
c+||
                      nc = nc + 1
                      ii = nc*(nc-1)/2*n_size-(nc-1)*(n_size-n_rank-1)
                      IF ( (n_size.EQ.1).OR.(kii+1.LE.nv(1)) ) THEN    ! <1||1> part 
                                                                       ! stored in aahlp
                        DO kj = 1,ki
                           aahlp(ii+kj) = aahlp(ii+kj) + cmplx(
     +                          ar(ki,lmp,iintsp)*axr(kj) -
     -                          ai(ki,lmp,iintsp)*axi(kj) +
     +                          br(ki,lmp,iintsp)*bxr(kj) - 
     _                          bi(ki,lmp,iintsp)*bxi(kj) ,
     +                          ar(ki,lmp,iintsp)*axi(kj) +
     +                          ai(ki,lmp,iintsp)*axr(kj) +
     +                          br(ki,lmp,iintsp)*bxi(kj) +
     +                          bi(ki,lmp,iintsp)*bxr(kj) )
                        ENDDO

                      ELSE                    ! components for <2||2> block stored in 
                                                               ! unused part of bbhlp 
                        DO kj = 1,ki
                           bbhlp(kj) = cmplx(
     +                          ar(ki,lmp,iintsp)*axr(kj) -
     -                          ai(ki,lmp,iintsp)*axi(kj) +
     +                          br(ki,lmp,iintsp)*bxr(kj) -
     _                          bi(ki,lmp,iintsp)*bxi(kj) ,
     +                          ar(ki,lmp,iintsp)*axi(kj) +
     +                          ai(ki,lmp,iintsp)*axr(kj) +
     +                          br(ki,lmp,iintsp)*bxi(kj) +
     +                          bi(ki,lmp,iintsp)*bxr(kj) )
                        ENDDO  
                         
c--->                   spin-down spin-down part
                        ij = ii + nv(1)
                        CALL CPP_BLAS_caxpy(ki,chi22,bbhlp,1,aa(ij+1),1)
c--->                   spin-down spin-up part, lower triangle
                        ij =  ii 
                        CALL CPP_BLAS_caxpy(ki,chi21,bbhlp,1,aa(ij+1),1)
                      ENDIF
c-||
                     ELSEIF ( l_noco .AND. l_ss ) THEN
                        IF ( iintsp.EQ.1 .AND. jintsp.EQ.1 ) THEN
c--->                      spin-up spin-up part
                           kjmax = ki
                           chihlp = chi11
                           ii = (ki-1)*(ki)/2
                        ELSEIF ( iintsp.EQ.2 .AND. jintsp.EQ.2 ) THEN
c--->                      spin-down spin-down part
                           kjmax = ki
                           chihlp = chi22
                           ii = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2+ 
     +                           nv(1)+nlotot
                        ELSE
c--->                      spin-down spin-up part
                           kjmax = nv(1)
                           chihlp = chi21
                           ii = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2
                        ENDIF
                            
                        DO kj = 1,kjmax
                           aa(ii+kj) = aa(ii+kj) + chihlp*cmplx(
     +                          ar(ki,lmp,iintsp)*axr(kj) -
     -                          ai(ki,lmp,iintsp)*axi(kj) +
     +                          br(ki,lmp,iintsp)*bxr(kj) - 
     _                          bi(ki,lmp,iintsp)*bxi(kj) ,
     +                          ar(ki,lmp,iintsp)*axi(kj) +
     +                          ai(ki,lmp,iintsp)*axr(kj) +
     +                          br(ki,lmp,iintsp)*bxi(kj) +
     +                          bi(ki,lmp,iintsp)*bxr(kj) )
                        ENDDO
c--->                pk non-collinear
#endif
                     ELSE
                      nc = nc + 1
                      ii = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
                        DO kj = 1,ki
                           ii = ii + 1
#ifdef CPP_INVERSION
                           aa(ii) = aa(ii) + 
     +                          ar(ki,lmp,iintsp)*axr(kj) -
     -                          ai(ki,lmp,iintsp)*axi(kj) +
     +                          br(ki,lmp,iintsp)*bxr(kj) - 
     _                          bi(ki,lmp,iintsp)*bxi(kj)
c                           IF (.not.apw(l)) THEN
c                           ENDIF
#else
                           aa(ii) = aa(ii) + cmplx(
     +                          ar(ki,lmp,iintsp)*axr(kj) -
     -                          ai(ki,lmp,iintsp)*axi(kj) ,
     +                          ar(ki,lmp,iintsp)*axi(kj) +
     +                          ai(ki,lmp,iintsp)*axr(kj) )
c                           IF (.not.apw(l)) THEN
                            aa(ii) = aa(ii) + cmplx(
     +                          br(ki,lmp,iintsp)*bxr(kj) -
     _                          bi(ki,lmp,iintsp)*bxi(kj) ,
     +                          br(ki,lmp,iintsp)*bxi(kj) +
     +                          bi(ki,lmp,iintsp)*bxr(kj) )
c                           ENDIF
#endif
                        ENDDO
                     ENDIF
c--->             end loop over ki
   90             CONTINUE
c--->          end loop over m'
  100          CONTINUE
c--->       end loop over l'
  110       CONTINUE
c--->       end loops over interstitial spin
            ENDDO
            ENDDO
         ENDIF                                ! invsat(na) = 0 or 1
c--->    end loop over equivalent atoms
  120    CONTINUE
         mlo = mlo + nlo(n)
         mlolo = mlolo + nlo(n)*(nlo(n)+1)/2

#ifndef CPP_INVERSION
         IF ( l_noco .AND. (.NOT. l_ss) ) THEN
c--->    pk non-collinear
          IF (n_size.EQ.1) THEN
            DO ki =  1, nv(1)+nlotot
c--->        spin-up spin-up part
               ii = (ki-1)*(ki)/2
               ij = (ki-1)*(ki)/2
               CALL CPP_BLAS_caxpy(ki,chi11,aahlp(ii+1),1,aa(ij+1),1)
               CALL CPP_BLAS_caxpy(ki,chi11,bbhlp(ii+1),1,bb(ij+1),1)
c--->        spin-down spin-down part
               ij = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2+nv(1)+nlotot
               CALL CPP_BLAS_caxpy(ki,chi22,aahlp(ii+1),1,aa(ij+1),1)
               CALL CPP_BLAS_caxpy(ki,chi22,bbhlp(ii+1),1,bb(ij+1),1)
c--->        spin-down spin-up part, lower triangle
               ij = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2
               CALL CPP_BLAS_caxpy(ki,chi21,aahlp(ii+1),1,aa(ij+1),1)
               CALL CPP_BLAS_caxpy(ki,chi21,bbhlp(ii+1),1,bb(ij+1),1)
c--->        spin-down spin-up part, upper triangle.
               DO kj = 1,ki-1
                  ij = (nv(1)+nlotot+kj-1)*(nv(1)+nlotot+kj)/2 + ki
                  aa(ij) = aa(ij) + conjg(aahlp(ii+kj))*chi21
                  bb(ij) = bb(ij) + conjg(bbhlp(ii+kj))*chi21
               ENDDO
            ENDDO
          ELSE 
            DO ii = 1,hlpmsize
                   aa(ii) = aa(ii)+aahlp(ii)*chi11
                aahlp(ii) = conjg(aahlp(ii))*chi21
            ENDDO
#ifdef CPP_MPI
#if ( defined(CPP_T90) || defined(CPP_HDF) )
            CALL mingeselle(
     >                      SUB_COMM,n_size,n_rank,nv,
     >                      hlpmsize,matsize,aahlp,
     X                      aa)
#endif
#endif
          ENDIF
c--->    pk non-collinear
         ENDIF            
#endif            
c--->    label 121: IF (lnonsph(n).LT.0) GOTO 121
  121    CONTINUE
c---> end loop over atom types (ntype)
  130 CONTINUE
      CALL cpu_time(time2)
      cp_nsph = cp_nsph + time2 - time1
      nmat = nmat + nkvecprevats

 77   CONTINUE        ! if secvar

c---> end loop over local spins (isp)
      ENDDO

      DEALLOCATE ( ar,ai,br,bi,axr,axi,bxr,bxi,fj,gj,plegend,
     +             rph,cph,gk,vk,alo,blo,clo,nkvec )
      IF (l_constr) THEN
        DEALLOCATE ( uun21, udn21, dun21, ddn21 )
      ENDIF
      IF ( l_noco .AND. (.NOT. l_ss) ) DEALLOCATE ( aahlp,bbhlp )
      IF (l_socfirst) DEALLOCATE(dplegend,cross_k)

      RETURN
      END SUBROUTINE hssphn
      END MODULE m_hssphn
