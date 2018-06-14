      MODULE m_pwden
      CONTAINS
      SUBROUTINE pwden(
     >     nvd,nbasfcn,neigd,n3d,nkptd,jspd,nobd,k1d,k2d,k3d,
     >     kq1d,kq2d,kq3d,mx1,mx2,mx3,ig_old,rgphs,jspins,isize,
     >     ikpt,jspin,l_noco,l_ss,kw1,kw2,kw3,ne,nv,nstr,kq1,kq2,kq3,
     >     nq3_fft,kmxq_fft,igq_fft,igfft,pgfft,omtil,we,z,
     >     ntypd,natd,nop,ntype,z1,vol,volint,l_f,nlotot,invtab,
     >     symor,cdinf,dos,vacdos,kv3,mrot,neq,ig,ig2,bbmat,
     >     taual,volmts,tau,bmat,rmt,sk3,eig,bkpt,lmax,pos,lmaxd,odi,
     X     qpw,cdom,
     <     qis,force,f_b8)
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c     In this subroutine the star function expansion coefficients of
c     the plane wave charge density is determined.
c
c     This subroutine is called for each k-point and each spin.
c
c
C     Two methods are implemented to calculate the charge density
C     1) which uses the FFT. The effort in calculating the charge
C        density is proportional to M * N * log(N) , M being number of
C        states and N being number of plane waves. This is the method
C        which we use for production runs
C     2) the traditional method for calculating the charge density
C        using the double summation. In this case the effort scales as
C        M * N * N. The method is only used for test purposes or for
C        special cases.
c
c
C     INPUT:    eigen vectors
C               reciprocal lattice information
C               Brillouine zone sampling
C               FFT information
C
C     OUTPUT:   qpw(s)
C               1) using FFT
C
C                2) traditional method
C
C                             -1             ef
C                qpw  (g) = vol * sum{ sum{ sum{ sum{ w(k) * f(nu) *
C                                  sp   k    nu   g'
C                                     *
C                                    c(g'-g,nu,k) * c(g',nu,k) } } } }
C                or :
C                             -1             ef
C                qpw  (g) = vol * sum{ sum{ sum{ sum{ w(k) * f(nu) *
C                                  sp   k    nu   g'
C                                     *
C                                    c(g',nu,k) * c(g'+g,nu,k) } } } }
C
c                qpw(g) are actuall 
C 
c                the weights w(k) are normalized: sum{w(k)} = 1
C                                                  k                -6
C                         a) 1                           for kT < 10
C                f(nu) = {                           -1             -6
C                         b){ 1+exp(e(k,nu) -ef)/kt) }   for kt >=10
C
C
C                                      Stefan Bl"ugel, JRCAT, Feb. 1997
C                                      Gustav Bihlmayer, UniWien       
C
c     In non-collinear calculations the density becomes a hermitian 2x2
c     matrix. This subroutine generates this density matrix in the 
c     interstitial region. The diagonal elements of this matrix 
c     (n_11 & n_22) are stored in qpw, while the real and imaginary part
c     of the off-diagonal element are store in cdom. 
c
c     Philipp Kurz 99/07
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c
#include"cpp_double.h"
      USE m_constants, ONLY : pimach
      USE m_forceb8
      USE m_pwint
      USE m_pwdenold
      USE m_fleurenv
      USE m_rfft
      USE m_set
      USE m_od_types, ONLY : od_inp

      IMPLICIT NONE
c
c-----> parameters and info for pwden_old:
c
      INTEGER, INTENT (IN) :: nvd,nbasfcn,neigd,n3d,nkptd,jspd,nobd
      INTEGER, INTENT (IN) :: kq1d,kq2d,kq3d,k1d,k2d,k3d
      INTEGER, INTENT (IN) :: mx1,mx2,mx3,nlotot
      INTEGER, INTENT (IN) :: ig_old(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (IN) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      LOGICAL, INTENT (IN) :: l_noco,l_ss
C
C----->  FFT  INFORMATION
C
      INTEGER :: ifftq2d,ifftq3d
      INTEGER, INTENT (IN) :: nq3_fft,kmxq_fft,kq1,kq2,kq3
      INTEGER, INTENT (IN) :: igq_fft(0:kq1d*kq2d*kq3d-1)
      INTEGER, INTENT (IN) :: igfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1,1)
      REAL,    INTENT (IN) :: pgfft(0:(2*k1d+1)*(2*k2d+1)*(2*k3d+1)-1)
C
C----->  INFORMATION FOR PWINT
C
      INTEGER, INTENT (IN) :: ntypd,natd,nop
      INTEGER, INTENT (IN) :: ntype,lmax(ntypd),lmaxd
      REAL,    INTENT (IN) :: z1,vol,volint
      LOGICAL, INTENT (IN) :: symor,cdinf,dos,vacdos,l_f

      INTEGER, INTENT (IN) :: kv3(3,n3d),mrot(3,3,nop),neq(ntypd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      REAL,    INTENT (IN) :: taual(3,natd),volmts(ntypd),tau(3,nop)
      REAL,    INTENT (IN) :: bmat(3,3),rmt(ntypd),sk3(n3d)
      REAL,    INTENT (IN) :: bbmat(3,3)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
c+odim
C
C----->  BRILLOUINE ZONE INTEGRATION
C
      REAL we(nobd), omtil
      REAL    eig(neigd)
      REAL bkpt(3),xk(3)
C
C----->  BASIS FUNCTION INFORMATION
C
      INTEGER ne,nv(jspd),kw1(nvd,jspd),kw2(nvd,jspd),kw3(nvd,jspd)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#else
      REAL,    INTENT (IN) :: z(nbasfcn,neigd)
#endif
C
C----->  SYMMETRY INFORMATION
C
      INTEGER, INTENT (IN) :: nstr(n3d),invtab(nop)
      REAL,    INTENT (IN) :: pos(3,natd)
C
C----->  CHARGE DENSITY INFORMATION
C
      INTEGER  isize,ikpt,jspin,jspins
      COMPLEX  qpw(n3d,jspd),cdom(n3d)
      REAL     qis(neigd,nkptd,jspd)
      COMPLEX, INTENT (INOUT) ::  f_b8(3,ntypd)
      REAL,    INTENT (INOUT) :: force(3,ntypd,jspd)
C
C-----> LOCAL VARIABLES
C
      LOGICAL  trdchg
      INTEGER  isn,nu,iv,ir,ik,il,im,in,istr,nw1,nw2,nw3,i,j
      INTEGER  ifftq1,ifftq2,ifftq3
      INTEGER  ist(-1:1)
      INTEGER  idens,ndens,ispin,jkpt,jsp_start,jsp_end
      REAL     q0,q0_11,q0_22,scale,tpi
      REAL     zero,  tol_3, s
      COMPLEX  x
      PARAMETER ( trdchg = .false. )
      PARAMETER ( zero   = 0.00,  tol_3=1.0e-3 )
C
C----->  WORK FILES
C
c-pk
      INTEGER  iv1d(nvd,jspd)
      REAL wtf(nobd),wsave(kq3d+15)
      REAL, ALLOCATABLE :: psir(:),psii(:),rhon(:)
      REAL, ALLOCATABLE :: psi1r(:),psi1i(:),psi2r(:),psi2i(:)
      REAL, ALLOCATABLE :: rhomat(:,:)
      REAL, ALLOCATABLE :: kpsir(:),kpsii(:)
      REAL, ALLOCATABLE :: ekin(:)
      COMPLEX  cwk(n3d),ecwk(n3d)
C
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL     CPP_BLAS_sdot
      EXTERNAL CPP_BLAS_sdot
#else
      COMPLEX  CPP_BLAS_cdotc
      EXTERNAL CPP_BLAS_cdotc
#endif
C
C----->  Intrinsic Functions
C
      intrinsic isign,real,cmplx,aimag
C
c-----> data
c
      DATA ist/1,0,0/

C------->          ABBREVIATIONS
C
C     rhon  : charge density in real space
c     ne    : number of occupied states
c     nv    : number of g-components in eigenstate
c     cv=z  : wavefunction in g-space (reciprocal space)
c     psir   : wavefunction in r-space (real-space)
c     cwk   : complex work array: charge density in g-space (as stars)
c     qpw   : charge density stored as stars
c     trdchg: logical key, determines the mode of charge density
c             calculation: false (default) : fft
c                          true            : double sum over stars
c     we    : weights for the BZ-integration for a particular k-point
c     omtil : volume (slab) unit cell, between -.5*D_tilde and +.5*D_tilde
c     kw1   : reciprocal lattice vectors G=G(k1,k2,k3) for wavefunction
c     kw2   :                             =kw1*a_1 + kw2*a_2 + kw3*a_3
c     kw3   : where a_i= Bravais lattice vectors in reciprocal space
c             kwi, depend on k-point.                            
c     kq1d  : dimension of the charge density FFT box in the pos. domain
c     kq2d  : defined in dimens.f program (subroutine apws). 1,2,3 indicate
c     kq3d  ; a_1, a_2, a_3 directions.
c     kq(i) : i=1,2,3 actual length of the fft-box for which FFT is done.
c     nstr  : number of members (arms) of reciprocal lattice (g) vector
c             of each star.
c     nq3_fft: number of stars in the  charge density  FFT-box
c     ng3   : number of 3 dim. stars in the charge density sphere defined
c             by gmax
c     kmxq_fft: number of g-vectors forming the nq3_fft stars in the
c               charge density sphere 
c     kimax : number of g-vectors forming the ng3 stars in the gmax-sphere
c     iv1d  : maps vector (kw1,kw2,kw3) of wave function into one
c             dimensional vector of cdn-fft box in positive domain.
c     ifftq3d: elements (g-vectors) in the charge density  FFT-box
c     igfft : pointer from the g-sphere (stored as stars) to fft-grid 
c             and     from fft-grid to g-sphere (stored as stars)
c     pgfft : contains the phases of the g-vectors of sph.     
c     isn   : isn = +1, FFT transform for g-space to r-space
c             isn = -1, vice versa
C
      tpi = 2 * pimach()

      IF (l_noco) THEN
         ALLOCATE ( psi1r(0:kq1d*kq2d*kq3d-1),
     +              psi1i(0:kq1d*kq2d*kq3d-1),
     +              psi2r(0:kq1d*kq2d*kq3d-1),
     +              psi2i(0:kq1d*kq2d*kq3d-1),
     +              rhomat(0:kq1d*kq2d*kq3d-1,4) )
      ELSE
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
         ALLOCATE ( psir(-kq1d*kq2d:2*kq1d*kq2d*(kq3d+1)-1),
     +              psii(1),
     +              rhon(-kq1d*kq2d:kq1d*kq2d*(kq3d+1)-1) )
         IF (l_f) ALLOCATE ( kpsii(1),
     +              kpsir(-kq1d*kq2d:2*kq1d*kq2d*(kq3d+1)-1),
     +               ekin(-kq1d*kq2d:2*kq1d*kq2d*(kq3d+1)-1))
#else
         ALLOCATE ( psir(0:kq1d*kq2d*kq3d-1),
     +              psii(0:kq1d*kq2d*kq3d-1),
     +              rhon(0:kq1d*kq2d*kq3d-1) )
         IF (l_f) ALLOCATE ( kpsir(0:kq1d*kq2d*kq3d-1),
     +                       kpsii(0:kq1d*kq2d*kq3d-1),
     +                        ekin(0:kq1d*kq2d*kq3d-1) )
#endif
      ENDIF
c
C========> CALCULATE CHARGE DENSITY IN TRADITIONAL WAY
C
      IF ( trdchg ) THEN
C
C--------> warning when calculating charge density in traditional way
C
            WRITE (6,'('' ===========> warning <========== '',/,
     +               '' calculating charge density in traditional way'',
     +               '' use fft-method to save cpu-time'')')
            IF (l_noco) THEN
               WRITE(6,*)'The interstitial charge density cannot be'
               WRITE(6,*)'calculated using the tradidional method'
               WRITE(6,*)'(without FFT) in a non-collinear calculation.'
               CALL fleur_err(
     +              'pwden: trdchg .AND. l_noco not possible')
            ENDIF
C
C--------> calculate charge density in reciprocal space
C
            CALL pwden_old( 
     >                     nvd,nbasfcn,nobd,neigd,n3d,jspd,
     >                     k1d,k2d,k3d,
     >                     omtil,mx1,mx2,mx3,nstr,rgphs,ig_old,
     >                     we,jspin,z,ne,nv(jspin),
     >                     kw1(1,jspin),kw2(1,jspin),kw3(1,jspin),
     <                     qpw)
C
C<=======  FINISH CALCULATION OF CHARE DENSITY IN TRADITIONAL WAY
C
      ELSE IF ( .not.trdchg ) THEN
C
C=======>  CALCULATE CHARGE DENSITY USING FFT
C 
C
C------> setup FFT
C
            ifftq1  = kq1
            ifftq2  = kq1*kq2
            ifftq3  = kq1*kq2*kq3
            ifftq3d = kq1d*kq2d*kq3d
            ifftq2d = kq1d*kq2d
c
            nw1=NINT(kq1/4.+0.3)
            nw2=NINT(kq2/4.+0.3)
            nw3=NINT(kq3/4.+0.3)
c
c------> g=0 star: calculate the charge for this k-point and spin
c                  analytically to test the quality of FFT
c
            q0 = zero
            q0_11 = zero
            q0_22 = zero
            IF (l_noco) THEN
               q0_11 = zero
               q0_22 = zero
               DO nu = 1 , ne
#if ( !defined(CPP_INVERSION) )
                  q0_11 = q0_11 + we(nu) *
     +             CPP_BLAS_cdotc(nv(1),z(1,nu),1,z(1,nu),1)
                  q0_22 = q0_22 + we(nu) *
     +             CPP_BLAS_cdotc(nv(2),z(nv(1)+nlotot+1,nu),1,
     +                                  z(nv(1)+nlotot+1,nu),1)
#endif
               ENDDO
               q0_11 = q0_11/omtil
               q0_22 = q0_22/omtil
            ELSE
               DO nu = 1 , ne
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
               q0=q0+we(nu)*CPP_BLAS_sdot(nv(jspin),z(1,nu),1,z(1,nu),1)
#else
               q0=q0+we(nu)
     +            *real(CPP_BLAS_cdotc(nv(jspin),z(1,nu),1,z(1,nu),1))
#endif
               ENDDO
               q0 = q0/omtil
            ENDIF
C
C--------> initialize charge density with zero
C
            IF (l_noco) THEN
               DO idens = 1,4
                  DO ir = 0,ifftq3d-1
                     rhomat(ir,idens) = 0.0
                  ENDDO
               ENDDO
               IF (ikpt.LE.isize) THEN
                  DO ispin = 1,jspd
                     DO jkpt = 1,nkptd
                        DO nu = 1,neigd
                           qis(nu,jkpt,ispin) = 0.0
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
               CALL dset(ifftq3d+2*ifftq2d,zero,rhon,1)
               IF (l_f) CALL dset(ifftq3d+2*ifftq2d,zero,ekin,1)
#else
               CALL dset(ifftq3d,zero,rhon,1)
               IF (l_f) CALL dset(ifftq3d,zero,ekin,1)
#endif
            ENDIF
C
C------> calculate:  wtf(nu,k) =  w(k)*f(nu,k)/vol
C
            DO nu = 1,ne
               wtf(nu) = we(nu)/omtil
            ENDDO
c
c------> prepare mapping from wave function box to cdn FFT box
c
            IF (l_ss) THEN
               jsp_start = 1
               jsp_end   = 2
            ELSE
               jsp_start = jspin
               jsp_end   = jspin
            ENDIF
            DO ispin = jsp_start,jsp_end
               DO iv = 1 , nv(ispin)
C                                              -kw1d <= L <= kw1d
C                                              -kw2d <= M <= kw2d
C                                              -kw3d <= N <= kw3d
                  il = kw1(iv,ispin)
                  im = kw2(iv,ispin)
                  in = kw3(iv,ispin)
C
C------>  L,M,N LATTICE POINTS OF G-VECTOR IN POSITIVE DOMAIN
c         (since charge density box = two times charge density box
c          wrap arround error should not occur )
C                                           0<= L <=2*kw1-1 = kq1-1
C                                           0<= M <=2*kw2-1 = kq2-1
C                                           0<= N <=2*kw3-1 = kq3-1
C
                  il = il  +  kq1 * ist( isign(1,il) )
                  im = im  +  kq2 * ist( isign(1,im) )
                  in = in  +  kq3 * ist( isign(1,in) )
C
                  iv1d(iv,ispin) =  in*ifftq2 + im*ifftq1 + il
               ENDDO    
            ENDDO       
          
c
c------------> LOOP OVER OCCUPIED STATES
c
            DO 40 nu = 1 , ne
C
C---> FFT transform c_nu,k(g) --> psi_nu,k(r), for each k-point
c                                              and each nu-state
               IF (l_noco) THEN
                  CALL dset(ifftq3d,zero,psi1r,1)
                  CALL dset(ifftq3d,zero,psi1i,1)
                  CALL dset(ifftq3d,zero,psi2r,1)
                  CALL dset(ifftq3d,zero,psi2i,1)
c------> map WF into FFTbox
#if ( !defined(CPP_INVERSION) )
c the compiler complains about the aimag if z is real, though these
c lines are never reached in a collinear calculation
                  IF (l_ss) THEN
                    DO iv = 1 , nv(1)
                      psi1r( iv1d(iv,1) )   = real( z(iv,nu) )
                      psi1i( iv1d(iv,1) )   = aimag( z(iv,nu) )
                    ENDDO
                    DO iv = 1 , nv(2)
                      psi2r( iv1d(iv,2) ) =  real(z(nv(1)+nlotot+iv,nu))
                      psi2i( iv1d(iv,2) ) = aimag(z(nv(1)+nlotot+iv,nu))
                    ENDDO
                  ELSE
                    DO iv = 1 , nv(jspin)
                      psi1r( iv1d(iv,jspin) ) = real( z(iv,nu) )
                      psi1i( iv1d(iv,jspin) ) = aimag( z(iv,nu) )
                      psi2r(iv1d(iv,jspin))=real( z(nv(1)+nlotot+iv,nu))
                      psi2i(iv1d(iv,jspin))=aimag(z(nv(1)+nlotot+iv,nu))
                     ENDDO
                  ENDIF
#endif
               ELSE
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
                  CALL dset(ifftq3d+ifftq2d,zero,psir(-ifftq2d),1)
                  CALL dset(ifftq3d+2*ifftq2d,zero,psir(ifftq3d),1)
#else
                  CALL dset(ifftq3d,zero,psir,1)
                  CALL dset(ifftq3d,zero,psii,1)
#endif
c------> map WF into FFTbox
                  DO iv = 1 , nv(jspin)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
                     psir( iv1d(iv,jspin) ) = z(iv,nu)
#else
                     psir( iv1d(iv,jspin) ) =  real(z(iv,nu))
                     psii( iv1d(iv,jspin) ) = aimag(z(iv,nu))
#endif
                  ENDDO
               ENDIF
c
c------> do (real) inverse FFT; notice that the array psir is filled from
c        0 to ifftq3-1, but starts at -ifftq2 to give work space for rfft
c
               IF (l_noco) THEN
                  isn = 1

                  CALL cfft(psi1r,psi1i,ifftq3,kq1,ifftq1,isn)
                  CALL cfft(psi1r,psi1i,ifftq3,kq2,ifftq2,isn)
                  CALL cfft(psi1r,psi1i,ifftq3,kq3,ifftq3,isn)
                  
                  CALL cfft(psi2r,psi2i,ifftq3,kq1,ifftq1,isn)
                  CALL cfft(psi2r,psi2i,ifftq3,kq2,ifftq2,isn)
                  CALL cfft(psi2r,psi2i,ifftq3,kq3,ifftq3,isn)
               ELSE
                  isn = 1
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
                  CALL rfft(isn,kq1,kq2,kq3+1,kq1,kq2,kq3,
     <                      nw1,nw2,nw3,wsave,psir(ifftq3d),
     X                      psir(-ifftq2))

c GM force part
                  IF (l_f) THEN
                     DO in=-1,kq3,2
                        DO im=0,ifftq2-1
                           ir = ifftq2 * in + im
                           ekin(ir) = ekin(ir) - wtf(nu) * eig(nu) *
     +                             (psir(ir)**2 + psir(ir+ifftq2)**2)
                        ENDDO
                     ENDDO

                     DO j = 1,3
                      CALL dset(ifftq3d+2*ifftq2d,zero,kpsir(ifftq3d),1)
                       CALL dset(ifftq3d+ifftq2d,zero,kpsir(-ifftq2d),1)
                        DO iv = 1 , nv(jspin)
                           xk(1)=kw1(iv,jspin)+bkpt(1)
                           xk(2)=kw2(iv,jspin)+bkpt(2)
                           xk(3)=kw3(iv,jspin)+bkpt(3)
                           s = 0.0
                           DO i = 1,3
                              s = s + xk(i)*bmat(i,j)
                           ENDDO
                           kpsir( iv1d(iv,jspin) ) = s * z(iv,nu)
                        ENDDO
                        CALL rfft(isn,kq1,kq2,kq3+1,kq1,kq2,kq3,
     <                            nw1,nw2,nw3,wsave,kpsir(ifftq3d),
     X                            kpsir(-ifftq2))
                        DO in=-1,kq3,2
                           DO im=0,ifftq2-1
                              ir = ifftq2 * in + im
                              ekin(ir) = ekin(ir) + wtf(nu) * 0.5 * 
     +                             (kpsir(ir)**2 + kpsir(ir+ifftq2)**2)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDIF
#else
                  CALL cfft(psir,psii,ifftq3,kq1,ifftq1,isn)
                  CALL cfft(psir,psii,ifftq3,kq2,ifftq2,isn)
                  CALL cfft(psir,psii,ifftq3,kq3,ifftq3,isn)
c GM force part
                  IF (l_f) THEN
                     DO ir = 0,ifftq3d-1
                        ekin(ir) = ekin(ir) - wtf(nu)*eig(nu)*
     +                                        (psir(ir)**2+psii(ir)**2)
                     ENDDO

                     DO j = 1,3
                        CALL dset(ifftq3d,zero,kpsir,1)
                        CALL dset(ifftq3d,zero,kpsii,1)
                        DO iv = 1 , nv(jspin)
                           xk(1)=kw1(iv,jspin)+bkpt(1)
                           xk(2)=kw2(iv,jspin)+bkpt(2)
                           xk(3)=kw3(iv,jspin)+bkpt(3)
                           s = 0.0
                           DO i = 1,3
                              s = s + xk(i)*bmat(i,j)
                           ENDDO
                           kpsir( iv1d(iv,jspin) ) = s *  real(z(iv,nu))
                           kpsii( iv1d(iv,jspin) ) = s * aimag(z(iv,nu))
                        ENDDO

                        CALL cfft(kpsir,kpsii,ifftq3,kq1,ifftq1,isn)
                        CALL cfft(kpsir,kpsii,ifftq3,kq2,ifftq2,isn)
                        CALL cfft(kpsir,kpsii,ifftq3,kq3,ifftq3,isn)
   
                           DO ir = 0,ifftq3d-1
                           ekin(ir) = ekin(ir) + wtf(nu) * 0.5 *
     +                                       (kpsir(ir)**2+kpsii(ir)**2)
                        ENDDO
                     ENDDO
                  ENDIF
#endif
               ENDIF
C----> calculate rho(r) = sum w(k)*f(nu)*conjg(psi_nu,k(r))*psi_nu,k(r)
C                         k,nu
C      again, we fill rhon() from -ifftq2 to ifftq3-1 for the rfft
C
               IF (l_noco) THEN
c--->             in the non-collinear case rho becomes a hermitian 2x2
c--->             matrix (rhomat).
                  DO ir = 0,ifftq3d-1
                     rhomat(ir,1) = rhomat(ir,1) 
     +                        + wtf(nu)*( psi1r(ir)**2 + psi1i(ir)**2 )
                     rhomat(ir,2) = rhomat(ir,2) 
     +                        + wtf(nu)*( psi2r(ir)**2 + psi2i(ir)**2 )
                     rhomat(ir,3) = rhomat(ir,3) + wtf(nu)*
     *                        (psi2r(ir)*psi1r(ir)+psi2i(ir)*psi1i(ir))
                     rhomat(ir,4) = rhomat(ir,4) + wtf(nu)*
     *                        (psi2r(ir)*psi1i(ir)-psi2i(ir)*psi1r(ir))
                  ENDDO
c--->             in a non-collinear calculation the interstitial charge
c--->             cannot be calculated by a simple substraction if the
c--->             muffin-tin (and vacuum) charge is know, because the
c--->             total charge does not need to be one in each spin-
c--->             channel. Thus it has to be calculated explicitly, if
c--->             it is needed.
                  IF (dos .OR. vacdos .OR. cdinf) THEN
                     DO ir = 0,ifftq3d-1
                        psi1r(ir) = (psi1r(ir)**2 + psi1i(ir)**2)
                        psi2r(ir) = (psi2r(ir)**2 + psi2i(ir)**2)
                     ENDDO
                     isn = -1
                     CALL dset(ifftq3d,zero,psi1i,1)
                     CALL cfft(psi1r,psi1i,ifftq3,kq1,ifftq1,isn)
                     CALL cfft(psi1r,psi1i,ifftq3,kq2,ifftq2,isn)
                     CALL cfft(psi1r,psi1i,ifftq3,kq3,ifftq3,isn)
                     CALL dset(ifftq3d,zero,psi2i,1)
                     CALL cfft(psi2r,psi2i,ifftq3,kq1,ifftq1,isn)
                     CALL cfft(psi2r,psi2i,ifftq3,kq2,ifftq2,isn)
                     CALL cfft(psi2r,psi2i,ifftq3,kq3,ifftq3,isn)
  
                     CALL cset(nq3_fft,(0.,0.),cwk,1)
                     DO ik = 0 , kmxq_fft - 1
                        cwk(igfft(ik,1))=cwk(igfft(ik,1))+pgfft(ik)*
     +                    cmplx(psi1r(igq_fft(ik)),psi1i(igq_fft(ik)))
                     ENDDO
                     DO istr = 1,nq3_fft
                        CALL pwint(
     >                            k1d,k2d,k3d,n3d,ntypd,natd,nop,invtab,
     >                         odi,ntype,neq,volmts,taual,z1,vol,volint,
     >                           symor,tpi,tau,mrot,rmt,sk3,bmat,ig2,ig,
     >                             kv3(1,istr),
     <                             x)
                        qis(nu,ikpt,1) = qis(nu,ikpt,1)
     +                            + real(cwk(istr)*x)/omtil/real(ifftq3)
                     ENDDO

                     CALL cset(nq3_fft,(0.,0.),cwk,1)
                     DO ik = 0 , kmxq_fft - 1
                        cwk(igfft(ik,1))=cwk(igfft(ik,1))+pgfft(ik)*
     +                    cmplx(psi2r(igq_fft(ik)),psi2i(igq_fft(ik)))
                     ENDDO
                     DO istr = 1,nq3_fft
                        CALL pwint(
     >                            k1d,k2d,k3d,n3d,ntypd,natd,nop,invtab,
     >                         odi,ntype,neq,volmts,taual,z1,vol,volint,
     >                           symor,tpi,tau,mrot,rmt,sk3,bmat,ig2,ig,
     >                             kv3(1,istr),
     <                             x)
                        qis(nu,ikpt,jspins) = qis(nu,ikpt,jspins)
     +                            + real(cwk(istr)*x)/omtil/real(ifftq3)
                     ENDDO
                  ENDIF
               ELSE
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
                  DO in=-1,kq3,2
                     DO im=0,ifftq2-1
                        ir = ifftq2 * in + im
                        rhon(ir) = rhon(ir) + wtf(nu) * ( psir(ir)**2 +
     +                                              psir(ir+ifftq2)**2 )
                     ENDDO
                  ENDDO
#else
                  DO ir = 0,ifftq3d-1
                     rhon(ir)=rhon(ir)+wtf(nu)*(psir(ir)**2+psii(ir)**2)
                  ENDDO
#endif
               ENDIF
c              DO ir = -ifftq2 , ifftq3-1
c     +                      + wtf(nu)*(psi(ir+ifftq3d) * psi(ir+ifftq3d)
c     +                               + psi(ir  ) * psi(ir  )
c     +                                 )
c              ENDDO

 40         CONTINUE
c
C<<<<<<<<<<<<<< END OUTER LOOP OVER STATES NU  >>>>>>>>>>>>>>>>>>
c
C
C----> perform back  FFT transform: rho(r) --> chgn(star)
c        ( do direct FFT)                    = cwk(star)

c--->  In a collinear calculation pwden is calles once per spin.
c--->  However in a non-collinear calculation pwden is only called once
c--->  and all four components of the density matrix (n_11 n_22 n_12
c--->  n_21) have to be calculated at once.
            ndens = 1
            IF (l_noco) ndens = 4
            DO idens = 1,ndens
               IF (l_noco) THEN
                  CALL dset(ifftq3d,zero,psi1r,1)
              
                  isn = -1
                  CALL cfft(rhomat(0,idens),psi1r,ifftq3,kq1,ifftq1,isn)
                  CALL cfft(rhomat(0,idens),psi1r,ifftq3,kq2,ifftq2,isn)
                  CALL cfft(rhomat(0,idens),psi1r,ifftq3,kq3,ifftq3,isn)
               ELSE
c--->  psir is used here as work array, charge is real ,but fft complex
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
                  CALL dset(ifftq3d,zero,psir(ifftq3d),1)
                  IF (l_f) CALL dset(ifftq3d,zero,kpsir(ifftq3d),1)
#else
                  CALL dset(ifftq3d,zero,psir,1)
                  CALL dset(ifftq3d,zero,psii,1)
                  IF (l_f) CALL dset(ifftq3d,zero,kpsir,1)
                  IF (l_f) CALL dset(ifftq3d,zero,kpsii,1)
#endif

                  isn = -1
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
                  CALL rfft(isn,kq1,kq2,kq3+1,kq1,kq2,kq3,
     <                  kq1,kq2,kq3,wsave,psir(ifftq3d),
     X                  rhon(-ifftq2))
                  IF (l_f) CALL rfft(isn,kq1,kq2,kq3+1,kq1,kq2,kq3,
     <                            kq1,kq2,kq3,wsave,kpsir(ifftq3d),
     X                            ekin(-ifftq2))
#else
                  CALL cfft(rhon,psir,ifftq3,kq1,ifftq1,isn)
                  CALL cfft(rhon,psir,ifftq3,kq2,ifftq2,isn)
                  CALL cfft(rhon,psir,ifftq3,kq3,ifftq3,isn)
c+apw
                  IF (l_f) THEN 
                     CALL cfft(ekin,psii,ifftq3,kq1,ifftq1,isn)
                     CALL cfft(ekin,psii,ifftq3,kq2,ifftq2,isn)
                     CALL cfft(ekin,psii,ifftq3,kq3,ifftq3,isn)
                  ENDIF
#endif
               ENDIF
c  ---> collect stars from the fft-grid
c
               CALL cset(nq3_fft,(0.,0.),cwk,1)
               CALL cset(nq3_fft,(0.,0.),ecwk,1)
               IF (l_noco) THEN
                  DO ik = 0 , kmxq_fft - 1
                     cwk(igfft(ik,1))=cwk(igfft(ik,1))+pgfft(ik)*
     +               cmplx(rhomat(igq_fft(ik),idens),psi1r(igq_fft(ik)))
                  ENDDO
               ELSE
                  DO ik = 0 , kmxq_fft - 1
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
                     cwk(igfft(ik,1))=cwk(igfft(ik,1))+pgfft(ik)*
     +                                     cmplx(rhon(igq_fft(ik)),zero)
#else
                     cwk(igfft(ik,1))=cwk(igfft(ik,1))+pgfft(ik)*
     +                        cmplx(rhon(igq_fft(ik)),psir(igq_fft(ik)))
#endif
                  ENDDO
c+apw
                  IF (l_f) THEN 
                    DO ik = 0 , kmxq_fft - 1
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
                     ecwk(igfft(ik,1))=ecwk(igfft(ik,1))+pgfft(ik)*
     +                                     cmplx(ekin(igq_fft(ik)),zero)
#else
                     ecwk(igfft(ik,1))=ecwk(igfft(ik,1))+pgfft(ik)*
     +                        cmplx(ekin(igq_fft(ik)),psii(igq_fft(ik)))
#endif
                    ENDDO
                  ENDIF
c-apw
               ENDIF
c
               scale=1.0/ifftq3
               DO istr = 1 , nq3_fft
                  cwk(istr) = scale * cwk(istr) / real( nstr(istr) )
               ENDDO
#ifdef CPP_APW
               IF (l_f) THEN
                 DO istr = 1 , nq3_fft
                  ecwk(istr) = scale * ecwk(istr) / real( nstr(istr) )
                 ENDDO
                 CALL force_b8(
     >                         ntype,ntypd,natd,jspd,neq,ecwk,nq3_fft,
     >                         tau,rmt,kv3,n3d,bmat,lmax,lmaxd,taual,
     >                         nop,mrot,bbmat,jspin,
     X                         force,f_b8)
               ENDIF
#endif
C
c---> check charge neutralilty
c
               IF ((idens.EQ.1).OR.(idens.EQ.2)) THEN
                  IF (l_noco) THEN
                     IF (idens.EQ.1) THEN
                        q0 = q0_11
                     ELSE
                        q0 = q0_22
                     ENDIF
                  ENDIF
                  IF ( abs( q0 ) .GT. 1.0e-9) THEN
                  IF ( abs( q0 - real(cwk(1)) )/q0 .gt. tol_3 ) then
                     WRITE ( 6,'(''bad quality of charge density'',
     +                    2f13.8)')q0, real( cwk(1) )
                     WRITE (16,'(''bad quality of charge density'',
     +                    2f13.8)')q0, real( cwk(1) )
                     CALL fleur_err('pwden: bad quality of charge')
                  ENDIF
                  ENDIF
               ENDIF
c
c---> add charge density to existing one
c
               IF(idens.LE.2) THEN
c--->       add to spin-up or -down density (collinear & non-collinear)
                  ispin = jspin
                  IF (l_noco) ispin = idens
                  DO istr = 1 , nq3_fft
                     qpw(istr,ispin) = qpw(istr,ispin) + cwk(istr)
                  ENDDO
               ELSE IF (idens.EQ.3) THEN
c--->       add to off-diag. part of density matrix (only non-collinear)
                  DO istr = 1 , nq3_fft
                     cdom(istr) = cdom(istr) + cwk(istr)
                  ENDDO
               ELSE
c--->       add to off-diag. part of density matrix (only non-collinear)
                  DO istr = 1 , nq3_fft
                     cdom(istr) = cdom(istr) + cmplx(0.0,1.0)*cwk(istr)
                  ENDDO
               ENDIF

            ENDDO
C<<<<<<< END LOOP OVER DIFFERENT PARTS OF THE DENSITY MATRIX >>>>>>>
      ELSE
           WRITE (16,'('' error in selecting trdchg '')')
           CALL fleur_err('pwden: error selecting trdch')
      ENDIF

      IF (l_noco) THEN
         DEALLOCATE ( psi1r,psi1i,psi2r,psi2i,rhomat )
      ELSE
         DEALLOCATE ( psir,psii,rhon )
         IF (l_f) DEALLOCATE ( kpsir,kpsii,ekin)
      ENDIF

      END SUBROUTINE pwden
      END MODULE m_pwden
