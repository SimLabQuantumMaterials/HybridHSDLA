      MODULE m_cdnovlp
      CONTAINS
      SUBROUTINE cdnovlp(irank,isize,
     >                   memd,nlhd,ntypsd,k1d,k2d,k3d,n2d,ntypd,nop,
     >                   n3d,natd,jmtd,lmaxd,jspd,nmzxyd,nmzd,msh,
     >                   ntype,mx3,nq2,nq3,ntypsy,lmax,nmzxy,
     >                   nmz,nvac,neq,bmat,taual,kv3,sk3,ig,kv2,nstr2,
     >                   nstr,film,zrfs,z1,delz,omtil,rmt,dx,rmsh,
     >                   jri,mlh,clnu,llh,nlh,nmem,mrot,tau,symor,
     >                   jspin,rh,invs,ncst,invtab,
     >                   odi,ods,amat,ig2,sk2,phi2,vol,
     X                   qpw,rhtxy,rho,rht)

c*****************************************************************
c     calculates the overlapping core tail density and adds
c     its contribution to the corresponging components of
c     valence density.      
c
c     OLD VERSION:
c     The previous version to calculate the overlap of the
c     core tails was done in real space:
c     A three dimensional real space grid was set up. On each
c     of these grid points the charge density of all atoms inside
c     the unit cell and neigboring unit cells was calculated.
c     This calculation required a lagrange fit since the core
c     charge is on a radial mesh. The same was done in the vacuum
c     region, except on each z-plane we worked on a two dimension
c     grid. After the charge density was generated on a equidistant
c     grid the charge density was FFTed into G-space. The set up
c     of the charge density on the real space grid is rather time
c     consuming. 3-loops are required for the 3D grid
c                1-loop over all atoms in the unit cells
c                Larange interpolation
c                3D and 2D FFTs
c     In order to save time the non-spherical contributions inside
c     the sphere had been ignored. It turns out that the later
c     approximation is pure in the context of force calculations.
c     
c     PRESENT VERSION:
c     The present version is written from scratch. It is based on the
c     idea that the FFT of an overlap of spherically symmetric
c     charges can be expressed by the product of
c
c     sum_natype{ F(G,ntype) * sum_atom(atype) {S(\vec{G},atom)}}
c
c     of form factor F and structure factor S. The Form factor
c     depends only G while the structure factor depends on \vec{G}
c     and can build up in G-space. F of a gaussian chargedensity can
c     be calculated analytically.

c     The core-tails to the vacuum region are described by an
c     exponentially decaying function into the vacuum:

c     rho(r_||,z,ivac)= sum_n{ rho(n,ivac) * exp(-kappa*(z-z_v))
c                                          * exp(iG(n)r_||) }

c     And the plane waves are expanded into lattice harmonics
c     up to a l_cutoff. Tests of the accuracy inside the sphere
c     have shown that a reduction of the l_cutoff inside the 
c     in order to save time leads to sizable errors and should 
c     be omitted.

c     rho_L(r) =  4 \pi i^l \sum_{g =|= 0}  \rho_int(g) r_i^{2} \times
c                              j_{l} (gr_i) \exp{ig\xi_i} Y^*_{lm} (g)

c     Tests habe shown that the present version is about 10 times
c     faster than the previous one also all nonspherical terms are
c     included up to l=8 and the previous one included only l=0.
c     For l=16 it is still faster by factor 2.

c     coded                  Stefan Bl"ugel, IFF Nov. 1997
c     tested                 RObert Abt    , IFF Dez. 1997
c*****************************************************************
c
      USE m_intgr, ONLY : intgr3,intgz0
      USE m_constants, ONLY : pimach
      USE m_spgrot
      USE m_qpwtonmt
      USE m_cylbes
      USE m_dcylbs
      USE m_rcerf
      USE m_od_cylbes
      USE m_od_chirot
      USE m_od_types, ONLY : od_inp, od_sym
#ifdef CPP_MPI
      USE m_mpi_bc_st
#endif
      IMPLICIT NONE
C
C     .. Parameters ..
      INTEGER    method1, method2
      PARAMETER (method1 = 2, method2 = 1)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT (IN) :: irank,isize
      INTEGER,INTENT (IN) :: memd,nlhd,ntypsd,k1d,k2d,k3d,n2d,ntypd,nop
      INTEGER,INTENT (IN) :: n3d,natd,jmtd,lmaxd,jspd,nmzxyd,nmzd,msh
      INTEGER,INTENT (IN) :: jspin,ntype,mx3,nq2,nq3,nmzxy,nmz,nvac
      REAL,   INTENT (IN) :: z1,delz,omtil
      LOGICAL,INTENT (IN) :: film,zrfs,symor,invs
C     ..
C     .. Array Arguments ..
      COMPLEX,INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER,INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),kv2(2,n2d)
      INTEGER,INTENT (IN) :: nstr2(n2d),ntypsy(natd),mrot(3,3,nop)
      INTEGER,INTENT (IN) :: neq(ntypd),kv3(3,n3d),nstr(n3d),lmax(ntypd)
      INTEGER,INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER,INTENT (IN) :: llh(0:nlhd,ntypsd),nlh(ntypsd),jri(ntypd)
      INTEGER,INTENT (IN) :: ncst(ntypd),invtab(nop)
      REAL,   INTENT (IN) :: bmat(3,3),dx(ntypd),sk3(n3d)
      REAL,   INTENT (IN) :: rmt(ntypd),taual(3,natd),rmsh(jmtd,ntypd)
      REAL,   INTENT (IN) :: tau(3,nop)
      COMPLEX,INTENT (INOUT) :: qpw(n3d,jspd)
      COMPLEX,INTENT (INOUT) :: rhtxy(nmzxyd,odi%n2d-1,2,jspd)
      REAL,   INTENT (INOUT) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL,   INTENT (INOUT) :: rht(nmzd,2,jspd)
      REAL,   INTENT (INOUT) :: rh(msh,ntypd)
c-odim
      REAL,   INTENT (IN) :: vol
      INTEGER,INTENT (IN) :: ig2(n3d)
      REAL,   INTENT (IN) :: amat(3,3),sk2(n2d),phi2(n2d)
c+odim
C     ..
C     .. Local Scalars ..
      COMPLEX czero,sf,carg,value,slope,ci
      REAL    fpi,tpi,pi,sfp
      REAL    ai,ar,a4,dif,dxx,f11,f12,g,gr,gz,qfin,qfout,dtildh,
     +        rkappa,sign,signz,tol_14,x,z,zero,zvac,alpha3,
     +        g2,phi,gamma,qq
      INTEGER ig3,imz,ir,ivac,j,j1,k,kz,k1,k2,l_cutoff,m0,
     +        n,nz,nrz,nat1,nat2,nzvac,n_out_p,nat,
     +        irec2,irec3,irec1,m,gzi
      LOGICAL tail
C     ..
C     .. Local Arrays ..
      COMPLEX qpwc(n3d)
      REAL    acoff(ntypd),alpha(ntypd),qf(n3d),rho_out(2)
      REAL    phas(nop),rhohelp(msh),rat(msh,ntypd)
      INTEGER kr(3,nop),mshc(ntypd)
      REAL    kro(3,ods%nop),fJ(-odi%M:odi%M),dfJ(-odi%M:odi%M)
      COMPLEX phaso(ods%nop)
C     ..
C     .. External Functions ..
      REAL     diflgr
      EXTERNAL diflgr
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cabs,cmplx,int,log,exp,max,min,real,sin,sqrt
C     ..
      DATA  czero /(0.0,0.0)/, zero /0.0/, tol_14 /1.0e-14/
c
c----> Abbreviation
c
c     l_cutoff : cuts off the l-expansion of the nonspherical charge
c                density due to the core-tails of neigboring atoms
c     mshc     : maximal radial meshpoint for which the radial coretail
c                density is larger than tol_14
c     method1  : two different ways to calculate the derivative of the
c                charge density at the sphere boundary.
c                (1) use subroutine diflgr based on lagrange interpol.
c                (2) use two point formular in real space, 
c                    see notes of SB.
c                Tests have shown that (2) is more accurate.
c     method2  : tow different integration routines to calculate form
c                factor of coretails outside the sphere.
c                (1) use subroutine intgrz to integrate the tails from
c                    outside to inside.
c                (2) use subroutine intgr3 to integrate the tails from
c                    muffin-tin radius to outside and include correction
c                    for start up.
c                Tests have shown that (1) is more accurate.
c           
c
      pi  = pimach()
      tpi = 2 * pi
      fpi = 4 * pi
      sfp = sqrt( fpi )
      ci = cmplx(0.0,1.0)
      IF (irank ==0) THEN
c
c----> prepare local array to store pw-expansion of pseudo core charge
c
      DO k = 1 , nq3    
         qpwc(k) = czero
      ENDDO
C 
c----> (1) set up radial mesh beyond muffin-tin radius
c      (2) cut_off core tails from noise 
c
      DO 45 n = 1 , ntype
         IF (ncst(n).GT.0) THEN
           DO  j = 1 , jri(n)
              rat(j,n) = rmsh(j,n)
           ENDDO
           dxx = exp(dx(n))
           DO j = jri(n) + 1 , msh
              rat(j,n) = rat(j-1,n)*dxx
           ENDDO
           DO j = jri(n) - 1 , msh
              rh(j,n) = rh(j,n)/ (fpi*rat(j,n)*rat(j,n))
           ENDDO
           DO j = msh , jri(n) , -1
              IF ( rh(j,n) .gt. tol_14 ) then
                 mshc(n) = j
                 GOTO 45
              END IF
           ENDDO
           mshc(n) = jri(n)
         ENDIF
  45  CONTINUE
c
c-----> the core density inside the spheres is replaced by a
c       gaussian-like pseudo density : n(r) = acoff*exp(-alpha*r*r)
c       acoff and alpha determined to obtain a continous and 
c       differentiable density at the sphere boundary.
c       IF mshc = jri  either core tail too small or no core (i.e. H)
c
      DO 60 n = 1,ntype
         IF ((mshc(n).GT.jri(n)).AND.(ncst(n).GT.0)) THEN

            j1 = jri(n) - 1
            IF ( method1 .EQ. 1) THEN
               dif = diflgr(rat(j1,n),rh(j1,n))
               WRITE (6,FMT=8000) n,rh(jri(n),n),dif
               alpha(n) = -0.5 * dif / ( rh(jri(n),n)*rmt(n) )
            ELSEIF ( method1 .EQ. 2) THEN
               alpha(n) = log( rh(j1,n) / rh(jri(n),n) )
               alpha(n) = alpha(n)
     +                   / ( rmt(n)*rmt(n)*( 1.0-exp( -2.0*dx(n) ) ) )
            ELSE
               WRITE (6,'('' error in choice of method1 in cdnovlp '')')
               STOP 'cdnovlp'
            ENDIF
            acoff(n) = rh(jri(n),n) * exp( alpha(n)*rmt(n)*rmt(n) )
            WRITE (6,FMT=8010) alpha(n),acoff(n)
            DO j = 1,jri(n) - 1
               rh(j,n) = acoff(n) * exp( -alpha(n)*rat(j,n)**2 )
            ENDDO

         ELSE
            alpha(n) = 0.0
         ENDIF
   60 ENDDO
c
 8000    FORMAT (/,10x,'core density and its first derivative ',
     +                 'at sph. bound. for atom type',
     +             i2,' is',3x,2e15.7)
 8010    FORMAT (/,10x,'alpha=',f10.5,5x,'acoff=',f10.5)
c
c=====> calculate the fourier transform of the core-pseudocharge
c
c     qpw(\vec{G}) = Sum_{n} [ F(G,n) * Sum_{atm{n}} S(\vec{G},atm) ]
c                  n = atom_type
c                  F = Formfactor = F_in_sphere + F_outsphere
c                  S = Structure factor
c
      tail = .false.
c
c*****> start loop over the atom type
c
      nat1 = 1
      DO 100 n = 1,ntype
         IF ( ( mshc(n) .GT. jri(n) ).AND.
     +        ( alpha(n) .GT. tol_14 ) )    THEN
c
         n_out_p = mshc(n)-jri(n)+1
c
c (1) Form factor for each atom type
c
         f11 = tpi * rmt(n) * rh(jri(n),n) / alpha(n)
         f12 = acoff(n) * ( pi/alpha(n) ) * sqrt( pi/alpha(n) )
         ar  = sqrt( alpha(n) ) * rmt(n) 
c
         DO 110 k = 1,nq3
            g = sk3(k)
c    first G=0
            IF ( k.eq.1 ) THEN
               ai = zero
c
c ---->     calculate form factor inside the mt-sphere
c           (use analytic integration of gaussian)
c
               qfin = - f11 + f12 * rcerf(ar,ai)
c
c ---->     calculate form factor outside the mt-sphere
c           (do numerical integration of tails)
c
               IF ( method2 .eq. 1) THEN

                  DO ir = 1 , n_out_p
                     j  = jri(n)+ir-1
                     rhohelp(mshc(n)+1-j) =  rat(j,n) * rat(j,n) 
     +                                     * rat(j,n) *  rh(j,n)
                  END DO
                  CALL intgz0(rhohelp,dx(n),n_out_p,qfout,tail)

               ELSE

                  DO ir = 1 , n_out_p
                     j  = jri(n)+ir-1
                     rhohelp(ir) = rat(j,n) * rat(j,n) * rh(j,n)
                  END DO
                  CALL intgr3(rhohelp,rat(jri(n),n),dx(n),
     +                        n_out_p,qfout)
c--->             have to remove the small r-correction from intgr3
                  qfout=qfout-rmt(n)*rhohelp(1)

               END IF

                  qfout = fpi * qfout
c
            ELSE 
c    then  G>0
               ai = 0.5*g/sqrt(alpha(n))
               gr = g*rmt(n)
               a4 = 0.25/alpha(n)
c
c ---->     calculate form factor inside the mt-sphere
c           (use analytic integration of gaussian)
c
               qfin = - f11 * sin(gr)/gr 
     +                + f12 * rcerf(ar,ai) * exp(-a4*g*g) 
c
c ---->     calculate form factor outside the mt-sphere
c           (do numerical integration of tails)
             
               IF ( method2 .eq. 1) THEN

                  DO ir = 1 , n_out_p 
                     j  = jri(n)+ir-1
                     rhohelp(mshc(n)-jri(n)+2-ir) =  rat(j,n)*rat(j,n) 
     +                                     * rh(j,n) * sin( g*rat(j,n) )
                  END DO
c
c--->       note we use here the integration routine for vacuum. Because 
c           the vacuum integration is made for an inwards integration 
c           from outside to inside. Outside the starting value will be 
c           nearly zero since the core density is small. if the tail 
c           correction (tail=.true.) is included for the integrals, the 
c           integrand is from infinity inward. This might not be 
c           necessary. Further the integration routine is made for 
c           equidistant meshpoints, therefore the term r(i) of
c           dr/di = dx*r(i) is included in rhohelp


                  CALL intgz0(rhohelp,dx(n),n_out_p,qfout,tail)

               ELSE

                  DO ir = 1 , n_out_p
                     j  = jri(n)+ir-1
                     rhohelp(ir) = rat(j,n) * rh(j,n) * sin(g*rat(j,n))
                  END DO
                  CALL intgr3(rhohelp,rat(jri(n),n),dx(n),
     +                        n_out_p,qfout)
c--->             have to remove the small r-correction from intgr3
croa...correction.from.intgr3.......................
                  if (rhohelp(1)*rhohelp(2).GT.zero) THEN
                  alpha3 = 1.0 + log(rhohelp(2)/rhohelp(1))/dx(n)
                  IF (alpha3.GT.zero)
     +                 qfout = qfout - rat(jri(n),n)*rhohelp(1)/alpha3
                  endif
croa...end.correction...............................


               END IF

               qfout = fpi * qfout / g
c
            END IF
c
            qf(k)    = (qfin + qfout)/omtil
  110    continue
c
c (2) structure constant for each atom of atom type
c
c
c    first G=0
c
         k=1
                 qpw(k,jspin) = qpw(k,jspin) + neq(n) * qf(k)
                 qpwc(k)      = qpwc(k)      + neq(n) * qf(k)
c
c    then  G>0
c
         DO 120 k = 2,nq3
           IF (.NOT.odi%d1) THEN
                 CALL spgrot(
     >                       nop,symor,tpi,mrot,tau,invtab,
     >                       kv3(1,k),
     <                       kr,phas)
c
c ----> start loop over equivalent atoms
c
                 nat2 = nat1 + neq(n) - 1
                 DO 140 nat = nat1,nat2
                    sf = czero
                    DO 160 j = 1,nop
                       x = -tpi* ( kr(1,j) * taual(1,nat)
     +                           + kr(2,j) * taual(2,nat)
     +                           + kr(3,j) * taual(3,nat))
                       sf = sf + cmplx(cos(x),sin(x))*phas(j)
  160               CONTINUE
                    sf = sf / real( nop )
                    qpw(k,jspin) = qpw(k,jspin) + sf * qf(k)
                    qpwc(k)      = qpwc(k)      + sf * qf(k)
  140            CONTINUE
             ELSE
c-odim
                CALL od_chirot(odi,ods,bmat,kv3(1,k),kro,phaso)
                nat2 = nat1 + neq(n) - 1
                DO 141 nat = nat1,nat2
c                  sf = cmplx(1.,0.)
                   sf = czero
                   DO 161 j = 1,ods%nop
                      x = -tpi* ( kro(1,j)*taual(1,nat)
     +                          + kro(2,j)*taual(2,nat)
     +                          + kro(3,j)*taual(3,nat))
                      sf = sf + cmplx(cos(x),sin(x))*phaso(j)
 161               CONTINUE
                   sf = sf / real( ods%nop )
                   qpw(k,jspin) = qpw(k,jspin) + sf * qf(k)
                   qpwc(k)      = qpwc(k)      + sf * qf(k)
 141            CONTINUE
c+odim
             END IF           
  120      CONTINUE
         END IF
         nat1 = nat1 + neq(n)
  100 CONTINUE
c
c=====> calculate core-tails to the vacuum region                
c       Coretails expanded in exponentially decaying functions.
c       Describe vacuum by: 
c       rho(r_||,z,ivac)= sum_n{ rho(n,ivac) * exp(-kappa*(z-z_v)) 
c                                           * exp(iG(n)r_||) }
      IF (film .AND. .NOT.odi%d1) THEN
c+gu
         dtildh = 0.5 * tpi / bmat(3,3)
         IF (nvac.EQ.1) THEN
           rho_out(1) = qpwc(1)*z1
           DO k = 2,nq3
             IF ((kv3(1,k).EQ.0).AND.(kv3(2,k).EQ.0)) THEN
               nz = 1
               IF (invs.OR.zrfs) nz = 2
               g = kv3(3,k) * bmat(3,3)
               rho_out(1) = rho_out(1) + nz*qpwc(k)*sin(g*z1)/g
             ENDIF
           ENDDO
           rho_out(1) =  qpwc(1) * dtildh - rho_out(1)
         ELSE
           DO ivac = 1, nvac
             carg = cmplx(0.0,0.0)
             DO k = 2,nq3
               IF ((kv3(1,k).EQ.0).AND.(kv3(2,k).EQ.0)) THEN
                 g = kv3(3,k) * bmat(3,3) * (3. - 2.*ivac)
                 carg = carg -qpwc(k)*(exp(ci*g*dtildh)-exp(ci*g*z1))/g
               ENDIF
             ENDDO
             rho_out(ivac) = qpwc(1) * ( dtildh-z1 ) - aimag(carg)
           ENDDO
         ENDIF
c-gu
c        nzvac = min(50,nmz)
         m0 = -mx3
         IF (zrfs) m0 = 0
c
c---> loop over 2D stars
c
         DO 280 k = 1,nq2
            k1 = kv2(1,k)
            k2 = kv2(2,k)
            DO 270 ivac = 1,nvac
               value = czero
               slope = czero
               sign = 3. - 2.*ivac
c
c ---> sum over gz-stars
               DO 250 kz = m0,mx3
                 ig3 = ig(k1,k2,kz)
c        ----> use only stars within the g_max sphere (oct.97 shz)
                 IF (ig3.ne.0) THEN
                   nz = 1
                   IF (zrfs) nz = nstr(ig3)/nstr2(k)
                   gz = kz*bmat(3,3)
                   DO 240 nrz = 1,nz
                      signz = 3. - 2.*nrz
                      carg = ci*sign*signz*gz
                      value = value + qpwc(ig3)* exp(carg*z1)
                      slope = slope + carg*qpwc(ig3)* exp(carg*z1)
  240              ENDDO
                 END IF
  250          ENDDO
c roa work-around
               IF (  abs(real(value)).GT.zero ) THEN
c roa work-around
c gb works also around
               rkappa = - real( slope/value )
               IF (k.EQ.1) rkappa = value/rho_out(ivac)
c               rkappa = - sign * real( slope/value )
               IF (rkappa.LE.zero) rkappa=min(rkappa,-tol_14)
               IF (rkappa.GT.zero) rkappa=max(rkappa,tol_14)
c gb works also around
               zvac   = - log( tol_14/cabs(value) ) / rkappa
               nzvac  = int( zvac/delz ) + 1
c               IF ( rkappa.GT.zero .AND. real(value).GT.zero ) THEN
               IF ( rkappa.GT.zero ) THEN
                  z = 0. 
                  IF ( k.EQ.1 ) THEN
                     DO imz = 1 , min( nzvac,nmz )
                        rht(imz,ivac,jspin) = 
     +                  rht(imz,ivac,jspin) + value*exp(-rkappa*z)
                        z = z + delz
  220                ENDDO
                  ELSE
                     DO imz = 1 , min( nzvac,nmzxy )
                        rhtxy(imz,k-1,ivac,jspin) = 
     +                  rhtxy(imz,k-1,ivac,jspin) + value*exp(-rkappa*z)
                        z = z + delz
  230                ENDDO
                  END IF
               END IF              
c roa work-around
               END IF              
c roa work-around
  270       ENDDO
  280    ENDDO
      ELSEIF (odi%d1) THEN
c-odim
c--->  rho_out is the charge lost between the interstitial and the
c--->  rectangular unit cell

         rho_out(1) = vol*qpwc(1)

         DO k = 2,nq3
            IF (kv3(3,k).EQ.0) THEN
               irec3 = ig(kv3(1,k),kv3(2,k),kv3(3,k))
               IF (irec3.NE.0) THEN
                  irec2 = ig2(irec3)
                  IF (irec2.NE.1) THEN
                     g2 = sk2(irec2)
                     CALL cylbes(odi%M,g2*z1,fJ)
                     rho_out(1) = rho_out(1) +
     +                    qpwc(k)*2.*vol*fJ(1)/(z1*g2)
                  END IF
               END IF
            END IF
         END DO

         rho_out(1) = bmat(3,3)*(qpwc(1)*omtil - rho_out(1))/(tpi*tpi)
c     then we are constructing our radial components of the vacuum
c     charge density so that the so that they are continuous in the
c     value and slope on the boundary
c     value = sum{gpar}[qpw(gz,gpar)exp{-i*m*phi(gpar))J_m(gpar*z1)]
c     slope = .......the same...................*gpar*dJ_m(gpar*z1)]
c     for determining the rht we need only continuity, because the
c     second condition is the charge neutrality of the system
c                         Y.Mokrousov
         DO irec1 = 1,odi%nq2
            value = czero
            slope = czero
            gzi = odi%kv(1,irec1)
            m = odi%kv(2,irec1)
            DO irec2 = 1,nq2
               irec3 = ig(kv2(1,irec2),kv2(2,irec2),gzi)
               IF (irec3.NE.0) THEN
                  g2 = sk2(irec2)
                  phi = phi2(irec2)
                  CALL cylbes(odi%M,g2*z1,fJ)
                  CALL dcylbs(odi%M,g2*z1,fJ,dfJ)
                  value = value + (ci**m)*qpwc(irec3)*
     *                 exp(cmplx(0.,-m*phi))*fJ(m)
                  slope = slope + (ci**m)*g2*qpwc(irec3)*
     *                 exp(cmplx(0.,-m*phi))*dfJ(m)
               END IF
            END DO

            IF (abs(real(value)).GT.zero) THEN
               IF (irec1.EQ.1) THEN
                  qq = real(value/rho_out(1))
                  rkappa = (z1*qq+sqrt(z1*z1*qq*qq + 4*qq))/2.
                  gamma = rkappa
               ELSE
                   rkappa = gamma
c                  rkappa = -real(slope/value)
               END IF
               IF (rkappa.LE.zero) rkappa=min(rkappa,-tol_14)
               IF (rkappa.GT.zero) rkappa=max(rkappa,tol_14)
               zvac=-log(tol_14/cabs(value))/rkappa
               nzvac=int(zvac/delz)+1
               IF (rkappa.GT.zero .AND. real(value).GT.zero) THEN
                  z = 0.0
                  IF (irec1.EQ.1) THEN
                     DO imz = 1,min(nzvac,nmz)
                        rht(imz,nvac,jspin)=rht(imz,nvac,jspin)+
     +                       value*exp(-rkappa*z)
                        z = z + delz
                     END DO
                  ELSE
                     DO imz = 1,min(nzvac,nmzxy)
                           rhtxy(imz,irec1-1,nvac,jspin)=
     +                          rhtxy(imz,irec1-1,nvac,jspin)+
     +                          value*exp(-rkappa*z)
                        z = z + delz
                     END DO
                  END IF
               END IF
            END IF
         END DO
c+odim

      END IF
c
c=====> update density inside the spheres                        
c
c ----> (1) subtract on-site contribution, because 
c           they are contained in the plane wave part 
c
      DO n = 1,ntype
         IF  ((mshc(n).GT.jri(n)).AND.(ncst(n).GT.0)) THEN
            DO j = 1,jri(n)
               rho(j,0,n,jspin) = rho(j,0,n,jspin)
     +                          - sfp*rat(j,n)*rat(j,n)*rh(j,n)
            ENDDO
         ENDIF
      ENDDO
c
c ----> (2) add the plane wave contribution of (core tails + on-site 
c           contribution) to the m.t. density, include full nonspherical 
c           components
c
      l_cutoff = 0
      DO n = 1,ntype
         l_cutoff = max( l_cutoff, lmax(n) )
      END DO
      ENDIF ! irank ==0

#ifdef CPP_MPI
      IF ( isize > 1 ) THEN
      CALL mpi_bc_st(
     >               irank,memd,nlhd,ntypsd,jmtd,ntypd,n3d,
     X               jspin,l_cutoff,nq3,ntype,nop,natd,symor,
     X               clnu,qpwc,lmax,ntypsy,jri,nmem,nlh,mlh,nstr,
     X               neq,kv3,mrot,invtab,llh,bmat,tau,taual,dx,
     X               rmsh,sk3)
      ENDIF
#endif

      CALL qpw_to_nmt(
     >                memd,nlhd,ntypsd,jmtd,ntypd,n3d,jspd,lmaxd,
     >                lmax,ntypsy,jri,dx,rmsh,ntype,nop,natd,
     >                symor,bmat,tau,taual,neq,kv3,mrot,invtab,
     >                nmem,nlh,mlh,llh,clnu,odi,ods,irank,isize,
     >                jspin,l_cutoff,nq3,nstr,sk3,fpi,sfp,qpwc,
     X                rho)

#ifdef CPP_MPI
      IF ( isize > 1 ) THEN
       CALL mpi_col_st(                        ! Collect rho()
     >                 jmtd,nlhd,ntypd,irank,
     X                 rho(1,0,1,jspin))
      ENDIF
#endif
 
      END SUBROUTINE cdnovlp
      END MODULE m_cdnovlp
