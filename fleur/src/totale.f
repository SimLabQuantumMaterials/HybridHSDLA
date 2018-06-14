      MODULE m_totale
      CONTAINS
      SUBROUTINE totale(
     >                  ntypd,memd,ntypsd,n3d,jspd,nmzxyd,n2d, 
     >                  natd,jmtd,nmzd,layerd,nop,nwdd,nlhd,nlod,
     >                  it,seigscv,seigc,ts,vr0,te_vcoul,te_veff,te_exc,
     >                  invs,invs2,film,nq2,nq3,l_f,nvac,jspins,ntype,
     >                  neq,mrot,llh,mlh,nlh,nmem,ngopr,jri,ntypsy,
     >                  dx,tau,pos,zatom,rmt,rmsh,clnu,l_noco,e_u_c,n_u,
     >                  odi,ods,odd,amat,bmat,
     X                  invtab,l_geo,force_old,force)
cf
c     ***************************************************
c     subroutine calculates the total energy of the slab
c                                  c.l.fu
c     ***************************************************
c     single particle energies
c     SEIGC  sum of the eigenvalues of the core states
c            calculated in cdngen.f
c     SEIGSCV  sum of the eigenvalues of the semicore and valence states
c              calculated in fermie.f 
C     TS         : entropy contribution to the free energy
c     SEIGC,SEIGSCV, TS are calculated in fermie.f
c     ***************************************************
c     TE_VCOUL  :   charge density-coulomb potential integral
c     TE_VEFF:   charge density-effective potential integral
c     TE_EXC :   charge density-ex-corr.energy density integral
c                 exchange-correlation energy
c     TE_VCOUL,TE_VEFF,TE_EXC are calculated in vgen.f
c     VMD :   Madelung term
c     ***************************************************
c     TOTE    :   total energy due to all electrons
c     TOTE = SEIGC + SEIGSCV + TE_VCOUL/2 -TE_VEFF + TE_EXC + VMD
c     ***************************************************
c     FREE ENRGY: F = TOTE - TS
c     total electron energy extrapolated for T->0
c     E0 = TOTE - TS/2
c     ***************************************************
c
      USE m_intgr, ONLY : intgr3
      USE m_constants, ONLY : pimach
      USE m_force_a4
      USE m_force_a3
      USE m_forcew
      USE m_loddop
      USE m_od_types, ONLY : od_inp, od_sym, od_dim
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT (IN) :: ntypd,memd,ntypsd,n3d,jspd,nmzxyd,n2d
      INTEGER,INTENT (IN) :: natd,jmtd,nmzd,layerd,nop,nwdd,nlhd,nlod
      INTEGER,INTENT (IN) :: it,nq2,nq3,nvac,jspins,ntype,n_u
      LOGICAL,INTENT (IN) :: invs,invs2,film,l_f,l_noco
      REAL,   INTENT (IN) :: seigc,seigscv,te_vcoul,te_veff,te_exc,ts
      REAL,   INTENT (IN) :: e_u_c
C     ..
C     .. Array Arguments ..
      INTEGER,INTENT (IN) :: mrot(3,3,nop),jri(ntypd),ntypsy(natd)
      INTEGER,INTENT (IN) :: llh(0:nlhd,ntypsd),mlh(memd,0:nlhd,ntypsd)
      INTEGER,INTENT (IN) :: nlh(ntypsd),nmem(0:nlhd,ntypsd),ngopr(natd)
      INTEGER,INTENT (IN) :: neq(ntypd),invtab(nop)
      COMPLEX,INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      REAL,INTENT (IN)    :: dx(ntypd),pos(3,natd),tau(3,nop)
      REAL,INTENT (IN)    :: rmt(ntypd),rmsh(jmtd,ntypd)
      REAL,INTENT (IN)    :: vr0(ntypd),zatom(ntypd)
      REAL,INTENT (INOUT) :: force(3,ntypd,jspd),force_old(3,ntypd)
      LOGICAL,INTENT (IN) :: l_geo(ntypd)
      REAL   ,INTENT(IN) :: amat(3,3),bmat(3,3)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
      TYPE (od_dim), INTENT (IN) :: odd
c+odim
C     ..
C     .. Local Scalars ..
      REAL rhs,tote,totz,vmd,zintn_r,sfp
      INTEGER n,j,nt,iter,i
      CHARACTER*8 iop,dop,name(10)
C     ..
C     .. Local Arrays ..
      REAL dpj(jmtd)
c.....density
      REAL,    ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:)
      COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
c.....potential
      REAL,    ALLOCATABLE :: vr(:,:,:,:),vz(:,:,:)
      COMPLEX, ALLOCATABLE :: vpw(:,:),vxy(:,:,:,:)
C     ..
      sfp = 2. * sqrt( pimach() )
c
      ALLOCATE ( rho(jmtd,0:nlhd,ntypd,jspd),rht(nmzd,2,jspd),
     +           qpw(n3d,jspd),rhtxy(nmzxyd,odi%n2d-1,2,jspd),
     +           vr(jmtd,0:nlhd,ntypd,jspd),vz(nmzd,2,jspd),
     +           vpw(n3d,jspd),vxy(nmzxyd,odi%n2d-1,2,jspd) )
c
      WRITE (6,FMT=8000)
      WRITE (16,FMT=8000)
 8000 FORMAT (/,/,/,5x,'t o t a l  e n e r g y')
c
c      ---> sum of eigenvalues (core, semicore and valence states)
c 
      tote = seigscv + seigc
      WRITE (6,FMT=8010) tote
      WRITE (16,FMT=8010) tote
 8010 FORMAT (/,10x,'sum of eigenvalues =',t40,f20.10)
c
c      ---> add contribution of coulomb potential
c
      tote = tote + 0.5e0*te_vcoul
      WRITE (6,FMT=8020) te_vcoul
      WRITE (16,FMT=8020) te_vcoul
 8020 FORMAT (/,10x,'density-coulomb potential integral =',t40,f20.10)
c
c      ---> add contribution of effective potential
c
      tote = tote - te_veff
      WRITE (6,FMT=8030) te_veff
      WRITE (16,FMT=8030) te_veff
 8030 FORMAT (/,10x,'density-effective potential integral =',t40,f20.10)
c
c      ---> add contribution of exchange-correlation energy
c
      tote = tote + te_exc
      WRITE (6,FMT=8040) te_exc
      WRITE (16,FMT=8040) te_exc
 8040 FORMAT (/,10x,'charge density-ex.-corr.energy density integral=',
     +     t40,f20.10)

c     ----> VM terms
c     ---> reload the density
c
      IF (l_noco) THEN
        nt = 70
        OPEN (nt,file='cdn',form='unformatted',status='old')
      ELSE
        nt = 71
        OPEN (nt,file='cdn1',form='unformatted',status='old')
      ENDIF
      CALL loddop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,nt,natd,neq,
     <            iop,dop,iter,rho,qpw,rht,rhtxy,name)
      CLOSE (nt)
c+for
c     ---> reload the COULOMB potential
c
      OPEN (11,file='potcoul',form='unformatted',status='old')
      REWIND 11
      CALL loddop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >            jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >            nlh,jri,ntypsd,ntypsy,11,natd,neq,
     <            iop,dop,iter,vr,vpw,vz,vxy,name)
      CLOSE (11)
C
C     CLASSICAL HELLMAN-FEYNMAN FORCE
C
      CALL force_a3(
     >              jmtd,memd,nlhd,ntypd,jspd,ntypsd,natd,
     >              jri,mlh,nlh,nmem,ntype,jspins,ntypsy,
     >              zatom,dx,rmt,rmsh,clnu,llh,rho,vr,neq,l_geo,
     X              force)
C
      IF (l_f) THEN
c
C       core contribution to force: needs TOTAL POTENTIAL and core charge
        OPEN (8,file='pottot',form='unformatted',status='old')
        REWIND 8
        CALL loddop(
     >              jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >              jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >              nlh,jri,ntypsd,ntypsy,8,natd,neq,
     <              iop,dop,iter,vr,vpw,vz,vxy,name)
        CLOSE (8)
c
        CALL force_a4(
     >                jmtd,memd,nlhd,ntypd,jspd,ntypsd,natd,
     >                jri,mlh,nlh,nmem,ntype,jspins,ntypsy,
     >                dx,rmsh,clnu,llh,vr,neq,l_geo,
     X                force)
c
      ENDIF
C

c-for
c     ---> add spin-up and spin-down charge density for lh=0
c
      IF (jspins.EQ.2) THEN
         DO 40 n = 1,ntype
            DO 20 i = 1,jri(n)
               rho(i,0,n,1) = rho(i,0,n,1) + rho(i,0,n,jspins)
   20       CONTINUE
   40    CONTINUE
      END IF   
c
c ----> coulomb interaction between electrons and nuclei of different m.t.s
c
       DO 60 n = 1,ntype
         DO 50 j = 1,jri(n)
            dpj(j) = rho(j,0,n,1)/rmsh(j,n)
  50     CONTINUE
         CALL intgr3(dpj,rmsh(1,n),dx(n),jri(n),rhs)
c
         tote = tote - neq(n)*zatom(n)*sfp*rhs/2.
c
         zintn_r = neq(n)*zatom(n)*sfp*rhs/2.
         WRITE (6,FMT=8045) zintn_r
         WRITE (16,FMT=8045) zintn_r
         CALL intgr3(rho(1,0,n,1),rmsh(1,n),dx(n),jri(n),totz)
         vmd = rmt(n)*vr0(n)/sfp + zatom(n) - totz*sfp
         vmd = -neq(n)*zatom(n)*vmd/ (2.*rmt(n))
         WRITE (6,FMT=8050) n,vmd
         WRITE (16,FMT=8050) n,vmd
         tote = tote + vmd
  60  CONTINUE
      IF (n_u.GT.0) THEN
        WRITE ( 6,FMT=8090) e_u_c
        WRITE (16,FMT=8090) e_u_c
        tote = tote - e_u_c             ! gu test
      ENDIF
      WRITE ( 6,FMT=8060) tote
      WRITE (16,FMT=8060) tote

      CALL force_w(
     >             film,ntypd,jspd,natd,nwdd,nop,layerd,nlod,
     >             ntype,jspins,neq,pos,tau,mrot,ngopr,l_f,
     >             invtab,l_geo,force,force_old,tote,amat,bmat,
     >             odi,ods,odd)
c
c     ---> calculate the free energy and the ground state energy,
c          extrapolated for T->0
c
      WRITE ( 6,FMT=8065) ts
      WRITE (16,FMT=8065) ts
      WRITE ( 6,FMT=8070) tote-ts
      WRITE (16,FMT=8070) tote-ts
      WRITE ( 6,FMT=8080) tote-0.5e0*ts
      WRITE (16,FMT=8080) tote-0.5e0*ts
 8060 FORMAT (/,/,' ----> total energy=',t40,f20.10,' htr')
 8050    FORMAT (/,10x,'Madelung term for atom type:',i3,t40,f20.10)
 8045    FORMAT (/,10x,'el.-nucl. inter. diff. m.t.',t40,f20.10)
 8065 FORMAT (/,/,' ---->(tkb*entropy) TS=',t40,f20.10,' htr')
 8070 FORMAT (/,/,' ----> free energy=',t40,f20.10,' htr')
 8080 FORMAT (/,/,'      extrapolation for T->0',
     +      /,' ----> total electron energy=',t40,f20.10,' htr')
 8090 FORMAT (/,/,' ----> correction for lda+U =',t40,f20.10,' htr')

      DEALLOCATE ( rho,rht,qpw,rhtxy,vr,vz,vpw,vxy )
 
      END SUBROUTINE totale
      END MODULE m_totale
