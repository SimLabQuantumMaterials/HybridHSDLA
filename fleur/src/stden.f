      MODULE m_stden
c     ************************************************************
c     generate flapw starting density by superposition of
c     atomic densities. the non-spherical terms inside
c     the spheres are obtained by a least squares fit
c     and the interstitial and vacuum warping terms are calculated
c     by a fast fourier transform.
c     e. wimmer   nov. 1984       c.l.fu diagonized 1987
c     *************************************************************
      CONTAINS
      SUBROUTINE stden(irank,isize,
     >                 memd,nlhd,ntypsd,k1d,k2d,k3d,n2d,ntypd,nop,
     >                 n3d,natd,jmtd,lmaxd,jspd,nmzxyd,nmzd,msh,
     >                 nlod,nvac,nspd,nstd,film,zrfs,invs,invs2,
     >                 ntype,mx3,nq2,nq3,nmzxy,nmz,z1,delz,omtil,
     >                 rmsh,taual,rmt,sk3,dx,bmat,amat,zatom,pos,
     >                 jri,nlh,llh,nmem,mlh,lmax,nstr,neq,ntypsy,
     >                 nstr2,ig,clnu,tau,mrot,nop2,symor,kv2,sk2,
     >                 namat,name,icorr,total,krla,kv3,ngopr,phi2,
     >                 ig2,vol,volint,volmts,area,jspins,bmu,nwd,
     >                 igrd,ndvgrd,idsprs,isprsv,sprsv,nlo,llo,
     >                 l_dulo,nwdd,ellow,elup,sigma,lepr,invtab,odi,ods)

      USE m_constants, ONLY : pimach
      USE m_enpara,    ONLY : w_enpara
      USE m_xcall,     ONLY : vxcall
      USE m_qsf
      USE m_checkdop
      USE m_cdnovlp
      USE m_wrtdop
      USE m_qfix
      USE m_atom2
      USE m_od_types, ONLY : od_inp, od_sym
      USE m_cylpts

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT (IN) :: irank,isize
      INTEGER,INTENT (IN) :: memd,nlhd,ntypsd,k1d,k2d,k3d,n2d,ntypd,nop
      INTEGER,INTENT (IN) :: n3d,natd,jmtd,lmaxd,jspd,nmzxyd,nmzd,msh
      INTEGER,INTENT (IN) :: icorr,krla,nop2,nvac,nspd,nstd,nwdd
      INTEGER,INTENT (IN) :: igrd,ndvgrd,idsprs,isprsv,jspins,lepr
      INTEGER,INTENT (IN) :: ntype,mx3,nq2,nq3,nmzxy,nmz,nwd,nlod
      LOGICAL,INTENT (IN) :: total,symor,film,zrfs,invs,invs2
      REAL,   INTENT (IN) :: sprsv,z1,delz,omtil,vol,volint,area,sigma
C     ..
C     .. Array Arguments ..
      COMPLEX,INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER,INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER,INTENT (IN) :: nstr2(n2d),ntypsy(natd),ngopr(natd)
      INTEGER,INTENT (IN) :: neq(ntypd),nstr(n3d),lmax(ntypd)
      INTEGER,INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER,INTENT (IN) :: llh(0:nlhd,ntypsd),nlh(ntypsd),jri(ntypd)
      INTEGER,INTENT (IN) :: mrot(3,3,nop),kv3(3,n3d),kv2(2,n2d)
      INTEGER,INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd),invtab(nop)
      REAL,   INTENT (IN) :: tau(3,nop)
      REAL,   INTENT (IN) :: pos(3,natd),zatom(ntypd),volmts(ntypd)
      REAL,   INTENT (IN) :: amat(3,3),bmat(3,3),dx(ntypd),sk3(n3d)
      REAL,   INTENT (IN) :: rmt(ntypd),taual(3,natd),rmsh(jmtd,ntypd)
      REAL,   INTENT (IN) :: bmu(ntypd),ellow(nwdd),elup(nwdd)
      LOGICAL,INTENT (IN) :: l_dulo(nlod,ntypd)
      CHARACTER*2, INTENT (IN) :: namat(0:103)
      CHARACTER*8, INTENT (IN) :: name(10)
c-odim
      REAL,   INTENT (IN) :: sk2(n2d),phi2(n2d)
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      REAL d,del,fix,h,r,rnot,sign,z,sfp,fpi,bm,qdel
      REAL denz1(1),vacpot(1),enmix
      INTEGER i,iter,ivac,iza,j,jr,k,n,n1,npd,ispin,nat
      INTEGER nw,ilo,natot,icorr_dummy,m
      COMPLEX czero
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
      REAL,    ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:),vbar(:,:)
      REAL,    ALLOCATABLE :: xp(:,:),rat(:,:),eig(:,:,:),sigm(:)
      REAL,    ALLOCATABLE :: rh(:,:,:),rh1(:,:,:),rhoss(:,:)
      REAL,    ALLOCATABLE :: ello0(:,:),energy(:,:),vacpar(:)
      INTEGER lnum(nstd,ntypd),nst(ntypd),skiplo(ntypd)
      INTEGER dummy(ntypd),jrc(ntypd)
      LOGICAL lchange(0:lmaxd,ntypd),llochg(nlod,ntypd)
      LOGICAL l_found(0:3),llo_found(nlod),l_enpara
      CHARACTER*8 name_l(10)
C     ..
C     .. External Subroutines ..
      EXTERNAL points,sphpts
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp
C     ..
C     .. Data statements ..
      DATA del/1.e-6/
      PARAMETER (czero=(0.0,0.0))
C     ..
c
      fpi = 4 * pimach()
      sfp = sqrt( fpi )
      IF (jspins > jspd) STOP 'stden: jspins > jspd' 
      ALLOCATE ( rho(jmtd,0:nlhd,ntypd,jspd) )
      IF (irank == 0) THEN

      ALLOCATE ( qpw(n3d,jspd),rhtxy(nmzxyd,odi%n2d-1,2,jspd) )
      ALLOCATE ( xp(3,nspd),rat(msh,ntypd),eig(nstd,jspd,ntypd) )
      ALLOCATE ( rh(msh,ntypd,jspd),rh1(msh,ntypd,jspd) )
      ALLOCATE ( ello0(nlod,ntypd),energy(0:lmaxd,ntypd),vacpar(2) )
      ALLOCATE ( rht(nmzd,2,jspd),vbar(2,ntypd),sigm(nmz) )
      ALLOCATE ( rhoss(msh,jspd) )

c--->    if sigma is not 0.0, then divide this charge among all atoms
      IF ( abs(sigma).lt. 1.e-6) THEN
         qdel = 0.0
      ELSE
         natot = 0
         DO n=1,ntype
            IF (zatom(n).GE.1.0) natot = natot + neq(n)
         ENDDO
         qdel = 2.*sigma/natot
      ENDIF
c
      WRITE (6,FMT=8000)
 8000 FORMAT (/,/,/,' superposition of atomic densities',/,/,
     +       ' original atomic densities:',/)
      DO n = 1,ntype
         r = rmsh(1,n)
         d = exp(dx(n))
!         DO i = 1, msh
!            rat(i,n) = r
!            r = r*d
!         ENDDO
         jrc(n) = 0
         DO WHILE (r < rmt(n) + 20.0) 
            IF ( jrc(n) > msh )  STOP 'increase msh in fl7para!'
            jrc(n) = jrc(n) + 1
            rat(jrc(n),n) = r
            r = r*d
         ENDDO 
         dummy(n) = 1
      ENDDO
c
c Generate the atomic charge densities
c
      DO n = 1,ntype
         z = zatom(n)
         r = rmt(n)
         h = dx(n)
         jr = jri(n)
         IF (jspins.EQ.2) THEN
           bm = bmu(n)
         ELSE
           bm = 0.
         ENDIF 
c--->    check whether this atom has been done already
         DO n1 = 1,n - 1
            IF (abs(z-zatom(n1)).GT.del) GO TO 40
            IF (abs(r-rmt(n1)).GT.del) GO TO 40
            IF (abs(h-dx(n1)).GT.del) GO TO 40
            IF (abs(bm-bmu(n1)).GT.del) GO TO 40
            IF (jr.NE.jri(n1)) GO TO 40
            DO ispin = 1, jspins
              DO i = 1,jrc(n) ! msh
                 rh(i,n,ispin) = rh(i,n1,ispin)
              ENDDO 
            ENDDO
            nst(n) = nst(n1)
            vbar(1,n) = vbar(1,n1)
            vbar(jspins,n) = vbar(jspins,n1)
            DO i = 1, nst(n1)
              lnum(i,n)  = lnum(i,n1)
              eig(i,1,n) = eig(i,1,n1)
              eig(i,jspins,n) = eig(i,jspins,n1)
            ENDDO
            GO TO 70
  40        CONTINUE
         ENDDO
c--->    new atom
         rnot = rmsh(1,n)
         IF (z.LT.1.0) THEN
            DO ispin = 1, jspins
              DO i = 1,jrc(n) ! msh
                 rh(i,n,ispin) = 1.e-10
              ENDDO
            ENDDO
         ELSE
            CALL atom2(
     >                 jrc(n),msh,nstd,z,rnot,h,icorr,total,krla,
     >                 igrd,ndvgrd,idsprs,isprsv,sprsv,
     >                 jspins,jspd,bmu(n),qdel,jri(n),
     <                 rhoss,nst(n),lnum(1,n),eig(1,1,n),vbar(1,n))
            DO ispin = 1, jspins
              DO i = 1, jrc(n) ! msh
                rh(i,n,ispin) = rhoss(i,ispin)
              ENDDO       
            ENDDO       
         END IF
c--->    list atomic density
         iza = zatom(n) + 0.0001
         WRITE (6,FMT=8030) namat(iza)
 8030    FORMAT (/,/,' atom: ',a2,/)
 8040    FORMAT (4 (3x,i5,f8.5,f12.6))
   70    CONTINUE
      ENDDO 

croa+
c..use cdnovlp to generate total density out of atom densities...
c
      DO ispin = 1, jspins
        nat = 1
        DO  n = 1,ntype
           DO  i = 1, jrc(n)
              rh1(i,n,ispin) = rh(i,n,ispin)*fpi*rat(i,n)**2
           ENDDO
           rh1(jrc(n):msh,n,ispin) = 0.0
c..prepare spherical mt charge
           DO i = 1,jri(n)
              rho(i,0,n,ispin) = rh(i,n,ispin)*sfp*rmsh(i,n)**2
           ENDDO
c..reset nonspherical mt charge
           DO k = 1,nlh(ntypsy(nat))
              DO j = 1,jri(n)
                 rho(j,k,n,ispin) = 0.e0
              ENDDO
           ENDDO
           nat = nat + neq(n)
        ENDDO
      ENDDO ! ispin

      qpw(:,:) = czero ; rht(:,:,:) = 0.e0 ; rhtxy(:,:,:,:) = czero

      ENDIF ! irank == 0
      DO ispin = 1, jspins
        CALL cdnovlp(irank,isize,
     >               memd,nlhd,ntypsd,k1d,k2d,k3d,n2d,ntypd,nop,
     >               n3d,natd,jmtd,lmaxd,jspd,nmzxyd,nmzd,msh,
     >               ntype,mx3,nq2,nq3,ntypsy,lmax,nmzxy,
     >               nmz,nvac,neq,bmat,taual,kv3,sk3,ig,kv2,nstr2,
     >               nstr,film,zrfs,z1,delz,omtil,rmt,dx,rmsh,
     >               jri,mlh,clnu,llh,nlh,nmem,mrot,tau,symor,
     >               ispin,rh1(1,1,ispin),invs,dummy,invtab,
     >               odi,ods,amat,ig2,sk2,phi2,vol,
     X               qpw,rhtxy,rho,rht)
croa-
c-spinloop
      ENDDO
      IF ( irank == 0 ) THEN
c
c Check the normalization of total density
c
      CALL qfix(
     >          k1d,k2d,k3d,n3d,ntypd,natd,nop,jspd,jmtd,nmzxyd,
     >          nlhd,nmzd,nmz,jspins,film,nvac,area,nq3,nmzxy,n2d,
     >          ntype,neq,volmts,taual,z1,vol,volint,nq2,invtab,
     >          symor,tau,mrot,rmt,sk3,bmat,ig2,ig,nlh,ntypsd,
     >          nstr,kv3,delz,jri,dx,rmsh,zatom,ntypsy,sigma,
     >          qpw,rhtxy,rho,rht,odi,
     <          fix)
c
c Write superposed density onto density file
c
      iter = 1
      name_l(:) = name(:)
      name_l(10) = 'ordered*'    ! always create ordered start density
      OPEN (71,file='cdn1',form='unformatted',status='new')
      CALL wrtdop(
     >            jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,natd,
     >            jspins,nq3,odi%nq2,nmzxy,nmz,nvac,ntype,neq,
     >            invs,invs2,film,delz,z1,dx,rmt,zatom,
     >            nlh,jri,ntypsd,ntypsy,namat,71,
     >            'input   ','density ',iter,rho,qpw,rht,rhtxy,name_l)
      CLOSE (71)
c
c Check continuity
c
      DO ispin = 1,jspins
        WRITE (6,'(a8,i2)') 'spin No.',ispin
        IF (film .AND. .NOT.odi%d1) THEN
c         ---> vacuum boundaries
           npd = min(nspd,25)
           CALL points(xp,npd)
           DO ivac = 1,nvac
              sign = 3. - 2.*ivac
              DO j = 1,npd
                 xp(3,j) = sign*z1/amat(3,3)
              ENDDO
              CALL checkdop(
     >                    xp,npd,0,0,ivac,1,ispin,.true.,nspd,jmtd,
     >                    memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >                    jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >                    nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >                    tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >                    qpw,rho,rhtxy,rht,odi,ods)
           ENDDO 
        ELSEIF (odi%d1) THEN
c-odim
           npd = min(nspd,25)
           CALL cylpts(xp,npd,z1)
           CALL checkdop(
     >          xp,npd,0,0,nvac,1,ispin,.true.,nspd,jmtd,
     >          memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >          jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >          nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >          tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >          qpw,rho,rhtxy,rht,odi,ods)
c+odim
        END IF
c         ---> m.t. boundaries
        nat = 1
        DO n = 1,ntype
           CALL sphpts(xp,nspd,rmt(n),pos(1,nat))
           CALL checkdop(
     >                   xp,nspd,n,nat,0,-1,ispin,.true.,nspd,jmtd,
     >                   memd,nlhd,ntypsd,n2d,n3d,ntypd,lmaxd,invtab,
     >                   jspd,natd,nmzd,nmzxyd,symor,lmax,ntypsy,nq2,
     >                   nq3,rmt,pos,amat,bmat,kv2,kv3,nop,nop2,ngopr,
     >                   tau,mrot,mlh,nlh,llh,nmem,clnu,jri,nstr,nstr2,
     >                   qpw,rho,rhtxy,rht,odi,ods)
           nat = nat + neq(n)
        ENDDO
      ENDDO

      l_enpara = .false.
      INQUIRE (file='enpara',exist=l_enpara)
c
c set up parameters for enpara-file
c
      IF (.not.l_enpara) THEN
      OPEN (40,file='enpara',form='formatted',status='unknown')
      
      DO n = 1,ntype
        DO i = 0,lmaxd
          lchange(i,n) = .true.
        ENDDO
        DO i = 1, nlo(n)
          llochg(i,n) = .true.
        ENDDO
      ENDDO

      DO nw = 1,nwd
        DO ispin = 1, jspins
c
c vacpar is taken as highest occupied level from atomic eigenvalues
c vacpar (+0.3)  serves as fermi level for film or bulk
c
          vacpar(1) = -999.9
          DO n = 1,ntype
            vacpar(1) = max( vacpar(1),eig(nst(n),ispin,n) )
          ENDDO
          IF (.not.film) vacpar(1) = vacpar(1) + 0.4
          vacpar(2) = vacpar(1)
          

          DO n = 1,ntype
            skiplo(n) = 0
            DO i = 0, 3
              l_found(i) = .false.
              energy(i,n) = vacpar(1)
            ENDDO
            DO i = 1, nlo(n)
              llo_found(i) = .false.
              ello0(i,n) = vacpar(1) - 1.5
            ENDDO
c
c take energy-parameters from atomic calculation
c
            DO i = nst(n), 1,-1 
               IF (.not.film) eig(i,ispin,n) = eig(i,ispin,n) + 0.4

              IF (.not.l_found(lnum(i,n)).AND.(lnum(i,n).LE.3)) THEN
                 energy(lnum(i,n),n) = eig(i,ispin,n)
                 IF (energy(lnum(i,n),n).LT.ellow(nw)) THEN
                   energy(lnum(i,n),n) = vacpar(1)
                   l_found(lnum(i,n))  = .true.
                 ENDIF
                 IF (energy(lnum(i,n),n).LT.elup(nw)) THEN
                   l_found(lnum(i,n))  = .true.
                 ENDIF
              ELSE
                IF (l_found(lnum(i,n)).AND.(nlo(n).GT.0)) THEN
                  DO ilo = 1,nlo(n) 
                    IF (llo(ilo,n).EQ.lnum(i,n)) THEN
                      IF ( .not.llo_found(ilo) ) THEN
                        ello0(ilo,n) = eig(i,ispin,n)
                        IF (( ello0(ilo,n).GT.elup(nw)).OR.
     +                      ( ello0(ilo,n).LT.ellow(nw))) THEN
                          ello0(ilo,n)= vacpar(1)
                        ELSEIF (l_dulo(ilo,n)) THEN
                          ello0(ilo,n)= energy(llo(ilo,n),n)
                        ELSE
                          skiplo(n) = skiplo(n) + 2*llo(ilo,n)+1
                        ENDIF
                        llo_found(ilo) = .true.
                        IF ((energy(llo(ilo,n),n)-ello0(ilo,n).LT.0.5).
     +                    AND.(llo(ilo,n).GE.0)) THEN
                          energy(llo(ilo,n),n) = vacpar(1)
                          IF (l_dulo(ilo,n)) 
     +                              ello0(ilo,n)= energy(llo(ilo,n),n)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF 
              ENDIF

            ENDDO
            IF (lepr.EQ.1) THEN
              DO i = 0, 3
                energy(i,n) = energy(i,n) - vbar(ispin,n)
              ENDDO 
              DO ilo = 1,nlo(n)
                ello0(ilo,n) = ello0(ilo,n) - vbar(ispin,n)
              ENDDO
            ENDIF
          ENDDO  ! atom types

          IF (film) THEN
! get guess for vacuum parameters
! YM : in 1D case should be modified, but also not that bad like this
!
!           generate coulomb potential by integrating inward to z1
!
            DO ivac = 1, nvac
              icorr_dummy = min(max(icorr,0),5)  ! use LDA for this purpose
              DO i=1,nmz
                   sigm(i) = (i-1)*delz*rht(i,ivac,ispin)
              ENDDO
              CALL qsf(delz,sigm,vacpar(ivac),nmz,0)
              denz1(1) = rht(1,ivac,ispin)          ! get estimate for potential at
              CALL  vxcall(6,icorr_dummy,krla,1,    !               vacuum boundary
     >                     1,1,denz1,
     <                     vacpot)
! seems to be the best choice for 1D not to substract vacpar
              IF (.NOT.odi%d1) THEN
                 vacpot(1) = vacpot(1) - fpi*vacpar(ivac)
              END IF
              IF (lepr.EQ.1) THEN
                vacpar(ivac) = -0.2 - vacpot(1)
                WRITE (6,'(" vacuum",i2," reference energy =",f12.6)')
     +                 ivac,vacpot
              ELSE
                vacpar(ivac) = vacpot(1)
              ENDIF
            ENDDO
            IF (nvac.EQ.1) vacpar(2) = vacpar(1)

          ENDIF
          IF (lepr.EQ.1) THEN
            enmix = 0.3
          ELSE
            enmix = 1.0
          ENDIF
c
c write enpara-file
c
          CALL w_enpara(
     >                  lmaxd,nlod,ntype,nw,ispin,film,nlo,skiplo,
     >          ello0,energy,vacpar,lchange,llochg,.true.,enmix,16)

        ENDDO  ! ispin
      ENDDO    ! nw

      CLOSE (40)
      ENDIF
      DEALLOCATE ( qpw,rhtxy,rht,xp,rat,eig,rh,rh1 )
      DEALLOCATE ( rhoss,ello0,energy,vacpar,vbar,sigm )
      ENDIF ! irank == 0
      DEALLOCATE ( rho )
c
      END SUBROUTINE stden
      END MODULE m_stden
