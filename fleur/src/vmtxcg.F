      MODULE m_vmtxcg
c.....------------------------------------------------------------------
c     fit spherical-harmonics expansion of exchange-correlation
c     potential inside muffint-tin spheres and add it to coulomb
c     potential
c                                     c.l.fu and r.podloucky           *
c     for the gradient correction. t.a. 1996.
c.....------------------------------------------------------------------
c     instead of vmtxcor.f: the different exchange-correlation
c     potentials defined through the key icorr are called through
c     the driver subroutine vxcallg.f, subroutines vectorized
c     ** r.pentcheva 22.01.96
c     *********************************************************
c     angular mesh calculated on speacial gauss-legendre points
c     in order to use orthogonality of lattice harmonics and
c     avoid a least square fit
c     ** r.pentcheva 04.03.96
c     *********************************************************
c
      CONTAINS
      SUBROUTINE vmtxcg(
     >                jspd,memd,nlhd,ntypsd,jmtd,ntypd,lmaxd,nspd,
     >                clnu,mlh,nmem,llh,nlh,rmsh,ntypsy,jri,dx,natd,neq,
     >                rho,icorr,total,krla,ntype,jspins,nsymt,
     >                igrd,ndvgrd,isprsv,idsprs0,idsprsl,chng,sprsv,
     X                vr,idsprs,rhmn,ichsmrg,
     <                excr)

      USE m_lhglptg
      USE m_grdchlh
      USE m_mkgylm
      USE m_gaussp
      USE m_xcallg, ONLY : vxcallg,excallg

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: jspd,memd,nlhd,ntypsd,jmtd,ntypd,lmaxd,nspd
      INTEGER, INTENT(IN) :: icorr,krla,ntype,jspins,nsymt,igrd,ndvgrd
      INTEGER, INTENT(IN) :: isprsv,idsprs0,idsprsl,natd
      LOGICAL, INTENT(IN) :: total
      REAL,    INTENT(IN) :: sprsv,chng
      INTEGER, INTENT(INOUT):: idsprs,ichsmrg
      REAL,    INTENT(INOUT):: rhmn
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: jri(ntypd),ntypsy(natd),neq(ntypd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),dx(ntypd)
      REAL,    INTENT (IN) :: rho(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (INOUT):: vr(jmtd,0:nlhd,ntypd,jspd)
      REAL,    INTENT (OUT)  :: excr(jmtd,0:nlhd,ntypd)
c     ..
c     .. Local Scalars ..
      INTEGER jr,js,k,lh,n,nd,ist,nsp,ixpm,idsprssv,nat,i
      REAL    rhmnm,d_15,elh,vlh
      LOGICAL lwbc              ! if true, white-bird trick
c     ..
c     .. Local Arrays ..
      REAL vxc(nspd,jspd),exc(nspd),rx(3,nspd)
      REAL vxcl(nspd,jspd),excl(nspd),divi
      REAL wt(nspd),rr2(jmtd),thet(nspd)
      REAL agr(nspd),agru(nspd),agrd(nspd),g2r(nspd),g2ru(nspd)
      REAL g2rd(nspd),gggr(nspd),gggru(nspd),gggrd(nspd)
      REAL grgru(nspd),grgrd(nspd),gzgr(nspd)
      REAL, ALLOCATABLE :: ylh(:,:,:),ylht(:,:,:),ylhtt(:,:,:)
      REAL, ALLOCATABLE :: ylhf(:,:,:),ylhff(:,:,:),ylhtf(:,:,:)
      REAL, ALLOCATABLE :: chlh(:,:,:),chlhdr(:,:,:),chlhdrr(:,:,:)
      REAL, ALLOCATABLE :: ch(:,:),chdr(:,:),chdt(:,:),chdf(:,:)
      REAL, ALLOCATABLE :: chdrr(:,:),chdtt(:,:),chdff(:,:),chdtf(:,:)
      REAL, ALLOCATABLE :: chdrt(:,:),chdrf(:,:) 
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC max,mod,min
cta+
c.....------------------------------------------------------------------
      idsprssv = idsprs
      idsprs   = min(idsprs0,idsprsl)

      WRITE (6,'(/'' 999idsprs in starting of nmtxad='',i2)') idsprs

      lwbc=.false.

      d_15 = 1.e-15
      ixpm = 1
      ist  = 1

c     generates nspd points on a sherical shell with radius 1.0
c     angular mesh equidistant in phi,
c     theta are zeros of the legendre polynomials
c
      CALL gaussp(
     >            lmaxd,
     <            rx,wt)

      nsp = nspd
c
c     generates the lattice harmonics on the angular mesh
c
      ALLOCATE ( ylh(nspd,0:nlhd,ntypsd),ylht(nspd,0:nlhd,ntypsd),
     +           ylhtt(nspd,0:nlhd,ntypsd),ylhf(nspd,0:nlhd,ntypsd),
     +           ylhff(nspd,0:nlhd,ntypsd),ylhtf(nspd,0:nlhd,ntypsd) )

      CALL lhglptg(
     >             nspd,nlhd,ntypsd,lmaxd,memd,
     >             rx,nsp,igrd,nsymt,
     >             clnu,nmem,mlh,nlh,llh,
     <             ylh,thet,ylht,ylhtt,ylhf,ylhff,ylhtf)
c
c
c     loop over topologically non-equivalent atoms
c

      ALLOCATE ( ch(nspd,jspd),chdr(nspd,jspd),chdt(nspd,jspd),
     +           chdf(nspd,jspd),chdrr(nspd,jspd),chdtt(nspd,jspd),
     +           chdff(nspd,jspd),chdtf(nspd,jspd),chdrt(nspd,jspd),
     +           chdrf(nspd,jspd),chlh(jmtd,0:nlhd,jspd),
     +           chlhdr(jmtd,0:nlhd,jspd),chlhdrr(jmtd,0:nlhd,jspd))

      nat = 1
      DO 200 n = 1,ntype

        nd = ntypsy(nat)

        DO jr = 1,jri(n)
          rr2(jr) = 1.e0/ (rmsh(jr,n)*rmsh(jr,n))
        ENDDO

        DO lh = 0,nlh(nd)

c         calculates gradients of radial charge densities of l=> 0.
c         rho*ylh/r**2 is charge density. chlh=rho/r**2.
c         charge density=sum(chlh*ylh).
c         chlhdr=d(chlh)/dr, chlhdrr=dd(chlh)/drr.

          DO js = 1,jspins

            DO jr = 1,jri(n)
              chlh(jr,lh,js) = rho(jr,lh,n,js)*rr2(jr)
            ENDDO

            IF (igrd.GT.0) THEN 
              CALL grdchlh(
     >                     ixpm,ist,jri(n),dx(n),rmsh(1,n),
     >                     chlh(1,lh,js),ndvgrd,
     <                     chlhdr(1,lh,js),chlhdrr(1,lh,js))
            ENDIF

          ENDDO ! js
        ENDDO   ! lh
!
!-->    loop over radial mesh 
!
        DO 190 jr = 1,jri(n)
!
!         following are at points on jr-th sphere.
!         initialize arrays
!
          ch(:,:)    = 0.0     ! charge density
          chdr(:,:)  = 0.0     ! d(ch)/dr
          chdt(:,:)  = 0.0     ! d(ch)/dtheta
          chdf(:,:)  = 0.0     ! d(ch)/dfai
          chdrr(:,:) = 0.0     ! dd(ch)/drr
          chdtt(:,:) = 0.0     ! dd(ch)/dtt
          chdff(:,:) = 0.0     ! dd(ch)/dff
          chdtf(:,:) = 0.0     ! dd(ch)/dtf
          chdrt(:,:) = 0.0     ! d(d(ch)/dr)dt
          chdrf(:,:) = 0.0     ! d(d(ch)/dr)df

          DO js = 1,jspins
            DO lh = 0,nlh(nd)

c             generate the densities on an angular mesh

              DO k = 1,nsp
                ch(k,js) = ch(k,js) + ylh(k,lh,nd)*chlh(jr,lh,js)
              ENDDO

              IF (igrd.GT.0) THEN
! 
                IF ((idsprs0.EQ.0 .AND. lh.EQ.0) .OR.
     &              (idsprsl.EQ.0 .AND. lh.GT.0)) THEN
                  DO k = 1,nsp
                  chdr(k,js) =chdr(k,js)+ ylh(k,lh,nd)*chlhdr(jr,lh,js)
                  chdrr(k,js)=chdrr(k,js)+ylh(k,lh,nd)*chlhdrr(jr,lh,js)
                  ENDDO
                ENDIF

                IF (idsprsl.EQ.0) THEN
                  DO k = 1,nsp
                  chdrt(k,js)=chdrt(k,js)+ylht(k,lh,nd)*chlhdr(jr,lh,js)
                  chdrf(k,js)=chdrf(k,js)+ylhf(k,lh,nd)*chlhdr(jr,lh,js)
                  chdt(k,js) =chdt(k,js) +ylht(k,lh,nd)*chlh(jr,lh,js)
                  chdf(k,js) =chdf(k,js) +ylhf(k,lh,nd)*chlh(jr,lh,js)
                  chdtt(k,js)=chdtt(k,js)+ylhtt(k,lh,nd)*chlh(jr,lh,js)
                  chdff(k,js)=chdff(k,js)+ylhff(k,lh,nd)*chlh(jr,lh,js)
                  chdtf(k,js)=chdtf(k,js)+ylhtf(k,lh,nd)*chlh(jr,lh,js)
                  ENDDO
                ENDIF

              ENDIF

            ENDDO ! lh
          ENDDO   ! js

          IF (igrd.GT.0) THEN
            CALL mkgylm(
     >                  jspins,rmsh(jr,n),thet,nsp,nspd,jspd,ch,chdr,
     >                  chdt,chdf,chdrr,chdtt,chdff,chdtf,chdrt,chdrf,
     <                  agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,
     <                  grgru,grgrd,gzgr)
          ELSE
            agr(:)   = 0.0 ; agru(:)  = 0.0 ; agrd(:)  = 0.0 
            g2r(:)   = 0.0 ; g2ru(:)  = 0.0 ; g2rd(:)  = 0.0 
            gggr(:)  = 0.0 ; gggru(:) = 0.0 ; gggrd(:) = 0.0 
            grgru(:) = 0.0 ; grgrd(:) = 0.0 ; gzgr(:)  = 0.0 
          ENDIF

c
c         rhmnm: rho_minimum_muffin-tin..

          rhmnm=10.e+10

          DO js=1,jspins
            DO i=1,nsp
              ch(i,js) = max(ch(i,js),d_15)
              rhmnm = min(rhmnm,ch(i,js))
            ENDDO
          ENDDO

          IF (rhmnm.LT.rhmn) THEN
            rhmn = rhmnm
            ichsmrg = 1
          ENDIF

          IF (rhmn.LT.chng) THEN
            WRITE (6,'(/'' rhmn.lt.chng in vmtxc. rhmn,chng='',
     +        2d9.2)') rhmn,chng
c            STOP 'vmtxcg: rhmn.lt.chng'
          ENDIF

c
c         calculate the ex.-cor. potential

          IF (mod(jr,1000).eq.0)
     &      WRITE (6,'(/'' 9999ic,kr,ig,js,nsp,iwb='',5i5,l3/
     &            '' ch''/(10d15.7))') icorr,krla,igrd,jspins,
     &            nsp,lwbc,((ch(k,js),k=1,nsp),js=1,jspins)

#ifdef CPP_600
#endif
          CALL vxcallg(
     >                 icorr,lwbc,jspins,nspd,nsp,ch,agr,agru,agrd,
     >                 g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <                 vxc,
     >                 idsprs,isprsv,sprsv)
#ifdef CPP_600
          IF (jri(n) < 600) THEN
          ELSE
          CALL vxcallg(
     >                 -1,lwbc,jspins,nspd,nsp,ch,agr,agru,agrd,
     >                 g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <                 vxcl,
     >                 idsprs,isprsv,sprsv)
          divi = 1.0 / (rmsh(jri(n),n) - rmsh(1,n))
          vxc(:,:) = ( vxcl(:,:) * ( rmsh(jri(n),n) - rmsh(jr,n) ) +
     +                  vxc(:,:) * ( rmsh(jr,n) - rmsh(1,n) ) ) * divi
          ENDIF
#endif

          IF (mod(jr,1000).eq.0)
     &              WRITE (6,'(/'' 999vxc''/(10d15.7))')
     &                    ((vxc(k,js),k=1,nsp),js=1,jspins)


c         now determine the corresponding potential number

          DO js = 1,jspins
c
c           multiplicate vxc with the weights of the k-points
c
            DO k = 1,nsp
              vxc(k,js) = vxc(k,js)*wt(k)
            ENDDO

            IF (mod(jr,1500).EQ.0)
     &        WRITE (6,'('' 999wt''/(10d15.7))') (wt(k),k=1,nsp)
            IF (mod(jr,1500).EQ.0)
     &        WRITE (6,'('' 999vxc''/(10d15.7))') (vxc(k,js),k=1,nsp)

            DO lh = 0,nlh(nd)
              vlh = 0
c
c --->        determine the corresponding potential number
cc            through gauss integration
c
              DO k = 1,nsp
                vlh = vlh + vxc(k,js)*ylh(k,lh,nd)
              ENDDO

c --->        add to the given potential

              vr(jr,lh,n,js) = vr(jr,lh,n,js) + vlh

              IF (mod(jr,1500).EQ.0) THEN
                WRITE(6,'('' 999lh,js,jr,n,vlh='',4i5,d15.7)')
     &            lh,js,jr,n,vlh
                WRITE(6,'('' 9ylh''/(10d15.7))') (ylh(k,lh,nd),k=1,nsp)
              ENDIF

            ENDDO ! lh
          ENDDO   ! js


          IF (total) then
c
c           calculate the ex.-cor energy density
c
            IF (mod(jr,2500).EQ.0)
     +        WRITE (6,'(/'' 9999ic,kr,ig,js,nsp='',5i5/
     +               '' ch''/(10d15.7))') icorr,krla,igrd,jspins,
     +               nsp,((ch(k,js),k=1,nsp),js=1,jspins)

#ifdef CPP_600
#endif
            CALL excallg(icorr,lwbc,jspins,nspd,nsp,
     >                   ch,agr,agru,agrd,g2r,g2ru,g2rd,
     >                   gggr,gggru,gggrd,gzgr,
     <                   exc,
     >                   idsprs,isprsv,sprsv)
#ifdef CPP_600
            IF (jri(n) < 600) THEN
            ELSE
            CALL excallg(-1,lwbc,jspins,nspd,nsp,
     >                   ch,agr,agru,agrd,g2r,g2ru,g2rd,
     >                   gggr,gggru,gggrd,gzgr,
     <                   excl,
     >                   idsprs,isprsv,sprsv)

            divi = 1.0 / ( rmsh(jri(n),n) - rmsh(1,n))
            exc(:) = ( excl(:) * ( rmsh(jri(n),n) - rmsh(jr,n) ) +
     +                exc(:) * ( rmsh(jr,n) - rmsh(1,n) ) ) * divi
            ENDIF
#endif


            IF (mod(jr,10000).EQ.0)
     &        WRITE (6,'(/'' 999exc''/(10d15.7))') (exc(k),k=1,nsp)

          ENDIF

c         now determine the corresponding energy density number
c
c         multiplicate exc with the weights of the k-points
c
          DO k = 1,nsp
            exc(k) = exc(k)*wt(k)
          ENDDO 

          DO lh = 0,nlh(nd)
c
c           determine the corresponding potential number through gauss
c           integration
c
            elh = 0.0
            DO k = 1,nsp
              elh = elh + exc(k)*ylh(k,lh,nd)
            ENDDO
            excr(jr,lh,n) =  elh

          ENDDO   

  190   ENDDO   

c        WRITE(6,'(/'' n='',i3/'' 9999vr''/(10d15.7))') n,
c     &   (((vr(jr,lh,n,js),jr=1,jri(n),100),lh=0,ntypsy(nat)),js=1,jspins)
c        WRITE(6,'(/'' 9999excr''/(10d15.7))')
c     &   ((excr(jr,lh,n),jr=1,jri(n),100),lh=0,ntypsy(nat))

        nat = nat + neq(n)
  200 ENDDO
      DEALLOCATE (ch,chdr,chdt,chdf,chdrr,chdtt,chdff,chdtf,chdrt,chdrf)
      DEALLOCATE (chlh,chlhdr,chlhdrr,ylh,ylht,ylhtt,ylhf,ylhff,ylhtf)
c
      idsprs = idsprssv
      RETURN
      END SUBROUTINE vmtxcg
      END MODULE m_vmtxcg
