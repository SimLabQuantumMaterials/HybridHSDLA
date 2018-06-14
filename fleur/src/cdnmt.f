      MODULE m_cdnmt
c***********************************************************************
c     This subroutine calculates the spherical and non-spherical charge-
c     density and the orbital moment inside the muffin-tin spheres.
c     Philipp Kurz 2000-02-03
c***********************************************************************
      CONTAINS
      SUBROUTINE cdnmt(
     >     jspd,ntypd,natd,lmaxd,ntypsd,llpd,nlhd,jmtd,nlod,llod,
     >     l_soc,l_mperp,l_fmpl,ntype,jsp_start,jsp_end,sfp,
     >     lmax,nlh,llh,neq,ntypsy,jri,nlo,llo,l_dulo,ulo_der,
     >     dx,rmsh,epar,ello,vr,uu,du,dd,uunmt,udnmt,dunmt,ddnmt,
     >     uulon,dulon,uloulopn,aclo,bclo,cclo,acnmt,bcnmt,ccnmt,
     >     orb,orbl,orblo,mt21,lo21,uloulopn21,uloulop21,
     >     uunmt21,ddnmt21,udnmt21,dunmt21,
     <     chmom,clmom,
     X     qa21,cp_mtsum,rho)

      USE m_types, ONLY : t_orb,t_orbl,t_orblo,t_mt21,t_lo21
      USE m_rhosphnlo
      USE m_radfun
      USE m_orbmom2
      USE m_cputime

      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,ntypd,natd,lmaxd,ntypsd,llpd,nlhd
      INTEGER, INTENT (IN) :: jmtd,nlod,llod
      INTEGER, INTENT (IN) :: ntype,jsp_start,jsp_end
      LOGICAL, INTENT (IN) :: l_soc,l_mperp,l_fmpl
      REAL, INTENT    (IN) :: sfp
      REAL, INTENT (INOUT) :: cp_mtsum
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nlh(ntypsd),llh(0:nlhd,ntypsd),lmax(ntypd)
      INTEGER, INTENT (IN) :: neq(ntypd),ntypsy(natd),jri(ntypd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
      REAL, INTENT    (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL, INTENT    (IN) :: epar(0:lmaxd,ntypd,jspd)
      REAL, INTENT    (IN) :: vr(jmtd,ntypd,jspd)
      REAL, INTENT    (IN) :: ello(nlod,ntypd,jspd)
      REAL, INTENT    (IN) :: uulon(nlod,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: dulon(nlod,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) ::uloulopn(nlod,nlod,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: aclo(nlod,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: bclo(nlod,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: cclo(nlod,nlod,ntypd,jsp_start:jsp_end)
      REAL,INTENT (IN)::acnmt(0:lmaxd,nlod,nlhd,ntypd,jsp_start:jsp_end)
      REAL,INTENT (IN)::bcnmt(0:lmaxd,nlod,nlhd,ntypd,jsp_start:jsp_end)
      REAL,INTENT (IN)::ccnmt(nlod,nlod,nlhd,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: uu(0:lmaxd,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: du(0:lmaxd,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: dd(0:lmaxd,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: uunmt(0:llpd,nlhd,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: udnmt(0:llpd,nlhd,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: dunmt(0:llpd,nlhd,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: ddnmt(0:llpd,nlhd,ntypd,jsp_start:jsp_end)
      REAL, INTENT    (IN) :: uloulopn21(nlod,nlod,ntypd)
      COMPLEX, INTENT (IN) :: uloulop21(nlod,nlod,ntypd)
      COMPLEX, INTENT (IN) :: ddnmt21((lmaxd+1)**2,nlhd,ntypd)
      COMPLEX, INTENT (IN) :: dunmt21((lmaxd+1)**2,nlhd,ntypd)
      COMPLEX, INTENT (IN) :: udnmt21((lmaxd+1)**2,nlhd,ntypd)
      COMPLEX, INTENT (IN) :: uunmt21((lmaxd+1)**2,nlhd,ntypd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      REAL, INTENT   (OUT) :: chmom(ntypd,jspd),clmom(3,ntypd,jspd)
      REAL, INTENT (INOUT) :: rho(jmtd,0:nlhd,ntypd,jspd)
      COMPLEX, INTENT(INOUT) :: qa21(ntypd)
      TYPE (t_orb),  INTENT (IN) :: 
     &              orb(0:lmaxd,-lmaxd:lmaxd,ntypd,jsp_start:jsp_end)
      TYPE (t_orbl), INTENT (IN) :: 
     &              orbl(nlod,-llod:llod,ntypd,jsp_start:jsp_end)
      TYPE (t_orblo),INTENT (IN) :: 
     &              orblo(nlod,nlod,-llod:llod,ntypd,jsp_start:jsp_end)
      TYPE (t_mt21), INTENT (IN) :: mt21(0:lmaxd,ntypd)
      TYPE (t_lo21), INTENT (IN) :: lo21(nlod,ntypd)
C     ..
C     .. Local Scalars ..
      INTEGER itype,na,nd,l,lp,llp,m,lh,j,ispin,noded,nodeu
      INTEGER ilo,ilop
      REAL time1,time2,s,wronk,sumlm,qmtt
      COMPLEX cs
C     ..
C     .. Local Arrays ..
      REAL qmtl(0:lmaxd),qmtllo(0:lmaxd)

C     ..
C     .. Allocatable Arrays ..
      REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:)
      REAL,    ALLOCATABLE :: us(:,:,:),uds(:,:,:)
      REAL,    ALLOCATABLE :: dus(:,:,:),duds(:,:,:)
      REAL,    ALLOCATABLE :: ddn(:,:,:)
      COMPLEX, ALLOCATABLE :: rho21(:,:,:)
c
      IF (l_mperp) THEN
         ALLOCATE ( f(jmtd,2,0:lmaxd,jspd),g(jmtd,2,0:lmaxd,jspd) )
         ALLOCATE ( us(0:lmaxd,ntypd,jspd),uds(0:lmaxd,ntypd,jspd) )
         ALLOCATE ( dus(0:lmaxd,ntypd,jspd),duds(0:lmaxd,ntypd,jspd) )
         ALLOCATE ( ddn(0:lmaxd,ntypd,jspd) )
         IF (l_fmpl) THEN
            ALLOCATE ( rho21(jmtd,0:nlhd,ntypd) )
            rho21(:,:,:) = cmplx(0.0,0.0)
         ENDIF
      ELSE
         ALLOCATE ( f(jmtd,2,0:lmaxd,jsp_start:jsp_end) )
         ALLOCATE ( g(jmtd,2,0:lmaxd,jsp_start:jsp_end) )
         ALLOCATE (   us(0:lmaxd,ntypd,jsp_start:jsp_end) )
         ALLOCATE (  uds(0:lmaxd,ntypd,jsp_start:jsp_end) )
         ALLOCATE (  dus(0:lmaxd,ntypd,jsp_start:jsp_end) )
         ALLOCATE ( duds(0:lmaxd,ntypd,jsp_start:jsp_end) )
         ALLOCATE (  ddn(0:lmaxd,ntypd,jsp_start:jsp_end) )
      ENDIF
      WRITE (6,FMT=8000)
      WRITE (16,FMT=8000)
 8000 FORMAT (/,5x,'l-like charge',/,t6,'atom',t15,'s',t24,'p',
     +     t33,'d',t42,'f',t51,'total')
      CALL cpu_time(time1)

      na = 1
      DO itype = 1,ntype
c--->    spherical component
         DO ispin = jsp_start,jsp_end
            DO l = 0,lmax(itype)
               CALL radfun(
     >              l,epar(l,itype,ispin),vr(1,itype,ispin),jri(itype),
     >              rmsh(1,itype),dx(itype),jmtd,
     <              f(1,1,l,ispin),g(1,1,l,ispin),us(l,itype,ispin),
     <              dus(l,itype,ispin),uds(l,itype,ispin),
     <              duds(l,itype,ispin),ddn(l,itype,ispin),
     <              nodeu,noded,wronk)
               DO j = 1,jri(itype)
                  s = uu(l,itype,ispin)*( f(j,1,l,ispin)*f(j,1,l,ispin)
     +                                +f(j,2,l,ispin)*f(j,2,l,ispin) )
     +             +   dd(l,itype,ispin)*( g(j,1,l,ispin)*g(j,1,l,ispin)
     +                                +g(j,2,l,ispin)*g(j,2,l,ispin) )
     +             + 2*du(l,itype,ispin)*( f(j,1,l,ispin)*g(j,1,l,ispin)
     +                                +f(j,2,l,ispin)*g(j,2,l,ispin) )
                  rho(j,0,itype,ispin) = rho(j,0,itype,ispin)
     +                                + s/(neq(itype)*sfp)
               ENDDO
            ENDDO

c--->       add the contribution of the local orbitals and flapw - lo
c--->       cross-terms to rho, qmtl. the latter are stored in
c--->       qmtllo. initialize qmtllo
            DO l = 0,lmaxd
               qmtllo(l) = 0.0
            END DO

            
            CALL rhosphnlo(
     >           nlod,lmaxd,nlhd,jmtd,ntypsd,natd,ntypsy,nlh,
     >           uloulopn(1,1,itype,ispin),dulon(1,itype,ispin),
     >           uulon(1,itype,ispin),llo(1,itype),nlo(itype),na,
     >           neq(itype),lmax(itype),ello(1,itype,ispin),
     >           vr(1,itype,ispin),jri(itype),rmsh(1,itype),dx(itype),
     >           sfp,aclo(1,itype,ispin),bclo(1,itype,ispin),
     >           cclo(1,1,itype,ispin),acnmt(0,1,1,itype,ispin),
     >           bcnmt(0,1,1,itype,ispin),ccnmt(1,1,1,itype,ispin),
     >           f(1,1,0,ispin),g(1,1,0,ispin),l_dulo(1,itype),
     >           ulo_der(1,itype),
     X           rho(1,0,itype,ispin),qmtllo)
           

c--->       l-decomposed density for each atom type
            qmtt = 0.
            DO l = 0,lmax(itype)
               qmtl(l) = ( uu(l,itype,ispin)+dd(l,itype,ispin)
     *              *ddn(l,itype,ispin) )/neq(itype) + qmtllo(l)
               qmtt = qmtt + qmtl(l)
            END DO
            chmom(itype,ispin) = qmtt
            WRITE (6,FMT=8100) itype, (qmtl(l),l=0,3),qmtt
            WRITE (16,FMT=8100) itype, (qmtl(l),l=0,3),qmtt
 8100       FORMAT (' -->',i2,2x,4f9.5,2x,f9.5)

c+soc
c--->       spherical angular component
            IF (l_soc) THEN
              CALL orbmom2(
     >                     lmaxd,lmax(itype),neq(itype),itype,nlod,llod,
     >                     nlo(itype),llo(1,itype),ddn(0,itype,ispin),
     >                   orb(0,-lmaxd,itype,ispin),uulon(1,itype,ispin),
     >                   dulon(1,itype,ispin),uloulopn(1,1,itype,ispin),
     >           orbl(1,-llod,itype,ispin),orblo(1,1,-llod,itype,ispin),
     <                     clmom(1,itype,ispin))
            ENDIF
c-soc
c--->       non-spherical components
            nd = ntypsy(na)
            DO lh = 1,nlh(nd)
               DO l = 0,lmax(itype)
                  DO lp = 0,l
                     llp = (l* (l+1))/2 + lp
                     DO j = 1,jri(itype)
                        s = uunmt(llp,lh,itype,ispin)*( 
     +                           f(j,1,l,ispin)*f(j,1,lp,ispin)
     +                         + f(j,2,l,ispin)*f(j,2,lp,ispin) )
     +                       + ddnmt(llp,lh,itype,ispin)*(
     +                           g(j,1,l,ispin)*g(j,1,lp,ispin)
     +                         + g(j,2,l,ispin)*g(j,2,lp,ispin) )
     +                       + udnmt(llp,lh,itype,ispin)*(
     +                           f(j,1,l,ispin)*g(j,1,lp,ispin)
     +                         + f(j,2,l,ispin)*g(j,2,lp,ispin) )
     +                       + dunmt(llp,lh,itype,ispin)*(
     +                           g(j,1,l,ispin)*f(j,1,lp,ispin)
     +                         + g(j,2,l,ispin)*f(j,2,lp,ispin) )
                        rho(j,lh,itype,ispin) = rho(j,lh,itype,ispin)
     +                                        + s/neq(itype)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO ! end of spin loop (ispin = jsp_start,jsp_end)

         IF (l_mperp) THEN

c--->      calculate off-diagonal integrated density
           DO l = 0,lmax(itype)
              qa21(itype) = qa21(itype) + conjg(
     +                mt21(l,itype)%uu * mt21(l,itype)%uun +
     +                mt21(l,itype)%ud * mt21(l,itype)%udn +
     +                mt21(l,itype)%du * mt21(l,itype)%dun +
     +                mt21(l,itype)%dd * mt21(l,itype)%ddn )/neq(itype)
           ENDDO
           DO ilo = 1, nlo(itype)
             qa21(itype) = qa21(itype) + conjg(
     +            lo21(ilo,itype)%ulou * lo21(ilo,itype)%uloun +
     +           lo21(ilo,itype)%ulod * lo21(ilo,itype)%ulodn +
     +            lo21(ilo,itype)%uulo * lo21(ilo,itype)%uulon +
     +            lo21(ilo,itype)%dulo * lo21(ilo,itype)%dulon )/
     /                                                 neq(itype)
             DO ilop = 1, nlo(itype)
               qa21(itype) = qa21(itype) + conjg(
     +                      uloulop21(ilo,ilop,itype) *
     +                      uloulopn21(ilo,ilop,itype) )/neq(itype)
             ENDDO
           ENDDO

           IF (l_fmpl) THEN
c--->        the following part can be used to calculate the full magnet.
c--->        density without the atomic sphere approximation for the
c--->        magnet. density, e.g. for plotting.
c--->        calculate off-diagonal part of the density matrix
c--->        spherical component
             DO l = 0,lmax(itype)
               DO j = 1,jri(itype)
                 cs = mt21(l,itype)%uu*( f(j,1,l,2)*f(j,1,l,1) +
     +                                   f(j,2,l,2)*f(j,2,l,1) )
     +              + mt21(l,itype)%ud*( f(j,1,l,2)*g(j,1,l,1) +
     +                                   f(j,2,l,2)*g(j,2,l,1) )
     +              + mt21(l,itype)%du*( g(j,1,l,2)*f(j,1,l,1) +
     +                                   g(j,2,l,2)*f(j,2,l,1) )
     +              + mt21(l,itype)%dd*( g(j,1,l,2)*g(j,1,l,1) +
     +                                   g(j,2,l,2)*g(j,2,l,1) )
                 rho21(j,0,itype) = rho21(j,0,itype)
     +                            + conjg(cs)/(neq(itype)*sfp)
               ENDDO
             ENDDO

c--->        non-spherical components
             nd = ntypsy(na)
             DO lh = 1,nlh(nd)
               DO l = 0,lmax(itype)
                 DO lp = 0,lmax(itype)
                   llp = lp*(lmax(itype)+1)+l+1
                   DO j = 1,jri(itype)
                      cs = uunmt21(llp,lh,itype)*(
     +                         f(j,1,lp,2)*f(j,1,l,1)
     +                       + f(j,2,lp,2)*f(j,2,l,1) )
     +                   + udnmt21(llp,lh,itype)*(
     +                         f(j,1,lp,2)*g(j,1,l,1)
     +                       + f(j,2,lp,2)*g(j,2,l,1) )
     +                   + dunmt21(llp,lh,itype)*(
     +                         g(j,1,lp,2)*f(j,1,l,1)
     +                       + g(j,2,lp,2)*f(j,2,l,1) )
     +                   + ddnmt21(llp,lh,itype)*(
     +                         g(j,1,lp,2)*g(j,1,l,1)
     +                       + g(j,2,lp,2)*g(j,2,l,1) )
                      rho21(j,lh,itype)= rho21(j,lh,itype)
     +                                 + conjg(cs)/neq(itype)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO

           ENDIF ! l_fmpl
         ENDIF ! l_mperp

         na = na + neq(itype)
      ENDDO ! end of loop over atom types
      CALL cpu_time(time2)
      cp_mtsum = cp_mtsum + time2 - time1

c---> for testing: to plot the offdiag. part of the density matrix it
c---> is written to the file rhomt21. This file can read in pldngen.
      IF (l_fmpl) THEN
         OPEN (26,file='rhomt21',form='unformatted',status='unknown')
         WRITE (26) rho21
         CLOSE (26)
         DEALLOCATE ( rho21 )
      ENDIF
c---> end of test output

      DEALLOCATE ( f,g,us,dus,uds,duds,ddn )

      END SUBROUTINE cdnmt
      END MODULE m_cdnmt
