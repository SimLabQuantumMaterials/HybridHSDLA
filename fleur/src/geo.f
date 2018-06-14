      MODULE m_geo
      CONTAINS
      SUBROUTINE geo(
     >               ntypd,natd,nwdd,nop,layerd,nlod,
     >               tau,mrot,ngopr,invtab,odi,ods,odd,tote,
     X               forcetot)

*********************************************************************
* calculates the new atomic positions after the force calculation   *
* subroutine is based on a BFGS method implemented by M. Weinert    *
*                                [cf. PRB 52 (9) p. 6313 (1995)]    * 
*                                                                   *
* as a first step we read in the file 'inp' with some additional    *
* information (subroutine rw_inp)                                   *
* then recover the old geometry optimisation information from file  *
* 'forces.dat' (subroutine bfsg0)                                   *
* this input together with the new forces (forcetot) are now used   *
* to calculate the new atomic positions (subroutine bfsg)           *
* finally the new 'inp' file is written (subroutine rw_inp)         *
*                                                           Gustav  *
c
c input: 
c        ntype .... total number of atom types
c        thetad ... approx. debye temperature
c        zat(ntype) mass number of the atom (or atomic number)
c        xa ....... mixing factor 
c        epsdisp .. limit for displacement to be converged
c        epsforce . the same for force
c        istepnow . steps to be done in this run
c
*********************************************************************

      USE m_cotra,     ONLY : cotra0,cotra1
      USE m_constants, ONLY : pimach
      USE m_rwinp
      USE m_types,     ONLY : t_utype
      USE m_inv3
      USE m_bfgs
      USE m_bfgs0
      USE m_od_types, ONLY : od_inp, od_sym, od_dim

      IMPLICIT NONE
C ..
C ..  Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,natd,nwdd,nop,layerd,nlod
      REAL,    INTENT (IN) :: tote
C ..
C ..  Array Arguments ..
      INTEGER, INTENT (IN) :: ngopr(natd),invtab(nop)
      INTEGER, INTENT (IN) :: mrot(3,3,nop)
      REAL,    INTENT (IN) :: tau(3,nop)
      REAL,    INTENT (INOUT) :: forcetot(3,ntypd)
c-odim
      TYPE (od_inp) odi
      TYPE (od_sym), INTENT (IN) :: ods
      TYPE (od_dim) odd
c+odim
C ..
C ..  Local Scalars ..
      REAL thetad,xa,epsdisp,epsforce,tpi,tworkf
      REAL dvac,dtild,scale,scpos,gmax,delgau,zc,tkb,alpha,spinf
      REAL e1s,e2s,chng,gmaxxc,omtil,theta,phi,e1_dos,e2_dos,sig_dos
      INTEGER i,j,na,ntype,istep0,istep,itype,jop
      INTEGER isec1,ndir,jspins,lpr,nwd,lepr,kcrel,itmax
      INTEGER imix,layers,kk,nnne,ieq,maxiter,nstars,nstm,gw,gw_neigd
      INTEGER igrd,ndvgrd,idsprs
      CHARACTER*3 latnam
      CHARACTER*4 namex,namgrp
      CHARACTER*12 relcor
      LOGICAL lconv,l_u2f,l_f2u,l_bmt,l_soc,starcoeff,l_noco,l_J
      LOGICAL strho,film,dos,secvar,invs,zrfs,invs2,vchk,cdinf,pot8
      LOGICAL form66,l_f,eonly,gauss,tria,frcor,slice,ctail,disp
      LOGICAL swsp,lflip,vacdos,integ,iplot,score,plpot,pallst,lwb
C ..
C ..  Local Arrays ..
      REAL a1(3),a2(3),a3(3),rmt(ntypd),dx(ntypd),amat(3,3),bmat(3,3)
      REAL ellow(nwdd),elup(nwdd),zelec(nwdd),bmu(ntypd),rkm(nwdd)
      REAL xold(3*ntypd),y(3*ntypd),h(3*ntypd,3*ntypd),zat(ntypd)
      REAL tau0(3,ntypd),tau0_i(3,ntypd),taual(3,natd),locx(2),locy(2)
      INTEGER neq(ntypd),nz(ntypd),relax(3,ntypd),llo(nlod,ntypd)
      INTEGER ncst(ntypd),jri(ntypd),lmax(ntypd)
      INTEGER lnonsph(ntypd),nflip(ntypd),izlay(layerd,2),nlo(ntypd)
      CHARACTER*3 noel(ntypd)
      CHARACTER*8 name(10)
      LOGICAL l_geo(ntypd),eig66(2),soc_opt(ntypd+2)
c+lda+u
      TYPE (t_utype) lda_u(ntypd)
c-lda+u

c---> define constant 2*Pi
      tpi = 2 * pimach()

      DO i = 1,3
         a1(i) = 0.0
         a2(i) = 0.0
         a3(i) = 0.0
      ENDDO

      gw=0
      CALL rw_inp(
     >            16,ntypd,natd,nwdd,layerd,nlod,5,'R',
     <  dvac,dtild,scale,scpos,gmax,rkm,delgau,zc,tkb,alpha,spinf,
     <  e1s,e2s,isec1,ndir,jspins,lpr,nwd,lepr,kcrel,itmax,imix,layers,
     <  l_u2f,l_f2u,l_bmt,kk,nnne,maxiter,latnam,noel,namex,namgrp,
     <  relcor,
     <  strho,film,dos,secvar,invs,zrfs,invs2,vchk,cdinf,pot8,form66,
     <  l_f,eonly,eig66,gauss,tria,frcor,slice,ctail,disp,swsp,lflip,
     <  vacdos,integ,iplot,score,plpot,pallst,a1,a2,a3,rmt,dx,ellow,
     <  elup,zelec,bmu,ncst,jri,lmax,lnonsph,nflip,izlay,
     <  name,igrd,ndvgrd,idsprs,lwb,chng,gmaxxc,l_soc,soc_opt,theta,phi,
     <  taual,
     <  ntype,neq,nz,xa,relax,thetad,epsdisp,epsforce,nlo,llo,tworkf,
     <  nstars,nstm,starcoeff,locx,locy,l_noco,l_J,l_geo,e1_dos,e2_dos,
     <  sig_dos,lda_u,gw,gw_neigd,odd)
c
      DO i = 1,3
         amat(i,1) = a1(i)*scale
         amat(i,2) = a2(i)*scale
         amat(i,3) = a3(i)*scale
      ENDDO
      CALL inv3(amat,bmat,omtil)
      DO j = 1,3
         DO i = 1,3
            bmat(i,j) = tpi*bmat(i,j)
         ENDDO
      ENDDO
c
      na = 1
      DO i = 1,ntype
        zat(i)=real(nz(i))
        IF (film) taual(3,na) = taual(3,na)/a3(3)
        DO j = 1,3
          tau0_i(j,i) = taual(j,na)
        ENDDO
        CALL cotra0(tau0_i(1,i),tau0(1,i),amat)
        na = na + neq(i)
      ENDDO

      CALL bfgs0(
     >           ntype,
     <           istep0,xold,y,h)
             
      DO itype=1,ntype
        IF (l_geo(itype)) THEN
          WRITE (6,'(6f10.5)') (tau0(j,itype),j=1,3),
     +                         (forcetot(i,itype),i=1,3)
          DO i = 1,3
            forcetot(i,itype)=forcetot(i,itype)*real(relax(i,itype))
          ENDDO
          WRITE (6,'(6f10.5,a,3i2)') (tau0(j,itype),j=1,3),
     +             (forcetot(i,itype),i=1,3),' relax: ',
     +             (relax(i,itype),i=1,3)
        ELSE
          DO i = 1,3
            forcetot(i,itype)=0.0
          ENDDO
        ENDIF
      ENDDO

      istep = 1
      CALL bfgs(
     >          ntype,istep,istep0,forcetot,
     >          zat,xa,thetad,epsdisp,epsforce,tote,
     X          xold,y,h,tau0,
     <          lconv)

      IF (lconv) THEN

        WRITE (6,'(a)') "Des woars!"
        STOP  ' GEO Des woars '

      ELSE

        na = 0
        DO itype=1,ntype
          CALL cotra1(tau0(1,itype),tau0_i(1,itype),bmat)
          DO ieq = 1,neq(itype)
            na = na + 1
            jop = invtab(ngopr(na))
            IF (odi%d1) jop = ods%ngopr(na)
            DO i = 1,3
              taual(i,na) = 0.0
              DO j = 1,3
                 IF (.NOT.odi%d1) THEN
                    taual(i,na) = taual(i,na) +
     +                   mrot(i,j,jop) * tau0_i(j,itype)
                 ELSE
                    taual(i,na) = taual(i,na) +
     +                   ods%mrot(i,j,jop) * tau0_i(j,itype)
                 END IF
              ENDDO
              IF (odi%d1) THEN
                 taual(i,na) = taual(i,na) +
     +              ods%tau(i,jop)/amat(3,3)
              ELSE
                 taual(i,na) = taual(i,na) + tau(i,jop)
              END IF
            ENDDO
          ENDDO
        ENDDO

        l_f = .false.
        CALL rw_inp(
     >              16,ntypd,natd,nwdd,layerd,nlod,45,'W',
     >  dvac,dtild,scale,scpos,gmax,rkm,delgau,zc,tkb,alpha,spinf,
     >  e1s,e2s,isec1,ndir,jspins,lpr,nwd,lepr,kcrel,itmax,imix,layers,
     >  l_u2f,l_f2u,l_bmt,kk,nnne,maxiter,latnam,noel,namex,namgrp,
     >  relcor,
     >  strho,film,dos,secvar,invs,zrfs,invs2,vchk,cdinf,pot8,form66,
     >  l_f,eonly,eig66,gauss,tria,frcor,slice,ctail,disp,swsp,lflip,
     >  vacdos,integ,iplot,score,plpot,pallst,a1,a2,a3,rmt,dx,ellow,
     >  elup,zelec,bmu,ncst,jri,lmax,lnonsph,nflip,izlay,
     >  name,igrd,ndvgrd,idsprs,lwb,chng,gmaxxc,l_soc,soc_opt,theta,phi,
     >  taual,
     >  ntype,neq,nz,xa,relax,thetad,epsdisp,epsforce,nlo,llo,tworkf,
     >  nstars,nstm,starcoeff,locx,locy,l_noco,l_J,l_geo,e1_dos,e2_dos,
     >  sig_dos,lda_u,gw,gw_neigd,odd)

      ENDIF

      RETURN
      END SUBROUTINE geo
      END MODULE m_geo
