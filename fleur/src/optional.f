      MODULE m_optional
      CONTAINS
      SUBROUTINE optional(
     > irank,isize,
     > ntypd,ntypsd,nmzxyd,nmzd,nlhd,jmtd,msh,memd,
     > jspd,n3d,n2d,ntype,jri,jspins,lmaxd,nop,natd,
     > k1d,k2d,k3d,nspd,nstd,mx3,omtil,zrfs,ig,taual,
     > iplot,strho,swsp,lflip,l_f2u,l_u2f,l_bmt,symor,sk3,
     > plpot,idsprs,ndvgrd,igrd,krla,icorr,lepr,
     > name,namat,nq3,nq2,nmzxy,nmz,delz,z1,volint,
     > invs2,invs,film,nflip,bmat,amat,pos,invtab,
     > neq,nstr2,nstr,lmax,mlh,llh,nmem,clnu,ngopr,
     > ntypsy,nlh,slice,bmu,sprsv,isprsv,mrot,tau,
     > nvac,zatom,dx,rmt,rmsh,kv2,kv3,nop2,
     > ig2,vol,volmts,area,score,l_noco,nwd,nlod,nlo,
     > llo,l_dulo,nwdd,ellow,elup,sigma,lda_u,n_u,
     > sk2,phi2,odi,ods)
c
c----------------------------------------
c this routine is called by: fleur.F
c
c optional --+-- plot -+- loddop
c            |         +- outcdn -+- cotra0
c            |                    +- cotra1
c            |                    +- starf2 -- spgrot
c            |                    +- starf3
c            |                    +- ylm3
c            +-- stden -+- atom2 -+- setcor
c            |          |         +- stpot1
c            |          |         +- differ -+- inwint
c            |          |         |          +- outint
c            |          |         +- vxcall (-> see vgen.F) or:
c            |          |         +- potl0 -+- grdchlh
c            |          |         |         +- mkgl0
c            |          |         |         +- vxcallg (-> see vgen.F)
c            |          |         +- intgr1
c            |          +- cdnovlp -+- spgrot
c            |          |           +- rcerf --wofz
c            |          |           +- diflgr
c            |          |           +- qpw_to_nmt -+- phasy1 -+- spgrot
c            |          |                          |          +- ylm3
c            |          |                          +- sphbes
c            |          +- qfix -- cdntot -+- intgr3
c            |          |                  +- qsf
c            |          |                  +- pwint -- spgrot
c            |          +- wrtdop
c            |          +- points -- qranf
c            |          +- sphpts -- qranf
c            |          +- checkdop -+- starf3
c            |                       +- cotra0
c            |                       +- starf2 -- spgrot
c            |                       +- fitchk
c            |                       +- cotra1
c            |                       +- ylm3
c            +-- cdnsp -+- loddop
c            |          +- wrtdop
c            |          +- intgr3
c            +-- flipcdn -+- loddop
c            |            +- wrtdop
c            +-- f2u -- wrtdop
c            +-- u2f -- loddop
c            +-- bmt -+- loddop
c                     +- wrtdop
c----------------------------------------

      USE m_plotdop
      USE m_stden
      USE m_cdnsp
      USE m_flipcdn
      USE m_f2u
      USE m_u2f
      USE m_cputime
      USE m_outtime
      USE m_od_types, ONLY : od_inp, od_sym

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: irank,isize
      INTEGER, INTENT (IN) :: ntypd,ntypsd,nmzxyd,nmzd,nlhd,jmtd,memd
      INTEGER, INTENT (IN) :: jspd,n3d,n2d,ntype,jspins,msh,nop,natd
      INTEGER, INTENT (IN) :: nq3,nq2,nmzxy,nmz,nvac,lmaxd,nop2,lepr
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,nspd,nstd,mx3,nwd,nlod,nwdd
      INTEGER, INTENT (IN) :: igrd,krla,icorr,idsprs,ndvgrd,isprsv,n_u
      REAL,    INTENT (IN) :: sprsv,delz,z1,vol,volint,area,omtil,sigma
      LOGICAL, INTENT (IN) :: iplot,strho,swsp,lflip,l_f2u,l_u2f,l_bmt
      LOGICAL, INTENT (IN) :: l_noco
      LOGICAL, INTENT (IN) :: plpot,invs,invs2,film,slice,symor,zrfs
      LOGICAL, INTENT (INOUT) :: score
C     ..
C     .. Array Arguments ..      
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (IN) :: nflip(ntypd),ntypsy(natd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: jri(ntypd),mrot(3,3,nop),ngopr(natd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),lmax(ntypd)
      INTEGER, INTENT (IN) :: nstr(n3d),nstr2(n2d),neq(ntypd)
      INTEGER, INTENT (IN) :: kv2(2,n2d),kv3(3,n3d),invtab(nop)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntypd),tau(3,nop),pos(3,natd)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3),bmu(ntypd),sk3(n3d)
      REAL,    INTENT (IN) :: taual(3,natd)
      REAL,    INTENT (IN) :: rmt(ntype),dx(ntype),zatom(ntype)
      REAL,    INTENT (IN) :: volmts(ntype),ellow(nwdd),elup(nwdd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
      CHARACTER*8 ,INTENT (IN) :: name(10)
      CHARACTER*2 ,INTENT (IN) :: namat(0:103)
      INTEGER     ,INTENT (IN) :: lda_u(ntype)
c-odim
      REAL        ,INTENT (IN) :: sk2(n2d),phi2(n2d)
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      REAL    :: time1,time2
      INTEGER :: it
      LOGICAL :: total
      CHARACTER*10 :: cdnfname
C     ..
      it = 1
      CALL cpu_time(time1)
      IF (irank == 0) THEN
   10 IF (plpot) score = .false.
      IF (iplot) THEN
         IF (strho) STOP 'strho=T and iplot=T'
         IF (l_noco) THEN
            cdnfname = 'cdn'
            CALL plotdop(
     >           odi,ods,n3d,1,nmzxyd,n2d,memd,nlhd,ntypsd,ntype,
     >           lmaxd,jmtd,ntypd,natd,nmzd,1,neq,
     >           nq3,nvac,nmz,nmzxy,nq2,nop,nop2,volint,film,
     >           plpot,score,slice,symor,invs,invs2,z1,delz,
     >           clnu,llh,nmem,mlh,nlh,ngopr,ntypsy,jri,pos,zatom,
     >           nstr,nstr2,lmax,kv2,kv3,mrot,tau,rmsh,invtab,
     >           amat,bmat,l_noco,cdnfname)
            cdnfname = 'mdnx'
            CALL plotdop(
     >           odi,ods,n3d,1,nmzxyd,n2d,memd,nlhd,ntypsd,ntype,
     >           lmaxd,jmtd,ntypd,natd,nmzd,1,neq,
     >           nq3,nvac,nmz,nmzxy,nq2,nop,nop2,volint,film,
     >           plpot,score,slice,symor,invs,invs2,z1,delz,
     >           clnu,llh,nmem,mlh,nlh,ngopr,ntypsy,jri,pos,zatom,
     >           nstr,nstr2,lmax,kv2,kv3,mrot,tau,rmsh,invtab,
     >           amat,bmat,l_noco,cdnfname)
            cdnfname = 'mdny'
            CALL plotdop(
     >           odi,ods,n3d,1,nmzxyd,n2d,memd,nlhd,ntypsd,ntype,
     >           lmaxd,jmtd,ntypd,natd,nmzd,1,neq,
     >           nq3,nvac,nmz,nmzxy,nq2,nop,nop2,volint,film,
     >           plpot,score,slice,symor,invs,invs2,z1,delz,
     >           clnu,llh,nmem,mlh,nlh,ngopr,ntypsy,jri,pos,zatom,
     >           nstr,nstr2,lmax,kv2,kv3,mrot,tau,rmsh,invtab,
     >           amat,bmat,l_noco,cdnfname)
            cdnfname = 'mdnz'
            CALL plotdop(
     >           odi,ods,n3d,1,nmzxyd,n2d,memd,nlhd,ntypsd,ntype,
     >           lmaxd,jmtd,ntypd,natd,nmzd,1,neq,
     >           nq3,nvac,nmz,nmzxy,nq2,nop,nop2,volint,film,
     >           plpot,score,slice,symor,invs,invs2,z1,delz,
     >           clnu,llh,nmem,mlh,nlh,ngopr,ntypsy,jri,pos,zatom,
     >           nstr,nstr2,lmax,kv2,kv3,mrot,tau,rmsh,invtab,
     >           amat,bmat,l_noco,cdnfname)
         ELSE
            IF (slice) THEN
               cdnfname = 'cdn_slice'
            ELSE
               cdnfname = 'cdn1'
            ENDIF
            CALL plotdop(
     >           odi,ods,n3d,jspd,nmzxyd,n2d,memd,nlhd,ntypsd,ntype,
     >           lmaxd,jmtd,ntypd,natd,nmzd,jspins,neq,
     >           nq3,nvac,nmz,nmzxy,nq2,nop,nop2,volint,film,
     >           plpot,score,slice,symor,invs,invs2,z1,delz,
     >           clnu,llh,nmem,mlh,nlh,ngopr,ntypsy,jri,pos,zatom,
     >           nstr,nstr2,lmax,kv2,kv3,mrot,tau,rmsh,invtab,
     >           amat,bmat,l_noco,cdnfname)
         ENDIF
      END IF
      ENDIF ! irank == 0
c
c     --->generate starting charge density
c
      IF (strho) THEN
         total = .false.
c
         CALL stden(irank,isize,
     >              memd,nlhd,ntypsd,k1d,k2d,k3d,n2d,ntypd,nop,
     >              n3d,natd,jmtd,lmaxd,jspd,nmzxyd,nmzd,msh,
     >              nlod,nvac,nspd,nstd,film,zrfs,invs,invs2,
     >              ntype,mx3,nq2,nq3,nmzxy,nmz,z1,delz,omtil,
     >              rmsh,taual,rmt,sk3,dx,bmat,amat,zatom,pos,
     >              jri,nlh,llh,nmem,mlh,lmax,nstr,neq,ntypsy,
     >              nstr2,ig,clnu,tau,mrot,nop2,symor,kv2,sk2,
     >              namat,name,icorr,total,krla,kv3,ngopr,phi2,
     >              ig2,vol,volint,volmts,area,jspins,bmu,nwd,
     >              igrd,ndvgrd,idsprs,isprsv,sprsv,nlo,llo,
     >              l_dulo,nwdd,ellow,elup,sigma,lepr,invtab,odi,ods)
c
         CALL cpu_time(time2)
         IF (irank == 0)
     &   CALL outtime('generation of start-density:',time2-time1)
      END IF
      IF (irank == 0) THEN
c
c     --->generate spin polarized charge density
c

      IF (swsp) THEN
         CALL cdnsp(
     >              ntype,jspins,nmz,nmzxy,ntypsd,ntypsy,
     >              nvac,odi%nq2,nq3,jmtd,natd,neq,
     >              jri,nlh,zatom,rmt,dx,film,invs,invs2,
     >              z1,delz,rmsh,msh,bmu,namat,n_u)
c
         CALL cpu_time(time2)
         CALL outtime('optional: spin polarized density:',time2-time1)
      END IF
c
c     --->flip magnetic moments
c
      IF (lflip) THEN
         CALL flipcdn(
     >                ntype,jspins,nmz,nmzxy,ntypsd,ntypsy,
     >                nvac,odi%nq2,nq3,film,invs,invs2,l_noco,z1,delz,
     >                jri,nlh,zatom,rmt,dx,natd,neq,nwd,
     >                nflip,lda_u,n_u,bmu,namat,nlo)
c
         CALL cpu_time(time2)
         CALL outtime('optional: flip magnetic moments:',time2-time1)
      END IF

      IF (l_u2f) THEN
        CALL u2f(
     >           nq3,odi%nq2,jspins,l_noco,jri,nlh,ntype,nmz,nmzxy,
     >           ntypsd,ntypsy,z1,delz,invs,invs2,film,
     >           natd,neq,nvac,zatom,namat,dx,rmt,nwd)
c
        CALL cpu_time(time2)
        CALL outtime('optional: conversion to formatted:',time2-time1)
      ENDIF

      IF (l_f2u) THEN
        CALL f2u(
     >           nq3,odi%nq2,jspins,l_noco,jri,nlh,ntype,nmz,nmzxy,
     >           ntypsd,ntypsy,z1,delz,invs,invs2,film,
     >           natd,neq,nvac,zatom,namat,dx,rmt,nwd)
c
        CALL cpu_time(time2)
        CALL outtime('optional: conversion to unformatted:',time2-time1)
      ENDIF

      IF (l_bmt) THEN
        CALL bmt(
     >           nq3,odi%nq2,jspins,l_noco,jri,nlh,ntype,nmz,nmzxy,
     >           ntypsd,ntypsy,z1,delz,invs,invs2,film,
     >           natd,neq,nvac,zatom,namat,dx,rmt)
      ENDIF

      ENDIF ! irank == 0

      END SUBROUTINE optional
      END MODULE m_optional
