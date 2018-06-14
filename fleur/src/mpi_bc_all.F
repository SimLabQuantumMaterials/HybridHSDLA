      MODULE m_mpi_bc_all
      CONTAINS
      SUBROUTINE mpi_bc_all(
     > n3d,memd,nlhd,ntypsd,ntypd,nwdd,natd,n2d,nop,k1d,k2d,k3d,
     > nkptd,nqptd,jspd,lmaxd,jmtd,nlod,irecl,irecl0,irecl1,l_zref,
     X nop2,llod,ntype,lpr,form66,eig66,eonly,secvar,isec1,l_f,
     X cdinf,dos,vacdos,integ,pallst,nn3d,iplot,strho,swsp,lflip,l_f2u,
     X l_u2f,l_bmt,slice,zrfs,film,invs,invs2,jspins,nvac,neq,lmax,
     X bk,wtkpt,rkmax,ellow,elup,nwd,nkpt,nsymt,jri,rmt,dx,ntypsy,
     X nlh,nq2,nq3,nmz,nmzxy,ustep,bbmat,bmat,invsat,ncv,pot8,
     X invsatnr,zatom,nstr,nstr2,kv2,kv3,sk3,igfft,pgfft,kmxq_fft,
     X nq3_fft,kq1_fft,kq2_fft,kq3_fft,zelec,izlay,ello0,ig2,ig,
     X rgphs,area,delz,z1,omtil,volint,vol,volmts,e1s,e2s,kk,layers,
     X nnne,ndir,mx1,mx2,mx3,layerd,lnonsph,el0,evac0,lepr,irank,
     X itmax,alpha,nmem,mlh,llh,rmsh,clnu,ngopr,mrot,tau,taual,invarop,
     X invarind,multab,invtab,amat,ncst,d_wgn,l_soc,soc_opt,theta,
     X phi,nlo,llo,nlol,lo1l,l_dulo,ulo_der,ctail,skiplo,tworkf,
     X nstars,nstm,starcoeff,locx,locy,l_noco,l_ss,l_mperp,l_constr,
     X alph,beta,b_con,l_relax,qss,lda_u,igrd,lapw_l,mix_b,kmxq2_fft,
     X nq2_fft,l_geo,neigd,sigma,zsigma,sig_b,pos,sk2,phi2,odd,MM,
     X nn2d_1,n2d_1,nop_1,igfft1,pgfft1,ig1,kv1,nstr1,ngopr1,tau1,
     X invtab1,multab1,mrot1,l_J,l_disp,sso_opt,qj,nmagn,mtypes,gw,
     X gw_neigd,l_magn,thetaJ,magtype,nmagtype,TYP_REAL,TYP_COMPLEX)
c
c**********************************************************************
      USE m_types, ONLY: t_utype
      USE m_od_types, ONLY : od_dim
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      TYPE (od_dim), INTENT(INOUT) :: odd
      INTEGER,INTENT(IN) :: irank,layerd,nn3d,irecl0,llod,nlod,neigd
      INTEGER,INTENT(IN) :: n3d,memd,nlhd,ntypsd,ntypd,nwdd,natd,n2d
      INTEGER,INTENT(IN) :: nop,k1d,k2d,k3d,nkptd,jspd,lmaxd,jmtd,nqptd
      INTEGER,INTENT(IN) :: n2d_1,nn2d_1,nop_1,MM
c+odim
      INTEGER,INTENT(INOUT) :: ig1(-odd%k3:odd%k3,-odd%M:odd%M)
      INTEGER,INTENT(INOUT) :: kv1(2,odd%n2d),nstr1(odd%n2d)
      INTEGER,INTENT(INOUT) :: ngopr1(natd)
      INTEGER,INTENT(INOUT) :: invtab1(odd%nop),multab1(odd%nop,odd%nop)
      REAL,   INTENT(INOUT) :: mrot1(3,3,odd%nop),tau1(3,odd%nop)
      INTEGER,INTENT(INOUT) :: igfft1(0:odd%nn2d-1,2)
      REAL,   INTENT(INOUT) :: pgfft1(0:odd%nn2d-1)
c-odim
      INTEGER,INTENT(INOUT) :: lpr,ntype,nwd,isec1,itmax,irecl1,nsymt
      INTEGER,INTENT(INOUT) :: nq2,nq3,nmz,nmzxy,lepr,jspins,nvac,irecl
      INTEGER,INTENT(INOUT) :: kk,layers,nnne,ndir,mx1,mx2,mx3,nop2
      INTEGER,INTENT(INOUT) :: nstars,nstm,igrd,gw,gw_neigd
      REAL,   INTENT(INOUT) :: omtil,area,delz,z1,e1s,e2s,theta,phi,vol
      REAL,   INTENT(INOUT) :: sigma,zsigma,alpha,mix_b,volint,tworkf
      LOGICAL,INTENT(INOUT) :: eonly,form66,secvar,zrfs,film,invs,invs2
      LOGICAL,INTENT(INOUT) :: l_f,cdinf,dos,vacdos,integ,l_zref
      LOGICAL,INTENT(INOUT) :: iplot,strho,swsp,lflip,l_f2u,l_u2f,l_bmt
      LOGICAL,INTENT(INOUT) :: slice,l_soc,l_noco,l_ss,l_mperp,l_constr
      LOGICAL,INTENT(INOUT) :: pallst,starcoeff,ctail,pot8
      INTEGER,INTENT(OUT) :: TYP_REAL,TYP_COMPLEX
c--- J<
      INTEGER,INTENT(INOUT) :: nmagn,mtypes
      INTEGER,INTENT(INOUT) :: magtype(ntypd),nmagtype(ntypd)
      LOGICAL,INTENT(INOUT) :: l_J,l_disp
      REAL,   INTENT(INOUT) :: thetaJ
      REAL,   INTENT(INOUT) :: qj(3,nqptd)
      LOGICAL, INTENT (INOUT) :: l_magn(ntypd)
c--- J>
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (INOUT) :: ustep(n3d)
      COMPLEX, INTENT (INOUT) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (INOUT) :: d_wgn(-3:3,-3:3,3,nop)
      INTEGER, INTENT (INOUT) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (INOUT) :: nmem(0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      INTEGER, INTENT (INOUT) :: jri(ntypd),ntypsy(natd),neq(ntypd)
      INTEGER, INTENT (INOUT) :: lnonsph(ntypd),nkpt(nwdd),lmax(ntypd)
      INTEGER, INTENT (INOUT) :: invsat(natd),ngopr(natd),mrot(3,3,nop)
      INTEGER, INTENT (INOUT) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      INTEGER, INTENT (INOUT) :: invsatnr(natd),izlay(layerd,2)
      INTEGER, INTENT (INOUT) :: nstr(n3d),skiplo(ntypd,jspd),ncv(ntypd)
      INTEGER, INTENT (INOUT) :: invarop(natd,nop),invarind(natd)
      INTEGER, INTENT (INOUT) :: multab(nop,nop),invtab(nop),kv2(2,n2d)
      INTEGER, INTENT (INOUT) :: kq2_fft(nwdd),kq3_fft(nwdd),nstr2(n2d)
      INTEGER, INTENT (INOUT) :: kv3(3,n3d),nq3_fft(nwdd),nq2_fft(nwdd)
      INTEGER, INTENT (INOUT) :: kmxq_fft(nwdd),kmxq2_fft(nwdd)
      INTEGER, INTENT (INOUT) :: igfft(0:nn3d-1,2),kq1_fft(nwdd)
      INTEGER, INTENT (INOUT) :: nlo(ntypd),llo(nlod,ntypd),ncst(ntypd)
      INTEGER, INTENT (INOUT) :: lo1l(0:llod,ntypd),nlol(0:llod,ntypd)
      INTEGER, INTENT (INOUT) :: lapw_l(ntypd),ulo_der(nlod,ntypd)
      REAL,    INTENT (INOUT) :: zatom(ntypd),sig_b(2)
      REAL,    INTENT (INOUT) :: pgfft(0:nn3d-1),zelec(nwdd)
      REAL,    INTENT (INOUT) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      REAL,    INTENT (INOUT) :: bk(3,nkptd,nwdd),wtkpt(nkptd,nwdd)
      REAL,    INTENT (INOUT) :: ellow(nwdd),elup(nwdd),rkmax(nwdd)
      REAL,    INTENT (INOUT) :: rmt(ntypd),dx(ntypd),evac0(2,jspd,nwdd)
      REAL,    INTENT (INOUT) :: sk3(n3d),volmts(ntypd),tau(3,nop)
      REAL,    INTENT (INOUT) :: amat(3,3),pos(3,natd)
      REAL,    INTENT (INOUT) :: bmat(3,3),bbmat(3,3),taual(3,natd)
      REAL,    INTENT (INOUT) :: el0(0:lmaxd,ntypd,jspd,nwdd)
      REAL,    INTENT (INOUT) :: ello0(nlod,ntypd,jspd)
      REAL,    INTENT (INOUT) :: rmsh(jmtd,ntypd)
      REAL,    INTENT (INOUT) :: alph(ntypd),beta(ntypd),qss(3)
      REAL,    INTENT (INOUT) :: b_con(2,ntypd),locx(2),locy(2)
      LOGICAL, INTENT (INOUT) :: l_geo(ntypd),l_relax(ntypd)
      LOGICAL, INTENT (INOUT) :: eig66(2),soc_opt(ntypd+2),sso_opt(2) 
      LOGICAL, INTENT (INOUT) :: l_dulo(nlod,ntypd)
      TYPE (t_utype), INTENT (INOUT) :: lda_u(ntypd)
      REAL,    INTENT (INOUT) :: sk2(n2d),phi2(n2d)
C     ..
C     .. Local Scalars ..
      INTEGER n,n_p,reclength
      CHARACTER*80    :: io_string
      INTEGER, PARAMETER :: n2max=13 ! some STM stuff in cdnval
C     ..
C     .. Local Arrays ..
      INTEGER i(31),ierr(3)
      REAL    r(19)
      LOGICAL l(37)
C     ..
C     .. External Subroutines.. 
      EXTERNAL MPI_BCAST,MPI_TYPE_CONTIGUOUS,MPI_TYPE_COMMIT
#ifdef CPP_T90
      EXTERNAL asnunit
#endif
C     ..
      INCLUDE 'mpif.h'
c
#ifdef CPP_T90
      CALL MPI_TYPE_CONTIGUOUS(1,MPI_REAL,TYP_REAL,ierr)
      CALL MPI_TYPE_CONTIGUOUS(1,MPI_COMPLEX,TYP_COMPLEX,ierr)
#else
      CALL MPI_TYPE_CONTIGUOUS(1,MPI_DOUBLE_PRECISION,TYP_REAL,ierr)
      CALL MPI_TYPE_CONTIGUOUS(1,MPI_DOUBLE_COMPLEX,TYP_COMPLEX,ierr)
#endif
      CALL MPI_TYPE_COMMIT(TYP_REAL,ierr)
      CALL MPI_TYPE_COMMIT(TYP_COMPLEX,ierr)

      IF (irank.EQ.0) THEN
      i(1)=irecl1 ; i(2)=lpr ; i(3)=ntype ; i(5)=nwd ; i(6)=isec1
      i(7)=nq2 ; i(8)=nq3 ; i(9)=nmz ; i(10)=nmzxy ; i(11)=lepr 
      i(12)=jspins ; i(13)=nvac ; i(14)=itmax ; i(15)=kk ; i(16)=layers
      i(17)=nnne ; i(18)=ndir ; i(19)=mx1 ; i(20)=mx2 ; i(21)=mx3
      i(22)=igrd ; i(23) = nop2 ; i(24) = nsymt ; i(25) = irecl
      i(26)=nstars ; i(27)=nstm ; i(28)=odd%nq2 ; i(29)=odd%nop
      i(30)=gw ; i(31)=gw_neigd
      r(1)=omtil ; r(2)=area ; r(3)=delz ; r(4)=z1 ; r(5)=alpha
      r(6)=e1s ; r(7)=e2s ; r(8)=theta ; r(9)=phi ; r(10)=tworkf 
      r(11)=locx(1) ; r(12)=locx(2); r(13)=locy(1) ; r(14)=locy(2)
      r(15)=sigma ; r(16)=zsigma ; r(17)=mix_b; r(18)=vol; r(19)=volint
      l(1)=eonly ; l(2)=form66 ; l(3)=secvar ; l(4)=zrfs ; l(5)=film
      l(6)=invs ; l(7)=invs2 ; l(8)=l_bmt ; l(9)=l_f ; l(10)=cdinf
      l(11)=dos ;l(12)=l_zref ; l(13)=vacdos ; l(14)=integ ; l(15)=iplot
      l(16)=strho ; l(17)=swsp ; l(18)=lflip ; l(19)=l_f2u ; l(20)=l_u2f
      l(21)=pallst ; l(22)=slice ; l(23)=l_soc ; l(24)=starcoeff
      l(25)=l_noco ; l(26)=l_ss; l(27)=l_mperp; l(28)=l_constr
      l(29)=odd%d1 ; l(30)=l_J ; l(31)=l_disp ; l(32)=ctail
      l(33)=eig66(1) ; l(34)=eig66(2) ; l(35)=sso_opt(1)
      l(36)=sso_opt(2) ; l(37)=pot8
      ENDIF
c
      CALL MPI_BCAST(i,31,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      gw=i(30) ; gw_neigd=i(31)
      nstars=i(26) ; nstm=i(27) ; odd%nq2=i(28) ; odd%nop=i(29)
      igrd=i(22) ; nop2=i(23) ; nsymt = i(24) ; irecl = i(25)
      nnne=i(17) ; ndir=i(18) ; mx1=i(19) ; mx2=i(20) ; mx3=i(21)
      jspins=i(12) ; nvac=i(13) ; itmax=i(14) ; kk=i(15) ; layers=i(16)
      nq2=i(7) ; nq3=i(8) ; nmz=i(9) ; nmzxy=i(10) ; lepr=i(11)
      irecl1=i(1) ; lpr=i(2) ; ntype=i(3) ; nwd=i(5) ; isec1=i(6)
c
      CALL MPI_BCAST(r,19,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      sigma=r(15) ; zsigma=r(16); mix_b=r(17); vol=r(18); volint= r(19)
      locx(1)=r(11); locx(2)=r(12); locy(1)=r(13); locy(2)=r(14)
      e1s=r(6) ; e2s=r(7) ; theta=r(8) ; phi=r(9) ; tworkf=r(10)
      omtil=r(1) ; area=r(2) ; delz=r(3) ; z1=r(4) ; alpha=r(5)
c
      CALL MPI_BCAST(l,37,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      sso_opt(2)=l(36) ; pot8=l(37)
      eig66(1)=l(33) ; eig66(2)=l(34) ; sso_opt(1)=l(35)
      odd%d1=l(29) ; l_J=l(30) ; l_disp=l(31) ; ctail=l(32)
      l_noco=l(25) ; l_ss=l(26) ; l_mperp=l(27) ; l_constr=l(28)
      pallst=l(21) ; slice=l(22) ; l_soc=l(23) ; starcoeff=l(24)
      strho=l(16) ; swsp=l(17) ; lflip=l(18) ; l_f2u=l(19) ; l_u2f=l(20)
      dos=l(11) ;l_zref=l(12) ; vacdos=l(13) ; integ=l(14) ; iplot=l(15)
      invs=l(6) ; invs2=l(7) ; l_bmt=l(8) ; l_f=l(9) ; cdinf=l(10)
      eonly=l(1) ; form66=l(2) ; secvar=l(3) ; zrfs=l(4) ; film=l(5)
c
c -> Broadcast the arrays:
 
      CALL MPI_BCAST(ustep,n3d,TYP_COMPLEX,0,MPI_COMM_WORLD,ierr)
      n = memd*(nlhd+1)*ntypsd
      CALL MPI_BCAST(clnu,n,TYP_COMPLEX,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mlh,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nlh,ntypsd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      n = (nlhd+1)*ntypsd
      CALL MPI_BCAST(nmem,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(llh,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(jri,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ncv,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ntypsy,natd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(neq,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lnonsph,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmax,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(invsat,natd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(invsatnr,natd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ngopr,natd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mrot,9*nop,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ig2,n3d,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      n = (2*k1d+1)*(2*k2d+1)*(2*k3d+1)
      CALL MPI_BCAST(ig,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rgphs,n,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ellow,nwdd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(elup,nwdd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rkmax,nwdd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rmt,ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(volmts,ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dx,ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sk3,n3d,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(evac0,2*jspd*nwdd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(amat,9,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bmat,9,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bbmat,9,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(taual,3*natd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(pos,3*natd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tau,3*nop,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      n = (lmaxd+1)*ntypd*jspd*nwdd
      CALL MPI_BCAST(el0,n,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      n = nlod*ntypd*jspd
      CALL MPI_BCAST(ello0,n,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rmsh,jmtd*ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
c
      CALL MPI_BCAST(nkpt,nwdd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bk,3*nkptd*nwdd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(wtkpt,nkptd*nwdd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
c
      n = natd*nop
      CALL MPI_BCAST(invarop,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(multab,nop**2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(invarind,natd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(invtab,nop,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(invsatnr,natd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kq2_fft,nwdd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kq3_fft,nwdd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kv2,2*n2d,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kv3,3*n3d,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nq3_fft,nwdd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nq2_fft,nwdd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kmxq_fft,nwdd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kmxq2_fft,nwdd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(izlay,layerd*2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nstr,n3d,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nstr2,n2d,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(igfft,nn3d*2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kq1_fft,nwdd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(pgfft,nn3d,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zatom,ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sig_b,2,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(zelec,nwdd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ncst,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nlo,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      n =  nlod*ntypd
      CALL MPI_BCAST(llo,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ulo_der,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(l_dulo,n,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      n = (llod+1)*ntypd
      CALL MPI_BCAST(lo1l,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nlol,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      n = jspd*ntypd
      CALL MPI_BCAST(skiplo,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(alph,ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(beta,ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(b_con,ntypd*2,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(l_relax,ntypd,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(l_geo,ntypd,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(soc_opt,ntypd+2,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(qss,3,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lda_u(:)%l,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lda_u(:)%u,ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lda_u(:)%j,ntypd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lda_u(:)%l_amf,ntypd,MPI_LOGICAL,0,
     +                                           MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lapw_l,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      n = 7*7*3*nop
      CALL MPI_BCAST(d_wgn,n,TYP_COMPLEX,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nstr1,n2d_1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sk2,n2d,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      IF (odd%d1) THEN
        CALL MPI_BCAST(phi2,n2d,TYP_REAL,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(tau1,3*nop_1,TYP_REAL,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(mrot1,9*nop_1,TYP_REAL,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(kv1,2*n2d_1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(ngopr1,natd,MPI_INTEGER,0,
     &       MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(igfft1,nn2d_1*2,MPI_INTEGER,0,
     &       MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(pgfft1,nn2d_1,TYP_REAL,0,MPI_COMM_WORLD,ierr)
        n = (2*k3d + 1)*(2*MM +1)
        CALL MPI_BCAST(ig1,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(invtab1,nop_1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(multab1,2*nop_1,MPI_INTEGER,
     &                                      0,MPI_COMM_WORLD,ierr)
      ENDIF
c--- J<
      IF (l_J) THEN
      CALL MPI_BCAST(nmagn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(mtypes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(thetaJ,1,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(qj,3*nqptd,TYP_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(l_magn,ntypd,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(magtype,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nmagtype,ntypd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      ENDIF
c--- J>
c
c
#ifdef CPP_T90
c
c-->  on cray assign the partitions of the global direct access files:
c
      CALL asnrm(ierr)
      WRITE (io_string,
     +     '("-O -n ",I5.5," -F global.privpos.page_size=",
     +     I5.5,".num_pages=1")')
     +     (irecl0*nkptd)/4096+1,irecl0/(8*512)+1
      CALL asnunit(28,io_string,ierr)
      IF (jspins.EQ.2) CALL asnunit(38,io_string,ierr)

      IF (irank.EQ.0) WRITE (*,*) 'for unit 28:',io_string
      IF (irank.EQ.0) WRITE (*,*) 'unit 28 assign with err=',ierr(1)

      WRITE (io_string,
     +     '("-O -n ",I5.5," -F global.privpos.page_size=",
     +     I5.5,".num_pages=1")')
     +     (irecl1*nkptd)/4096+1,irecl1/(8*512)+1
      CALL asnunit(29,io_string,ierr)

      IF (irank.EQ.0) WRITE (*,*) 'for unit 29:',io_string
      IF (irank.EQ.0) WRITE (*,*) 'unit 29 assign with err=',ierr(1)

      IF (.not.form66) THEN
        n_p = irecl/(128*4096)+1 
        WRITE (io_string,
     +       '("-O -n ",I9.9," -F global.privpos.page_size=",
     +       I5.5,".num_pages=",I5.5)')
     +       (irecl*nkptd*jspd)/4096+1,128,n_p
        CALL asnunit(66,io_string,ierr)

        IF (irank.EQ.0) WRITE (*,*) 'for unit 66:',io_string
        IF (irank.EQ.0) WRITE (*,*) 'unit 66 assign with err=',ierr(1)
      ENDIF

      IF (dos.or.vacdos.or.cdinf) THEN
c
c      calculate record length & assign tmp_dos-file
c
         reclength = neigd*6 + 4*ntypd*neigd + neigd*layerd*2*jspd
         reclength = 8*(reclength + 80*neigd*layerd +5)
c
         WRITE (io_string,
     +     '(" -O -n ",I9.9," -F global.privpos.page_size=",
     +     I5.5,".num_pages=1")')
     +     (reclength*nkptd*jspins)/4096+1,reclength/(8*512)+1
         CALL asnunit(84,io_string,ierr)
         IF (irank.EQ.0) WRITE (*,*) 'for unit 84:',io_string
         IF (irank.EQ.0) WRITE (*,*) 'unit 84 assign with err=',ierr(1)
c
      ENDIF
      IF (nstm.EQ.3) THEN
c
c+dw  if on T3E assign to tmp_vacwave file
c
         reclength = 8*( neigd + 4*n2max*neigd + 4 )
         WRITE (io_string,
     +     '(" -O -n ",I9.9," -F global.privpos.page_size=",
     +     I5.5,".num_pages=1")')
     +     (reclength*nkptd)/4096+1,reclength/(8*512)+1
         CALL asnunit(88,io_string,ierr)
         IF (irank.EQ.0) WRITE (*,*) 'for unit 88:',io_string
         IF (irank.EQ.0) WRITE (*,*) 'unit 88 assign with err=',ierr(1)
      ENDIF
#endif

      RETURN
      END SUBROUTINE mpi_bc_all
      END MODULE m_mpi_bc_all
