      MODULE m_vxcl91
c.....-----------------------------------------------------------------
c.....(local) pw91 exchange-correlation potential in hartree.
c.....------------------------------------------------------------------
      CONTAINS
      SUBROUTINE vxcl91(
     >                   jspins,mirm,irmx,rh,agr,agru,agrd,
     >                   g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,
     <                   vxc,
     >                   isprsv,sprsv)

      USE m_corl91
      USE m_corg91
      USE m_xch91
      USE m_constants, ONLY: pimach

      IMPLICIT NONE

      ! .. Arguments ..
      INTEGER, INTENT (IN) :: jspins,irmx,mirm,isprsv
      REAL,    INTENT (IN) :: sprsv

      REAL,    INTENT (IN) :: rh(mirm,jspins)
      REAL,    INTENT (IN) :: agr(mirm),agru(mirm),agrd(mirm)
      REAL,    INTENT (IN) :: g2r(mirm),g2ru(mirm),g2rd(mirm)
      REAL,    INTENT (IN) :: gggr(mirm),gggru(mirm)
      REAL,    INTENT (IN) :: gggrd(mirm),gzgr(mirm)
      REAL,    INTENT (OUT) :: vxc(mirm,jspins)

      ! .. local variables
      REAL pi,ro,zta,alf,alfc,c13,c23,c43,c53,
     +     cedg,cedl,dbrod,dbrou,dsprs,ec,ecrs,eczta,fk,gz,
     +     ro13,ro2,rod,rod3,rod43,rod53,rou,rou3,rou43,rou53,
     +     rs,sd,sk,su,tc,td,tksg,tu,uc,ud,uu,vc,
     +     vcgd,vcgu,vcld,vclu,vxgd,vxgu,vxld,
     +     vxlu,wc,xced,xedg,xedgd,xedgu,xedl,xedld,xedlu,
     +     rou13,rod13,xcptu,xcptd

      INTEGER i
      REAL, PARAMETER :: sml = 1.e-14

      pi = pimach()

      DO 300 i = 1,irmx

        IF (jspins.eq.1) THEN
          rou=rh(i,1)/2
          rod=rou
        ELSE
          rou=rh(i,1)
          rod=rh(i,jspins)
        ENDIF
          ro=rou+rod
          zta=(rou-rod)/ro

          if(zta.gt.1.e0-sml) zta = 1.e0 - sml
          if(zta.lt.-1.e0-sml) zta = -1.e0 + sml

c.....
c       vxlu,vxld,vxgu,vxgd: exchange potential in ry.(local,grad),
cc        (up,dw).
c       vclu,vcld,vcgu,vcgd: correl. potential in ry.(local,grad),
cc        (up,dw).
c       all later in hartree.
c.....
          vxlu = 0.0e0
          vclu = 0.0e0
          vxld = 0.0e0
          vcld = 0.0e0
          vxgu = 0.0e0
          vcgu = 0.0e0
          vxgd = 0.0e0
          vcgd = 0.0e0

c.....
          if(ro.lt.sml) then
              ro = sml
              zta = 0.e0
              go to 200
          endif
c.....
          c13 = 1.e0/3.e0
          c23 = 2.e0/3.e0
          c43 = 4.e0/3.e0
          c53 = 5.e0/3.e0
c.....  alf=-3*(3/4*pai)**(1/3).
          alf = -1.861051473e0
c.....
          ro2 = ro*ro
          ro13 = ro**c13
c.....
          rou3 = rou**3
          rou13 = rou**c13
          rou43 = rou**c43
c.....
          rod3 = rod**3
          rod13 = rod**c13
          rod43 = rod**c43
c.....
c       gr2=drr*drr
c       gr2u=drru**2
c       drrd=drr-drru
c       gr2d=drrd**2
c       ddrrd=ddrr-ddrru
c.....
c.....  gz: for wang-perdew ssf.
c.....
          gz = ((1.e0+zta)**c23+ (1.e0-zta)**c23)/2.e0
c.....
          rs = 0.620350491e0/ro13
c.....
c.....  exchange-potential, vxp,vxlu,vxld: v-exchange-(para,up,dw).
          vxlu = c43*alf*rou13
          vxld = c43*alf*rod13

c     stop911
c.....
c....  .gradient correction.
c.....
c         if(abs(agr(i)).lt.sml) go to 200
c.....

c.....
c         dsprs = 1.e0
c         if(idsprs.eq.1) dsprs = 1.e-19
          dsprs = 1.e-19
          if(isprsv.eq.1) dsprs = dsprs*sprsv

c.....
c      agr,agru,agrd: abs(grad(rho)), for all, up, and down.
cc     gr2,gr2u,gr2d: grad(rho_all)**2, grad(rho_up)**2, grad(rho_d)**2.
c      g2r,g2ru,g2rd: laplacian rho_all, _up and _down.
c      gggru,-d: grad(rho)*grad(abs(grad(rho))) for all,up and down.
c      grgru,-d: grad(rho_all)*grad(rhor_up) and for down.

c         g2r=ddrr+2*drr/rv
c.....
          rou53 = rou**c53
c.....
c.....  edrru: d(abs(d(rou)/dr))/dr, edrrd for down.
c         edrru=ddrru
c         if(drru.lt.0.) edrru=-ddrru
c.....
c         agr,agbru,-d: abs(grad(rho)),for rou, rod.
c         gggru,-d: grad(rho)*grad(abs(grad(rho))) for up and down.
c.....  su:at ro=2*rou. 1/(2(3*pai**2)**(1/3))*|grad(rou)|/rou**(4/3).
          su = 0.128278244e0*agru(i)/rou43
c         if(su.gt.huges) go to 200
c         g2ru=ddrru+2*drru/rv
          tu = .016455307e0*g2ru(i)/rou53
          uu = 0.002110857e0*gggru(i)/rou3

          dbrou = rou*2

          call xch91(dbrou,su,uu,tu,xedlu,xedgu,vxlu,vxgu)

          vxgu = dsprs*vxgu


c.....
          rod53 = rod**c53
c         edrrd=ddrrd
c         if(drrd.lt.0.) edrrd=-ddrrd

          sd = 0.128278244e0*agrd(i)/rod43
c         if(sd.gt.huges) go to 200

c         g2rd=ddrrd+2*drrd/rv

          td = .016455307e0*g2rd(i)/rod53
          ud = 0.002110857e0*gggrd(i)/rod3


          dbrod = rod*2

          call xch91(dbrod,sd,ud,td,xedld,xedgd,vxld,vxgd)

          vxgd = dsprs*vxgd


c....   cro: c(n) of (6),phys.rev..b33,8822('86). in ry.
c....   dcdr: d(cro)/d(ro).
c.....  0.001625816=1.745*f(=0.11)*cro(rs=0).


c     stop912
c         pw91

          call corl91(rs,zta,ec,vclu,vcld,ecrs,eczta,alfc)

          vclu = vclu*2.e0
          vcld = vcld*2.e0

          fk = 1.91915829e0/rs
          sk = sqrt(4.e0*fk/pi)
          tksg = 2.e0*sk*gz
          tc = agr(i)/ (ro*tksg)
          uc = gggr(i)/ (ro2*tksg**3)
          vc = g2r(i)/ (ro*tksg**2)
          wc = gzgr(i)/ (ro*tksg**2)

          call corg91(fk,sk,gz,ec,ecrs,eczta,rs,zta,tc,uc,vc,wc,cedg,
     +                vcgu,vcgd)

          vcgu = dsprs*vcgu*2.e0
          vcgd = dsprs*vcgd*2.e0
c     stop914

  200     continue

          xcptu = vxlu + vclu + vxgu + vcgu
          xcptd = vxld + vcld + vxgd + vcgd


c     change to hartree.
c     stop915
          vxc(i,1)      = xcptu
          vxc(i,jspins) = xcptd

c     stop913
  300 continue

      END SUBROUTINE vxcl91
      END MODULE m_vxcl91
