      MODULE m_exchpbe
c----------------------------------------------------------------------
c  pbe exchange for a spin-unpolarized electronic system
c  k burke's modification of pw91 codes, may 14, 1996
c  modified again by k. burke, june 29, 1996, with simpler fx(s)
c----------------------------------------------------------------------
c references:
c [a]j.p.~perdew, k.~burke, and m.~ernzerhof, submiited to prl, may96
c [b]j.p. perdew and y. wang, phys. rev.  b {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (e).
c [c] B.~Hammer, L.~B.~Hansen and J.~K.~Norskov PRB 59 7413 (1999)
c----------------------------------------------------------------------
      CONTAINS
      SUBROUTINE exchpbe(
     >                   icorr,rho,s,u,v,lgga,lpot,
     <                   ex,vx)

      IMPLICIT NONE

      ! .. Arguments
      INTEGER, INTENT (IN) :: icorr
      INTEGER, INTENT (IN) :: lgga ! =0=>don't put in gradient corrections, just lda
      INTEGER, INTENT (IN) :: lpot ! =0=>don't get potential and don't need u and v
      REAL,    INTENT (IN) :: rho  ! density
      REAL,    INTENT (IN) :: s    ! abs(grad rho)/(2*kf*rho), where kf=(3 pi^2 rho)^(1/3)
      REAL,    INTENT (IN) :: u    ! (grad rho)*grad(abs(grad rho))/(rho**2 * (2*kf)**3)
      REAL,    INTENT (IN) :: v    ! (laplacian rho)/(rho*(2*kf)**2) [pw86(24)]
      REAL,    INTENT (OUT) :: ex,vx ! exchange energy per electron (ex) and potential (vx)

      ! .. local variables ..
      REAL :: um,uk,ul,exunif,fs,fss,fxpbe,p0,s2
      REAL :: xwu,css,dxwu,ddx               ! for wu-cohen
      REAL, PARAMETER :: teo = 10.e0/81.e0   ! for wu-cohen
      REAL, PARAMETER :: cwu = 0.0079325     ! for wu-cohen
      REAL, PARAMETER :: thrd=1.e0/3.e0
      REAL, PARAMETER :: thrd4=4.e0/3.e0
      REAL, PARAMETER :: ax=-0.738558766382022405884230032680836e0  ! -0.75*(3/pi)^(1/3)

c----------------------------------------------------------------------
c uk, ul defined after [a](13)  (for icorr==7)
c----------------------------------------------------------------------
      IF ((icorr == 7).OR.(icorr == 10).OR.(icorr == 11)) THEN
        uk=0.8040
      ELSEIF (icorr.EQ.8) THEN
        uk=1.2450
      ELSEIF (icorr.EQ.9) THEN    ! changed to [c]
        uk=0.8040
      ELSE
        WRITE (6,'(//'' icorr is not correctly transferred. icorr='',i5)
     &    ') icorr
        STOP
      ENDIF
      IF (icorr == 11) THEN  ! pbe_sol
        um=0.123456790123456e0
      ELSE
        um=0.2195149727645171e0
      ENDIF

      ul=um/uk
c----------------------------------------------------------------------
c construct lda exchange energy density: 
                          ! e_x[unif] = -0.75 * (3/pi)^(1/3) * rho^(4/3)  
c----------------------------------------------------------------------
      exunif = ax*rho**thrd
      IF (lgga.EQ.0) THEN
          ex = exunif
          vx = ex*thrd4
          RETURN
      ENDIF
c----------------------------------------------------------------------
c construct pbe enhancement factor
c       e_x[pbe]=e_x[unif]*fxpbe(s)
c       fxpbe(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c----------------------------------------------------------------------
      s2 = s*s
c+gu
      IF (icorr.EQ.7 .OR. icorr.EQ.8 .OR. icorr.EQ.11) THEN
        p0 = 1.e0 + ul*s2
        fxpbe = 1e0 + uk - uk/p0
      ELSEIF (icorr.EQ.9) THEN
        p0 = exp( - ul*s2 )
        fxpbe = 1e0 + uk * ( 1e0 - p0 )
      ELSEIF (icorr.EQ.10) THEN
        css = 1+cwu*s2*s2
        xwu = teo*s2 + (um-teo)*s2*exp(-s2) + log(css)
        p0 = 1.e0 + xwu/uk
        fxpbe = 1e0 + uk - uk/p0
      ENDIF
c-gu
      ex = exunif*fxpbe
      IF (lpot.EQ.0) RETURN
c----------------------------------------------------------------------
c  energy done. now the potential:
c  find first and second derivatives of fx w.r.t s.
c  fs=(1/s)*d fxpbe/ ds
c  fss=d fs/ds
c----------------------------------------------------------------------
c+gu
      IF (icorr.EQ.7 .OR. icorr.EQ.8 .OR. icorr.EQ.11) THEN
        fs = 2.e0*uk*ul/ (p0*p0)
        fss = -4.e0*ul*s*fs/p0
      ELSEIF (icorr.EQ.9) THEN
        fs = 2.e0*ul*p0
        fss = -2.e0*ul*s*fs
      ELSEIF (icorr.EQ.10) THEN
        dxwu = 2*teo + 2*(um-teo)*exp(-s2)*(1-s2) + 4*cwu*s2/css
        fs = dxwu / (p0*p0)
        ddx = 4*s*((um-teo)*exp(-s2)*(s2-2)+2*cwu*(1-cwu*s2*s2)/css**2)
        fss = ( ddx - 2*s*dxwu*dxwu/(p0*uk) ) / (p0*p0)
      ENDIF
c-gu
c----------------------------------------------------------------------
c calculate potential from [b](24)
c----------------------------------------------------------------------
      vx = exunif* (thrd4*fxpbe- (u-thrd4*s2*s)*fss-v*fs)

      RETURN
      END SUBROUTINE exchpbe
      END MODULE m_exchpbe
