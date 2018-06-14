      MODULE m_vintcz
c     *************************************************************
c     z-dependent part of coulomb potential in interstitial       *
c     [-d/2,d/2] region            c.l.fu, r.podloucky            *
c     *************************************************************
!     modified for thick films to avoid underflows gb`06
!---------------------------------------------------------------
      CONTAINS
      COMPLEX FUNCTION vintcz(
     >                        n3d,nmzxyd,n2d,jspd,k1d,k2d,k3d,nmzd,
     >                        z,nrec2,sk3,sk2,nstr2,kv2,nstr,bmat,
     >                        z1,invs,zrfs,delz,nmzxy,mx3,ig,
     >                        psq,vxy,vz,rhobar,sig1dh,vz1dh,alphm)

      USE m_constants, ONLY : pimach
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n3d,nmzxyd,n2d,jspd,k1d,k2d,k3d,nmzd
      INTEGER, INTENT (IN) :: nmzxy,mx3,nrec2
      COMPLEX, INTENT (IN) :: rhobar
      REAL,    INTENT (IN) :: sig1dh,vz1dh,z1,delz,z
      LOGICAL, INTENT (IN) :: invs,zrfs
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: psq(n3d),vxy(nmzxyd,n2d-1,2,jspd)
      COMPLEX, INTENT (IN) :: alphm(n2d,2)
      INTEGER, INTENT (IN) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      INTEGER, INTENT (IN) :: kv2(2,n2d),nstr(n3d),nstr2(n2d)
      REAL,    INTENT (IN) :: vz(nmzd,2,jspd),sk2(n2d),sk3(n3d)
      REAL,    INTENT (IN) :: bmat(3,3)
C     ..
C     .. Local Scalars ..
      COMPLEX argr,sumrr,vcons1,ci,test
      REAL bj0,dh,fit,g,g3,q,qdh,signz,vcons2,zf,tpi,fpi
      REAL maxexp,minexp,e_m,e_p,cos_q,sin_q
      INTEGER ig3n,im,iq,ivac,k1,k2,m0,nrz,nz
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,cmplx,conjg,cos,exp,sin
C     ..
c
      maxexp = log(2.0)*MAXEXPONENT(tpi)
      minexp = log(2.0)*MINEXPONENT(tpi)

      tpi = 2 * pimach()
      fpi = 2 * tpi
      dh = z1
      ci = (0.0,1.0)
      sumrr = (0.,0.)
      vintcz = (0.,0.)
c--->    if z is in the vacuum, use vacuum representations (m.w.)
      IF (abs(z).GE.z1) THEN
         ivac = 1
         IF (z.LT.0.0) THEN
            ivac = 2
            IF (invs .OR. zrfs) ivac = 1
         END IF
         zf = (abs(z)-z1)/delz + 1.0
         im = zf
         q = zf - im
         IF (nrec2.EQ.1) THEN
            fit = 0.5* (q-1.)* (q-2.)*vz(im,ivac,1) -
     +            q* (q-2.)*vz(im+1,ivac,1) +
     +            0.5*q* (q-1.)*vz(im+2,ivac,1)
            vintcz = cmplx(fit,0.0)
         ELSE IF (im+2.LE.nmzxy) THEN
            vintcz = 0.5* (q-1.)* (q-2.)*vxy(im,nrec2-1,ivac,1) -
     +               q* (q-2.)*vxy(im+1,nrec2-1,ivac,1) +
     +               0.5*q* (q-1.)*vxy(im+2,nrec2-1,ivac,1)
            IF ((invs.AND. (.NOT.zrfs)) .AND.
     +          z.LT.0) vintcz = conjg(vintcz)
         END IF
         RETURN
      END IF
c
      IF (nrec2.EQ.1) THEN
         m0 = -mx3
         IF (zrfs .OR. invs) m0 = 0
c     ---->    g=0 coefficient
c           -----> v1(z)
         DO 20 iq = m0,mx3
            IF (iq.EQ.0) GO TO 20
            ig3n = ig(0,0,iq)
c     ----> use only stars within the g_max sphere (oct.97 shz)
            IF (ig3n.ne.0) THEN
            q = iq*bmat(3,3)
            nz = nstr(ig3n)
            sumrr = (0.0,0.0)
c      --> sum over gz-stars
            DO 10 nrz = 1,nz
               signz = 3. - 2.*nrz
               q = signz*q
               qdh = q*dh
               bj0 = sin(qdh)/qdh
               argr = ci*q*z
               sumrr = sumrr + (exp(argr)-exp(ci*qdh))/ (q*q) +
     +                 ci*cos(qdh)* (dh-z)/q + bj0* (z*z-dh*dh)/2.
   10       CONTINUE
            vintcz = vintcz + fpi*psq(ig3n)*sumrr
            ENDIF
   20    CONTINUE
c           -----> v2(z)
         vintcz = vintcz + vz1dh - fpi* (dh-z)*
     +            (sig1dh-rhobar/2.* (dh-z))
c     ---->    (g.ne.0)  coefficients
      ELSE
         m0 = -mx3
         IF (zrfs) m0 = 0
         k1 = kv2(1,nrec2)
         k2 = kv2(2,nrec2)
         DO 40 iq = m0,mx3
            ig3n = ig(k1,k2,iq)
c     ----> use only stars within the g_max sphere (oct.97 shz)
            IF (ig3n.ne.0) THEN
c           -----> v3(z)
               q = iq*bmat(3,3)
               g = sk2(nrec2)
               g3 = sk3(ig3n)
               vcons1 = fpi*psq(ig3n)/ (g3*g3)
               vcons2 = - 1.0 / (2.*g)
! underflow
               IF (-g*(z+dh) >= minexp ) THEN
                 e_m = exp( -g*(z+dh) )
               ELSE
                 e_m = exp( minexp )
               ENDIF
               IF ( g*(z-dh) >= minexp ) THEN
                 e_p = exp( g*(z-dh) )
               ELSE
                 e_p = exp( minexp )
               ENDIF
! underflow
               nz = 1
               IF (zrfs) nz = nstr(ig3n)/nstr2(nrec2)
               sumrr = (0.0,0.0)
c                --> sum over gz-stars
               vacua: DO nrz = 1,nz
                  signz = 3. - 2.*nrz
                  q = signz*q
                  cos_q = cos(q*dh)
                  sin_q = sin(q*dh)
                  sumrr = sumrr + cmplx(cos(q*z),sin(q*z)) + vcons2 *
     +                      ( (g + ci*q) * e_p * (cos_q + ci*sin_q) +
     +                        (g - ci*q) * e_m * (cos_q - ci*sin_q) )
               ENDDO vacua
               vintcz = vintcz + vcons1*sumrr
            ENDIF 
   40    CONTINUE
c  ----> v4(z)
         IF (-g*z >= minexp ) THEN
           e_m = exp( -g*z  )
         ELSE
           e_m = exp( minexp )
         ENDIF
         IF ( g*z  <= maxexp ) THEN
           e_p = exp( g*z  )
         ELSE
           e_p = exp( maxexp )
         ENDIF
         test = e_m*alphm(nrec2-1,2) + e_p*alphm(nrec2-1,1)
         IF ( 2.0 * test == test ) test = cmplx(0.0,0.0)
         vintcz = vintcz + tpi/g* test
      ENDIF

      END FUNCTION vintcz
      END MODULE m_vintcz
