      MODULE m_hnonmuff
c*********************************************************************
c     updates hamiltonian by adding non-spherical matrix elements in
c     the second-variation scheme. usage of tlmplm-nonmuff required
c                r. p  1995
c*********************************************************************
      CONTAINS
      SUBROUTINE h_nonmuff(
     >                     nvd,lmaxd,ntypd,neigd,natd,nop,
     >                     lmd,lmplmd,ntype,neq,lmax,rmt,taual,
     >                     bmat,mrot,ngopr,omtil,invsat,
     >                     nb0,jsp,z,ne,nbasfcn,
     >                     us,dus,uds,duds,
     >                     bkpt,rk,k1,k2,k3,nv,
     X                     h)

      USE m_constants, ONLY : pimach
      USE m_sphbes
      USE m_dsphbs
      USE m_ylm

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nvd,lmaxd,ntypd,neigd,natd,nop
      INTEGER, INTENT (IN) :: jsp,nb0,ne,nv,ntype,lmd,lmplmd,nbasfcn
      REAL,    INTENT (IN) :: omtil
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: k1(nvd),k2(nvd),k3(nvd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: lmax(ntypd),neq(ntypd),ngopr(natd)
      INTEGER, INTENT (IN) :: invsat(natd)
      REAL,    INTENT (IN) :: bkpt(3),rk(nvd),rmt(ntypd),taual(3,natd)
      REAL,    INTENT (IN) :: bmat(3,3)
      REAL,    INTENT (IN) :: us(0:lmaxd,ntypd),uds(0:lmaxd,ntypd)
      REAL,    INTENT (IN) :: dus(0:lmaxd,ntypd),duds(0:lmaxd,ntypd)
      REAL,    INTENT (INOUT) :: h(ne*(ne+1)/2)
#ifdef CPP_INVERSION
      REAL,    INTENT (IN) :: z(nbasfcn,ne)
#else
      COMPLEX, INTENT (IN) :: z(nbasfcn,ne)
#endif
C     ..
C     .. Local Scalars ..
      COMPLEX dtd,dtu,hij,phase,sij,utd,utu
      REAL con1,ff,gg,gs,th,ws,tpi,fpi
      INTEGER l,l1,ll1,lm,lmp,lwn,invsfct
      INTEGER i,im,in,j,k,ke,m,m1,n,na,nn,np,nrec,ii,ij
C     ..
C     .. Local Arrays ..
      COMPLEX,DIMENSION (0:lmplmd) :: tdd,tdu,tud,tuu
      COMPLEX a(neigd,0:lmd),ax(neigd)
      COMPLEX b(neigd,0:lmd),bx(neigd), ylm( (lmaxd+1)**2 )
      INTEGER ind(0:lmd,0:lmd)
      REAL vmult(3),vsmult(3),f(0:lmaxd,nvd),g(0:lmaxd,nvd)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cexp,cmplx,conjg,real,sqrt
C     ..
      tpi = 2 * pimach()
      fpi = 2 * tpi
      con1 = fpi/sqrt(omtil)
c--->    loop over each atom type
      na = 0
      DO 110 n = 1,ntype
         lwn = lmax(n)
c--->    set up wronskians for the matching conditions for each ntype
         DO k = 1,nv
            gs = rk(k)*rmt(n)
            CALL sphbes(
     >                  lwn,gs,
     <                  f(0,k))
            CALL dsphbs(
     >                  lwn,gs,f(0,k),
     <                  g(0,k))
         ENDDO
         DO 30 l = 0,lwn
            ws = uds(l,n)*dus(l,n) - us(l,n)*duds(l,n)
            DO 20 k = 1,nv
               ff = f(l,k)
               gg = rk(k)*g(l,k)
               f(l,k) = con1* (uds(l,n)*gg-ff*duds(l,n))/ws
               g(l,k) = con1* (dus(l,n)*ff-gg*us(l,n))/ws
   20       CONTINUE
   30    CONTINUE
c--->    loop over equivalent atoms
         DO 100 nn = 1,neq(n)
           na = na + 1
c
c--->      load the tuu, etc. matrices for this atom 
c
           nrec = nb0 + na
           IF (jsp.EQ.1) READ (28,rec=nrec) tuu,tdd,tud,tdu,ind
           IF (jsp.EQ.2) READ (38,rec=nrec) tuu,tdd,tud,tdu,ind
c+inv
           IF ((invsat(na).EQ.0) .OR. (invsat(na).EQ.1)) THEN
              IF (invsat(na).EQ.0) invsfct = 1
              IF (invsat(na).EQ.1) invsfct = 2
              np = ngopr(na)
c---> a and b
            DO lm = 0,lmd
               DO ke = 1,ne
                  a(ke,lm) = cmplx(0.0,0.0)
                  b(ke,lm) = cmplx(0.0,0.0)
               END DO
            END DO
            DO 50 k = 1,nv
               vmult(1) = bkpt(1) + k1(k)
               vmult(2) = bkpt(2) + k2(k)
               vmult(3) = bkpt(3) + k3(k)
               th = tpi* (vmult(1)*taual(1,na)+vmult(2)*taual(2,na)+
     +              vmult(3)*taual(3,na))
!              phase = cexp(cmplx(0.0,-th))
               phase = cmplx(cos(th),-sin(th))
c-->     apply the rotation that brings this atom into the
c-->     representative for hamiltonian (this is the definition
c-->     of ngopr(na)) and transform to cartesian coordinates
               DO 40 j = 1,3
                  vsmult(j) = vmult(1)*mrot(1,j,np) +
     +                        vmult(2)*mrot(2,j,np) +
     +                        vmult(3)*mrot(3,j,np)
   40          CONTINUE
               vmult(1) = bmat(1,1)*vsmult(1) + bmat(2,1)*vsmult(2) +
     +                    bmat(3,1)*vsmult(3)
               vmult(2) = bmat(1,2)*vsmult(1) + bmat(2,2)*vsmult(2) +
     +                    bmat(3,2)*vsmult(3)
               vmult(3) = bmat(1,3)*vsmult(1) + bmat(2,3)*vsmult(2) +
     +                    bmat(3,3)*vsmult(3)
               CALL ylm4(
     >                   lwn,vmult,
     <                   ylm)
c-->     synthesize the complex conjugates of a and b
               DO l = 0,lwn
                  ll1 = l* (l+1)
                  DO m = -l,l
                     lm = ll1 + m
                     hij = f(l,k) * ( phase * ylm(lm+1) )
                     sij = g(l,k) * ( phase * ylm(lm+1) )
                     DO ke = 1,ne
                        a(ke,lm) = a(ke,lm) + hij*z(k,ke)
                        b(ke,lm) = b(ke,lm) + sij*z(k,ke)
                     END DO
                  END DO
               END DO
   50       CONTINUE
            DO 90 l = 0,lwn
               DO 80 m = -l,l
                  lmp = l* (l+1) + m
c--->    initialize ax and bx
                  DO k = 1,ne
                     ax(k) = cmplx(0.0,0.0)
                     bx(k) = cmplx(0.0,0.0)
                  END DO
c--->    loop over l,m
                  DO 70 l1 = 0,lwn
                     DO 60 m1 = -l1,l1
                        lm = l1* (l1+1) + m1
                        in = ind(lmp,lm)
                        IF (in.NE.-9999) THEN
                           IF (in.GE.0) THEN
                              utu = conjg(tuu(in))*invsfct
                              dtu = conjg(tdu(in))*invsfct
                              utd = conjg(tud(in))*invsfct
                              dtd = conjg(tdd(in))*invsfct
                           ELSE
                              im = -in
                              utu = tuu(im)*invsfct
                              dtd = tdd(im)*invsfct
                              utd = tdu(im)*invsfct
                              dtu = tud(im)*invsfct
                           END IF
c--->    update ax, bx
                           DO k = 1,ne
                              ax(k) = ax(k) + conjg(utu*a(k,lm)+
     +                                utd*b(k,lm))
                              bx(k) = bx(k) + conjg(dtu*a(k,lm)+
     +                                dtd*b(k,lm))
                           END DO
                        END IF
   60                CONTINUE
   70             CONTINUE
c--->    update hamiltonian in upper packed storage mode
                  DO i = 1,ne
                     ii = (i-1)*i/2
                     DO j = 1,i - 1
                        ij = ii + j
                        hij = a(i,lmp)*ax(j) + b(i,lmp)*bx(j)
                        h(ij) = h(ij) + real(hij)
                     END DO
                     h(ii+i) = h(ii+i) + real(a(i,lmp)*ax(i)+
     +                        b(i,lmp)*bx(i))
                  END DO
   80          CONTINUE
   90       CONTINUE
            ENDIF
c-inv
  100    CONTINUE
  110 CONTINUE

      END SUBROUTINE h_nonmuff
      END MODULE m_hnonmuff
