      MODULE m_hssrwu
c
c*********************************************************************
c     updates the hamiltonian and overlap matrices with the
c     contributions from the spheres, both spherical and non-
c     spherical, for step forward approach
c                r. wu  1992
c*********************************************************************
      CONTAINS
      SUBROUTINE hssr_wu(
     >                   nvd,lmaxd,ntypd,neigd,natd,nop,jspd,nwdd,
     >                   lmd,lmplmd,ntype,neq,lmax,nbasfcn,invsat,
     >                   nlod,llod,jsp,nw,el,lnonsph,ne,ddn,nv,
     >                   tuu,tud,tdu,tdd,tuulo,tdulo,tuloulo,mlotot,
     >                   mlolotot,ind,invtab,acof,bcof,ccof,
     X                   cp_sph,cp_nsph,h,s)
c
      USE m_constants, ONLY : pimach
      USE m_sphbes
      USE m_dsphbs
      USE m_cputime
      USE m_ylm

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nvd,lmaxd,ntypd,neigd,natd,nop,jspd,nwdd
      INTEGER, INTENT (IN) :: jsp,nw,ne,nv,ntype,lmd,lmplmd,nbasfcn
      INTEGER, INTENT (IN) :: nlod,llod,mlotot,mlolotot
      REAL,    INTENT (INOUT) :: cp_sph,cp_nsph
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: lnonsph(ntypd),invsat(natd)
      INTEGER, INTENT (IN) :: lmax(ntypd),neq(ntypd)
      REAL,    INTENT (IN) :: el(0:lmaxd,ntypd,jspd,nwdd)
      REAL,    INTENT (IN) :: ddn(0:lmaxd,ntypd)
      COMPLEX, INTENT (IN) :: acof(neigd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: bcof(neigd,0:lmd,natd)
      COMPLEX, INTENT (IN) :: ccof(-llod:llod,neigd,nlod,natd)
      COMPLEX, INTENT (IN) :: tdulo(0:lmd,-llod:llod,mlotot)
      COMPLEX, INTENT (IN) :: tuulo(0:lmd,-llod:llod,mlotot)
      COMPLEX, INTENT (IN) :: tuloulo(-llod:llod,-llod:llod,mlolotot)
      COMPLEX, INTENT (IN) :: tdd(0:lmplmd,ntypd),tdu(0:lmplmd,ntypd)
      COMPLEX, INTENT (IN) :: tud(0:lmplmd,ntypd),tuu(0:lmplmd,ntypd)
      INTEGER, INTENT (IN) :: ind(0:lmd,0:lmd,ntypd),invtab(nop)
#ifdef CPP_INVERSION
      REAL,    INTENT (INOUT) :: h(neigd,neigd),s(neigd,neigd)
#else
      COMPLEX, INTENT (INOUT) :: h(neigd,neigd),s(neigd,neigd)
#endif  

C     ..
C     .. Local Scalars ..
      COMPLEX dtd,dtu,hij,sij,utd,utu
      REAL invsfct,time1,time2
      INTEGER i,im,in,j,k,ke,l,l1,ll1,lm,lmp,lwn,m,m1,n,na,nn,lmplm
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: a(:,:),b(:,:),ax(:),bx(:)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,conjg,exp,real,sqrt
C     ..

      ALLOCATE ( a(neigd,0:lmd),ax(neigd) )
      ALLOCATE ( b(neigd,0:lmd),bx(neigd) )
      na = 0
      DO n = 1,ntype        ! loop over atom-types
        lwn = lmax(n)
        DO nn = 1,neq(n)    ! loop over atoms
          na = na + 1
c+inv
          IF ((invsat(na).EQ.0) .OR. (invsat(na).EQ.1)) THEN
            CALL cpu_time(time1)
            IF (invsat(na).EQ.0) invsfct = 1.0
            IF (invsat(na).EQ.1) invsfct = sqrt(2.0)
            DO lm = 0, lmd
              DO ke = 1, ne
                a(ke,lm) = invsfct*acof(ke,lm,na)
                b(ke,lm) = invsfct*bcof(ke,lm,na)
              ENDDO
            ENDDO

            DO l = 0,lwn                    ! l loop
               DO m = -l,l                  ! m loop
                  lmp = l* (l+1) + m
                  DO i = 1,ne               ! matrix update
                     DO j = 1,i - 1
                        sij = a(i,lmp)*conjg(a(j,lmp)) +
     +                        b(i,lmp)*conjg(b(j,lmp))*ddn(l,n)
                        hij = sij * el(l,n,jsp,nw) +
     +                      0.5 * ( a(i,lmp)*conjg(b(j,lmp)) +
     +                              b(i,lmp)*conjg(a(j,lmp)) )
#ifdef CPP_INVERSION
                        s(i,j) = s(i,j) + real(sij)
                        h(i,j) = h(i,j) + real(hij)
#else
                        s(i,j) = s(i,j) + sij
                        h(i,j) = h(i,j) + hij
#endif
                     ENDDO
                  ENDDO
                  DO i = 1,ne
                     sij = a(i,lmp)*conjg(a(i,lmp)) +
     +                     b(i,lmp)*conjg(b(i,lmp))*ddn(l,n)
                     hij = sij * el(l,n,jsp,nw) +
     +                  0.5 * ( a(i,lmp)*conjg(b(i,lmp)) +
     +                          b(i,lmp)*conjg(a(i,lmp)) )
#ifdef CPP_INVERSION
                     s(i,i) = s(i,i) + real(sij)
                     h(i,i) = h(i,i) + real(hij)
#else
                     s(i,i) = s(i,i) + sij
                     h(i,i) = h(i,i) + hij
#endif
                  ENDDO
               ENDDO        ! m
            ENDDO           ! l
            CALL cpu_time(time2)
            cp_sph = cp_sph + time2 - time1

            IF (lnonsph(n) >= 0 ) THEN
              DO l = 0,lnonsph(n)
                DO m = -l,l

                  lmp = l* (l+1) + m
                  ax(:) = cmplx(0.0,0.0)
                  bx(:) = cmplx(0.0,0.0)

                  DO l1 = 0,lnonsph(n)         ! l', m' loop
                     DO m1 = -l1,l1
                        lm = l1* (l1+1) + m1
                        in = ind(lmp,lm,n)
                        IF (in.NE.-9999) THEN

                           IF (in.GE.0) THEN
                              utu = conjg(tuu(in,n))
                              dtu = conjg(tdu(in,n))
                              utd = conjg(tud(in,n))
                              dtd = conjg(tdd(in,n))
                           ELSE
                              im = -in
                              utu = tuu(im,n)
                              dtd = tdd(im,n)
                              utd = tdu(im,n)
                              dtu = tud(im,n)
                           END IF
c--->    update ax, bx
                           DO k = 1,ne
                              ax(k) = ax(k) + utu*conjg(a(k,lm))+
     +                                utd*conjg(b(k,lm))
                              bx(k) = bx(k) + dtu*conjg(a(k,lm))+
     +                                dtd*conjg(b(k,lm))
                           ENDDO

                        ENDIF ! in =/= -9999
                     ENDDO    ! m1
                  ENDDO       ! l1
c
c
c--->    update hamiltonian
                 DO i = 1,ne
                    DO j = 1,i - 1
                        hij = a(i,lmp)*ax(j) + b(i,lmp)*bx(j)
#ifdef CPP_INVERSION
                        h(i,j) = h(i,j) + real(hij)
#else
                        h(i,j) = h(i,j) + hij
#endif
                     ENDDO
                  ENDDO
                  DO i = 1,ne
#ifdef CPP_INVERSION
                     h(i,i) = h(i,i) + real(a(i,lmp)*ax(i)+
     +                                      b(i,lmp)*bx(i))
#else
                     h(i,i) = h(i,i) + a(i,lmp)*ax(i)+
     +                                 b(i,lmp)*bx(i)
#endif
                  ENDDO

                ENDDO ! m
              ENDDO   ! l
            ENDIF     ! lnonsph >=0
            CALL cpu_time(time1)
            cp_nsph = cp_nsph + time1 - time2

c-inv
          ENDIF ! invsatom = 0 or 1
        ENDDO   ! loop over atoms
      ENDDO     ! loop over atom-types

      DEALLOCATE ( a, b, ax, bx )
      END SUBROUTINE hssr_wu
      END MODULE m_hssrwu
