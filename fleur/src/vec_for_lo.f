      MODULE m_vecforlo
      CONTAINS
      SUBROUTINE vec_for_lo(
     >                      llod,nlod,ntypd,natd,nvd,jspd,nintsp,nop,na,
     >                      n,np,lmaxd,con1,tpi,eps,l_ss,nlo,llo,invsat,
     >                      lnonsph,nv,mrot,k1,k2,k3,qss,taual,rmt,bmat,
     >                      rk,gk,vk,
     <                      nkvec,kvec)

      USE m_orthoglo
      USE m_ylm
      USE m_matmul   , ONLY : matmul3
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: llod,nlod,ntypd,natd,nvd,jspd
      INTEGER, INTENT (IN) :: nintsp,nop,na,n,np,lmaxd
      REAL,    INTENT (IN) :: con1,tpi,eps
      LOGICAL, INTENT (IN) :: l_ss
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd),invsat(natd)
      INTEGER, INTENT (IN) :: lnonsph(ntypd),nv(jspd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      REAL,    INTENT (IN) :: qss(3),taual(3,natd),rmt(ntypd)
      REAL,    INTENT (IN) :: gk(nvd,3,nintsp),vk(nvd,3,nintsp)
      REAL,    INTENT (IN) :: bmat(3,3),rk(nvd,jspd)
      INTEGER, INTENT (OUT):: kvec(2*(2*llod+1),nlod),nkvec(nlod,nintsp)
C     ..
C     .. Local Scalars ..
      COMPLEX term1 
      REAL t1nn,t2nn,t3nn,th
      INTEGER l,lo,m,mind,ll1,lm,iintsp,k,nkmin,ntyp
      LOGICAL linind,enough,l_lo1
C     ..
C     .. Local Arrays ..
      REAL qssbti(3),bmrot(3,3),v(3),gkrot(nvd,3,nintsp)
      REAL rph(nvd,nintsp),cph(nvd,nintsp),vmult(3)
      COMPLEX ylm( (lmaxd+1)**2 )
      COMPLEX cwork(-2*llod:2*llod+1,2*(2*llod+1),nlod,nintsp)
C     ..
C     .. Data statements ..
      REAL, PARAMETER :: linindq = 1.0e-7

      ntyp = n
      DO iintsp = 1,nintsp
         IF (iintsp.EQ.1) THEN
            qssbti(1) = - qss(1)/2
            qssbti(2) = - qss(2)/2
            qssbti(3) = - qss(3)/2
         ELSE
            qssbti(1) = + qss(1)/2
            qssbti(2) = + qss(2)/2
            qssbti(3) = + qss(3)/2
         ENDIF

c--->    set up phase factors
         t1nn =  tpi*taual(1,na)
         t2nn =  tpi*taual(2,na)
         t3nn =  tpi*taual(3,na)
         DO k = 1,nv(iintsp)
            th = ( k1(k,iintsp) + qssbti(1) )*t1nn +
     +           ( k2(k,iintsp) + qssbti(2) )*t2nn +
     +           ( k3(k,iintsp) + qssbti(3) )*t3nn
            rph(k,iintsp) = cos(th)
            cph(k,iintsp) = -sin(th)
         END DO

         IF (np.EQ.1) THEN
            DO k = 1,nv(iintsp)
               gkrot(k,1,iintsp) = gk(k,1,iintsp)
               gkrot(k,2,iintsp) = gk(k,2,iintsp)
               gkrot(k,3,iintsp) = gk(k,3,iintsp)
            END DO
         ELSE
            CALL matmul3(mrot(1,1,np),bmat,bmrot)
            DO k = 1,nv(iintsp)
c-->           apply the rotation that brings this atom into the
c-->           representative (this is the definition of ngopr(na))
c-->           and transform to cartesian coordinates
               v(1) = vk(k,1,iintsp)
               v(2) = vk(k,2,iintsp)
               v(3) = vk(k,3,iintsp)
               gkrot(k,1,iintsp) = bmrot(1,1)*v(1) + 
     +                             bmrot(2,1)*v(2) +
     +                             bmrot(3,1)*v(3)
               gkrot(k,2,iintsp) = bmrot(1,2)*v(1) +
     +                             bmrot(2,2)*v(2) +
     +                             bmrot(3,2)*v(3)
               gkrot(k,3,iintsp) = bmrot(1,3)*v(1) + 
     +                             bmrot(2,3)*v(2) +
     +                             bmrot(3,3)*v(3)
            END DO
         END IF
c--->   end loop over interstitial spin
      ENDDO

      nkvec(:,:) = 0
      cwork(:,:,:,:) = cmplx(0.0,0.0)
      enough=.false.
      DO k = 1,min(nv(1),nv(nintsp))

        IF (.NOT.enough) THEN
          DO iintsp = 1,nintsp

c-->        generate spherical harmonics
            vmult(1) =  gkrot(k,1,iintsp)
            vmult(2) =  gkrot(k,2,iintsp)
            vmult(3) =  gkrot(k,3,iintsp)
            CALL ylm4(
     >                lnonsph(ntyp),vmult,
     <                ylm)
            l_lo1=.false.
            IF ((rk(k,iintsp).LT.eps).AND.(.not.l_ss)) THEN
                l_lo1=.true.
            ELSE
                l_lo1=.false.
            ENDIF
! --> here comes a part of abccoflo() 
            IF ( l_lo1) THEN
              DO lo = 1,nlo(ntyp)
               IF ((nkvec(lo,iintsp).EQ.0).AND.(llo(lo,ntyp).EQ.0)) THEN
                 enough = .false.
                 nkvec(lo,iintsp) = 1
                 kvec(nkvec(lo,iintsp),lo) = k
                 term1 = con1* ((rmt(ntyp)**2)/2)
                 cwork(0,1,lo,iintsp) = term1 / sqrt(2*tpi)
                 IF((invsat(na).EQ.1).OR.(invsat(na).EQ.2)) THEN
                    cwork(1,1,lo,iintsp) = conjg(term1) / sqrt(2*tpi)
                 ENDIF
               ENDIF
              ENDDO
            ELSE
              enough = .true.
              term1 = con1* ((rmt(ntyp)**2)/2)*
     *                cmplx(rph(k,iintsp),cph(k,iintsp))
              DO lo = 1,nlo(ntyp)
                IF (invsat(na).EQ.0) THEN
                  IF ((nkvec(lo,iintsp)).LT. (2*llo(lo,ntyp)+1)) THEN
                      enough = .false.
                      nkvec(lo,iintsp) = nkvec(lo,iintsp) + 1
                      l = llo(lo,ntyp)
                      ll1 = l*(l+1) + 1
                      DO m = -l,l
                        lm = ll1 + m
                        cwork(m,nkvec(lo,iintsp),lo,iintsp) = 
     +                                                   term1*ylm(lm)
                      END DO
                      CALL orthoglo(
     >                  llod,nlod,nkvec(lo,iintsp),lo,l,linindq,.false.,
     <                                 cwork(-2*llod,1,1,iintsp),linind)
                      IF (linind) THEN
                        kvec(nkvec(lo,iintsp),lo) = k
!                        write(*,*) nkvec(lo,iintsp),k,' <- '
                      ELSE
                        nkvec(lo,iintsp) = nkvec(lo,iintsp) - 1
                      ENDIF
                  ENDIF
                ELSE
                  IF ((invsat(na).EQ.1) .OR. (invsat(na).EQ.2)) THEN
                    IF (nkvec(lo,iintsp).LT.2*(2*llo(lo,ntyp)+1)) THEN
                        enough = .false.
                        nkvec(lo,iintsp) = nkvec(lo,iintsp) + 1
                        l = llo(lo,ntyp)
                        ll1 = l*(l+1) + 1
                        DO m = -l,l
                           lm = ll1 + m
                           mind = -l + m
                           cwork(mind,nkvec(lo,iintsp),lo,iintsp) = 
     +                                                   term1*ylm(lm)
                           mind = l + 1 + m
                           cwork(mind,nkvec(lo,iintsp),lo,iintsp) = 
     +                              ((-1)** (l+m))*conjg(term1*ylm(lm))
                        END DO
                        CALL orthoglo(
     >                  llod,nlod,nkvec(lo,iintsp),lo,l,linindq,.true.,
     <                                cwork(-2*llod,1,1,iintsp),linind)
                        IF (linind) THEN
                           kvec(nkvec(lo,iintsp),lo) = k
                        ELSE
                           nkvec(lo,iintsp) = nkvec(lo,iintsp) - 1
                        END IF
                      END IF
                    END IF
                  END IF
                END DO
                IF ((k.EQ.nv(iintsp)) .AND. (.NOT.enough)) THEN
                   WRITE (6,FMT=*)
     +               'abccoflo did not find enough linearly independent'
                   WRITE (6,FMT=*)
     +               'clo coefficient-vectors. the linear independence'
                   WRITE (6,FMT=*) 'quality, linindq, is set: ',linindq
                   WRITE (6,FMT=*) 'this value might be to large.'
                   write(*,*) na,k,nv 
                   STOP 'vec_for_lo: not enough lin. indep. clo-vectors'
                END IF
              ENDIF
! -- >        end of abccoflo-part           
            ENDDO
         ENDIF

! -->    check whether we have already enough k-vecs
         enough=.true.
         DO lo = 1,nlo(ntyp)
           IF (nkvec(lo,1).EQ.nkvec(lo,nintsp)) THEN   ! k-vec accepted by both spin channels
             IF (invsat(na).EQ.0) THEN
                IF ( nkvec(lo,1).LT.(2*llo(lo,ntyp)+1) ) THEN 
                  enough=.false.
                ENDIF
             ELSE
                IF ( nkvec(lo,1).LT.(2*(2*llo(lo,ntyp)+1))  ) THEN
                  enough=.false.
                ENDIF
             ENDIF
           ELSE
             nkmin = min(nkvec(lo,1),nkvec(lo,nintsp)) ! try another k-vec
             nkvec(lo,1) = nkmin ; nkvec(lo,nintsp) = nkmin
             enough=.false.
           ENDIF
         ENDDO
         IF ( enough ) THEN
!           DO  lo = 1,nlo(ntyp)
!             DO l = 1, nkvec(lo,1)
!              write(*,*) lo,l,kvec(l,lo)
!             ENDDO
!           ENDDO
           RETURN
         ENDIF
      ENDDO

      END SUBROUTINE vec_for_lo
      END MODULE m_vecforlo
