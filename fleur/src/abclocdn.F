      MODULE m_abclocdn
c*********************************************************************
c Calculates the (upper case) A, B and C coefficients for the local
c orbitals. The difference to abccoflo is, that a summation over the
c Gs ist performed. The A, B and C coeff. are set up for each eigen-
c state.
c Philipp Kurz 99/04
c*********************************************************************
c*************** ABBREVIATIONS ***************************************
c nkvec   : stores the number of G-vectors that have been found and
c           accepted during the construction of the local orbitals.
! kvec    : k-vector used in hssphn to attach the local orbital 'lo'
!           of atom 'na' to it.
c*********************************************************************
      CONTAINS
      SUBROUTINE abclocdn(
     >                    nobd,natd,nlod,llod,lmaxd,neigd,ntypd,
     >                    lmd,nbasfcn,nlo,llo,invsat,invsatnr,
     >                    l_noco,ccchi,kspin,l_ss,iintsp,
     >                    con1,rmt,phase,ylm,ntyp,na,k,s,nv,ne,z,
     >                    nbasf0,alo1,blo1,clo1,kvec,
     <                    nkvec,enough,acof,bcof,ccof)
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nobd,natd,nlod,llod,lmaxd,neigd,ntypd
      INTEGER, INTENT (IN) :: lmd,nbasfcn,kspin,iintsp
      INTEGER, INTENT (IN) :: k,na,ne,ntyp,nv
      REAL,    INTENT (IN) :: con1,rmt,s
      COMPLEX, INTENT (IN) :: phase
      LOGICAL, INTENT (IN) :: l_noco,l_ss
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nbasf0(nlod,natd),invsat(natd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd),invsatnr(natd)
      REAL,    INTENT (IN) :: alo1(nlod,ntypd),blo1(nlod,ntypd)
      REAL,    INTENT (IN) :: clo1(nlod,ntypd)
      COMPLEX, INTENT (IN) :: ylm( (lmaxd+1)**2 )
      COMPLEX, INTENT (IN) :: ccchi(2)
      INTEGER, INTENT (IN) :: kvec(2*(2*llod+1),nlod)
      LOGICAL, INTENT (OUT) :: enough(natd)
      COMPLEX, INTENT (INOUT) :: acof(nobd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: bcof(nobd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: ccof(-llod:llod,nobd,nlod,natd)
      INTEGER, INTENT (INOUT) :: nkvec(nlod,natd)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#else
      REAL,    INTENT (IN) :: z(nbasfcn,neigd)
#endif
C     ..
C     .. Local Scalars ..
      COMPLEX ctmp,term1
      REAL linindq,eps
      INTEGER i,l,ll1,lm,lo,m,mind,nbasf,na2,lmp
      LOGICAL linind
C     ..
C     .. Local Arrays ..
      COMPLEX clotmp(-llod:llod)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg
C     ..
C     .. Data statements ..

      DATA linindq/1.0e-4/,eps/1.0e-30/
C     ..
      enough(na) = .true.
      term1 = con1 * ((rmt**2)/2) * phase
c---> the whole program is in hartree units, therefore 1/wronskian is
c---> (rmt**2)/2. the factor i**l, which usually appears in the a, b
c---> and c coefficients, is included in the t-matrices. thus, it does
c---> not show up in the formula above.
      DO lo = 1,nlo(ntyp)
         l = llo(lo,ntyp)
         IF (.NOT.((s.LE.eps).AND.(l.GE.1))) THEN
           IF (invsat(na).EQ.0) THEN

             IF ((nkvec(lo,na)).LT. (2*llo(lo,ntyp)+1)) THEN
               enough(na) = .false.
               nkvec(lo,na) = nkvec(lo,na) + 1
               nbasf = nbasf0(lo,na) + nkvec(lo,na)
               l = llo(lo,ntyp)
               ll1 = l* (l+1)
               DO m = -l,l
                  clotmp(m) = term1*conjg(ylm(ll1+m+1))
               END DO
               IF ( kvec(nkvec(lo,na),lo) == k ) THEN
                  DO i = 1,ne
                     DO m = -l,l
                        lm = ll1 + m
c+gu_con
                        IF (l_noco) THEN
                          IF (l_ss) THEN
                            ctmp = clotmp(m)*
     +                             ccchi(iintsp)*z(kspin+nbasf,i)
                          ELSE
                            ctmp = clotmp(m)*( ccchi(1)*z(nbasf,i)
     +                               +   ccchi(2)*z(kspin+nbasf,i) )
                          ENDIF
                        ELSE
                           ctmp = z(nbasf,i)*clotmp(m)
                        ENDIF
                        acof(i,lm,na) = acof(i,lm,na) +
     +                                   ctmp*alo1(lo,ntyp)
                        bcof(i,lm,na) = bcof(i,lm,na) +
     +                                   ctmp*blo1(lo,ntyp)
                        ccof(m,i,lo,na) = ccof(m,i,lo,na) +
     +                                   ctmp*clo1(lo,ntyp)
                     END DO
                  END DO
!                  write(6,9000) nbasf,k,lo,na,
!     +                          (clo1(lo,ntyp)*clotmp(m),m=-l,l)
! 9000             format(2i4,2i2,7(' (',e9.3,',',e9.3,')'))
               ELSE
                  nkvec(lo,na) = nkvec(lo,na) - 1
               ENDIF ! linind
             ENDIF   ! nkvec < 2*llo

           ELSEIF (invsat(na).EQ.1) THEN
             IF ((nkvec(lo,na)).LT. (2* (2*llo(lo,ntyp)+1))) THEN
                enough(na) = .false.
                nkvec(lo,na) = nkvec(lo,na) + 1
                nbasf = nbasf0(lo,na) + nkvec(lo,na)
                l = llo(lo,ntyp)
                ll1 = l* (l+1)
                DO m = -l,l
                   clotmp(m) = term1*conjg(ylm(ll1+m+1))
                END DO
                IF ( kvec(nkvec(lo,na),lo) == k ) THEN
c                  write(*,*)'k vector nr ',k,' has been accepted'
                   DO i = 1,ne
                      DO m = -l,l
                         lm = ll1 + m
!                        if(i.eq.1 .and. l.eq.1) then
!              write(*,*)'k=',k,' z=',z(nbasf,i),' clotmp=',clotmp(m)
!              write(*,*)'clo1=',clo1(lo,ntyp),' term1=',term1
!                         endif
c+gu_con
                         IF (l_noco) THEN
                            IF (l_ss) THEN
                               ctmp = clotmp(m)*
     +                              ccchi(iintsp)*z(kspin+nbasf,i)
                            ELSE
                               ctmp = clotmp(m)*( ccchi(1)*z(nbasf,i)
     +                                  + ccchi(2)*z(kspin+nbasf,i) )
                            ENDIF
                         ELSE
                            ctmp = z(nbasf,i)*clotmp(m)
                         ENDIF
                         acof(i,lm,na) = acof(i,lm,na) +
     +                                               ctmp*alo1(lo,ntyp)
                         bcof(i,lm,na) = bcof(i,lm,na) +
     +                                               ctmp*blo1(lo,ntyp)
                         ccof(m,i,lo,na) = ccof(m,i,lo,na) +
     +                                               ctmp*clo1(lo,ntyp)
#if ( defined(CPP_SOC) && defined(CPP_INVERSION) )
                         ctmp = z(nbasf,i)*conjg(clotmp(m))*(-1)**(l+m)
                         na2 = invsatnr(na)
                         lmp = ll1 - m
                         acof(i,lmp,na2) = acof(i,lmp,na2) +
     +                                               ctmp*alo1(lo,ntyp)
                         bcof(i,lmp,na2) = bcof(i,lmp,na2) +
     +                                               ctmp*blo1(lo,ntyp)
                         ccof(-m,i,lo,na2) = ccof(-m,i,lo,na2) +
     +                                               ctmp*clo1(lo,ntyp)
#endif
                      ENDDO  ! m
                   ENDDO     ! i = 1,ne
                ELSE       
                    nkvec(lo,na) = nkvec(lo,na) - 1
                ENDIF       ! linind
              ENDIF         ! nkvec < 2*llo
           ELSE
              STOP 'invsat =/= 0 or 1'
           ENDIF 
         ELSE
            enough(na) = .false.
         ENDIF  ! s > eps  & l >= 1
      END DO
      IF ((k.EQ.nv) .AND. (.NOT.enough(na))) THEN
         WRITE (6,FMT=*)
     +     'abclocdn did not find enough linearly independent'
         WRITE (6,FMT=*)
     +     'ccof coefficient-vectors. the linear independence'
         WRITE (6,FMT=*) 'quality, linindq, is set to: ',linindq,'.'
         WRITE (6,FMT=*) 'this value might be to large.'
         STOP 'abclocdn: did not find enough lin. ind. ccof-vectors'
      END IF

      END SUBROUTINE abclocdn
      END MODULE m_abclocdn
