      MODULE m_abccoflo
c*********************************************************************
c Calculates the (upper case) A, B and C coefficients for the local
c orbitals.
c Philipp Kurz 99/04
c*********************************************************************
      CONTAINS
      SUBROUTINE abccoflo(
     >                    llod,nlod,lmaxd,ntypd,natd,
     >                    con1,rmt,rph,cph,ylm,ntyp,na,k,nv,
     >                    l_lo1,sfp,llo,nlo,invsat,alo1,blo1,clo1,
     X                    nkvec,cwork,
     <                    enough,alo,blo,clo,kvec)
c
c*************** ABBREVIATIONS ***************************************
c kvec    : stores the number of the G-vectors, that have been used to
c           construct the local orbitals
c nkvec   : stores the number of G-vectors that have been found and
c           accepted during the construction of the local orbitals.
c enough  : enough is set to .true. when enough G-vectors have been
c           accepted.
c linindq : if the norm of that part of a local orbital (contructed 
c           with a trial G-vector) that is orthogonal to the previous
c           ones is larger than linindq, then this G-vector is 
c           accepted.
c*********************************************************************
c
      IMPLICIT NONE
C     .. 
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: llod,nlod,lmaxd,ntypd,natd
      REAL,    INTENT (IN) :: con1,cph,rmt,rph,sfp
      INTEGER, INTENT (IN) :: k,na,ntyp,nv
      LOGICAL, INTENT (IN) :: l_lo1
      LOGICAL, INTENT (OUT):: enough
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd),invsat(natd)
      INTEGER, INTENT (IN)::  kvec(2* (2*llod+1),nlod)
      REAL,    INTENT (IN) :: alo1(nlod),blo1(nlod),clo1(nlod)
      COMPLEX, INTENT (IN) :: ylm( (lmaxd+1)**2 )
      COMPLEX, INTENT (OUT):: alo(-llod:llod,2* (2*llod+1),nlod)
      COMPLEX, INTENT (OUT):: blo(-llod:llod,2* (2*llod+1),nlod)
      COMPLEX, INTENT (OUT):: clo(-llod:llod,2* (2*llod+1),nlod)
      INTEGER,INTENT (INOUT):: nkvec(nlod)
      COMPLEX,INTENT (INOUT):: cwork(-2*llod:2*llod+1,2*(2*llod+1),nlod)
C     ..
C     .. Local Scalars ..
      COMPLEX term1
      REAL linindq
      INTEGER l,lo,m,mind,ll1,lm
      LOGICAL linind
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,conjg,sqrt
C     ..
C     .. Data statements ..
      DATA linindq/1.0e-4/
C     ..
c
c---> the whole program is in hartree units, therefore 1/wronskian is
c---> (rmt**2)/2. the factor i**l, which usually appears in the a, b
c---> and c coefficients, is included in the t-matrices. thus, it does
c---> not show up in the formula above.
c
c-abccoflo1
      IF ( l_lo1) THEN
      DO lo = 1,nlo(ntyp)
         IF ( (nkvec(lo).EQ.0).AND.(llo(lo,ntyp).EQ.0) ) THEN
            enough = .false.
            nkvec(lo) = 1
            m = 0
            clo(m,nkvec(lo),lo) = con1* ((rmt**2)/2) / sfp
            alo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*alo1(lo)
            blo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*blo1(lo)
            clo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*clo1(lo)
            IF (kvec(nkvec(lo),lo).NE.k) STOP 'abccoflo:1'
            IF (invsat(na).EQ.1) THEN
               cwork(0,nkvec(lo),lo) = 1/sqrt(2.0)
               cwork(1,nkvec(lo),lo) = 1/sqrt(2.0)
            ENDIF
         ENDIF
      ENDDO
      ELSE
      enough = .true.
      term1 = con1* ((rmt**2)/2)*cmplx(rph,cph)
      DO lo = 1,nlo(ntyp)
         IF (invsat(na).EQ.0) THEN
            IF ((nkvec(lo)).LT. (2*llo(lo,ntyp)+1)) THEN
               enough = .false.
               nkvec(lo) = nkvec(lo) + 1
               l = llo(lo,ntyp)
               ll1 = l*(l+1) + 1
               DO m = -l,l
                  lm = ll1 + m
                  clo(m,nkvec(lo),lo) = term1*ylm(lm)
               END DO
               IF ( kvec(nkvec(lo),lo) == k ) THEN
                  DO m = -l,l
                     alo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*alo1(lo)
                     blo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*blo1(lo)
                     clo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*clo1(lo)
                  END DO
!                  WRITE(6,9000) nkvec(lo),k,lo,na,
!     +                          (clo(m,nkvec(lo),lo),m=-l,l)
! 9000             format(2i4,2i2,7(' (',e9.3,',',e9.3,')'))
               ELSE
                  nkvec(lo) = nkvec(lo) - 1
               ENDIF
            ENDIF
         ELSE
            IF ((invsat(na).EQ.1) .OR. (invsat(na).EQ.2)) THEN
c           only invsat=1 is needed invsat=2 for testing
               IF ((nkvec(lo)).LT. (2* (2*llo(lo,ntyp)+1))) THEN
                  enough = .false.
                  nkvec(lo) = nkvec(lo) + 1
                  l = llo(lo,ntyp)
                  ll1 = l*(l+1) + 1
                  DO m = -l,l
                     lm = ll1 + m
                     clo(m,nkvec(lo),lo) = term1*ylm(lm)
                  END DO
                   IF ( kvec(nkvec(lo),lo) == k ) THEN
                     DO m = -l,l
!                            if(l.eq.1) then
!               WRITE(*,*)'k=',k,' clotmp=',clo(m,nkvec(lo),lo)
!               WRITE(*,*)'clo1=',clo1(lo),' term1=',term1
!                            endif
                        alo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*
     +                                        alo1(lo)
                        blo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*
     +                                        blo1(lo)
                        clo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*
     +                                        clo1(lo)
!                        kvec(nkvec(lo),lo) = k
                     END DO
                  ELSE
                     nkvec(lo) = nkvec(lo) - 1
                  END IF
               END IF
            END IF
         END IF
      END DO
      IF ((k.EQ.nv) .AND. (.NOT.enough)) THEN
         WRITE (6,FMT=*)
     +     'abccoflo did not find enough linearly independent'
         WRITE (6,FMT=*)
     +     'clo coefficient-vectors. the linear independence'
         WRITE (6,FMT=*) 'quality, linindq, is set to: ',linindq,'.'
         WRITE (6,FMT=*) 'this value might be to large.'
         STOP 'abccoflo: did not find enough lin. ind. clo-vectors'
      END IF
      ENDIF  ! abccoflo1

      END SUBROUTINE abccoflo
      END MODULE m_abccoflo
