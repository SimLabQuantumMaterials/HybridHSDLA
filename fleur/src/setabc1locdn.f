      MODULE m_setabc1locdn
c***********************************************************************
c calculates the (lower case) a, b and c coefficients for the local
c orbitals. The radial function of the local orbital is a linear 
c combination of the apw radial function and its derivative and the
c extra radial funtion (a*u + b*udot + c*ulo). This function is zero
c and has zero derivative at the muffin tin boundary.
c In addition the the total number of basisfuntions (apw + lo) nbasf and
c the number of the first basisfunction of each local orbital nbasf0 is
c determined.
c Philipp Kurz 99/04
c***********************************************************************
      CONTAINS
      SUBROUTINE setabc1locdn(
     >                        ntypd,natd,nlod,llod,nobd,lmaxd,nv,
     >                        nmat,ne,ntype,neq,l_noco,iintsp,
     >                        nlo,llo,l_dulo,invsat,invsatnr,ddn,
     >                        us,dus,uds,duds,dulos,ulos,dulon,uulon,
     >                        nlotot,kveclo,
     <                        enough,nkvec,kvec,nbasf0,ccof,
     <                        alo1,blo1,clo1)
c
c*************** ABBREVIATIONS *****************************************
c nbasf   : total number of basisfunctions (apw + lo)
c nbasf0  : number of the first basisfunction of each local orbital
c nkvec   : stores the number of G-vectors that have been found and
c           accepted during the construction of the local orbitals.
c***********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,natd,nlod,llod,nobd,lmaxd
      INTEGER, INTENT (IN) :: ne,nmat,ntype,nv,iintsp,nlotot
      LOGICAL, INTENT (IN) :: l_noco
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: neq(ntypd),invsat(natd),invsatnr(natd)
      INTEGER, INTENT (IN)  :: nlo(ntypd),llo(nlod,ntypd)
      INTEGER, INTENT (IN)  :: kveclo(nlotot)
      REAL,    INTENT (IN)  ::  us(0:lmaxd,ntypd), dus(0:lmaxd,ntypd)
      REAL,    INTENT (IN)  :: uds(0:lmaxd,ntypd),duds(0:lmaxd,ntypd)
      REAL,    INTENT (IN)  ::  ulos(nlod,ntypd),dulos(nlod,ntypd)
      REAL,    INTENT (IN)  :: uulon(nlod,ntypd),dulon(nlod,ntypd)
      REAL,    INTENT (IN)  :: ddn(0:lmaxd,ntypd)
      LOGICAL, INTENT (IN)  :: l_dulo(nlod,ntypd)
      INTEGER, INTENT (OUT) :: nbasf0(nlod,natd),nkvec(nlod,natd)
      INTEGER, INTENT (OUT) :: kvec(2*(2*llod+1),nlod,natd)
      REAL,    INTENT (OUT) :: alo1(nlod,ntypd),blo1(nlod,ntypd)
      REAL,    INTENT (OUT) :: clo1(nlod,ntypd)
      COMPLEX, INTENT (INOUT) :: ccof(-llod:llod,nobd,nlod,natd)
      LOGICAL, INTENT (OUT) :: enough(natd)
C     ..
C     .. Local Scalars ..
      REAL ka,kb,ws
      INTEGER i,l,lo,m,natom,nbasf,nn,ntyp,lm
      LOGICAL apw_at
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,sqrt
C     ..
      DO natom = 1,natd
         enough(natom) = .true.
      END DO
      DO ntyp = 1,ntype
C     ..
c look, whether 'ntyp' is a APW atom; then set apw_at=.true.
c
         apw_at = .false.
         DO lo = 1,nlo(ntyp)
            IF (l_dulo(lo,ntyp)) apw_at = .true.
         ENDDO

         DO lo = 1,nlo(ntyp)
           l = llo(lo,ntyp)
           IF (apw_at) THEN
             IF (l_dulo(lo,ntyp)) THEN
c udot lo
               ka = sqrt( 1+(us(l,ntyp)/uds(l,ntyp))**2 * ddn(l,ntyp))
               alo1(lo,ntyp)=1.00 / ka
               blo1(lo,ntyp)=-us(l,ntyp)/ (uds(l,ntyp) * ka)
               clo1(lo,ntyp)=0.00
             ELSE
c u2 lo
               alo1(lo,ntyp)=1.00
               blo1(lo,ntyp)=0.00
               clo1(lo,ntyp)=-us(l,ntyp)/ulos(lo,ntyp)
             ENDIF
           ELSE
             ws = uds(l,ntyp)*dus(l,ntyp) - us(l,ntyp)*duds(l,ntyp)
             ka = 1.0/ws* (duds(l,ntyp)*ulos(lo,ntyp)-
     +            uds(l,ntyp)*dulos(lo,ntyp))
             kb = 1.0/ws* (us(l,ntyp)*dulos(lo,ntyp)-
     +            dus(l,ntyp)*ulos(lo,ntyp))
             clo1(lo,ntyp) = 1.0/sqrt(ka**2+ (kb**2)*ddn(l,ntyp)+1.0+
     +                   2.0*ka*uulon(lo,ntyp)+2.0*kb*dulon(lo,ntyp))
             alo1(lo,ntyp) = ka*clo1(lo,ntyp)
             blo1(lo,ntyp) = kb*clo1(lo,ntyp)
           ENDIF
         END DO
      END DO
c---> set up enough, nbasf0 and initialize nkvec
      natom = 0
      nbasf = nv
      DO ntyp = 1,ntype
         DO nn = 1,neq(ntyp)
            natom = natom + 1
            DO lo = 1,nlo(ntyp)
               enough(natom) = .false.
               nkvec(lo,natom) = 0
               l = llo(lo,ntyp)
               IF (invsat(natom).EQ.0) THEN
                  nbasf0(lo,natom) = nbasf
                  nbasf = nbasf + 2*l + 1
               END IF
               IF (invsat(natom).EQ.1) THEN
                  nbasf0(lo,natom) = nbasf
                  nbasf0(lo,invsatnr(natom)) = nbasf
                  nbasf = nbasf + 2* (2*l+1)
               END IF
c--->          initialize ccof
               IF (iintsp.NE.2) THEN
               DO i = 1,nobd ! ne
                  DO m = -llod,llod ! -l,l
                     ccof(m,i,lo,natom) = cmplx(0.0,0.0)
                  END DO
               END DO
               ENDIF
            END DO
         END DO
      END DO
c      write (*,*) 'in setabc1locdn: nmat = ',nmat,' nbasf = ',nbasf
c      write (*,*) 'array nbasf0 :'
c      do natom = 1,natd
c         write (*,fmt='(15i4)') (nbasf0(lo,natom),lo=1,nlod)
c      enddo
c      write (*,*)
      IF ( .NOT. l_noco ) THEN
        IF ((nmat).NE.nbasf) THEN
          write (*,*) 'in setabc1locdn: nmat = ',nmat,' nbasf = ',nbasf
          STOP 'setabc1locdn: number of bas.-fcn.'
        ENDIF
      ENDIF
!
!--> sort the k-vectors used for the LO's according to atom & lo:
!
      natom = 0
      lm = 0
      DO ntyp = 1, ntype
        DO nn = 1, neq(ntyp)
          natom = natom + 1
          IF ((invsat(natom).EQ.0) .OR. (invsat(natom).EQ.1)) THEN
            DO lo = 1,nlo(ntyp)
              m = ( invsat(natom) +1 ) * ( 2 * llo(lo,ntyp) + 1 )
              DO l = 1, m
                lm = lm + 1
                kvec(l,lo,natom) =  kveclo(lm)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE setabc1locdn
      END MODULE m_setabc1locdn
