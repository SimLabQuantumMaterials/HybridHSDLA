      MODULE m_eparas
c***********************************************************************
c Calculates qlo, enerlo and sqlo, which are needed to determine the 
c new energy parameters.
c Philipp Kurz 99/04
c***********************************************************************
c also the 'normal' energy parameters are now included...
c
c if (l_mcd) then mcd contains mcd spectrum: first index = polarization
c second = core level ; third = band index                  gb.2001
c corrected to work also for multiple LO's of same l at the same atom
c                                                           gb.2005
c*************** ABBREVIATIONS *****************************************
c qlo     : charge density of one local orbital at the current k-point
c sqlo    : qlo integrated over the Brillouin zone
c enerlo  : qlo*energy integrated over the Brillouin zone
c***********************************************************************
c
      CONTAINS
      SUBROUTINE eparas(
     >                  llod,noccbd,nlod,natd,neigd,ntypd,lmaxd,lmd,
     >                  isize,ikpt,ntype,neq,ne,we,eig,ccof,nlo,llo,
     >                  invsat,skip_t,l_evp,acof,bcof,ddn,
     >                  ncore,nstd,l_mcd,m_mcd,uulon,dulon,uloulopn,
     <                  enerlo,sqlo,ener,sqal,qal,mcd)

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: llod,noccbd,nlod,natd,neigd,ntypd,lmaxd
      INTEGER, INTENT (IN) :: lmd,ne,ikpt,ntype,isize,skip_t
      INTEGER, INTENT (IN) :: nstd
      LOGICAL, INTENT (IN) :: l_mcd,l_evp
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: neq(ntypd),invsat(natd),ncore(ntypd)
      INTEGER, INTENT (IN)  :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (IN)  :: eig(neigd),we(noccbd),ddn(0:lmaxd,ntypd)
      REAL,    INTENT (IN)  :: uloulopn(nlod,nlod,ntypd)
      REAL,    INTENT (IN)  :: uulon(nlod,ntypd),dulon(nlod,ntypd)
      COMPLEX, INTENT (IN)  :: ccof(-llod:llod,noccbd,nlod,natd)
      COMPLEX, INTENT (IN)  :: acof(noccbd,0:lmd,natd)
      COMPLEX, INTENT (IN)  :: bcof(noccbd,0:lmd,natd)
      COMPLEX, INTENT (IN)  :: m_mcd(nstd,(3+1)**2,3*ntypd,2)
      REAL,    INTENT (OUT) :: enerlo(nlod,ntypd),sqlo(nlod,ntypd)
      REAL,    INTENT (OUT) :: ener(0:3,ntypd),sqal(0:3,ntypd)
      REAL,    INTENT (OUT) :: qal(0:3,ntypd,neigd)
      REAL,    INTENT (OUT) :: mcd(3*ntypd,nstd,neigd)

C     ..
C     .. Local Scalars ..
      INTEGER i,l,lo,lop,m,natom,nn,ntyp
      INTEGER nt1,nt2,lm,n,ll1,ipol,icore,index
      REAL fac
      COMPLEX suma,sumb,sumab,sumba
C     ..
C     .. Local Arrays ..
      REAL qlo(noccbd,nlod,nlod,ntypd)
      REAL qaclo(noccbd,nlod,ntypd),qbclo(noccbd,nlod,ntypd)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg
c
c---> initialize ener, sqal, enerlo and sqlo on first call
c

      IF ((ikpt.LE.isize).AND..NOT.l_evp) THEN
         IF (l_mcd) THEN
           mcd(:,:,:) = 0.0
         ENDIF
         ener(:,:) = 0.0
         sqal(:,:) = 0.0
         qal(:,:,:) = 0.0
         enerlo(:,:) = 0.0
         sqlo(:,:) = 0.0
      END IF
c
c--->    l-decomposed density for each occupied state
c
!         DO 140 i = (skip_t+1),ne    ! this I need for all states
         DO 140 i = 1,ne              ! skip in next loop
            nt1 = 1
            DO 130 n = 1,ntype
               fac = 1./neq(n)
               nt2 = nt1 + neq(n) - 1
               DO 120 l = 0,3
                  suma = cmplx(0.,0.)
                  sumb = cmplx(0.,0.)
                  ll1 = l* (l+1)
                  DO 110 m = -l,l
                     lm = ll1 + m
                     IF ( .not.l_mcd ) THEN
                       DO natom = nt1,nt2
                         suma = suma + acof(i,lm,natom)*
     +                          conjg(acof(i,lm,natom))
                         sumb = sumb + bcof(i,lm,natom)*
     +                          conjg(bcof(i,lm,natom))
                       ENDDO
                     ELSE
                       suma = cmplx(0.,0.) ; sumab = cmplx(0.,0.) 
                       sumb = cmplx(0.,0.) ; sumba = cmplx(0.,0.)
                       DO natom = nt1,nt2
                         suma = suma + acof(i,lm,natom)*
     +                          conjg(acof(i,lm,natom))
                         sumb = sumb + bcof(i,lm,natom)*
     +                          conjg(bcof(i,lm,natom))
                         sumab= sumab + acof(i,lm,natom) *
     +                            conjg(bcof(i,lm,natom))
                         sumba= sumba + bcof(i,lm,natom) *
     +                            conjg(acof(i,lm,natom))
                       ENDDO
                       DO icore = 1, ncore(n)
                         DO ipol = 1, 3
                           index = 3*(n-1) + ipol
                           mcd(index,icore,i)=mcd(index,icore,i) + fac*(
     +                         suma * conjg(m_mcd(icore,lm+1,index,1))*
     +                                      m_mcd(icore,lm+1,index,1)  +
     +                         sumb * conjg(m_mcd(icore,lm+1,index,2))*
     +                                      m_mcd(icore,lm+1,index,2)  +
     +                         sumab* conjg(m_mcd(icore,lm+1,index,2))*
     +                                      m_mcd(icore,lm+1,index,1)  +
     +                         sumba* conjg(m_mcd(icore,lm+1,index,1))*
     +                                      m_mcd(icore,lm+1,index,2)  ) 
                         ENDDO
                       ENDDO 
                     ENDIF     ! end MCD
  110             CONTINUE
                  qal(l,n,i) = (suma+sumb*ddn(l,n))/neq(n)
  120          CONTINUE
               nt1 = nt1 + neq(n)
  130       CONTINUE
  140    CONTINUE
c
c--->    perform Brillouin zone integration and summation over the
c--->    bands in order to determine the energy parameters for each
c--->    atom and angular momentum
c
         DO l = 0,3
            DO n = 1,ntype
               DO i = (skip_t+1),noccbd
                  ener(l,n) = ener(l,n) + qal(l,n,i)*we(i)*eig(i)
                  sqal(l,n) = sqal(l,n) + qal(l,n,i)*we(i)
               ENDDO
            ENDDO
         ENDDO

c---> initialize qlo

      qlo(:,:,:,:) = 0.0
      qaclo(:,:,:) = 0.0
      qbclo(:,:,:) = 0.0

c---> density for each local orbital and occupied state

      natom = 0
      DO ntyp = 1,ntype
         DO nn = 1,neq(ntyp)
            natom = natom + 1
            DO lo = 1,nlo(ntyp)
               l = llo(lo,ntyp)
               ll1 = l* (l+1)
               DO m = -l,l
                  lm = ll1 + m
                  DO i = 1,ne
                     qbclo(i,lo,ntyp) = qbclo(i,lo,ntyp) +       real(
     +                      bcof(i,lm,natom)*conjg(ccof(m,i,lo,natom)) +
     +                      ccof(m,i,lo,natom)*conjg(bcof(i,lm,natom)) )
                     qaclo(i,lo,ntyp) = qaclo(i,lo,ntyp) +       real(
     +                      acof(i,lm,natom)*conjg(ccof(m,i,lo,natom)) +
     +                      ccof(m,i,lo,natom)*conjg(acof(i,lm,natom)) )
                  ENDDO
               ENDDO
               DO lop = 1,nlo(ntyp)
                 IF (llo(lop,ntyp).EQ.l) THEN
                   DO m = -l,l
                     DO i = 1,ne
                       qlo(i,lop,lo,ntyp) = qlo(i,lop,lo,ntyp) +  real(
     +                   conjg(ccof(m,i,lop,natom))*ccof(m,i,lo,natom))
                     ENDDO
                   ENDDO
                 ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

c---> perform brillouin zone integration and sum over bands

      DO ntyp = 1,ntype
         DO lo = 1,nlo(ntyp)
            l = llo(lo,ntyp)
            DO i = 1,ne
               qal(l,ntyp,i)= qal(l,ntyp,i)  + ( 1.0/neq(ntyp) )* (
     +                        qaclo(i,lo,ntyp)*uulon(lo,ntyp)     +
     +                        qbclo(i,lo,ntyp)*dulon(lo,ntyp)     )
            END DO
            DO lop = 1,nlo(ntyp)
               IF (llo(lop,ntyp).EQ.l) THEN
               DO i = 1,ne
                 enerlo(lo,ntyp) = enerlo(lo,ntyp) +
     +                             qlo(i,lop,lo,ntyp)*we(i)*eig(i)
                 sqlo(lo,ntyp) = sqlo(lo,ntyp) + 
     +                             qlo(i,lop,lo,ntyp)*we(i)
                 qal(l,ntyp,i)= qal(l,ntyp,i)  + ( 1.0/neq(ntyp) ) *
     +                      qlo(i,lop,lo,ntyp)*uloulopn(lop,lo,ntyp)
               ENDDO
               ENDIF
            ENDDO
         END DO
      END DO

      END SUBROUTINE eparas
      END MODULE m_eparas
