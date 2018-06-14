      MODULE m_qmtsl
      CONTAINS
c***********************************************************************
c Calculates the mt-spheres contribution to the layer charge for states 
c  {En} at the current k-point. 
c                                      Yury Koroteev 2003
c                     from eparas.F  by  Philipp Kurz 99/04
c
c***********************************************************************
c
      SUBROUTINE q_mt_sl(
     >                  llod,nobd,nlod,natd,neigd,ntypd,lmaxd,lmd,nsld,
     >                  ikpt,ntype,neq,ne,ccof,nlo,llo,
     >                  invsat,skip_t,noccbd,acof,bcof,ddn,
     >                  uulon,dulon,uloulopn,nmtsl,nsl,lmax,
     <                  qmtslk)

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: llod,nobd,nlod,natd,neigd,ntypd,lmaxd,lmd
      INTEGER, INTENT (IN) :: ne,ikpt,ntype,skip_t,noccbd
      INTEGER, INTENT (IN) :: nsl,nsld
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: neq(ntypd),invsat(natd)
      INTEGER, INTENT (IN)  :: nlo(ntypd),llo(nlod,ntypd),lmax(ntypd)
      REAL,    INTENT (IN)  :: ddn(0:lmaxd,ntypd)
      REAL,    INTENT (IN)  :: uloulopn(nlod,nlod,ntypd)
      REAL,    INTENT (IN)  :: uulon(nlod,ntypd),dulon(nlod,ntypd)
      COMPLEX, INTENT (IN)  :: ccof(-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (IN)  :: acof(nobd,0:lmd,natd)
      COMPLEX, INTENT (IN)  :: bcof(nobd,0:lmd,natd)
      INTEGER, INTENT (IN)  :: nmtsl(ntypd,natd)
      REAL,    INTENT (OUT) :: qmtslk(nsl,neigd)
C     ..
C     .. Local Scalars ..
      INTEGER i,l,lo,m,natom,nn,ntyp,nt1,nt2
      INTEGER lm,n,ll1,ipol,icore,index,nl
      REAL fac,sabd,ss,qq
      COMPLEX suma,sumb,sumab,sumba
C     ..
C     .. Local Arrays ..
      REAL, ALLOCATABLE :: qlo(:,:,:),qmt(:,:),qmtlo(:,:)
      REAL, ALLOCATABLE :: qaclo(:,:,:),qbclo(:,:,:),qmttot(:,:)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg,cmplx


      ALLOCATE ( qlo(nobd,nlod,ntypd),qmt(ntypd,neigd) )
      ALLOCATE ( qaclo(nobd,nlod,ntypd),qbclo(nobd,nlod,ntypd) )
      ALLOCATE ( qmttot(ntypd,neigd),qmtlo(ntypd,neigd) )
c
c--->    l-decomposed density for each valence state
c
!         DO 140 i = (skip_t+1),ne    ! this I need for all states
         DO 140 i = 1,ne              ! skip in next loop
            nt1 = 1
            DO 130 n = 1,ntype
               fac = 1./neq(n)
               nt2 = nt1 + neq(n) - 1
               sabd = 0.0
               DO 120 l = 0,lmax(n)
                  suma = cmplx(0.,0.)
                  sumb = cmplx(0.,0.)
                  ll1 = l* (l+1)
                  DO 110 m = -l,l
                     lm = ll1 + m
                     DO natom = nt1,nt2
                        suma = suma + acof(i,lm,natom)*
     +                          conjg(acof(i,lm,natom))
                        sumb = sumb + bcof(i,lm,natom)*
     +                          conjg(bcof(i,lm,natom))
                     ENDDO
  110             CONTINUE
                  ss = suma + sumb*ddn(l,n)
                  sabd = sabd + ss
  120          CONTINUE
               qmt(n,i) = sabd*fac
               nt1 = nt1 + neq(n)
  130       CONTINUE
  140    CONTINUE
c                  
c---> initialize qlo
c
      DO ntyp = 1,ntypd
         DO lo = 1,nlod
            DO i = 1,nobd
               qlo(i,lo,ntyp) = 0.0
               qaclo(i,lo,ntyp) = 0.0
               qbclo(i,lo,ntyp) = 0.0
            END DO
         END DO
      END DO
c
c---> density for each local orbital and valence state
c
      natom = 0
      DO ntyp = 1,ntype
         DO nn = 1,neq(ntyp)
            natom = natom + 1
            DO lo = 1,nlo(ntyp)
               l = llo(lo,ntyp)
               ll1 = l* (l+1)
               DO i = 1,ne
                  DO m = -l,l
                     lm = ll1 + m
                     qlo(i,lo,ntyp) = qlo(i,lo,ntyp) +
     +                      ccof(m,i,lo,natom)*conjg(ccof(m,i,lo,natom))
                     qbclo(i,lo,ntyp) = qbclo(i,lo,ntyp) +
     +                      bcof(i,lm,natom)*conjg(ccof(m,i,lo,natom)) +
     +                      ccof(m,i,lo,natom)*conjg(bcof(i,lm,natom))
                     qaclo(i,lo,ntyp) = qaclo(i,lo,ntyp) +
     +                      acof(i,lm,natom)*conjg(ccof(m,i,lo,natom)) +
     +                      ccof(m,i,lo,natom)*conjg(acof(i,lm,natom))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      natom = 1
      DO ntyp = 1,ntype
         IF (invsat(natom).EQ.1) THEN
           DO lo = 1,nlo(ntyp)
              DO i = 1,ne
                qlo(i,lo,ntyp) = 2*qlo(i,lo,ntyp)
              ENDDO
           ENDDO
         ENDIF
         natom = natom + neq(ntyp)
      ENDDO
c
c--->  l-decomposed density for each valence state
c--->      ( a contribution from local orbitals)
c--->                       and
c--->  total  l-decomposed density for each valence state
c
      DO i = 1,ne
         DO ntyp = 1,ntype
	      fac = 1.0/neq(ntyp)
	      qq = 0.0
              DO lo = 1,nlo(ntyp)
		 qq = qq + qlo(i,lo,ntyp)*uloulopn(lo,lo,ntyp) +
     +                     qaclo(i,lo,ntyp)*uulon(lo,ntyp)     +
     +                     qbclo(i,lo,ntyp)*dulon(lo,ntyp)    
              ENDDO
	      qmtlo(ntyp,i) = qq*fac
	      qmttot(ntyp,i) = qmt(ntyp,i) + qmtlo(ntyp,i) 
	   ENDDO
      ENDDO
c
	DO i = 1,ne
         DO nl = 1,nsl
	      qq = 0.0 	
              DO ntyp = 1,ntype
	         qq = qq + qmttot(ntyp,i)*nmtsl(ntyp,nl)
              ENDDO
              qmtslk(nl,i) = qq
	   ENDDO
	ENDDO
!        DO ntyp = 1,ntype
!        write(*,*) qmttot(ntyp,1)
!        write(*,*) (nmtsl(ntyp,nl),nl=1,nsl)
!	ENDDO
c
      DEALLOCATE ( qlo,qmt,qmtlo,qaclo,qbclo,qmttot )

      END SUBROUTINE q_mt_sl
      END MODULE m_qmtsl
