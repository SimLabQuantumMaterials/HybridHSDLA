      MODULE m_qal21 
!***********************************************************************
! Calculates qal21  needed to determine the off-diagonal parts of the 
! DOS
!***********************************************************************
c
      CONTAINS
      SUBROUTINE qal_21(
     >                  llod,nlod,natd,neigd,ntypd,lmaxd,lmd,
     >                  jspins,ntype,neq,noccbd,we,ccof,nlo,llo,
     >                  alph,beta,acof,bcof,mt21,lo21,uloulopn21,
     <                  qal,qmat)

      USE m_types, ONLY : t_mt21, t_lo21
      USE m_rotdenmat

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: llod,nlod,natd,neigd,ntypd,lmaxd,lmd
      INTEGER, INTENT (IN) :: ntype,noccbd,jspins
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: neq(ntypd)
      INTEGER, INTENT (IN)  :: nlo(ntypd),llo(nlod,ntypd)
      REAL,    INTENT (INout)  :: we(noccbd),qal(0:3,ntypd,neigd,jspins)
      REAL,    INTENT (IN)  :: uloulopn21(nlod,nlod,ntypd)
      REAL,    INTENT (IN)  :: alph(ntypd),beta(ntypd)
      COMPLEX, INTENT (IN)  :: ccof(-llod:llod,noccbd,nlod,natd,jspins)
      COMPLEX, INTENT (IN)  :: acof(noccbd,0:lmd,natd,jspins)
      COMPLEX, INTENT (IN)  :: bcof(noccbd,0:lmd,natd,jspins)
      REAL,    INTENT (OUT) :: qmat(0:3,ntypd,neigd,4)
      TYPE (t_mt21), INTENT (IN) :: mt21(0:lmaxd,ntypd)
      TYPE (t_lo21), INTENT (IN) :: lo21(0:lmaxd,ntypd)

C     ..
C     .. Local Scalars ..
      INTEGER i,l,lo,lop,m,natom,nn,ntyp
      INTEGER nt1,nt2,lm,n,ll1,ipol,icore,index
      REAL fac
      COMPLEX sumaa,sumbb,sumab,sumba
      COMPLEX, PARAMETER :: ci = (0.0,1.0)

C     ..
C     .. Local Arrays ..
      COMPLEX qlo(noccbd,nlod,nlod,ntypd)
      COMPLEX qaclo(noccbd,nlod,ntypd),qbclo(noccbd,nlod,ntypd)
      COMPLEX qcloa(noccbd,nlod,ntypd),qclob(noccbd,nlod,ntypd)
      COMPLEX qal21(0:3,ntypd,neigd)
      COMPLEX q_loc(2,2),q_hlp(2,2),chi(2,2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg
c
c--->    l-decomposed density for each occupied state
c
      states : DO i = 1, noccbd
         nt1 = 1
         types : DO n = 1,ntype
            nt2 = nt1 + neq(n) - 1
            ls : DO l = 0,3
               IF (i==1) THEN
               ENDIF
               sumaa = cmplx(0.,0.) ; sumab = cmplx(0.,0.) 
               sumbb = cmplx(0.,0.) ; sumba = cmplx(0.,0.)
               ll1 = l* (l+1)
               ms : DO m = -l,l
                  lm = ll1 + m
                  atoms : DO natom = nt1,nt2
                    sumaa = sumaa + acof(i,lm,natom,1)*
     +                        conjg(acof(i,lm,natom,jspins))
                    sumbb = sumbb + bcof(i,lm,natom,1)*
     +                        conjg(bcof(i,lm,natom,jspins))
                    sumba = sumba + acof(i,lm,natom,1) *
     +                        conjg(bcof(i,lm,natom,jspins))
                    sumab = sumab + bcof(i,lm,natom,1) *
     +                        conjg(acof(i,lm,natom,jspins))
                  ENDDO atoms
               ENDDO ms
               qal21(l,n,i) = sumaa * mt21(l,n)%uun +
     +                        sumbb * mt21(l,n)%ddn +
     +                        sumba * mt21(l,n)%dun +
     +                        sumab * mt21(l,n)%udn 
             ENDDO ls
            nt1 = nt1 + neq(n)
         ENDDO types
      ENDDO states

c---> initialize qlo

      qlo(:,:,:,:) = cmplx(0.,0.)
      qaclo(:,:,:) = cmplx(0.,0.)
      qcloa(:,:,:) = cmplx(0.,0.)
      qclob(:,:,:) = cmplx(0.,0.)
      qbclo(:,:,:) = cmplx(0.,0.)

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
                  DO i = 1, noccbd
                     qbclo(i,lo,ntyp) = qbclo(i,lo,ntyp) +      
     +               bcof(i,lm,natom,1)*conjg(ccof(m,i,lo,natom,jspins)) 
                     qbclo(i,lo,ntyp) = qbclo(i,lo,ntyp) +      
     +               ccof(m,i,lo,natom,1)*conjg(bcof(i,lm,natom,jspins)) 
                     qaclo(i,lo,ntyp) = qaclo(i,lo,ntyp) +       
     +               acof(i,lm,natom,1)*conjg(ccof(m,i,lo,natom,jspins)) 
                     qaclo(i,lo,ntyp) = qaclo(i,lo,ntyp) +       
     +               ccof(m,i,lo,natom,1)*conjg(acof(i,lm,natom,jspins)) 
                  ENDDO
               ENDDO
               DO lop = 1,nlo(ntyp)
                 IF (llo(lop,ntyp).EQ.l) THEN
                   DO m = -l,l
                     DO i = 1, noccbd
                       qlo(i,lop,lo,ntyp) = qlo(i,lop,lo,ntyp) +  
     +          conjg(ccof(m,i,lop,natom,jspins))*ccof(m,i,lo,natom,1) +
     +          conjg(ccof(m,i,lo,natom,jspins))*ccof(m,i,lop,natom,1)
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
            DO i = 1, noccbd
               qal21(l,ntyp,i)= qal21(l,ntyp,i)  + 
     +                        qaclo(i,lo,ntyp)*lo21(lo,ntyp)%uulon +
     +                        qcloa(i,lo,ntyp)*lo21(lo,ntyp)%uloun +
     +                        qclob(i,lo,ntyp)*lo21(lo,ntyp)%ulodn +
     +                        qbclo(i,lo,ntyp)*lo21(lo,ntyp)%dulon 
            END DO
            DO lop = 1,nlo(ntyp)
               IF (llo(lop,ntyp).EQ.l) THEN
               DO i = 1, noccbd
                 qal21(l,ntyp,i)= qal21(l,ntyp,i)  + 
     +                    qlo(i,lop,lo,ntyp)*uloulopn21(lop,lo,ntyp)
               ENDDO
               ENDIF
            ENDDO
         END DO
      END DO

      DO n = 1,ntype
        fac = 1./neq(n)
        qal21(:,n,:) = qal21(:,n,:) * fac
      ENDDO
!
! rotate into global frame
!
      type : DO n = 1,ntype
         chi(1,1) =  exp(-ci*alph(n)/2)*cos(beta(n)/2)
         chi(1,2) = -exp(-ci*alph(n)/2)*sin(beta(n)/2)
         chi(2,1) =  exp( ci*alph(n)/2)*sin(beta(n)/2)
         chi(2,2) =  exp( ci*alph(n)/2)*cos(beta(n)/2)
         state : DO i = 1, noccbd
            lls : DO l = 0,3
               CALL rot_den_mat(alph(n),beta(n),
     +                          qal(l,n,i,1),qal(l,n,i,2),qal21(l,n,i))
               IF (.false.) THEN
               IF (n==1) write(*,'(3i3,4f10.5)') l,n,i,qal21(l,n,i),
     +                                                  qal(l,n,i,:)
               q_loc(1,1) = qal(l,n,i,1); q_loc(2,2) = qal(l,n,i,2)
               q_loc(1,2) = qal21(l,n,i); q_loc(2,1) = conjg(q_loc(1,2))
               q_hlp = matmul( transpose( conjg(chi) ) ,q_loc)
               q_loc = matmul(q_hlp,chi)
               qmat(l,n,i,1) = real(q_loc(1,1))
               qmat(l,n,i,2) = real(q_loc(1,2))
               qmat(l,n,i,3) = aimag(q_loc(1,2))
               qmat(l,n,i,4) = real(q_loc(2,2))
               IF (n==1) write(*,'(3i3,4f10.5)') l,n,i,qmat(l,n,i,:)
               ENDIF
             ENDDO lls
         ENDDO state
      ENDDO type

      END SUBROUTINE qal_21 
      END MODULE m_qal21 
