      MODULE m_hsoham
c
c*********************************************************************
c set up spin-orbit contribution to hamiltonian
c*********************************************************************
c
      CONTAINS
      SUBROUTINE hsoham(
     > jspd,neigd,natd,lmaxd,ntypd,nlod,llod,
     > ntype,soc_opt,jspins,nsz,lmax,neq,nlo,llo,chelp,
     > rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop,
     > ahelp,bhelp,rsopp,rsoppd,rsopdp,rsopdpd,soangl,
     < hsomtx)

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,neigd,natd,lmaxd,ntypd,nlod,llod
      INTEGER, INTENT (IN) :: jspins,ntype
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nsz(jspd),lmax(ntypd),neq(ntypd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd)
      LOGICAL, INTENT (IN) :: soc_opt(ntypd+2) 
      REAL,    INTENT (IN) :: rsopp  (ntypd,lmaxd,2,2)
      REAL,    INTENT (IN) :: rsoppd (ntypd,lmaxd,2,2)
      REAL,    INTENT (IN) :: rsopdp (ntypd,lmaxd,2,2)
      REAL,    INTENT (IN) :: rsopdpd(ntypd,lmaxd,2,2)
      REAL,    INTENT (IN) :: rsoplop (ntypd,nlod,2,2)
      REAL,    INTENT (IN) :: rsoplopd(ntypd,nlod,2,2)
      REAL,    INTENT (IN) :: rsopdplo(ntypd,nlod,2,2)
      REAL,    INTENT (IN) :: rsopplo (ntypd,nlod,2,2)
      REAL,    INTENT (IN) :: rsoploplop(ntypd,nlod,nlod,2,2)
      COMPLEX, INTENT (IN) :: ahelp(-lmaxd:lmaxd,lmaxd,natd,neigd,jspd)
      COMPLEX, INTENT (IN) :: bhelp(-lmaxd:lmaxd,lmaxd,natd,neigd,jspd)
      COMPLEX, INTENT (IN) :: chelp(-llod :llod ,neigd,nlod,natd ,jspd)
      COMPLEX, INTENT (IN) :: soangl(lmaxd,-lmaxd:lmaxd,2,
     +                               lmaxd,-lmaxd:lmaxd,2)
      COMPLEX, INTENT (OUT):: hsomtx(2,2,neigd,neigd)
C     ..
C     .. Local Scalars ..
      COMPLEX c_1,c_2,c_3,c_4,c_5
      INTEGER i,j,jsp,jsp1,l,lwn,m,m1,n,na,nn,i1,j1,ilo,ilop
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: c_b(:,:,:),c_a(:,:,:),c_c(:,:,:)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,conjg
c
c---------------------------------------------------------------------
c  ss'  _
c H  = \  (xhelp(s,i,na,l,m) conjg(yhelp(s',j,na,l,m')*rsoxy(na,l,s,s')
c           *<slm|L*S|s'lm'>
c  ij  /_
c       na,l,m,m'
c                       x,y = a,b
c---------------------------------------------------------------------
c
c---> update hamiltonian matrices: upper triangle
c
      ALLOCATE ( c_b(-lmaxd:lmaxd,lmaxd,natd),
     +           c_a(-lmaxd:lmaxd,lmaxd,natd),
     +           c_c(-llod :llod ,nlod ,natd) )

      DO i1 = 1,2
        jsp = i1
        IF (jspins.EQ.1) jsp = 1
        DO j1 = 1,2
          jsp1 = j1
          IF (jspins.EQ.1) jsp1 = 1
          DO j = 1,nsz(jsp1)
c
c prepare \sum_m' conjg( xhelp(m',l,na,j,jsp1) ) * soangl(l,m,i1,l,m',j1)
c
            na = 0
            DO n = 1,ntype
              DO nn = 1, neq(n)
                na = na + 1
!--> regular part
                DO l = 1,lmax(n)
                  DO m = -l,l
                    c_a(m,l,na) = cmplx(0.,0.)
                    c_b(m,l,na) = cmplx(0.,0.)
                    DO m1 = -l,l
                      c_a(m,l,na) = c_a(m,l,na) + soangl(l,m,i1,l,m1,j1)
     *                                     *conjg(ahelp(m1,l,na,j,jsp1))
                      c_b(m,l,na) = c_b(m,l,na) + soangl(l,m,i1,l,m1,j1)
     *                                     *conjg(bhelp(m1,l,na,j,jsp1))
                    ENDDO
                  ENDDO
                ENDDO
!--> LO contribution
                DO ilo = 1,nlo(n)
                  l = llo(ilo,n)
                  IF (l.GT.0) THEN
                    DO m = -l,l
                      c_c(m,ilo,na) = cmplx(0.,0.)
                      DO m1 = -l,l
                        c_c(m,ilo,na) = c_c(m,ilo,na) + conjg(
     *                   chelp(m1,j,ilo,na,jsp1))*soangl(l,m,i1,l,m1,j1)
                      ENDDO
                    ENDDO
                  ENDIF
                ENDDO
! end lo's
              ENDDO
            ENDDO
c
c continue loop structure
c
            DO i = 1,nsz(jsp)
              hsomtx(i1,j1,i,j) = cmplx(0.,0.)
              na = 0
c
c--->    loop over each atom type
c
              DO 160 n = 1,ntype
                IF ( (.not. soc_opt(ntype+1)) .or. soc_opt(n) ) THEN 

                lwn = lmax(n)
c
c--->    loop over equivalent atoms
c
                DO 150 nn = 1,neq(n)
                  na = na + 1
                  DO l = 1,lwn
c 
                    DO m = -l,l
                      c_1 =   rsopp(n,l,i1,j1) * ahelp(m,l,na,i,jsp) +
     +                       rsopdp(n,l,i1,j1) * bhelp(m,l,na,i,jsp)
                      c_2 =  rsoppd(n,l,i1,j1) * ahelp(m,l,na,i,jsp) +
     +                      rsopdpd(n,l,i1,j1) * bhelp(m,l,na,i,jsp)
                     hsomtx(i1,j1,i,j) = hsomtx(i1,j1,i,j) +
     +                               c_1*c_a(m,l,na) + c_2*c_b(m,l,na)  
                    ENDDO
c 
                  ENDDO
!--> LO contribution
                 DO ilo = 1,nlo(n)
                   l = llo(ilo,n)
                   IF (l.GT.0) THEN
                     DO m = -l,l
                       c_3 = rsopplo(n,ilo,i1,j1) *ahelp(m,l,na,i,jsp) +
     +                      rsopdplo(n,ilo,i1,j1) *bhelp(m,l,na,i,jsp)
                       c_4 = rsoplop(n,ilo,i1,j1) *chelp(m,i,ilo,na,jsp)
                       c_5 =rsoplopd(n,ilo,i1,j1) *chelp(m,i,ilo,na,jsp)
                       hsomtx(i1,j1,i,j) = hsomtx(i1,j1,i,j) + 
     +                               c_4*c_a(m,l,na) + c_5*c_b(m,l,na) +
     +                               c_3*c_c(m,ilo,na)
                     ENDDO
                     DO ilop = 1,nlo(n)
                       IF (llo(ilop,n).EQ.l) THEN
                         DO m = -l,l
                           hsomtx(i1,j1,i,j) = hsomtx(i1,j1,i,j) + 
     +                                    rsoploplop(n,ilop,ilo,i1,j1) * 
     +                            chelp(m,i,ilop,na,jsp) * c_c(m,ilo,na)
                         ENDDO
                       ENDIF
                     ENDDO
                   ENDIF
                 ENDDO
! end lo's
  150           ENDDO

                ELSE

                  na = na + neq(n) 

                ENDIF 
  160         ENDDO
c
            ENDDO
c!i
          ENDDO
c!j
        ENDDO
c!jsp1
      ENDDO
c!jsp
      DEALLOCATE (c_a,c_b,c_c)
c
c---> update hamiltonian matrices: lower triangle
c
       DO i = 1,nsz(1)
        DO j = 1,nsz(jspins)
          hsomtx(2,1,j,i) = conjg(hsomtx(1,2,i,j))
        ENDDO
      ENDDO
c
      RETURN
      END SUBROUTINE hsoham
      END MODULE m_hsoham
