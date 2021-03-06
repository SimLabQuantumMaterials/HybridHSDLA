      MODULE m_umtx
*********************************************************************
* The calculation of the "U"-contribution to Hartree-Fock matrix.   *
*********************************************************************
      CONTAINS
      SUBROUTINE umtx(
     >                lmaxb,ntype,n_u,lda_u,f0,f2,f4,f6,
     <                u)

      USE m_constants, ONLY : pimach
      USE m_sgaunt
      IMPLICIT NONE

      INTEGER, PARAMETER   :: lmaxw=3,lmmaxw1=(2*lmaxw+2)**2
      INTEGER, INTENT (IN) :: n_u,lmaxb,ntype
      INTEGER, INTENT (IN) :: lda_u(ntype)
      REAL,    INTENT (IN) :: f0(n_u),f2(n_u),f4(n_u),f6(n_u)
      REAL,    INTENT (OUT) :: u(-lmaxb:lmaxb,-lmaxb:lmaxb,
     +                           -lmaxb:lmaxb,-lmaxb:lmaxb,n_u)
      
      INTEGER i,j,k,l,m,mk,nfk,n,itype
      INTEGER m1,m2,m3,m4,lm1,lm2,lm3,lm4,kf,l_l(n_u)
      REAL    uk,uq,avu,avj,cgk1,cgk2,pi,tol
      REAL    fk(lmaxb+1,n_u)
      REAL,   ALLOCATABLE :: c(:,:,:)
c
      pi = pimach()
      tol = 1.0e-14
C
C transformation to Hr-units:
C
      n = 0
      DO itype = 1,ntype
         IF (lda_u(itype).GE.0) THEN
           n = n + 1
           l_l(n) = lda_u(itype)
           fk(1,n) = f0(n) / 27.21
           fk(2,n) = f2(n) / 27.21
           fk(3,n) = f4(n) / 27.21
           IF (l_l(n).EQ.3) THEN
              fk(4,n) = f6(n) / 27.21
           ELSEIF (l_l(n).GT.3) THEN
              STOP ' LDA+U for p, d or f-states! '
           ENDIF
         ENDIF
      ENDDO
c
c evaluate Gaunt parameter
c
      ALLOCATE( c(0:2*lmaxw+1,lmmaxw1,lmmaxw1) )
      DO k = 1,lmmaxw1
        DO j = 1,lmmaxw1
          DO i = 0,2*lmaxw+1
            c(i,j,k) = 0.0
          ENDDO
        ENDDO
      ENDDO
      CALL sgaunt(lmaxw,lmmaxw1,lmaxb,pi,
     <            c)
c
c lda_u(n) is here only the 'l' for atom 'n'
c
      DO 150 n = 1,n_u                      !!! over d-atoms
         l = l_l(n)
         kf = 2*l
         DO m1 = -l,l
           lm1 = l*(l+1)+m1+1
           DO m2 = -l,l
             lm2 = l*(l+1)+m2+1
             DO m3 = -l,l
               lm3 = l*(l+1)+m3+1
               DO m4 = -l,l
                  lm4 = l*(l+1)+m4+1
                  uk = 0.0e0
                  DO k=0,kf,2
                    uq = 0.e0
                    DO mk=-k,k
                      IF (mk.NE.m1-m3)  GOTO 3
                      cgk1 = c(k/2,lm1,lm3)
                      IF (abs(cgk1).LT.tol) GOTO 3
                      IF (mk.NE.m4-m2)  GOTO 3
                      cgk2 = c(k/2,lm4,lm2)
                      IF (abs(cgk2).LT.tol) GOTO 3
                      uq = uq+cgk1*cgk2
 3                    CONTINUE
                    ENDDO                   ! mk
                    IF (abs(uq).LT.tol) GOTO 2
                    nfk=k/2+1
                    uk=uk+uq*fk(nfk,n)*4*pi/(2*k+1)
 2                  CONTINUE
                  ENDDO                     ! k
                  u(m1,m2,m3,m4,n)=uk
               ENDDO                       ! m4 etc.
             ENDDO
           ENDDO
         ENDDO
         avu=0.e0
         avj=0.e0
	       
         DO i = -l,l
            DO j = -l,l
               avu = avu+u(i,j,i,j,n)
               avj = avj+(u(i,j,i,j,n)-u(i,j,j,i,n))
            ENDDO
         ENDDO
         avu = avu/(2*l+1)/(2*l+1)
         avj = avj/(2*l+1)/(2*l)
         avj = avu-avJ
!        WRITE (6,*) 'U-matr:'
!        IF (l.eq.2) WRITE (6,111) ((u(i,j,i,j,n),i=-l,l),j=-l,l)
!        IF (l.eq.3) WRITE (6,211) ((u(i,j,i,j,n),i=-l,l),j=-l,l)
!      	 WRITE (6,*) 'J-matr:'
!        IF (l.eq.2) WRITE (6,111) ((u(i,j,j,i,n),i=-l,l),j=-l,l)
!        IF (l.eq.3) WRITE (6,211) ((u(i,j,j,i,n),i=-l,l),j=-l,l)
!         PRINT*,'U-av:',avu*27.21
!         PRINT*,'J-av:',avj*27.21
  111    FORMAT (5f8.4)
  211    FORMAT (7f8.4)
  112    FORMAT (10e20.10)
cc         WRITE (9,112) ((((u(i,j,k,m,n),i=-l,l),j=-l,l),k=-l,l),m=-l,l)
  150 ENDDO
      DEALLOCATE (c)

      END SUBROUTINE umtx
      END MODULE m_umtx
