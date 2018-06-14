      MODULE m_abcrot2
      PRIVATE
      PUBLIC :: abcrot2
      CONTAINS
      SUBROUTINE abcrot2(
     >                 ntypd,natd,neigd,lmaxd,lmd,llod,nlod,ntype,neq,
     >                 neig,lmax,nlo,llo,
     X                 acof,bcof,ccof)
      USE m_dwigner
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,natd,neigd,lmd,llod,nlod,ntype
      INTEGER, INTENT (IN) :: lmaxd,neig
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd),nlo(ntypd)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd)

      COMPLEX, INTENT (INOUT) :: acof(neigd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: bcof(neigd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: ccof(-llod:llod,neigd,nlod,natd)
C     ..
C     .. Local Scalars ..
      INTEGER itype,ineq,iatom,iop,ilo,i,l,m,lm,lmp,ifac
      REAL alpha,beta,gamma
      REAL amx(3,3,1),imx(3,3)
      COMPLEX d_wgn(-lmaxd:lmaxd,-lmaxd:lmaxd,1:lmaxd,1)

      OPEN (333,file='orbcomprot')
      READ (333,*) alpha
      READ (333,*) beta
      READ (333,*) gamma
      CLOSE (333)

      CALL euler(
     >           alpha,beta,gamma,
     <           amx)

      imx(:,:) = 0. ; imx(1,1) = 1. ; imx(2,2) = 1. ; imx(3,3) = 1.

      CALL d_wigner(
     >              1,amx,imx,lmaxd,
     <              d_wgn)

      iatom = 0
      iop = 1
      DO itype = 1, ntype
        DO ineq = 1, neq(itype)
          iatom = iatom + 1
          DO l = 1, lmax(itype)

            DO i = 1, neig
              acof(i,l**2:l*(l+2),iatom) = matmul(
     &                                 conjg(d_wgn(-l:l,-l:l,l,iop)),
     &                                 acof(i,l**2:l*(l+2),iatom))
              bcof(i,l**2:l*(l+2),iatom) = matmul(
     &                                 conjg(d_wgn(-l:l,-l:l,l,iop)),
     &                                 bcof(i,l**2:l*(l+2),iatom))
            ENDDO
          ENDDO
          DO ilo = 1, nlo(itype)
            l = llo(ilo,itype)
            IF (l.gt.0) THEN
              DO i = 1 ,neig
                ccof(-l:l,i,ilo,iatom) = matmul(
     &                               conjg(d_wgn(-l:l,-l:l,l,iop)),
     &                               ccof(-l:l,i,ilo,iatom))
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      END SUBROUTINE abcrot2

c********************************************************************
c********************************************************************
      SUBROUTINE euler(alpha,beta,gamma,amx)

      IMPLICIT NONE

      REAL,    INTENT (IN)  :: alpha,beta,gamma 
      REAL,    INTENT (OUT) :: amx(3,3,1)

      REAL alph,bet,gamm
      REAL bmx(3,3),cmx(3,3),dmx(3,3),hmx(3,3)
      INTEGER nwf,i,j,ii

c..define the D,C,B-matrices
      amx(:,:,:)=0.

      alph = alpha ; bet = beta ; gamm = gamma

      dmx(1,1) = cos(alph) ; dmx(1,2) = sin(alph) ; dmx(1,3) = 0. 
      dmx(2,1) =-sin(alph) ; dmx(2,2) = cos(alph) ; dmx(2,3) = 0. 
      dmx(3,1) = 0.        ; dmx(3,2) = 0.        ; dmx(3,3) = 1. 

      cmx(1,1) = 1.  ; cmx(1,2) = 0.        ; cmx(1,3) = 0. 
      cmx(2,1) = 0.  ; cmx(2,2) = cos(bet)  ; cmx(2,3) = sin(bet)
      cmx(3,1) = 0.  ; cmx(3,2) =-sin(bet)  ; cmx(3,3) = cos(bet)
  
      bmx(1,1) = cos(gamm) ; bmx(1,2) = sin(gamm) ; bmx(1,3) = 0. 
      bmx(2,1) =-sin(gamm) ; bmx(2,2) = cos(gamm) ; bmx(2,3) = 0. 
      bmx(3,1) = 0.        ; bmx(3,2) = 0.        ; bmx(3,3) = 1. 

      hmx(:,:) = 0. 
      DO i = 1,3
        DO j = 1,3
          DO ii = 1,3
             hmx(i,j) = hmx(i,j) + cmx(i,ii)*dmx(ii,j)
          ENDDO  
        ENDDO
      ENDDO 

      DO i = 1,3
        DO j = 1,3
          DO ii = 1,3
             amx(i,j,1) = amx(i,j,1) + bmx(i,ii)*hmx(ii,j)
          ENDDO  
        ENDDO
      ENDDO
 
      END SUBROUTINE euler  

      END MODULE m_abcrot2
