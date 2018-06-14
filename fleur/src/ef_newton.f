      MODULE m_efnewton
      CONTAINS
      SUBROUTINE ef_newton(
     >                     nkptd,neigd,jspd,n,
     >                     inkem,nocst,index,tkb,e,
     X                     w_near_ef,ef,we)

c***********************************************************************
c
c     This subroutine adjusts the Fermi-Energy to obtain charge
c     neutrality. This is done by solving
c
c     ( sum ff(e) * we ) - w_near_ef = 0
c        e
c
c     using the Newton-Method.
c     Here the sum is over the eigenvalues between ef-8kt and ef+8kt, 
c     ff is the Fermi-Function, we is the weight of the according 
c     k-point and w_near_ef is the weight required between ef-8kt 
c     and ef+8kt to obtain neutrality.
c
c***********************************************************************

      IMPLICIT NONE

C     .. Parameters ..
      REAL zero,one,eps
      PARAMETER (zero=0.0e0,one=1.0e0,eps=1.0e-10)
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN)    :: nkptd,neigd,jspd,n
      INTEGER, INTENT (IN)    :: inkem,nocst
      REAL,    INTENT (IN)    :: tkb
      REAL,    INTENT (INOUT) :: ef,w_near_ef
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)    :: index(n)
      REAL,    INTENT (IN)    :: e(nkptd*neigd*jspd)
      REAL,    INTENT (INOUT) :: we(nkptd*neigd*jspd)
C     ..
C     .. Local Scalars ..
      REAL dff,expo,sdff,sff
      INTEGER icnt,idim,itry
C     ..
C     .. Local Arrays ..
      REAL ff(neigd*nkptd*jspd)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp
C     ..
c***********************************************************************
c                  ABBREVIATIONS
c
c     e          : linear list of the eigenvalues within the highest
c                  energy-window
c     we         : list of weights of the eigenvalues in e
c     ef         : fermi energy determined to obtain charge neutrality
c     tkb        : value of temperature (kt) broadening around fermi
c                  energy in htr units
c     index      : index list, e(index(n)) is the list of the 
c                  eigenvalues in ascending order
c     ikem       : number of eigenvalues below ef-8kt
c     nocst      : number of eigenvalues below ef+8kt
c     w_near_ef  : weight (charge/spindg) between ef-8kt and ef+8kt
c                  needed to obtain charge neutrality
c
c***********************************************************************

      WRITE (6,FMT='(/,5x,''EF_NEWTON:  '',
     +               ''Adjust Fermi-Energy by Newton-Method.'',/)')
c
c --->    NEWTON CYCLE
c
      itry = 0
 20   CONTINUE
      DO icnt = 1,50
         sff = zero 
         sdff = zero
         DO idim = inkem + 1,nocst
c
c --->    COMPUTE FERMI-FUNCTION
c
            expo = exp((e(index(idim))-ef)/tkb)
            ff(idim) = one/ (expo+one)
c
c --->    COMPUTE THE DERIVATIVE
c
            dff = ff(idim)*ff(idim)*expo/tkb
c
c --->    MULTIPLY WITH THE WEIGHT
c
            ff(idim) = we(index(idim))*ff(idim)
            dff = we(index(idim))*dff
c
c --->    AND SUM IT UP
c
            sff = sff + ff(idim)
            sdff = sdff + dff
         END DO
         sff = sff - w_near_ef
         IF (abs(sff).LT.eps) GO TO 10
         IF (abs(sdff).LT.1e-29) THEN
            WRITE(*,*) 'Instability in determination of fermi-level'
            STOP 'ef_newton: change temperature broadening (tkb)'
         ENDIF

         ef = ef - sff/sdff
      END DO
c
c---> NOT CONVERGED AFTER 50 ITERATIONS
c
      WRITE (6,FMT=8000) icnt,sff,-sff/sdff
      IF (itry < 20 ) THEN
         ef = ef + 0.001
         itry = itry + 1
         GOTO 20
      ENDIF
      Write(*,*) 'Determination of fermi-level did not converge'
      STOP 'ef_newton: change temperature broad. or set gauss=T'

   10 CONTINUE

      w_near_ef = sff + w_near_ef
      WRITE (6,FMT=8010) icnt,sff,-sff/sdff

      DO idim = inkem + 1,nocst
         we(index(idim)) = ff(idim)
      END DO

 8000 FORMAT (15x,'ef_newton failed after      :',i3,'iterations.',/,
     +        15x,'The error in the weight is  : ',e12.5,/,
     +        15x,'The error in ef is          : ',e12.5,' htr',/)
 8010 FORMAT (15x,'Number of iterations needed  : ',i3,/,
     +        15x,'The error in the weight is   : ',e12.5,/,
     +        15x,'The error in ef is           : ',e12.5,' htr',/)

      END SUBROUTINE ef_newton
      END MODULE m_efnewton
