      MODULE m_ferhis
      CONTAINS
      SUBROUTINE ferhis(jspins,nkpt,tkb,index,idxeig,idxkpt,idxjsp,n,
     +              nstef,ws,weight,spindg,ef,seigv,ts,
     >              neigd,nkptd,jspd,e,ne,wtkpt,
     X              we,
     <              w,
     >              qss,l_J,l_disp,bmat)
c***********************************************************************
c
C     This subroutine determines the fermi energy and the sum of the
c     single particle eigenvalues by histogram method.
C
C
C     Theory :   zelec(nwd) = sum{ sum{ sum{ we * f(e) } } }
C                             sp    e    k
C
C
C                seigv = sum{ sum{ sum{ w(k) * f(e) * e }
C                         sp   e    k
C
c                the weights w(k) are normalized: sum{w(k)} = 1
C                                                  k                -6
C                         a) 1                           for kT < 10
c                we    = {                           -1             -6
c                         b){ 1+exp(e(k,nu) -ef)/kt) }   for kt >=10
C
C                in case a) we choose the Fermi energy the highest
C                           valence state
C                        b) we choose as Fermi energy the arithmetric
C                           mean between the highest occupied and lowest
C                           unoccupied state if this energy difference
C                           Delta E <= kT, otherwise as a).
C
c                                      stefan bl"ugel , kfa , oct 1991
C
C               free energy and extrapolation T -> 0  implemented
C                         (see M.J.Gillan, J.Phys.: Condens. Matter 1,
C                          (1989) 689-711 )
C
C                                      peter richard, kfa, jun. 1995
c
c               adapted to flapw7
c
c                                      philipp kurz, kfa, oct. 1995
c               entropy calculation changed
c                     
c                                      r.pentcheva, kfa, may  1996
c
c***********************************************************************
      USE m_efnewton

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: neigd,nkptd,jspd
      INTEGER jspins,n,nkpt,nstef
      REAL ef,seigv,spindg,tkb,ts,weight,ws
      LOGICAL, INTENT (IN) :: l_disp,l_J
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: idxeig(neigd*nkptd*jspd)
      INTEGER, INTENT (IN) :: idxjsp(neigd*nkptd*jspd)
      INTEGER, INTENT (IN) :: idxkpt(neigd*nkptd*jspd)
      INTEGER, INTENT (IN) ::  index(neigd*nkptd*jspd)
      INTEGER, INTENT (IN) ::     ne(nkptd,jspd)
      REAL,    INTENT (IN) ::  wtkpt(nkptd)
      REAL,    INTENT (IN) ::      e(nkptd*neigd*jspd)
      REAL,    INTENT (OUT)::      w(neigd,nkptd,jspd)
      REAL,    INTENT (INOUT) ::  we(nkptd*neigd*jspd)

c--- J constants
      REAL,    INTENT (IN) :: qss(3)
      REAL,    INTENT (IN) :: bmat(3,3)
c--- J constants

C     ..
C     .. Local Scalars ..
      REAL del,efermi,eight,emax,emin,entropy,fermikn,gap,ql,
     +     half,kb,one,wfermi,wvals,w_below_emin,w_near_ef,zero
      INTEGER ink,inkem,j,js,k,kpt,nocc,nocst

C     .. Local Arrays ..      
      REAL :: qc(3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC log
C     ..
c***********************************************************************
C------->          ABBREVIATIONS
C
c     eig        : array of eigenvalues within all energy-windows
c     wtkpt      : list of the weights of each k-point (from inp-file)
c     e          : linear list of the eigenvalues within the highest
c                  energy-window
c     we         : list of weights of the eigenvalues in e
c     w          : array of weights (output, needed to calculate the
c                  new charge-density)
c     zelec      : number of electrons in a window
c     spindg     : spindegeneracy (2 in nonmagnetic calculations)
c     seigv      : weighted sum of the occupied valence eigenvalues
c     seigsc     : weighted sum of the semi-core eigenvalues
c     seigscv    : sum of seigv and seigsc
C     ts         : entropy contribution to the free energy
c     kb         : boltzmann constant in htr/k
c     tkb        : value of temperature (kt) broadening around fermi
c                  energy in htr units
c     ef         : fermi energy determined to obtain charge neutrality
c     wfermi     : weight of the occupation number for states at the
c                  fermi energy.
c     fd         : fermi dirac distribution
c     fermikn    : fermi dirac distribution for the k-th point 
c                  and n-th state
C**********************************************************************
C     ..
C     .. Data statements ..
      DATA del/1.0e-6/,zero/0.0/,half/0.5/,one/1.0/,
     +     eight/8.0/,kb/3.0553e-6/
C     ..
      WRITE (6,FMT='(/)')
      WRITE (6,FMT='(''FERHIS:  Fermi-Energy by histogram:'')')

      efermi = ef
      IF (nstef.LT.n) THEN
         gap = e(index(nstef+1)) - ef
         WRITE (6,FMT=8050) gap
      END IF
      WRITE (6,FMT=8010) spindg* (ws-weight)
      WRITE (16,FMT=8010) spindg* (ws-weight)
C
C---> DETERMINE OCCUPATION AT THE FERMI LEVEL
C
      wfermi = ws - weight
C                                          -6
C======> DETERMINE FERMI ENERGY for kT >= 10
C
C
      IF (tkb.GE.del) THEN
c
c---> TEMPERATURE BROADENING
c
         IF (nstef.LT.n) THEN
c
c--->    STATES ABOVE EF AVAILABLE           
c
            ef = half* (e(index(nstef+1))+ef)
            emax = ef + eight*tkb
            emin = ef - eight*tkb
            w_near_ef = zero
            w_below_emin = zero
            inkem = 0
            DO 10 ink = 1,n

               IF (e(index(ink)).LT.emin) THEN
                  inkem = ink
                  w_below_emin = w_below_emin + we(index(ink))
               ELSE IF (e(index(ink)).GT.emax) THEN
                  GO TO 20
               END IF

   10       CONTINUE

            WRITE (6,*) 'CAUTION!!!  All calculated eigenvalues are ',
     +                  'below ef + 8kt.'
            WRITE (16,*) 'CAUTION!!!  All calculated eigenvalues are ',
     +                   'below ef + 8kt.'

   20       CONTINUE

            w_near_ef = weight - w_below_emin

            IF (w_near_ef.GT.del) THEN
C
c--->       STATES BETWEEN EF-8kt AND EF+8kt AVAILABLE:
C--->            ADJUST FERMI-ENERGY BY NEWTON-METHOD
C
               nocst = ink - 1
               CALL ef_newton(
     >                        nkptd,neigd,jspd,n,
     >                        inkem,nocst,index,tkb,e,
     X                        w_near_ef,ef,we)
c
               WRITE (16,FMT=8030) ef,spindg*weight,spindg*w_below_emin,
     +           spindg* (w_below_emin+w_near_ef)
               WRITE (6,FMT=8030) ef,spindg*weight,spindg*w_below_emin,
     +           spindg* (w_below_emin+w_near_ef)

            ELSE
c
c--->       NO STATES BETWEEN EF-8kt AND EF+8kt AVAILABLE
c
               WRITE (6,FMT=8020)
               nocst = nstef
               we(index(nocst)) = we(index(nocst)) - wfermi
               ef = efermi
               tkb = zero
            END IF
         ELSE
c
c--->    NO STATES ABOVE EF AVAILABLE
c
            tkb = zero
            nocst = nstef
            we(index(nocst)) = we(index(nocst)) - wfermi
         END IF

      ELSE
c
c---> NO TEMPERATURE BROADENING IF tkb < del
c
         nocst = nstef
         we(index(nocst)) = we(index(nocst)) - wfermi
      END IF
c
c      write(6,*) nocst,'    nocst in ferhis'
c      do  ink = 1,nocst
c         write(6,*) ink,index(ink),we(index(ink)),
c     +      '    ink,index(ink),we(index(ink)): weights for eigenvalues'
c      end do
C
C
C=======>   DETERMINE OCCUPATION NUMBER AND WEIGHT OF EIGENVALUES
C                     FOR EACH K_POINT
C
      DO 50 js = 1,jspins
         DO 40 k = 1,nkpt
            DO 30 j = 1,neigd ! ne(k,js)
               w(j,k,js) = zero
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE

      WRITE (6,FMT=8080) nocst
      DO 60 ink = 1,nocst
         w(idxeig(index(ink)),idxkpt(index(ink)),
     +     idxjsp(index(ink))) = we(index(ink))
c         WRITE (6,*) ink,we(index(ink))
   60 CONTINUE
c
c======>   CHECK SUM OF VALENCE WEIGHTS
c
      wvals = zero
      DO 90 js = 1,jspins
         DO 80 k = 1,nkpt
            DO 70 j = 1,ne(k,js)
               wvals = wvals + w(j,k,js)
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
      WRITE (6,FMT=8070) wvals
C
C
C=======>   DETERMINE ENTROPY
C
c
c --->   formula for entropy:
c
c        entropy = - two * sum wtkpt(kpt) * 
c                          kpt
c                       { sum ( f(e(kpt,nu))*log(f(e(kpt,nu)))
c                          n
c                              +(1-f(e(kpt,nu)))*log(1-f(e(kpt,nu))) )  }
c
c        here we have   w(n,kpt,js)= spindg*wghtkp(kpt)*f(e(kpt,n))
c
      entropy = zero
      DO 100 js = 1,jspins
         DO kpt = 1 , nkpt
c            nocc=1
            DO nocc=1,ne(kpt,js) 
c            DO WHILE ( ( nocc.LE.ne(kpt,js) ) .AND. 
c     +                 ( w(nocc,kpt,js).GT.0.0 ) )
               fermikn = w(nocc,kpt,js)/wtkpt(kpt)
               IF ( fermikn .GT. zero .AND. fermikn .LT. one )
     +              entropy = entropy + wtkpt(kpt) * 
     +                   ( fermikn * log( fermikn)
     +                 + ( one - fermikn) * log( one - fermikn) )
c               WRITE (6,FMT='(i5,i5,f10.6,f10.8)') nocc,kpt,entropy,
c     +         fermikn
c               nocc=nocc+1
            END DO
         END DO
 100  CONTINUE
      entropy = -spindg*entropy
      ts = tkb*entropy
      WRITE (6,FMT=8060) entropy,entropy*kb

C
C=======>   DETERMINE SINGLE PARTICLE ENERGY
C
C
      seigv = zero
      DO 110 ink = 1,nocst
         seigv = seigv + e(index(ink))*we(index(ink))
  110 CONTINUE
      seigv = spindg*seigv
      WRITE (6,FMT=8040) seigv

c--- J constants
      IF (l_J) THEN
        IF (l_disp) THEN
          DO j = 1,3
            qc(j) = qss(1)*bmat(1,j)+qss(2)*bmat(2,j)+qss(3)*bmat(3,j)
          ENDDO
          ql=sqrt(qc(1)**2+qc(2)**2+qc(3)**2)
          WRITE (114,FMT=1001) qss(1),qss(2),qss(3),ql,seigv
        ELSE
          WRITE (114,FMT=1002) qss(1),qss(2),qss(3),seigv
        ENDIF
      ENDIF
 1001 FORMAT (4(f14.10,1x),f20.10)
 1002 FORMAT (3(f14.10,1x),f20.10)
c--- J constants 

c
c 7.12.95 r.pentcheva   seigscv = seigsc + seigv   will be
c calculated in fermie
c
 8000 FORMAT (/,10x,'==>efrmhi: not enough wavefunctions.',i10,2e20.10)
 8010 FORMAT (10x,'charge neutrality (T=0)     :',f11.6,'    (zero if ',
     +       'the highest occ. eigenvalue is "entirely" occupied)')
 8020 FORMAT (/,10x,'no eigenvalues within 8 tkb of ef',
     +       ' reverts to the t=0 k method.')
 8030 FORMAT (/,5x,'-->  new fermi energy            :',f11.6,' htr',
     +       /,10x,'valence charge              :',f11.6,' e ',/,10x,
     +       'actual charge blw ef-8kt    :',f11.6,' e ',/,10x,
     +       'actual charge blw ef+8kt    :',f11.6,' e ')
 8040 FORMAT (/,10x,'sum of val. single particle energies: ',f20.10,
     +       ' htr',/)
 8050 FORMAT (/,10x,'bandgap                     :',f11.6,' htr')
 8060 FORMAT (10x,'entropy         :',f11.6,' *kb htr/K =',
     +       f10.5,' htr/K')
 8070 FORMAT (10x,'sum of the valence weights  :',f11.6)
 8080 FORMAT (10x,'number of occ. states       :',i10)

      END SUBROUTINE ferhis
      END MODULE m_ferhis
