      MODULE m_efield
      CONTAINS
      SUBROUTINE efield(
     >                  nwdd,ntypd,nstd,ncst,
     >                  ntype,nwd,area,neq,zelec,zatom,
     <                  zsigma,sigma,sig_b)
c
c*********************************************************************
c     sets the values of the sheets of charge for external electric
c     fields by requiring charge neutality.
c     the position of the sheets of charge relative to the vaccum
c     boundary is fixed (10 a.u.), but can be set to a different
c     value in the file apwefl.
c
c     modified and fixed 10-99 
c*********************************************************************

      USE m_setcor
      IMPLICIT NONE
c     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nwdd,ntypd,nstd
      INTEGER, INTENT (IN) :: ntype,nwd
      REAL,    INTENT (IN) :: area
      REAL,    INTENT (OUT):: zsigma,sigma,sig_b(2)
c     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),ncst(ntypd)
      REAL,    INTENT (IN) :: zelec(nwdd),zatom(ntypd)
C     ..
C     .. Local Scalars ..
      REAL qn,qe,sigm,efi,bmu
      INTEGER n,iwd,nst,nc,ios
      LOGICAL ext33
C     ..
C     .. Local Array ..
      INTEGER kappa(nstd),nprnc(nstd)
      REAL occ(nstd)
c
c--->    position of sheets relative to the vacuum boundary
c--->    read off of unit33 apwefl if it exists, or use default value
c
      INQUIRE (file='apwefl',exist=ext33)
      IF (ext33) THEN
         OPEN (33,file='apwefl',form='formatted',status='old')
         READ (33,*) zsigma
         READ (33,*,IOSTAT=ios) sig_b(:)
         IF (ios.ne.0) sig_b(:) = 0.0
         CLOSE(33)
      ELSE
         zsigma=10.
      ENDIF
c
c--->    obtain total nuclear charge
c
      qn=0.0
      DO n=1,ntype
         qn=qn+neq(n)*zatom(n)
      ENDDO
c
c--->    obtain total electronic charge (in electron units)
c
      qe=0.0
c---> core electrons 
      DO n = 1,ntype
         IF (zatom(n).GE.1.0) THEN
            bmu = 0.0
            CALL setcor(
     >                  zatom(n),nstd,1,1,bmu,
     <                  nst,kappa,nprnc,occ)
            DO nc=1,ncst(n)
              qe=qe+neq(n)*occ(nc)
            ENDDO
            WRITE (6,*) 'neq= ',neq(n),'  ncore= ',qe
         ENDIF
      ENDDO
c---> semi-core and valence electrons
      DO iwd=1,nwd
        qe=qe+zelec(iwd)
        WRITE (6,*) 'zelec=  ',zelec(iwd)
      ENDDO
c
c--->   this choice of sign means that sheet of charge is given in
c--->   units of nuclear charge; then enters as -sigma in Coulomb
c--->   potential. charge neutrality is qn + 2 sigma = qe
c
      sigma = 0.5*(qe-qn)
      sigm = sigma/area
c
c--->   sign of electric field: E>0 repels electrons, consistent with
c--->   conventional choice that current is opposite to electron flow.
c--->   with this choice, qe < qn corresponds to E>0
c
      IF ( maxval(abs(sig_b(:))) < 0.0000001 ) THEN  ! test
      efi = -(6.462e+10)*sigm
c
      WRITE (6,1000) qe,qn,sigma,sigm,efi
      WRITE (6,1001) zsigma
      IF (efi.GT.0.00001) THEN
         WRITE (16,1000) qe,qn,sigma,sigm,efi
         WRITE (16,1001) zsigma
      ENDIF
 1000 FORMAT (1x//,' parameters for external electric field:'/,
     +  3x,'total electronic charge   =',f12.5/,
     +  3x,'total nuclear charge      =',f12.5/,
     +  3x,'charge per external sheet =',f12.5,5x,'(surface density=',
     +  f12.5,'  e/a.u.**2)'/,
     +  3x,'external field            =',e12.5,5x,'V/cm')
 1001 FORMAT (3x,'z-z1 of external sheet    =',f12.5/)

      ELSE

      WRITE (6,1002) qe,qn
      WRITE (6,1001) zsigma
      efi = -(6.462e+10)*(sigm + sig_b(1)/area)
      WRITE (6,1003) 1,sigma+sig_b(1),(sigm+ sig_b(1)/area),1,efi
      efi = -(6.462e+10)*(sigm + sig_b(2)/area)
      WRITE (6,1003) 2,sigma+sig_b(2),(sigm+ sig_b(2)/area),2,efi

 1002 FORMAT (1x//,' parameters for external electric field:'/,
     +  3x,'total electronic charge   =',f12.5/,
     +  3x,'total nuclear charge      =',f12.5/)
 1003 FORMAT (3x,'charge on external sheet ',i1,' = ',f12.5,5x,
     +  '(surface density=',f12.5,'  e/a.u.**2)'/,
     +  3x,'external field on sheet ',i1,' = ',e12.5,5x,'V/cm')

      ENDIF

      IF (abs(sigma).GT.0.49) THEN
        WRITE ( 6,*) 'If you really want to calculate an e-field this'
        WRITE ( 6,*) 'big, uncomment STOP in efield.f !'
        WRITE (16,*) 'If you really want to calculate an e-field this' 
        WRITE (16,*) 'big, uncomment STOP in efield.f !'
        STOP 'E-field too big or No. of e- not correct'
      ENDIF

      END SUBROUTINE efield
      END MODULE m_efield
