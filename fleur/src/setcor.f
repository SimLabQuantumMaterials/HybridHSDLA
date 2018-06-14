      MODULE m_setcor
      CONTAINS
      SUBROUTINE setcor(
     >                  z,nstd,jspd,jspins,bmu,
     <                  nst,kappa,nprnc,occ)
c
c     *****************************************************
c     sets the values of kappa and occupation numbers of
c     the neutral atoms.
c         m. weinert  february 1982
c     *****************************************************

      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER,INTENT (IN) :: nstd,jspd,jspins
      REAL,INTENT (IN)    :: z
      REAL,INTENT (INOUT) :: bmu
C     ..
      INTEGER,INTENT (OUT) :: nst
      INTEGER,INTENT (OUT) :: kappa(nstd),nprnc(nstd)
      REAL,INTENT (OUT)    :: occ(nstd,jspd)
C     ..
C     .. Local Scalars ..
      INTEGER iz,jz,jz0,k,n,nmax,jspin,d1,d10,aoff
      INTEGER k_h(2),n_h(2)
      REAL fj,l,bmu_l,o_h(2), fac(2)
      LOGICAL l_clf
      CHARACTER(len=13) :: fname
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC min0
C     ..

      l_clf = .false.
      d1  = mod(nint(z),10)
      d10 = int( (nint(z) + 0.5)/10 )
      aoff = iachar('1')-1 
      fname = 'corelevels.'//achar(d10+aoff)//achar(d1+aoff)
      INQUIRE (file=fname, exist=l_clf)

      IF (l_clf) THEN
        OPEN (61,file=fname,form='formatted')
        READ (61,'(i3)') nst
        IF (bmu.LT.0.001) bmu = 999.
        IF (nst > nstd) STOP 'corelevels: nst > nstd'
        DO n = 1, nst
          fac(1) = 1.0 ; fac(2) = 1.0
          READ (61,'(4i3)') nprnc(n),kappa(n),n_h(1),n_h(2)
          IF ( n_h(1) < 0 )  fac(1) = -0.5
          IF ( n_h(2) < 0 )  fac(2) = -0.5
          IF (jspins.EQ.1) THEN
            occ(n,1) = fac(1) * n_h(1) + fac(2) * n_h(2) 
          ELSE
            occ(n,1) = fac(1) * n_h(1) ; occ(n,2) = fac(2) * n_h(2)
          ENDIF
!          write(*,*) nprnc(n),kappa(n),occ(n,1), occ(n,2)
        ENDDO
        CLOSE (61)
        RETURN
      ELSE
        jspin=1
      ENDIF

      IF (z.GT.92.01e0) STOP 'setcor: z > 92'

      jz0 = z + 0.01e0
      jz = jz0
      k = 0
      nmax = 7
      DO 50 n = 1,nmax
         IF (jz.LE.0) GO TO 60
c--->    s states
         k = k + 1
         nprnc(k) = n
         kappa(k) = -1
         occ(k,1) = 2
         jz = jz - 2
         IF (jz.LE.0) GO TO 60
c--->    p states
         IF (n.EQ.1) GO TO 50
         k = k + 1
         nprnc(k) = n
         kappa(k) = 1
         occ(k,1) = 2
         jz = jz - 2
         IF (jz.LE.0) GO TO 60
         k = k + 1
         nprnc(k) = n
         kappa(k) = -2
         occ(k,1) = 4
         jz = jz - 4
         IF (jz.LE.0) GO TO 60
c--->    d functions
         iz = 0
         IF (n.EQ.3 .AND. jz0.GT.20) iz = min0(jz0-20,4)
         IF (n.EQ.4 .AND. jz0.GT.38) iz = min0(jz0-38,4)
         IF (n.EQ.4 .AND. jz0.EQ.41) iz = 4
         IF (n.EQ.5 .AND. jz0.GT.70) iz = min0(jz0-70,4)
         IF (n.EQ.5 .AND. (jz0.EQ.57.OR.jz0.EQ.64)) iz = 1
         IF (n.EQ.6 .AND. jz0.GT.88) iz = 1
         IF (n.EQ.6 .AND. jz0.EQ.90) iz = 2
         IF (iz.EQ.0) GO TO 10
         k = k + 1
         nprnc(k) = n
         kappa(k) = 2
         occ(k,1) = iz
         jz = jz - iz
         IF (iz.LT.4) GO TO 10
         IF (n.EQ.6) GO TO 50
         IF (n.EQ.4 .AND. jz0.EQ.41) GO TO 10
         IF (n.EQ.5 .AND. jz0.EQ.74) GO TO 10
         iz = 1
         IF (n.EQ.3 .AND. jz0.GT.25) iz = min0(jz0-24,6)
         IF (n.EQ.3 .AND. jz0.EQ.29) iz = 6
         IF (n.EQ.4 .AND. jz0.GT.43) iz = jz0 - 41
         IF (n.EQ.4 .AND. jz0.GT.45) iz = 6
         IF (n.EQ.5 .AND. jz0.GT.75) iz = jz0 - 74
         IF (n.EQ.5 .AND. jz0.GT.77) iz = 6
         k = k + 1
         nprnc(k) = n
         kappa(k) = -3
         occ(k,1) = iz
         jz = jz - iz
c--->    f states
   10    IF (n.NE.4) GO TO 40
c+gu  IF (jz0.LE.57) GO TO 50
         IF (jz0.LE.62) GO TO 50
         k = k + 1
         nprnc(k) = n
         kappa(k) = 3
         iz = 6
         IF (jz0.GE.62) GO TO 20
         iz = jz0 - 56
         occ(k,1) = iz
         jz = jz - iz
         GO TO 50
   20    occ(k,1) = iz
         jz = jz - iz
         iz = 8
         k = k + 1
         nprnc(k) = n
         kappa(k) = -4
         IF (jz0.GE.70) GO TO 30
         iz = jz0 - 62
         IF (jz0.EQ.64) iz = 1
         occ(k,1) = iz
         jz = jz - iz
         GO TO 50
   30    occ(k,1) = iz
         jz = jz - iz
   40    IF (n.NE.5) GO TO 50
         IF (jz0.LE.90) GO TO 50
         k = k + 1
         nprnc(k) = n
         kappa(k) = 3
         iz = jz0 - 89
         occ(k,1) = iz
         jz = jz - iz
   50 CONTINUE
   60 nst = k
      IF (k.GE.1) occ(k,1) = occ(k,1) + jz
c
c add magnetic moments
c
      IF (jspins.EQ.2) THEN
        bmu_l = bmu
        DO k = 1,nst
          occ(k,jspins) = occ(k,1)/2.0 
          occ(k,1) = occ(k,jspins)  
        ENDDO
        DO k = nst,1,-1
          fj = iabs(kappa(k)) - 0.5e0
          l = fj + 0.5e0*isign(1,kappa(k)) + 0.01e0
c polarize (d,f) only
          IF (l.GT.1.99) THEN
            IF (2*occ(k,1).GE.abs(bmu_l)) THEN
              occ(k,1) = occ(k,1) + bmu_l/2.
              occ(k,jspins) = occ(k,jspins) - bmu_l/2.
              GOTO 70
            ELSE
              IF (bmu_l.GT.0) THEN
                occ(k,1) = 2.0*occ(k,1)
                occ(k,jspins) = 0.0
                bmu_l = bmu_l - occ(k,1)
              ELSE
                occ(k,jspins) = 2.0*occ(k,jspins)
                occ(k,1) = 0.0
                bmu_l = bmu_l + occ(k,jspins)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
  70    CONTINUE
      ENDIF
c
      IF (z.EQ.65) THEN
        k_h(1) = kappa(15) ; n_h(1) = nprnc(15) ; o_h(1) = occ(15,1)
        k_h(2) = kappa(16) ; n_h(2) = nprnc(16) ; o_h(2) = occ(16,1)
        kappa(15)= kappa(17) ; nprnc(15)=nprnc(17) ; occ(15,1)=occ(17,1)
        kappa(16)= kappa(18) ; nprnc(16)=nprnc(18) ; occ(16,1)=occ(18,1)
        kappa(17)= kappa(19) ; nprnc(17)=nprnc(19) ; occ(17,1)=occ(19,1)
        kappa(18) = k_h(1) ; nprnc(18) =  n_h(1)  ; occ(18,1)= o_h(1)
        kappa(19) = k_h(2) ; nprnc(19) =  n_h(2)  ; occ(19,1)= o_h(2)

        IF (jspins.EQ.2) THEN
          o_h(1) = occ(15,jspins) ; o_h(2) = occ(16,jspins)
          occ(15,jspins) = occ(17,jspins) 
          occ(16,jspins) = occ(18,jspins) 
          occ(17,jspins) = occ(19,jspins) 
          occ(18,jspins) = o_h(1) 
          occ(19,jspins) = o_h(2) 
        ENDIF
      ENDIF

      END SUBROUTINE setcor
      END MODULE m_setcor
