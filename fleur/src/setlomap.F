      SUBROUTINE setlomap(
     >                    ntypd,nlod,llod,
     >                    ntyp,nlo,llo,
     <                    lo1l,nlol,l_dulo,ulo_der)
c***********************************************************************
c sets up nlol and lo1l
c 
c nlo     : number of local orbitals for each atom type
c llo     : the l quantum numbers of the local orbitals
c nlol    : the of local orbitals with a certain l (for each atom type)
c lo1l    : the number of the first local orbital with that l
c l_dulo  : if .true., this is a local orbital formed with a $\dot u$
c
c p.kurz jul. 1996
c***********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,nlod,llod,ntyp
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)    :: nlo(ntypd)
      INTEGER, INTENT (INOUT) :: llo(nlod,ntypd)
      INTEGER, INTENT (OUT)   :: lo1l(0:llod,ntypd),nlol(0:llod,ntypd)
      INTEGER, INTENT (OUT)   :: ulo_der(nlod,ntypd)
      LOGICAL, INTENT (OUT)   :: l_dulo(nlod,ntypd)
C     ..
C     .. Local Scalars ..
      INTEGER ilo,l
C     ..
      WRITE (6,FMT=8000) nlo(ntyp), (llo(ilo,ntyp),ilo=1,nlo(ntyp))
 8000 FORMAT ('the number of local orbitals for this atom type is: ',i2,
     +       /,'the corresponding values of l are: ',20i2)
      IF (nlo(ntyp).GT.nlod) STOP 'setlomap: nlo > nlod!!!'

      DO l = 0,llod
         nlol(l,ntyp) = 0
         lo1l(l,ntyp) = 0
      END DO
      DO ilo = 1,nlod
        l_dulo(ilo,ntyp) = .false.
        ulo_der(ilo,ntyp)= 0
      ENDDO
      l = -1

      DO ilo = 1,nlo(ntyp)
c+gu
         IF (llo(ilo,ntyp).LT.0) THEN
#ifdef CPP_APW 
           l_dulo(ilo,ntyp) = .true.
#else
           l_dulo(ilo,ntyp) = .false.
           ulo_der(ilo,ntyp)= -llo(ilo,ntyp)/10+1
           llo(ilo,ntyp) = llo(ilo,ntyp)+10*(ulo_der(ilo,ntyp)-1)
#endif
           llo(ilo,ntyp) = abs( llo(ilo,ntyp) ) - 1
           WRITE(6,'(A,I2,A,I2)') 'I use',ulo_der(ilo,ntyp),
     &       '. derivative of l =',llo(ilo,ntyp)
         ELSE
           l_dulo(ilo,ntyp) = .false.
         ENDIF
c-gu
         IF (llo(ilo,ntyp).GT.llod) STOP 'setlomap: l > llod!!!'
         IF (llo(ilo,ntyp).LT.l) THEN
            WRITE (6,FMT=*)
     +        'setlomap: please specify the l quantum numbers ',
     +        'of the local orbitals is ascending order.'
            STOP 'setlomap'
         END IF
         IF (llo(ilo,ntyp).GT.l) THEN
            l = llo(ilo,ntyp)
            lo1l(l,ntyp) = ilo
         END IF
         nlol(l,ntyp) = nlol(l,ntyp) + 1
      END DO

      RETURN
      END
