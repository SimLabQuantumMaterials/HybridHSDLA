      MODULE m_bandstr1
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      CONTAINS
      SUBROUTINE bandstr1(
     >                    idsyst,idtype,bmat,nkpt,name,jspins,film)

      USE m_cotra,  ONLY : cotra3
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: idsyst,idtype,nkpt,jspins
      REAL,    INTENT (IN)  :: bmat(3,3)
      LOGICAL, INTENT (IN)  :: film
      CHARACTER(len=8), INTENT (IN) :: name(10)
      REAL,    POINTER              :: kpts(:,:)      

      REAL, POINTER             :: syp(:,:) 
      CHARACTER(len=1), POINTER :: ssy(:)
      REAL,    ALLOCATABLE      :: rsyp(:,:),del(:),d(:)
      INTEGER, ALLOCATABLE      :: nk(:)
      INTEGER nosyp,i,j,n,ntot
      REAL dk(3),syp1(3),rsyp1(3)
      
      CALL get_points(
     >                idsyst,idtype,
     <                nosyp,syp,ssy)
!
! --> transform to cartesian units
!
      ALLOCATE ( rsyp(3,nosyp),del(nosyp),d(nosyp),nk(nosyp) ) 
      DO i = 1,nosyp
        syp1(:) = syp(:,i)
        CALL cotra3(syp1,rsyp1,bmat)
        rsyp(:,i) = rsyp1(:)
      ENDDO
!
! --> calculate length between points
!
      d(1) = 0.0
      DO i = 2,nosyp
        del(i) = ( rsyp(1,i) - rsyp(1,i-1) )**2 +
     +           ( rsyp(2,i) - rsyp(2,i-1) )**2 +
     +           ( rsyp(3,i) - rsyp(3,i-1) )**2
        del(i) = sqrt(del(i))
        d(i) = d(i-1) + del(i)
      ENDDO
!
! --> distibute points evenly
!
      ntot = nosyp
      DO i = 2,nosyp
        nk(i) = NINT ( (nkpt-nosyp)*( del(i) / d(nosyp) ) )
        ntot = ntot + nk(i)
      ENDDO 
      ALLOCATE ( kpts(3,ntot) )
!
! --> generate k-points mesh 
!
      n = 1
      DO i = 2,nosyp
         dk(1) = ( syp(1,i) - syp(1,i-1) ) / (nk(i)+1)
         dk(2) = ( syp(2,i) - syp(2,i-1) ) / (nk(i)+1)
         dk(3) = ( syp(3,i) - syp(3,i-1) ) / (nk(i)+1)
         kpts(:,n) = syp(:,i-1)
         n = n + 1
         DO j = 1, nk(i)
            kpts(:,n) = kpts(:,n-1) + dk(:)
            n = n + 1
         ENDDO
      ENDDO
      kpts(:,n) = syp(:,nosyp)

      OPEN (41,file='kpts',form='formatted',status='new')
      IF (.NOT.film) THEN
        WRITE(41,'(i5,f20.10)') ntot,1.0
        DO n = 1,ntot
           WRITE (41,'(4f10.5)') kpts(:,n),1.0
        ENDDO
      ELSE
        WRITE(41,'(i5,f20.10,3x,l1)') ntot,1.0,.false.
        DO n = 1,ntot
           WRITE (41,'(3f10.5)') kpts(1:2,n),1.0
        ENDDO
      END IF
      CLOSE (41)

      CALL write_gnu(
     >               nosyp,d,ssy,name,jspins)

      DEALLOCATE ( rsyp,syp,del,nk,ssy,d )

      END SUBROUTINE bandstr1
!----------------------------------------------------------------------
! once the file "bands.1" and "bands.2" are created, activate with:
! gnuplot < band.gnu > band.ps
!----------------------------------------------------------------------
      SUBROUTINE write_gnu(
     >                     nosyp,d,ssy,name,jspins)
!
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: nosyp,jspins
      REAL,    INTENT (IN) :: d(nosyp)
      CHARACTER(len=1), INTENT (IN) :: ssy(nosyp)
      CHARACTER(len=8), INTENT (IN) :: name(10)
      
      INTEGER n,aoff,adel
      aoff = iachar('a')-1
      adel = iachar('a')-iachar('A')
      write(*,*) aoff,adel 

      OPEN (27,file='band.gnu',status='unknown')
      WRITE (27,900)
      WRITE (27,901)
      WRITE (27,902)
      WRITE (27,903)
      WRITE (27,904) name(:)
      DO n = 1, nosyp
        WRITE (27,905) d(n),d(n)
      ENDDO
      WRITE (27,906) d(1),d(nosyp)
!
! nomal labels
!
      IF (iachar(ssy(1)) < aoff ) THEN
        WRITE (27,907) ssy(1),d(1)
      ELSE
        WRITE (27,907) " ",d(1)
      ENDIF
      DO n = 2, nosyp-1
        IF (iachar(ssy(n)) < aoff ) THEN 
          WRITE (27,908) ssy(n),d(n)
        ELSE
          WRITE (27,908) " ",d(n)
        ENDIF
      ENDDO
      IF (iachar(ssy(nosyp)) < aoff ) THEN
        WRITE (27,909) ssy(nosyp),d(nosyp)
      ELSE
        WRITE (27,909) " ",d(nosyp)
      ENDIF
!
! greek labels
!
      DO n = 1, nosyp
        IF (iachar(ssy(n)) > aoff ) THEN
          WRITE (27,914) achar(iachar(ssy(n))-adel),d(n)
        ENDIF
      ENDDO
!
! now write the rest
!
      WRITE (27,910) 
      WRITE (27,911) d(nosyp)+0.00001
      IF (jspins == 2) WRITE (27,912) "'","'"
      WRITE (27,913) "'","'"
      CLOSE (27)

 900  FORMAT ('set terminal postscript enhanced "Times-Roman" 20')
 901  FORMAT ('set xlabel ""')
 902  FORMAT ('set ylabel "E - E_F (eV)"')
 903  FORMAT ('set nokey')
 904  FORMAT ('set title "',10a8,'"')
 905  FORMAT ('set arrow from',f9.5,', -9.0 to',f9.5,',  5.0 nohead')
 906  FORMAT ('set arrow from',f9.5,', 0.0 to',f9.5,', 0.0 nohead lt 3')
#ifdef CPP_AIX
 907  FORMAT ('set xtics ("',a1,'"',f9.5,', \\')
 908  FORMAT ('           "',a1,'"',f9.5,', \\')
#else
 907  FORMAT ('set xtics ("',a1,'"',f9.5,', \')
 908  FORMAT ('           "',a1,'"',f9.5,', \')
#endif
 909  FORMAT ('           "',a1,'"',f9.5,'  )')
 910  FORMAT ('set ytics -8,2,4')
#ifdef CPP_AIX
 911  FORMAT ('plot [0:',f9.5,'] [-9:5] \\')
 912  FORMAT ('"< awk ',a1,'{print  $1,$2+0.00}',a1,
     +        ' bands.2" w p pt 12 ps 0.5, \\')
#else
 911  FORMAT ('plot [0:',f9.5,'] [-9:5] \')
 912  FORMAT ('"< awk ',a1,'{print  $1,$2+0.00}',a1,
     +        ' bands.2" w p pt 12 ps 0.5, \')
#endif
 913  FORMAT ('"< awk ',a1,'{print  $1,$2+0.00}',a1,
     +        ' bands.1" w p pt  7 ps 0.5')
 914  FORMAT ('set label "',a1,'" at ',f9.5,
     +        ', -9.65 center font "Symbol,20"')

      END SUBROUTINE write_gnu
!----------------------------------------------------------------------
! given a bravais-lattice, determine <nosyp> symmetry points (syp)
! and their names (lowercase = greek)
!----------------------------------------------------------------------
      SUBROUTINE get_points(
     >                      idsyst,idtype,
     <                      nosyp,syp,ssy)

      IMPLICIT NONE
      INTEGER, INTENT (IN)  :: idsyst,idtype
      INTEGER, INTENT (OUT) :: nosyp
      REAL, POINTER         :: syp(:,:)  ! actually intent out
      CHARACTER(len=1), POINTER :: ssy(:)

      LOGICAL         :: band_inp_file
      INTEGER         :: n
      REAL, PARAMETER :: f12 = 1./2., f14 = 1./4., zro = 0.0
      REAL, PARAMETER :: f34 = 3./4., f38 = 3./8., one = 1.0
      REAL, PARAMETER :: f13 = 1./3., f23 = 2./3.

      !Check if band_inp file exists for lines
      INQUIRE(file ="band_inp",exist= band_inp_file)
      IF (band_inp_file) THEN
         OPEN(99,file ="band_inp",status ="old")
         nosyp=0
         !count the lines
         DO
            READ(99,*,END = 100)
            nosyp=nosyp+1
         ENDDO
 100     REWIND 99
         ALLOCATE(syp(3,nosyp),ssy(nosyp) )
         DO n = 1,nosyp
            READ(99,*,err = 110,END=110) ssy(n),syp(:,n)
         ENDDO
         CLOSE(99)
         RETURN
 110     WRITE(*,*) "Error reading band_inp file"
         STOP "Bandstr1"
      ENDIF

      !No band_inp file -> determine default lines for bandstructure
      write(*,*) idsyst,idtype
      nosyp = -1
      IF ( (idsyst == 1).AND.(idtype ==  3) ) THEN       ! fcc
        nosyp = 7
        ALLOCATE ( syp(3,nosyp),ssy(nosyp) )
        syp(:,1) = (/f12,f12,one/)  ; ssy(1) = "X"
        syp(:,2) = (/f38,f38,f34/)  ; ssy(2) = "K"
        syp(:,3) = (/zro,zro,zro/)  ; ssy(3) = "g"
        syp(:,4) = (/f12,f12,f12/)  ; ssy(4) = "L"
        syp(:,5) = (/f12,f14,f34/)  ; ssy(5) = "W"
        syp(:,6) = (/f12,zro,f12/)  ; ssy(6) = "X"
        syp(:,7) = (/zro,zro,zro/)  ; ssy(7) = "g"
      ENDIF
      IF ( (idsyst == 5).AND.(idtype ==  1) ) THEN       ! rhombohedric (trigonal)
        nosyp = 8
        ALLOCATE ( syp(3,nosyp),ssy(nosyp) )
        syp(:,1) = (/f12,f12, f12/)  ; ssy(1) = "Z"
        syp(:,2) = (/zro,zro, zro/)  ; ssy(2) = "g"
        syp(:,3) = (/f14,f14,-f14/)  ; ssy(3) = "K"
        syp(:,4) = (/f12,f12,-f12/)  ; ssy(4) = "Z"
        syp(:,5) = (/f14,f12,-f14/)  ; ssy(5) = "W"
        syp(:,6) = (/zro,f12, zro/)  ; ssy(6) = "L"
        syp(:,7) = (/zro,zro, zro/)  ; ssy(7) = "g"
        syp(:,8) = (/f12,f12, zro/)  ; ssy(8) = "F"
      ENDIF
      IF ( (idsyst == 4).AND.(idtype ==  1) ) THEN       ! hexagonal (angle = 60)
        nosyp = 8
        ALLOCATE ( syp(3,nosyp),ssy(nosyp) )
        syp(:,1) = (/zro,zro, zro/)  ; ssy(1) = "g"
        syp(:,2) = (/f12,f12, zro/)  ; ssy(2) = "M"
        syp(:,3) = (/f13,f23, zro/)  ; ssy(3) = "K"
        syp(:,4) = (/zro,zro, zro/)  ; ssy(4) = "g"
        syp(:,5) = (/zro,zro, f12/)  ; ssy(5) = "A"
        syp(:,6) = (/f12,f12, f12/)  ; ssy(6) = "L"
        syp(:,7) = (/f13,f23, f12/)  ; ssy(7) = "H"
        syp(:,8) = (/zro,zro, f12/)  ; ssy(8) = "A"
      ENDIF
      IF ( (idsyst == 1).AND.(idtype ==  1) ) THEN       ! simple cubic
        nosyp = 8
        ALLOCATE ( syp(3,nosyp),ssy(nosyp) )
        syp(:,1) = (/f12,f12, zro/)  ; ssy(1) = "M"
        syp(:,2) = (/zro,zro, zro/)  ; ssy(2) = "g"
        syp(:,3) = (/f12,zro, zro/)  ; ssy(3) = "X"
        syp(:,4) = (/f12,f12, zro/)  ; ssy(4) = "M"
        syp(:,5) = (/f12,f12, f12/)  ; ssy(5) = "R"
        syp(:,6) = (/f12,zro, zro/)  ; ssy(6) = "X"
        syp(:,7) = (/zro,zro, zro/)  ; ssy(7) = "g"
        syp(:,8) = (/f12,f12, f12/)  ; ssy(8) = "R"
      ENDIF
      IF ( (idsyst == 1).AND.(idtype ==  2) ) THEN       ! body centered cubic
        nosyp = 6
        ALLOCATE ( syp(3,nosyp),ssy(nosyp) )
        syp(:,1) = (/zro,zro, zro/)  ; ssy(1) = "g"
        syp(:,2) = (/f12,-f12,f12/)  ; ssy(2) = "H"
        syp(:,3) = (/zro,zro, f12/)  ; ssy(3) = "N"
        syp(:,4) = (/f14,f14, f14/)  ; ssy(4) = "P"
        syp(:,5) = (/zro,zro, zro/)  ; ssy(5) = "g"
        syp(:,6) = (/zro,zro, f12/)  ; ssy(6) = "N"
      ENDIF
      IF ( (idsyst == 2).AND.(idtype ==  2) ) THEN       ! body centered tetragonal (a > c)
        nosyp = 8
        ALLOCATE ( syp(3,nosyp),ssy(nosyp) )
        syp(:,1) = (/f12,f12,-f12/)  ; ssy(1) = "Z"    ! via Lambda and V
        syp(:,2) = (/zro,zro, zro/)  ; ssy(2) = "g"    ! via Sigma
        syp(:,3) = (/-f12,f12,f12/)  ; ssy(3) = "Z"    ! via Y
        syp(:,4) = (/zro,zro, f12/)  ; ssy(4) = "X"    ! via Delta
        syp(:,5) = (/zro,zro, zro/)  ; ssy(5) = "g"    
        syp(:,6) = (/zro,f12, zro/)  ; ssy(6) = "N"    ! via Q
        syp(:,7) = (/f14,f14, f14/)  ; ssy(7) = "P"    ! via W
        syp(:,8) = (/zro,zro, f12/)  ; ssy(8) = "X"
      ENDIF
      IF ( (idsyst == 2).AND.(idtype ==  2) ) THEN       ! body centered tetragonal (a < c)
        nosyp = 9
        ALLOCATE ( syp(3,nosyp),ssy(nosyp) )
        syp(:,1) = (/-f12,f12,f12/)  ; ssy(1) = "Z"    ! via F and Sigma
        syp(:,2) = (/zro,zro, zro/)  ; ssy(2) = "g"    ! via Delta
        syp(:,3) = (/zro,zro, f12/)  ; ssy(3) = "X"    ! via W
        syp(:,4) = (/f14,f14, f14/)  ; ssy(4) = "P"    ! via Q
        syp(:,5) = (/zro,f12, zro/)  ; ssy(5) = "N"     
        syp(:,6) = (/zro,zro, zro/)  ; ssy(6) = "g"    ! via Lambda
        syp(:,7) = (/f12,f12,-f12/)  ; ssy(7) = "Z"    ! via U and Y
        syp(:,8) = (/f12,f12, zro/)  ; ssy(8) = "X"
        syp(:,9) = (/f14,f14, f14/)  ; ssy(9) = "P"
      ENDIF
      IF ( (idsyst == 2).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
        nosyp = 8
        ALLOCATE ( syp(3,nosyp),ssy(nosyp) )
        syp(:,1) = (/zro,zro, zro/)  ; ssy(1) = "g"    ! via Delta
        syp(:,2) = (/f12,zro, zro/)  ; ssy(2) = "X"    ! via Y
        syp(:,3) = (/f12,f12, zro/)  ; ssy(3) = "M"    ! via Sigma
        syp(:,4) = (/zro,zro, zro/)  ; ssy(4) = "g"    ! via Lambda
        syp(:,5) = (/zro,zro, f12/)  ; ssy(5) = "Z"    ! via U
        syp(:,6) = (/f12,zro, f12/)  ; ssy(6) = "R"    ! via T
        syp(:,7) = (/f12,f12, f12/)  ; ssy(7) = "A"    ! via S
        syp(:,8) = (/zro,zro, f12/)  ; ssy(8) = "Z"
      ENDIF
      IF ( (idsyst == 3).AND.(idtype ==  1) ) THEN       ! primitive tetragonal (a < c)
        nosyp = 10
        ALLOCATE ( syp(3,nosyp),ssy(nosyp) )
        syp(:,1) = (/zro,zro, zro/)  ; ssy(1) = "g"    ! via Sigma
        syp(:,2) = (/f12,zro, zro/)  ; ssy(2) = "X"    ! via D
        syp(:,3) = (/f12,f12, zro/)  ; ssy(3) = "S"    ! via C
        syp(:,4) = (/zro,f12, zro/)  ; ssy(4) = "Y"    ! via Delta
        syp(:,5) = (/zro,zro, zro/)  ; ssy(5) = "g"    ! via Lambda
        syp(:,6) = (/zro,zro, f12/)  ; ssy(6) = "Z"    ! via A
        syp(:,7) = (/f12,zro, f12/)  ; ssy(7) = "U"    ! via P
        syp(:,8) = (/f12,f12, f12/)  ; ssy(8) = "R"    ! via E
        syp(:,9) = (/zro,f12, f12/)  ; ssy(8) = "T"    ! via B
        syp(:,10)= (/zro,zro, f12/)  ; ssy(8) = "Z"
      ENDIF

      IF ( nosyp == -1 ) STOP 'bandstr1: lattice not implemented'

      END SUBROUTINE get_points
!----------------------------------------------------------------------
      END MODULE m_bandstr1
