      MODULE m_chkmt
!---------------------------------------------------------------------
!  Check muffin tin radii and determine a reasonable choice for MTRs.
!  Derive also other parameters for the input file, to provide some
!  help in the out-file.                                        gb`02
!---------------------------------------------------------------------
      CONTAINS
      SUBROUTINE chkmt(
     >                 natd,ntype,neq,film,pos,dvac,rmt,amat,z,
     >                 l_gga,noel,l_test,odd,jspins,
     <                 nel,kmax,dtild,dvac1,
     <                 nlo,llo,ncst,lmax,jri,rmt1,dx)

      USE m_od_types, ONLY : od_dim
      USE m_enpara,   ONLY : w_enpara
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: natd,ntype,jspins
      REAL,    INTENT (IN) :: dvac
      LOGICAL, INTENT (IN) :: film,l_gga,l_test
      INTEGER, INTENT (OUT):: nel
      REAL,    INTENT (OUT):: kmax,dtild,dvac1
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntype),z(ntype)
      REAL,    INTENT (IN) :: pos(3,natd),rmt(ntype),amat(3,3)
      INTEGER, INTENT (OUT):: nlo(ntype),llo(2,ntype)
      INTEGER, INTENT (OUT):: ncst(ntype),lmax(ntype),jri(ntype)
      REAL,    INTENT (OUT):: rmt1(ntype),dx(ntype)
      CHARACTER*3, INTENT (IN) :: noel(ntype)
c-odim
      TYPE (od_dim), INTENT (IN) :: odd
c+odim
C     ..
C     .. Local Scalars ..
      INTEGER na,n,nna,nn,n1,nn1,k1,k2,k3
      INTEGER i,j,jri1,lmax1,ncst2,locore
      REAL    sum,sss,xmin,dx1,rkm,fac,sum_r,fac_1,fac_2
      LOGICAL error
C     ..
C     .. Local Arrays ..
      INTEGER minni(2),ncst1(0:103),nce(0:24)
      INTEGER nval(4,ntype),noflo(2,ntype),skiplo(ntype)
      REAL    el0(0:3,ntype),ello0(2,ntype),evac0(2)
      LOGICAL lchange(0:3,ntype),llochg(2,ntype)
      REAL    dist(ntype,ntype),dist1(ntype,ntype),t_rmt(0:103)
      
!
! number of core levels for each element; the INT(ncst1/100) number
! provides information about possible local orbitals: 1...(s,p)-LO
! 2...p-LO and 3...d-LO
!
      ncst1 =(/0,0,                                                0,  ! Va,H,He
     +     01, 01,                                  1, 1, 1, 1, 1, 1,  ! Li - Ne
     +     04, 04,                                  4, 4, 4, 4, 4, 4,  ! Na - Ar
     +    107,107,207,207, 7, 7, 7, 7, 7, 7, 7, 7,309, 9, 9, 9, 9, 9,  ! K - Kr
     +    112,112,212,212,12,12,12,12,12,12,12,12,314,14,14,14,14,14,  ! Rb - Xe
     +    117,117,217,217,17,17,17,17,17,17,17,17, 17,17,17,17,17,     ! Cs - Lu
     +                219,19,19,19,19,19,19,19,19,319,21,21,21,21,21,  ! Hf - Rn
     +    124,124,224,224,24,24,24,24,24,24,24,24, 24,24,24,24,24/)    ! Fr - Lw
!
! electrons associated with a given mumber of core-levels
!
      nce(0) = 0   ; nce(1) = 2   ; nce(4) = 10  ; nce(7) = 18 
      nce(9) = 28  ; nce(12) = 36 ; nce(14) = 46 ; nce(17) = 54
      nce(19) = 68 ; nce(21) = 78 ; nce(24) = 86
!
! typical muffin-tin radii
!
      t_rmt(0:103) = 2.3 ! default value
      t_rmt(1) = 1.0 ; t_rmt(5:9) = 1.3 ; t_rmt(16:17) = 1.8

      error=.false.
      dist(:,:) = 9.99e19
      na=0
      DO n=1,ntype
        DO n1=1,neq(n)
          na=na+1
!
! check distance to other atoms:
!
          nna=0
          DO nn=1,ntype
            DO nn1=1,neq(nn)
              nna=nna+1

              sss=9.9e19
              DO k1=-3,3
                DO k2=-3,3
                  DO k3=-3,3

                    IF (.NOT.((nn.eq.n).AND.(nn1.eq.n1).AND.(k1.eq.0)
     +                 .AND.(k2.eq.0).AND.(k3.eq.0))) THEN
                      sum=0.0e0
                      DO i=1,3
                        sum=sum+(pos(i,na)-pos(i,nna)+
     +                           k1*amat(i,1)+k2*amat(i,2)+
     +                           k3*amat(i,3))**2
                      ENDDO
                      IF ( sum.lt.sss ) sss=sum
                    ENDIF

                  ENDDO
                ENDDO
              ENDDO
              dist(n,nn) = min( dist(n,nn),sqrt(sss) ) 
              IF ( sss.LE.(rmt(nn)+rmt(n))**2 ) THEN
                 error=.true.
                 IF (l_test)
     +           WRITE(6,240) nn,nn1,(pos(i,nna),i=1,3),rmt(nn),
     +                        n ,n1 ,(pos(i,na),i=1,3),rmt(n )
              ENDIF

            ENDDO ! nn1
          ENDDO   ! nn
!
! distance to vacuum
!
          IF (film) THEN
             IF (odd%d1) THEN
                IF ((sqrt(pos(1,na)**2+pos(2,na)**2)+
     +               rmt(n)).GT.dvac/2.) THEN
                   error=.true.
                   WRITE(6,241) n ,n1
                   WRITE(6,*) sqrt(pos(1,na)**2+pos(2,na)**2),
     &                  rmt(n),dvac/2.
                END IF
             ELSE
                IF ( ( (pos(3,na)+rmt(n) ).GT. dvac/2.).OR.
     +               ( (pos(3,na)-rmt(n) ).LT.-dvac/2.) ) THEN
                   error=.true.
                   WRITE(6,241) n ,n1
                   WRITE(6,*) pos(3,na),rmt(n),dvac/2.
                ENDIF
             ENDIF
          END IF
        ENDDO
      ENDDO
       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
!     DO n = 1, ntype
!       WRITE (*,'(12f12.6)') dist(:,n)
!     ENDDO
      dist1 = dist
      rmt1(:) = 999.
      WRITE (6,*) '----------------------------------------------------'
      WRITE (6,*) 'Suggested values for input: '
      WRITE (6,*) 

      IF (film) THEN
        fac = 0.95
      ELSE
        fac = 0.975
      ENDIF

      minni = minloc(dist)       ! minni(1) and minni(2) are the indices of the closest atoms
      xmin  = minval(dist)       ! xmin is their distance

      DO WHILE ( xmin < 999.0 )

        sum_r = 1.0 / ( t_rmt(z(minni(1))) + t_rmt(z(minni(2))) )
        fac_1 = t_rmt(z(minni(1))) * sum_r
        fac_2 = t_rmt(z(minni(2))) * sum_r

        IF (rmt1(minni(1)) > 990.) THEN         ! if not set, determine MTR 
          IF (rmt1(minni(2)) > 990.) THEN       ! both not set, choose in between
            rmt1(minni(1)) = fac * xmin * fac_1 ! / 2
            rmt1(minni(2)) = fac * xmin * fac_2 ! / 2
          ELSE
            rmt1(minni(1)) = fac * ( xmin - rmt1(minni(2)) )
            IF (2*rmt1(minni(1)).GT.dist(minni(1),minni(1))) THEN
              rmt1(minni(1)) = fac * dist(minni(1),minni(1)) / 2
            ENDIF
          ENDIF
        ELSEIF (rmt1(minni(2)) > 990.) THEN
          rmt1(minni(2)) = fac * ( xmin - rmt1(minni(1)) )
          IF (2*rmt1(minni(2)).GT.dist(minni(2),minni(2))) THEN
            rmt1(minni(2)) = fac * dist(minni(2),minni(2)) / 2
          ENDIF
        ENDIF

        dist(minni(1),minni(1)) = 999.0
        dist(minni(2),minni(1)) = 999.0
        dist(minni(1),minni(2)) = 999.0
        dist(minni(2),minni(2)) = 999.0

        DO j = 1, 2
          DO n = 1, ntype
            IF (z(n) == z(minni(j))) THEN
              IF (rmt1(n) > 990.) THEN 
                rmt1(n) = rmt1(minni(j))
              ENDIF
            ENDIF
          ENDDO
        ENDDO

        minni = minloc(dist)
        xmin  = minval(dist)

      ENDDO

      dvac1 = 0.0
      rkm = 0.0
      na = 0
      nel = 0
      WRITE (6,230)
      DO n= 1,ntype
!
!--> determine M.T. radii 
!
        DO j= 1,ntype
           dist(j,n) = rmt1(n)+rmt1(j)
           IF ( dist1(j,n)-dist(j,n) < 0.0 ) THEN
             WRITE(*,*) j,n,dist1(j,n)-dist(j,n)
             rmt1(n) = fac * dist1(j,n) / 2.
             rmt1(j) = fac * dist1(j,n) / 2.
           ENDIF
        ENDDO
        IF (film) THEN
          DO nn = 1, neq(n)
            na = na + 1
            IF (odd%d1) THEN
               dvac1 = max( dvac1, sqrt(pos(1,na)**2+pos(2,na)**2)
     +              +rmt1(n) )
            ELSE
               dvac1 = max( dvac1, pos(3,na)+rmt1(n) )
            END IF
          ENDDO
        ENDIF 
!
!--> calculate jri, dx and lmax
!
        IF (l_gga) THEN
          jri1 = nint( 240*rmt1(n) ) 
        ELSE
          jri1 = nint( 160*rmt1(n) ) 
        ENDIF
        jri1 = nint( jri1*0.5 ) * 2 + 1
        dx1 =  log(3200*z(n)*rmt1(n))/(jri1-1)
        IF (rmt1(n).LT.1.8) THEN
          lmax1 = 6
        ELSEIF (rmt1(n).LT.2.4) THEN
          lmax1 = 8
        ELSEIF (rmt1(n).LT.3.0) THEN
          lmax1 = 10
        ELSEIF (rmt1(n).LT.3.2) THEN
          lmax1 = 12 
        ELSE
          WRITE (6,'("Atom Nr.",i3,"( ",a3,") has a M.T. radius of",
     +                                     f8.4)') n,noel(n),rmt1(n)
          WRITE (6,'("that was truncated to 3.0")')
          rmt1(n) = 3.0
          lmax1 = 12 
        ENDIF
        rkm = max( rkm , lmax1/rmt1(n) )
!
!--> determine core levels
!
        ncst2 = ncst1( z(n) )
        IF (ncst2.GT.300) THEN     ! should add d-LO
          ncst2 = ncst2 - 300 ; locore = 10
          nlo(n) = 1 ; llo(1,n) = 2
        ELSEIF (ncst2.GT.200) THEN ! should add p-LO
          ncst2 = ncst2 - 200 ; locore = 6
          nlo(n) = 1 ; llo(1,n) = 1
        ELSEIF (ncst2.GT.100) THEN ! should add (s,p)-LO
          ncst2 = ncst2 - 100 ; locore = 8
          nlo(n) = 2 ; llo(1,n) = 0 ; llo(2,n) = 1
        ELSE
          nlo(n) = 0 ; locore = 0
        ENDIF
        nel = nel + ( z(n) - nce(ncst2) + locore ) * neq(n)
        IF ((locore == 6).OR.(locore == 10)) ncst2 = ncst2 - 2
        IF (locore == 8) ncst2 = ncst2 - 3
           
        WRITE (6,9070) noel(n),z(n),ncst2,lmax1,jri1,rmt1(n),dx1
 9070   FORMAT (a3,i3,3i5,2f10.6)
        WRITE (6,9090) neq(n),.false.,nlo(n),(llo(i,n),i=1,nlo(n))
 9090   FORMAT (i2,',force =',l1,',nlo=',i2,',llo=',20i2)
        dx(n) = dx1 ; ncst(n) = ncst2 ; lmax(n) = lmax1 ; jri(n) = jri1
!
!--> determine valence states
!
        IF ( z(n) < 3 ) THEN
          nval(:,n) = (/1,2,3,4/)
        ELSEIF ( z(n) < 11 ) THEN
          nval(:,n) = (/2,2,3,4/)
        ELSEIF ( z(n) < 19 ) THEN
          nval(:,n) = (/3,3,3,4/)
        ELSEIF ( z(n) < 31 ) THEN
          nval(:,n) = (/4,4,3,4/)
        ELSEIF ( z(n) < 37 ) THEN
          nval(:,n) = (/4,4,4,4/)
        ELSEIF ( z(n) < 49 ) THEN
          nval(:,n) = (/5,5,4,4/)
        ELSEIF ( z(n) < 55 ) THEN
          nval(:,n) = (/5,5,5,4/)
        ELSEIF ( z(n) < 72 ) THEN
          nval(:,n) = (/6,6,5,4/)
        ELSEIF ( z(n) < 81 ) THEN
          nval(:,n) = (/6,6,5,5/)
        ELSEIF ( z(n) < 87 ) THEN
          nval(:,n) = (/6,6,6,5/)
        ELSE
          nval(:,n) = (/7,7,6,5/)
        ENDIF
        
        skiplo(n) = 0
        DO i = 1, nlo(n)
          noflo(i,n) = nval(llo(i,n)+1,n) - 1 
          skiplo(n) = skiplo(n) + (2*llo(i,n) + 1)
        ENDDO

      ENDDO ! loop over atom types
  230 FORMAT ('Atom Z  ncst lmax jri    rmt         dx')

      IF (film) THEN
        dvac1 = 2* (dvac1+0.3)
        dtild = dvac1 + 1.5 * maxval( rmt1(:) )
        WRITE (6,'("vacuum distance dvac =",f10.5)') dvac1
        WRITE (6,'("extra vac.dist. dtild=",f10.5)') dtild
      ENDIF
      WRITE (6,'("k_max =",f8.5)') rkm
      WRITE (6,'("G_max =",f8.5)') 3*rkm
      WRITE (6,'("Valence Electrons =",i5)') nel
      kmax = rkm

      WRITE (6,*) '----------------------------------------------------'
 
      IF ( error.AND.l_test ) STOP 'Error checking M.T. radii'
  240 FORMAT('   error in muffin tin radii  , pair  ',2(/,2i5,4f10.5))
  241 FORMAT('   error: atom ',i3,' # ',i3,'reaches out into vaccuum')

      IF (.not.l_test) THEN
        OPEN (40,file='enpara',form='formatted',status='unknown')
        DO j = 1, jspins
          ello0 = real(noflo)
          el0 = real(nval)
          evac0(:) = -0.1
          lchange = .false.
          llochg  = .false.

          DO n = 1,ntype ! chek for hydrogen, where this mode does not work
            IF (z(n) == 1) THEN
              IF (film) THEN
                el0(0:3,n) = -0.1
              ELSE
                el0(0:3,n) =  0.1
              ENDIF
              lchange(0:3,n) = .true.
            ENDIF
          ENDDO

          CALL w_enpara(
     >                   3,2,ntype,1,j,film,nlo,
     >                   skiplo,ello0,el0,evac0,lchange,
     >                   llochg,.true.,1.0,93)
        ENDDO
        CLOSE (40)
      ENDIF

      END SUBROUTINE chkmt
      END MODULE m_chkmt
