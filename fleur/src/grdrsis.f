      MODULE m_grdrsis
c.....-----------------------------------------------------------------
c     evaluates gradient of an interstitial function in real space 
c  
c     based on 'rhzgrd' coded by t.asada. june,1995.
c.....-----------------------------------------------------------------
      CONTAINS
      SUBROUTINE grdrsis(
     >                   ro,bmat,xmax1,xmax2,xmax3,ndvgrd,
     <                   dro)
c.....-----------------------------------------------------------------
c input: 
c  ro(0:xmax1*xmax2*xmax3-1)  
c   any quantity stored in usual interst. box (xmax1 x xmax2 x xmax3)
c  bmat
c   bravais matrix of reciprocal space 
c  ndvgrd
c   number of ponts used when calculating derivative (3 <= ndvgrd <= 6)
c       
c output:
c  dro(0:xmax1*xmax2*xmax3-1,3)
c   gradient of ro in non-internal coordinates 
c       
c.....-----------------------------------------------------------------

      USE m_constants, ONLY : pimach
      IMPLICIT NONE 
c     ..
c     .. Scalar arguments ..
      INTEGER, INTENT (IN) :: xmax1,xmax2,xmax3,ndvgrd
c     ..
c     .. Array arguments ..
      REAL,    INTENT(IN)  :: ro(0:xmax1*xmax2*xmax3-1)
      REAL,    INTENT(IN)  :: bmat(3,3) 
c     ..
c     .. Array output ..
      REAL,    INTENT(OUT) :: dro(0:xmax1*xmax2*xmax3-1,3)  
c     ..
c     .. Locals ..
      INTEGER :: xmax(3)
      INTEGER :: direction,xyz(3),x1,x2,x3,i,ii(-3:2) 
      REAL    :: pi  
      REAL    :: dx
      REAL    :: drointern(0:xmax1*xmax2*xmax3-1,3)  

      pi = pimach() 

      xmax(1)= xmax1
      xmax(2)= xmax2
      xmax(3)= xmax3
      DO i=1,3
        IF ( xmax(i) < 3 ) THEN
          STOP 'stopped in grdrsis:  grid to small' 
        END IF
      END DO 
      IF ( (ndvgrd < 3) .or. (ndvgrd > 6) ) THEN
        STOP 'stopped in grdrsis:  ndvgrd notin [3,6]' 
      ENDIF


      DO direction=1,3 

        dx= 1./REAL(xmax(direction)) 

        DO x1=0,xmax(1)-1
          DO x2=0,xmax(2)-1
            DO x3=0,xmax(3)-1 

              DO i= -3,2
                xyz(1)= x1
                xyz(2)= x2
                xyz(3)= x3 
                xyz(direction)= xyz(direction)+i 
                ! make use of periodic boundary cond. in interstitial: 
                IF ( xyz(direction) < 0 ) THEN
                  xyz(direction)= xyz(direction)+xmax(direction)
                END IF 
                IF ( xyz(direction) >= xmax(direction) ) THEN
                  xyz(direction)= xyz(direction)-xmax(direction) 
                END IF 
                ! find coordinates in 1-dim array ro:
                ii(i)= xyz(3)*xmax(1)*xmax(2) + xyz(2)*xmax(1) + xyz(1) 
              END DO

              IF (ndvgrd.EQ.3) THEN
                drointern(ii(0),direction)=  
     &            df3( ro(ii(-1)), 
     &                 ro(ii(0)),ro(ii(1)), dx)
              ELSEIF (ndvgrd.EQ.4) THEN
                drointern(ii(0),direction)= 
     &            df4( ro(ii(-1)),
     &                 ro(ii(0)),ro(ii(1)),ro(ii(2)), dx)
              ELSEIF (ndvgrd.EQ.5) THEN
                drointern(ii(0),direction)= 
     &            df5( ro(ii(-2)),ro(ii(-1)),
     &                 ro(ii(0)),ro(ii(1)),ro(ii(2)), dx)
              ELSEIF (ndvgrd.EQ.6) THEN
                drointern(ii(0),direction)= 
     &            df6( ro(ii(-3)),ro(ii(-2)),ro(ii(-1)),
     &                 ro(ii(0)),ro(ii(1)),ro(ii(2)), dx)
              ENDIF

            END DO
          END DO
        END DO 

      END DO

 
      DO i=0,xmax(1)*xmax(2)*xmax(3)-1 

        DO direction=1,3
          dro(i,direction)=   bmat(1,direction)*drointern(i,1) 
     &                      + bmat(2,direction)*drointern(i,2)
     &                      + bmat(3,direction)*drointern(i,3)
          dro(i,direction)= dro(i,direction)/(2.*pi) 
        END DO 

      END DO     

      END SUBROUTINE grdrsis
!--------------------------------------------------------------------
! Functions: formulae for 1st deriv.:
!
      REAL FUNCTION df3(g1,f0,f1,d)             ! three point formula 
        REAL g1,f0,f1,d
        df3 = (-1*g1-0*f0+f1)/ (2*d)
      END FUNCTION df3

      REAL FUNCTION df4(g1,f0,f1,f2,d)          ! four point formula
        REAL g1,f0,f1,f2,d
        df4 = (-2*g1-3*f0+6*f1-f2)/ (6*d)
      END FUNCTION df4

      REAL FUNCTION df5(g2,g1,f0,f1,f2,d)       ! five point formula
        REAL g2,g1,f0,f1,f2,d
        df5 = (2*g2-16*g1-0*f0+16*f1-2*f2)/ (24*d)
      END FUNCTION df5

      REAL FUNCTION df6(g3,g2,g1,f0,f1,f2,d)   ! six point formula 
        REAL g3,g2,g1,f0,f1,f2,d
        df6 = (-4*g3+30*g2-120*g1+40*f0+60*f1-6*f2)/ (120*d)
      END FUNCTION df6

!----------------------------------------------------------------------
      END MODULE m_grdrsis
