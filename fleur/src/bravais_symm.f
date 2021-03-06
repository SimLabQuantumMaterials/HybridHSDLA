      MODULE m_bravaissymm
!********************************************************************
!     determines the point group of the bravais lattice given the
!     lattice vectors. the idea is to determine all the lattice
!     vectors that have the same length as a_{1,2,3}, and then use
!     these to determine the possible rotation matrices.
!     these rotation matrices are in lattice coordinates.    mw 12-99
!********************************************************************
      CONTAINS
      SUBROUTINE bravais_symm(
     >                        as,bs,scale,nop48,neig12,
     <                        nops,mrot)

      IMPLICIT NONE

!==> Arguments
      INTEGER, INTENT (IN) :: nop48  ! max. number of operations (typically 48)
      INTEGER, INTENT (IN) :: neig12 ! max. number of lattice vectors with same length
                                     ! (max occurs for close-packed fcc: 12)
      REAL,    INTENT (IN) :: as(3,3),bs(3,3),scale(3) ! lattice information
      INTEGER, INTENT(OUT) :: nops, mrot(3,3,nop48)    ! point group operations

!==> Locals
      REAL    amet(3,3),b1,b2,b3,d1,d2,d3,dmax,dt
      INTEGER i,k,k1,k2,k3,m1,m2,m3,n1,n2,n3
      INTEGER irot(3,3)
      INTEGER lv1(3,neig12),lv2(3,neig12),lv3(3,neig12)

      REAL, PARAMETER :: eps=1.0e-7

!---> set up metric for distances
      amet = 0.0
      DO i = 1,3
         amet(1,1) = amet(1,1) + (scale(i)**2)*as(i,1)*as(i,1)
         amet(2,2) = amet(2,2) + (scale(i)**2)*as(i,2)*as(i,2)
         amet(3,3) = amet(3,3) + (scale(i)**2)*as(i,3)*as(i,3)
         amet(2,1) = amet(2,1) + (scale(i)**2)*as(i,1)*as(i,2)
         amet(3,2) = amet(3,2) + (scale(i)**2)*as(i,2)*as(i,3)
         amet(3,1) = amet(3,1) + (scale(i)**2)*as(i,3)*as(i,1)
      ENDDO

      amet(1,2) = amet(2,1)
      amet(2,3) = amet(3,2)
      amet(1,3) = amet(3,1)

!---> distances for the lattice vectors
      d1 = amet(1,1)
      d2 = amet(2,2)
      d3 = amet(3,3)
      b1 = ( bs(1,1)/scale(1) )**2 + ( bs(1,2)/scale(2) )**2             &
     &   + ( bs(1,3)/scale(3) )**2
      b2 = ( bs(2,1)/scale(1) )**2 + ( bs(2,2)/scale(2) )**2             &
     &   + ( bs(2,3)/scale(3) )**2
      b3 = ( bs(3,1)/scale(1) )**2 + ( bs(3,2)/scale(2) )**2             &
     &   + ( bs(3,3)/scale(3) )**2

!---> determine the cutoffs along each direction a_i:
      dmax = max( d1,d2,d3)

      m1 = nint( dmax * b1 )
      m2 = nint( dmax * b2 )
      m3 = nint( dmax * b3 )

!---->loop over all possible lattice vectors to find those with the
!---->length, i.e., ones that could be rotations
      n1 = 1
      n2 = 1
      n3 = 1

      lv1(1:3,1) = (/ 1,0,0 /)
      lv2(1:3,1) = (/ 0,1,0 /)
      lv3(1:3,1) = (/ 0,0,1 /)

      DO k3=-m3,m3
         DO k2=-m2,m2
            DO k1=-m1,m1

               dt = distance2(k1,k2,k3)

!---->    check if the same length
               IF ( abs( dt - d1 ) < eps ) THEN
                  IF (.not.( k1==1 .and. k2==0 .and. k3==0 ) ) THEN
                     n1 = n1+1
                     if(n1.gt.neig12) stop 'n1>neig12'
                     lv1(1,n1) = k1
                     lv1(2,n1) = k2
                     lv1(3,n1) = k3
                  ENDIF
               ENDIF

               IF ( abs( dt - d2 ) < eps ) THEN
                  IF (.not.( k1==0 .and. k2==1 .and. k3==0 ) ) THEN
                     n2 = n2+1
                     if(n2.gt.neig12) stop 'n2>neig12'
                     lv2(1,n2) = k1
                     lv2(2,n2) = k2
                     lv2(3,n2) = k3
                  ENDIF
               ENDIF

               IF ( abs( dt - d3 ) < eps ) THEN
                  IF (.not.( k1==0 .and. k2==0 .and. k3==1 ) ) THEN
                     n3 = n3+1
                     if(n3.gt.neig12) stop 'n3>neig12'
                     lv3(1,n3) = k1
                     lv3(2,n3) = k2
                     lv3(3,n3) = k3
                  ENDIF
               ENDIF

            ENDDO
         ENDDO
      ENDDO

!---> the possible rotation matrices are given by the matrix of
!---> column vectors of lv_{1,2,3}
      nops = 0
      DO k3 = 1,n3
         DO k2 = 1,n2
            DO k1 = 1,n1

!--->          check whether determinant is +/-1 (needs to be for rotation)
               IF ( abs(mdet(k1,k2,k3)) .NE. 1 ) CYCLE

!--->          check whether this maintains lengths correctly
!--->          if M is the metric, then must have R^T M R = M 
               irot = reshape( (/ lv1(:,k1),lv2(:,k2),lv3(:,k3) /) ,
     &                         (/ 3,3 /) )
               IF ( any( abs(
     &           matmul( transpose(irot), matmul(amet,irot) ) - amet 
     &                     ) > eps ) ) CYCLE

               nops = nops + 1
               IF ( nops > nop48 ) STOP 'nop > nop48'
               mrot(:,:,nops) = irot

            ENDDO
         ENDDO
      ENDDO

      WRITE (6,'(//," Point group of the Bravais lattice has ",i2,
     &        " operations")') nops

      RETURN

      CONTAINS   ! INTERNAL routines

      REAL FUNCTION distance2(l1,l2,l3)
!*********************************************************************
!     calculates the magnitude square for a vector (l1,l2,l3) given in
!     lattice units
!*********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l1,l2,l3

      distance2 = l1*(l1*amet(1,1) + 2*l2*amet(2,1))
     &          + l2*(l2*amet(2,2) + 2*l3*amet(3,2))
     &          + l3*(l3*amet(3,3) + 2*l1*amet(1,3))

      RETURN
      END FUNCTION distance2

      INTEGER FUNCTION mdet(k1,k2,k3)
!*********************************************************************
!     determines the determinant for possible rotation matrix
!     ( lv1(:,k1) ; lv2(:,k2) ; lv3(:,k3) )
!*********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: k1,k2,k3

      mdet = lv1(1,k1)*( lv2(2,k2)*lv3(3,k3) - lv2(3,k2)*lv3(2,k3) )
     &     + lv1(2,k1)*( lv2(3,k2)*lv3(1,k3) - lv2(1,k2)*lv3(3,k3) )
     &     + lv1(3,k1)*( lv2(1,k2)*lv3(2,k3) - lv2(2,k2)*lv3(1,k3) )

      RETURN
      END FUNCTION mdet

      END SUBROUTINE bravais_symm
      END MODULE m_bravaissymm
