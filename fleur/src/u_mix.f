      MODULE m_umix
!
! mix the old and new density matrix for the lda+U method
!                                                 gb.2001
      CONTAINS
      SUBROUTINE u_mix(
     >                 n_u,jspins,n_mmp_new)

      USE m_nmat_rot
!
! ... Arguments
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: n_u,jspins
      COMPLEX, INTENT (INOUT) :: n_mmp_new(-3:3,-3:3,n_u,jspins)
!
! ... Locals ...
      INTEGER n,j,k,iofl
      REAL gam,del,sum1,sum2,mix_u
      REAL alpha,spinf,theta,phi
      LOGICAL n_exist
      COMPLEX,ALLOCATABLE :: n_mmp(:,:,:,:),n_mmp_old(:,:,:,:)
!
! check for possible rotation of n_mmp
!
      INQUIRE (file='n_mmp_rot',exist=n_exist)
      IF (n_exist) THEN
         OPEN (68,file='n_mmp_rot',status='old',form='formatted')
         READ(68,*) theta,phi
         CLOSE (68)
        CALL nmat_rot(0.0,-theta,-phi,3,n_u,jspins,
     X                 n_mmp_new)
      ENDIF
!
! check for LDA+U and open density-matrix - file
!
      INQUIRE (file='n_mmp_mat',exist=n_exist)
      OPEN (69,file='n_mmp_mat',status='unknown',form='formatted')
      

      IF (n_exist) THEN
         ALLOCATE (  n_mmp_old(-3:3,-3:3,n_u,jspins) )
         ALLOCATE (      n_mmp(-3:3,-3:3,n_u,jspins) )
         READ (69,9000) n_mmp_old(:,:,:,:)

         READ (69,'(2(6x,f5.3))',IOSTAT=iofl) alpha,spinf
         IF ( iofl == 0 ) THEN
!
! mix here straight with given mixing factors 
!
          REWIND (69)
          sum1 = 0.0
          IF (jspins.EQ.1) THEN
             DO  n = 1,n_u
               DO j = -3,3
                  DO k = -3,3
                     sum1 = sum1 + abs( n_mmp_new(k,j,n,1) - 
     +                                  n_mmp_old(k,j,n,1) )
                     n_mmp(k,j,n,1) = alpha  * n_mmp_new(k,j,n,1) +
     +                          (1. - alpha) * n_mmp_old(k,j,n,1)
                  ENDDO
               ENDDO
             ENDDO
             WRITE (6,'(a16,f12.6)') 'n_mmp distance =',sum1
           ELSE
             sum2 = 0.0
             gam = 0.5 * alpha * ( 1.0 + spinf )
             del = 0.5 * alpha * ( 1.0 - spinf )
             DO  n = 1,n_u
               DO j = -3,3
                  DO k = -3,3
                     sum1 = sum1 + abs( n_mmp_new(k,j,n,1) - 
     +                                  n_mmp_old(k,j,n,1) )
                     sum2 = sum2 + abs( n_mmp_new(k,j,n,2) - 
     +                                  n_mmp_old(k,j,n,2) )
                     n_mmp(k,j,n,1) = gam  * n_mmp_new(k,j,n,1) +
     +                          (1. - gam) * n_mmp_old(k,j,n,1) +
     +                                del  * n_mmp_new(k,j,n,2) -
     +                                del  * n_mmp_old(k,j,n,2)
                     n_mmp(k,j,n,2) = gam  * n_mmp_new(k,j,n,2) +
     +                          (1. - gam) * n_mmp_old(k,j,n,2) +
     +                                del  * n_mmp_new(k,j,n,1) -
     +                                del  * n_mmp_old(k,j,n,1)
                  ENDDO
               ENDDO
             ENDDO
             WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 1 =',sum1
             WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 2 =',sum2
           ENDIF
           WRITE (69,9000) n_mmp
           WRITE (69,'(2(a6,f5.3))') 'alpha=',alpha,'spinf=',spinf

         ELSEIF (iofl > 0 ) THEN
!
! read error ; stop
!
           WRITE (6,*) 'ERROR READING mixing factors in n_mmp_mat'
           WRITE (6,'(2(a6,f5.3))') 'alpha=',alpha,'spinf=',spinf
           STOP 'ERROR READING n_mmp_mat'
         ELSE
!
! calculate distance and write new n_mmp to mix in broyden.F
!
           sum1 = 0.0
           DO  n = 1,n_u
             DO j = -3,3
                DO k = -3,3
                   sum1 = sum1 + abs( n_mmp_new(k,j,n,1) -
     +                                n_mmp_old(k,j,n,1) )
                ENDDO
             ENDDO
           ENDDO
           IF (jspins.EQ.1) THEN
             WRITE (6,'(a16,f12.6)') 'n_mmp distance =',sum1
           ELSE
             sum2 = 0.0
             WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 1 =',sum1
             DO  n = 1,n_u
               DO j = -3,3
                  DO k = -3,3
                     sum2 = sum2 + abs( n_mmp_new(k,j,n,2) -
     +                                  n_mmp_old(k,j,n,2) )
                  ENDDO
               ENDDO
             ENDDO
             do j=-3,3
              write(6,'(14f12.6)') (n_mmp_old(k,j,1,2),k=-3,3)
             enddo
             WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 2 =',sum2
             do j=-3,3
              write(6,'(14f12.6)') (n_mmp_new(k,j,1,2),k=-3,3)
             enddo
           ENDIF
           REWIND(69)
           WRITE (69,9000) n_mmp_old
           WRITE (69,9000) n_mmp_new
         ENDIF !  iofl == 0 

         DEALLOCATE ( n_mmp_old,n_mmp )
      ELSE
!
! first time with lda+u; write new n_mmp  
!
         WRITE (69,9000) n_mmp_new
         WRITE (69,'(2(a6,f5.3))') 'alpha=',0.05,'spinf=',1.0
      ENDIF
      
 9000 FORMAT(7f20.13)

      CLOSE (69)
      END SUBROUTINE u_mix
      END MODULE m_umix 
