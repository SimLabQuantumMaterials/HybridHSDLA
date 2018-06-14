      MODULE m_xyavden
      CONTAINS
      SUBROUTINE xy_av_den(
     >                    n3d,k3d,nq3,nmzd,nmz,dvac,delz,
     >                    area,ig2,kv3,amat,psq,rht)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n3d,k3d,nq3,nmzd,nmz
      REAL,    INTENT(IN) :: dvac,delz,area
      INTEGER, INTENT(IN) :: ig2(n3d),kv3(3,n3d)
      REAL,    INTENT(IN) :: amat(3,3),rht(nmzd,2)
      COMPLEX, INTENT(IN) :: psq(n3d)

      INTEGER  ivfft,i,j,k
      REAL     ani,z
      REAL,    ALLOCATABLE :: af1(:),bf1(:)
      EXTERNAL cfft

      ivfft =  3*k3d
      ALLOCATE (af1(ivfft),bf1(ivfft))

      af1(:) = 0.0 ; bf1(:) = 0.0
      DO i = 1, nq3
        IF (ig2(i) == 1) THEN
          k = kv3(3,i)
          IF ( k < 0 ) THEN
            k = ivfft + k + 1 
          ELSE
            k = k + 1
          ENDIF
          af1(k) = real(psq(i))
          bf1(k) = aimag(psq(i))
        ENDIF
      ENDDO

      CALL cfft(af1,bf1,ivfft,ivfft,ivfft,+1)

      OPEN(77,file='qws',status='unknown')
      j = 1
      k = 3 - 2*j
      DO i = nmz,1,-1
        z = (dvac/2 + (i-1)*delz) * k 
        WRITE(77,'(2f20.10)') z,rht(i,j)*area
      ENDDO
      ani = 1.0/real(ivfft)
      j = 0
      DO i = 0,ivfft - 1
        j = j + 1
        z = amat(3,3)*i*ani
        IF (z > amat(3,3)/2) z = z - amat(3,3)
        IF ( abs(z) < dvac/2 ) THEN
          WRITE(77,'(2f20.10)') z,af1(j)*area
        ENDIF
      ENDDO
      j = 1
      k = 3 - 2*j
      DO i = 1, nmz
        z = (dvac/2 + (i-1)*delz) * k 
        WRITE(77,'(2f20.10)') z,rht(i,j)*area
      ENDDO
      
      CLOSE(77)
      DEALLOCATE (af1,bf1)
      STOP

      END SUBROUTINE xy_av_den
      END MODULE m_xyavden

