      MODULE m_nmat_rot

! Calculate the Wigner rotation matrices for complex spherical
! harmonics for all space-group rotations and l=1,2,3. Needed 
! for the calculation of the density matrix in nmat.
!                 gb 2002

      CONTAINS
      SUBROUTINE nmat_rot(
     >                    alpha,beta,gamma,l_in,n_u,jspins,
     X                    n_mmp)

      USE m_constants, ONLY : pimach
      USE m_inv3

      IMPLICIT NONE

! .. arguments:
      INTEGER, INTENT(IN)  :: l_in,n_u,jspins
      REAL,    INTENT(IN)  :: alpha,beta,gamma
      COMPLEX, INTENT(INOUT) :: n_mmp(-3:3,-3:3,n_u,jspins)

! .. local variables:
      INTEGER ns,signum,ispin,n
      INTEGER i,j,k,l,m,mp,x_lo,x_up,x,e_c,e_s
      REAL fac_l_m,fac_l_mp,fac_lmpx,fac_lmx,fac_x,fac_xmpm
      REAL pi,co_bh,si_bh,zaehler,nenner,cp,sp
      REAL sina,sinb,sinc,cosa,cosb,cosc,determ,dt
      COMPLEX ci,phase_g,phase_a,bas,d(-l_in:l_in,-l_in:l_in)
      COMPLEX n_tmp(-l_in:l_in,-l_in:l_in)
      COMPLEX nr_tmp(-l_in:l_in,-l_in:l_in)
      LOGICAL, SAVE :: written = .false.

      REAL dmat(3,3),dmati(3,3)

      INTRINSIC sqrt,max,min

      ci = cmplx(0.0,1.0)
      pi = pimach()

        co_bh = cos(beta*0.5)
        si_bh = sin(beta*0.5)

      DO l = 1, l_in

        DO m = -l,l
          fac_l_m = fac(l+m) * fac(l-m)
          phase_g = exp( - ci * gamma * m )

          DO mp = -l,l
            fac_l_mp = fac(l+mp) * fac(l-mp)

            zaehler = sqrt( real(fac_l_m * fac_l_mp) )
            phase_a = exp( - ci * alpha * mp ) 
            x_lo = max(0, m-mp)
            x_up = min(l-mp, l+m)

            bas = zaehler * phase_a * phase_g 
            d(m,mp) = cmplx(0.0,0.0)
            DO x = x_lo,x_up
              fac_lmpx = fac(l-mp-x)
              fac_lmx  = fac(l+m-x)
              fac_x    = fac(x)
              fac_xmpm = fac(x+mp-m)
              nenner = fac_lmpx * fac_lmx * fac_x * fac_xmpm
              e_c = 2*l + m - mp - 2*x 
              e_s = 2*x + mp - m
              IF (e_c.EQ.0) THEN
                cp = 1.0
              ELSE
                cp = co_bh ** e_c
              ENDIF
              IF (e_s.EQ.0) THEN
                sp = 1.0
              ELSE
                sp = si_bh ** e_s
              ENDIF
              d(m,mp) = d(m,mp) + bas * (-1)**x * cp * sp / nenner
            ENDDO

          ENDDO ! loop over mp
        ENDDO   ! loop over m

      ENDDO

      DO ispin = 1, jspins
        DO n = 1, n_u
           n_tmp(:,:) = n_mmp(-l_in:l_in,-l_in:l_in,n,ispin)

           nr_tmp = matmul( transpose( conjg(d) ) , n_tmp)
           n_tmp =  matmul( nr_tmp, d )

           n_mmp(-l_in:l_in,-l_in:l_in,n,ispin) = n_tmp(:,:)
         ENDDO
      ENDDO
      write(*,'(3f10.5)') alpha,beta,gamma
      write(*,'(14f8.4)') n_mmp

      END SUBROUTINE nmat_rot

      ELEMENTAL REAL FUNCTION  fac(n)

      INTEGER, INTENT (IN) :: n
      INTEGER :: i
 
      fac = 0
      IF (n.LT.0) RETURN
      fac = 1
      IF (n.EQ.0) RETURN
      DO i = 2,n
        fac = fac * i
      ENDDO

      END FUNCTION  fac
      
      END MODULE m_nmat_rot
