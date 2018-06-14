      MODULE m_mrot2su
      CONTAINS
      SUBROUTINE mrot2su(mrot_k,amat,su)
      USE m_grp_k, ONLY : euler

      INTEGER, INTENT(IN)  :: mrot_k(:,:,:)
      REAL,    INTENT(IN)  :: amat(3,3)
      COMPLEX, INTENT(OUT) :: su(:,:,:)


      INTEGER n,nopk
      REAL    alpha,beta,gamma,cb,sb
      COMPLEX eia

      nopk = SIZE(mrot_k,3)
      DO n = 1, nopk
        CALL euler(mrot_k(:,:,n),amat,alpha,beta,gamma)
        cb = cos(beta/2) ; sb = sin(beta/2)
        eia = exp( cmplx( 0.0 , alpha ) )
        su(1,1,n) =  conjg(eia)*cb 
        su(2,1,n) = -conjg(eia)*sb
        su(1,2,n) = eia*sb 
        su(2,2,n) = eia*cb
      ENDDO


      END SUBROUTINE mrot2su
      END MODULE m_mrot2su
