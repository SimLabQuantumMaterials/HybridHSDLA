      MODULE m_od_types

! types for 1D calculations

        TYPE od_dim
          LOGICAL :: d1
          INTEGER :: mb,M,k3,m_cyl
          INTEGER :: chi,rot
          LOGICAL :: invs,zrfs
          INTEGER :: n2d,nq2,nn2d
          INTEGER :: kimax2
          INTEGER :: nop,nat      
        END TYPE od_dim

        TYPE od_inp
          LOGICAL :: d1
          INTEGER :: mb,M,k3,m_cyl
          INTEGER :: chi,rot
          LOGICAL :: invs,zrfs 
          INTEGER :: n2d,nq2,nn2d
          INTEGER :: kimax2
          INTEGER, POINTER :: ig(:,:)  !(-k3:k3,-M:M)
          INTEGER, POINTER :: kv(:,:)        !(2,n2d) 
          INTEGER, POINTER :: nst2(:)        !(n2d)
        END TYPE od_inp

        TYPE od_sym
          INTEGER :: nop,nat
          INTEGER, POINTER :: ngopr(:)     !(nat)
          REAL   , POINTER :: mrot(:,:,:)  !(3,3,nop)
	  REAL   , POINTER :: tau(:,:)     !(3,nop) 
          INTEGER, POINTER :: invtab(:)    !(nop)
          INTEGER, POINTER :: multab(:,:)  !(nop,nop)
        END TYPE od_sym 

        TYPE od_lda
          INTEGER :: nn2d
          INTEGER, POINTER :: igf(:,:)  !(0:nn2d-1,2)
          REAL   , POINTER :: pgf(:)    !(0:nn2d-1)
        END TYPE od_lda

        TYPE od_gga
          INTEGER          :: nn2d
          REAL, POINTER    :: pgfx(:)  ! (0:nn2d-1)
          REAL, POINTER    :: pgfy(:)
          REAL, POINTER    :: pgfxx(:)
          REAL, POINTER    :: pgfyy(:)
          REAL, POINTER    :: pgfxy(:)
        END TYPE od_gga

      END MODULE m_od_types
