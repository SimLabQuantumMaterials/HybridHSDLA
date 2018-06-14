      MODULE m_types
!
! Types for orbital moment calculation:
!
      TYPE t_orb                           ! 'normal' contributions
        SEQUENCE
        REAL :: uu,dd                      ! z   component
        COMPLEX :: uup,uum,ddp,ddm         ! +/- component
      END TYPE t_orb

      TYPE t_orbl                          ! local orbitals & (u,d)
        SEQUENCE
        REAL :: uulo,dulo
        COMPLEX :: uulop,uulom,dulop,dulom
      END TYPE t_orbl

      TYPE t_orblo                         ! lo,lo' contributions
        SEQUENCE
        REAL :: z
        COMPLEX :: p,m
      END TYPE t_orblo
!
! Types for spin-off-diagonal charge density:
!
      TYPE t_mt21                          ! 'normal' contributions
        SEQUENCE
        REAL ::  uun,udn,dun,ddn           ! normes of radial overlaps
        COMPLEX :: uu,ud,du,dd             ! values
      END TYPE t_mt21

      TYPE t_lo21                          ! ocal orbitals & (u,d)
        SEQUENCE
        REAL ::  uulon,dulon,uloun,ulodn   ! normes of radial overlaps
        COMPLEX :: uulo,dulo,ulou,ulod     ! values
      END TYPE t_lo21
!
! Type for LDA+U:
!
      TYPE t_utype
        SEQUENCE
        REAL u,j
        INTEGER l
        LOGICAL :: l_amf
      END TYPE t_utype

      END MODULE m_types
