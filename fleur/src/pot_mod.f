      MODULE m_potmod
      CONTAINS
      SUBROUTINE pot_mod(
     >                   jmtd,nlhd,ntypd,nmzd,nmzxyd,n2d,n3d,natd,
     >                   ntypsd,jspins,film,nq2,nq3,ntype,
     >                   neq,nvac,nmz,nmzxy,jri,ntypsy,nlh,
     X                   vr,vxy,vz,vpw,vpw_w)

      IMPLICIT NONE
      INTEGER, INTENT (IN) :: jmtd,nlhd,ntypd,nmzd,nmzxyd,n2d,n3d,natd
      LOGICAL, INTENT (IN) :: film
      INTEGER, INTENT (IN) :: nq2,nq3,jspins,ntype,nmz,nmzxy,nvac,ntypsd

      INTEGER, INTENT (IN) :: neq(ntype),jri(ntype)
      INTEGER, INTENT (IN) :: ntypsy(natd),nlh(ntypsd)
      REAL,    INTENT (INOUT) :: vr(jmtd,0:nlhd,ntypd,jspins)
      REAL,    INTENT (INOUT) :: vz(nmzd,2,jspins)
      COMPLEX, INTENT (INOUT) :: vpw(n3d,jspins),vpw_w(n3d,jspins)
      COMPLEX, INTENT (INOUT) :: vxy(nmzxyd,n2d-1,2,jspins)
      
      INTEGER i,j,n,ivac,nat,typmag,bxcflag,bxc_r,bxc_c

! --- modify mag.-pot.  >> ---------------------------------------------

      ! what do you want ?
      bxcflag= 0
      typmag= ntype-1
      !  0 : potential is not changed
      !  1 : B_xc is kept in MTs and set to zero elsewhere
      !  2 : B_xc is kept in first typmag MTs and set to zero elsewhere
      !  3 : B_xc is read from file
      !  4 : B_xc is written to file

      IF (bxcflag/=0) THEN

        IF ( (bxcflag<0) .or. (bxcflag>5) ) THEN
          STOP 'Stop in vgen:  bxcflag out of bounds'
        ENDIF
        IF ( (bxcflag==2) .and. ((typmag<0).or.(typmag>ntype)) ) THEN
          STOP 'Stop in vgen:  typmag out of bounds'
        ENDIF
        IF (jspins/=2) THEN
          STOP 'Stop in vgen:  no B-field as jspins==1'
        ENDIF

        IF (bxcflag/=4) THEN

          IF (bxcflag==3) THEN
            OPEN(201,file='bxc',form='unformatted',action='read')
            DO j= 1,ntype+2
              READ(201)
            ENDDO
          ENDIF
          IF (bxcflag/=1) THEN
            nat= 1
            DO n= 1,ntype
              IF ( (bxcflag==3) .or. (n>typmag) ) THEN
                DO j= 0,nlh(ntypsy(nat))
                  DO i= 1,jri(n)
                    vr(i,j,n,1)= ( vr(i,j,n,1)+vr(i,j,n,2) )/2.
                    vr(i,j,n,2)= vr(i,j,n,1)
                    IF (bxcflag==3) THEN
                      READ(201) bxc_r
                      vr(i,j,n,1)= vr(i,j,n,1) + bxc_r
                      vr(i,j,n,2)= vr(i,j,n,2) - bxc_r
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF
              nat= nat + neq(n)
            ENDDO
          ENDIF
          DO j= 1,nq3
            vpw(j,1)= ( vpw(j,1)+vpw(j,2) )/2.
            vpw(j,2)= vpw(j,1)
            vpw_w(j,1)= ( vpw_w(j,1)+vpw_w(j,2) )/2.
            vpw_w(j,2)= vpw_w(j,1)
            IF (bxcflag==3) THEN
              READ(201) bxc_c
              vpw(j,1)= vpw(j,1) + bxc_c
              vpw(j,2)= vpw(j,2) - bxc_c
              READ(201) bxc_c
              vpw_w(j,1)= vpw_w(j,1) + bxc_c
              vpw_w(j,2)= vpw_w(j,2) - bxc_c
            ENDIF
          ENDDO
          IF (film) THEN
            DO ivac= 1,nvac
              DO i= 1,nmz
                vz(i,ivac,1)= ( vz(i,ivac,1)+vz(i,ivac,2) )/2.
                vz(i,ivac,2)= vz(i,ivac,1)
                IF (bxcflag==3) THEN
                  READ(201) bxc_r
                  vz(i,ivac,1)= vz(i,ivac,1) + bxc_r
                  vz(i,ivac,2)= vz(i,ivac,2) - bxc_r
                ENDIF
              ENDDO
              DO i= 1,nmzxy
                DO j= 1,nq2-1
                  vxy(i,j,ivac,1)=
     &             ( vxy(i,j,ivac,1)+vxy(i,j,ivac,2) )/2.
                  vxy(i,j,ivac,2)= vxy(i,j,ivac,1)
                  IF (bxcflag==3) THEN
                    READ(201) bxc_c
                    vxy(i,j,ivac,1)= vxy(i,j,ivac,1) + bxc_c
                    vxy(i,j,ivac,2)= vxy(i,j,ivac,2) - bxc_c
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          WRITE(6,fmt='(1x)')
          SELECT CASE (bxcflag)
            CASE(1)
              WRITE(6,fmt='(A)') 'B_xc outside MTs is set to zero !!'
            CASE(2)
              WRITE(6,fmt='(A,i3,1x,A)')
     &         'B_xc outside the first',typmag,'MTs is set to zero !!'
            CASE(3)
              CLOSE(201)
              WRITE(6,fmt='(A)') 'B_xc is read from file "bxc" !!'
          END SELECT
          WRITE(6,fmt='(1x)')

        ELSE ! (bxc==4), write B_xc it file

          OPEN(201,file='bxc',form='unformatted',status='replace')
          WRITE(201) nq3, ntype, film
          nat= 1
          DO n= 1,ntype
            WRITE(201) nlh(ntypsy(nat)), jri(n), neq(n)
            nat= nat + neq(n)
          ENDDO
          IF (film) THEN
            WRITE(201) nvac, nmz, nmzxy, nq2
          ELSE
            WRITE(201) 0,1,1,1
          ENDIF
          nat= 1
          DO n= 1,ntype
            DO j= 0,nlh(ntypsy(nat))
              DO i= 1,jri(n)
                WRITE(201) ( vr(i,j,n,1)-vr(i,j,n,2) )/2.
              ENDDO
            ENDDO
            nat= nat + neq(n)
          ENDDO
          DO j= 1,nq3
            WRITE(201) ( vpw(j,1)-vpw(j,2) )/2.
            WRITE(201) ( vpw_w(j,1)-vpw_w(j,2) )/2.
          ENDDO
          IF (film) THEN
            DO ivac= 1,nvac
              DO i= 1,nmz
                WRITE(201) ( vz(i,ivac,1)-vz(i,ivac,2) )/2.
              ENDDO
              DO i= 1,nmzxy
                DO j= 1,nq2-1
                  WRITE(201) ( vxy(i,j,ivac,1)-vxy(i,j,ivac,2) )/2.
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          CLOSE(201)
          STOP 'Stop:  B_xc is written to "bxc" '

        ENDIF

      ENDIF ! bxcflag/=0


      END SUBROUTINE pot_mod
      END MODULE m_potmod
