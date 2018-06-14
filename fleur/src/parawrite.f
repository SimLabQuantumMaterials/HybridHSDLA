      SUBROUTINE parawrite(
     >                     nop,k1d,k2d,k3d,n3d,n2d,nlod,llod,
     >                     kq1d,kq2d,kq3d,kxc1d,kxc2d,kxc3d,
     >                     ntype,nat,jmtd,ntypsd,nlhd,memd,lmaxd,
     >                     jspins,nvac,nvd,nv2d,nwdd,msh,
     >                     nstd,nkptd,nobd,ncvd,layerd,odd)

      USE m_od_types, ONLY : od_dim
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: k1d,k2d,k3d,n3d,n2d,nlod,llod
      INTEGER, INTENT(IN) :: kq1d,kq2d,kq3d,kxc1d,kxc2d,kxc3d
      INTEGER, INTENT(IN) :: ntype,nat,jmtd,ntypsd,lmaxd
      INTEGER, INTENT(IN) :: jspins,nvac,nvd,nv2d,nwdd
      INTEGER, INTENT(IN) :: nkptd,nobd,ncvd,msh,nstd
      INTEGER, INTENT(IN) :: nop,nlhd,memd,layerd
c-odim
      TYPE (od_dim), INTENT(IN) :: odd
c+odim
   
      INTEGER nmz,nmzxy
c
c--->    write parameter deck
c
      OPEN (3,file='fl7para',form='formatted',status='unknown')
      REWIND 3
c
      WRITE (3,'(6x,''Symmetry elements and FFT parameters '')')
      WRITE (3,8080) nop,k1d,k2d,k3d,n3d,n2d
 8080 FORMAT (6x,'parameter (nop =',i2,',k1d=',i3,',k2d=',i3,',k3d=',i3,
     +       ',n3d=',i6,',n2d=',i4,')')

c+sb(cdn_fft;Feb.97)
      WRITE (3,'(6x,''FFT-parameters for charge density'')')
      WRITE (3,8090) kq1d,kq2d,kq3d
 8090 FORMAT (6x,'parameter (kq1d=',i3,',kq2d=',i3,',kq3d=',i3,')')
c-sb(cdn_fft;Feb.97)
c+sb(xc_fft;Oct.97)
      WRITE (3,'(6x,''FFT-parameters for XC-potential'')')
      WRITE (3,8100) kxc1d,kxc2d,kxc3d
 8100 FORMAT (6x,'parameter (kxc1d=',i3,',kxc2d=',i3,',kxc3d=',i3,')')
c-sb(xc_fft;Oct.97)

      WRITE (3,'(6x,''(Inequivalent) atoms and radial mesh'')')
      WRITE (3,8110) ntype,nat,jmtd
 8110 FORMAT (6x,'parameter (ntypd=',i3,',natd=',i3,',jmtd=',i4,')')

c
      WRITE (3,'(6x,''Different lattice harmonics components'')')
      WRITE (3,8120) ntypsd,nlhd,memd
 8120 FORMAT (6x,'parameter (ntypsd=',i3,',nlhd=',i3,',memd=',i2,')')

      WRITE (3,'(6x,''L-cutoff of potential, charge & wf '')')
      WRITE (3,8130) lmaxd
 8130 FORMAT (6x,'parameter (lmaxd=',i2,')')

      WRITE (3,'(6x,''Number of spins and vacua'')')
      WRITE (3,8140) jspins,nvac
 8140 FORMAT (6x,'parameter (jspd=',i1,',nvacd=',i1,')')

      nmz = 250
      nmzxy = 100
      WRITE (3,'(6x,''Vacuum layers for G=0 and G=/=0'')')
      WRITE (3,8150) nmz,nmzxy
 8150 FORMAT (6x,'parameter (nmzd=',i3,',nmzxyd=',i3,')')

c+gu
      WRITE (3,'(6x,''3 & 2D planewaves, windows, k-points'')')
      WRITE (3,8180) nvd,nv2d,nwdd,nkptd
 8180 FORMAT (6x,'parameter (nvd=',i5,',nv2d=',i4,',nwdd=',
     +       i1,',nkptd=',i5,')')

      WRITE (3,'(6x,''Number of (occupied) bands'')')
      WRITE (3,8190) nobd,nobd
 8190 FORMAT (6x,'parameter (nobd=',i4,',neigd=',i4,')')
c-gu

      WRITE (3,'(6x,''Nuclear mesh and core levels'')')
      WRITE (3,8200) msh,nstd
 8200 FORMAT(6x,'parameter (msh=',i4,',nstd=',i2,')')

      WRITE (3,'(6x,''Max. l-value for pseudocharge exp.'')')
      WRITE (3,8210) ncvd
 8210 FORMAT (6x,'parameter (ncvd=',i3,')')

      WRITE (3,'(6x,''Layers for vacuum DOS'')')
      WRITE (3,'(6x,''parameter(layerd='',i3,'')'')') layerd

      WRITE (3,'(6x,''Local Orbital Parameters'')')
      WRITE (3,8220) nlod,llod
 8220 FORMAT(6x,'parameter (nlod=',i3,',llod=',i3,')')

      IF (odd%d1) THEN
        WRITE (3,'(6x,''One-dimensional parameters'')')
        WRITE (3,8230) odd%mb,odd%M,odd%m_cyl,odd%chi,odd%rot,
     &                           odd%nop,odd%n2d,odd%d1
      ELSE
        WRITE (3,'(6x,''One-dimensional parameters'')')
        WRITE (3,8230) 0,0,0,0,0,nop,n2d,.false.
      END IF
 8230   FORMAT (6x,'parameter (vM=',i3,',MM=',i3,',m_cyl=',i3,
     &                       ',chi=',i3,
     &                       ',rot=',i3,',nop=',i3,',n2d=',i6,
     &                       ',d1=',l1,')')

      CLOSE (3)
      RETURN

      END
