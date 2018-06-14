      SUBROUTINE m_perp(
     >     jspd,ntypd,jmtd,itype,l_relax,l_constr,mix_b,
     >     jri,dx,rmsh,vr0,volmts,chmom,qa21,alphdiff,
     X     alph,beta,b_con)
c***********************************************************************
c calculates the perpendicular part of the local moment.
c if l_relax is true the angle of the output local moment is calculated
c and mixed with the input angles using mix_b as the mixing parameter
c if l_constr is true the output constraint b-field is calculated and
c mixed with the input contraint field using mix_b
c Philipp Kurz 2000-02-09
c***********************************************************************

      USE m_intgr, ONLY : intgr3
      USE m_constants, ONLY : pimach
      USE m_polangle
      USE m_rotdenmat
      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspd,ntypd,jmtd
      INTEGER, INTENT (IN) :: itype
      LOGICAL, INTENT (IN) :: l_constr
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: jri(ntypd)
      LOGICAL, INTENT (IN) :: l_relax(ntypd)
      REAL, INTENT    (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL, INTENT    (IN) :: vr0(jmtd,ntypd,jspd),volmts(ntypd)
      REAL, INTENT    (IN) :: mix_b
      REAL, INTENT    (IN) :: chmom(ntypd,jspd)
      REAL, INTENT    (IN) :: alphdiff(ntypd) 
      REAL, INTENT (INOUT) :: alph(ntypd),beta(ntypd),b_con(2,ntypd)
      COMPLEX, INTENT (IN) :: qa21(ntypd)
C     ..
C     .. Local Scalars ..
      INTEGER iri
      REAL b_xavh,scale,b_con_outx,b_con_outy,mx,my,mz,
     +     alphh,betah,fpi,mz_tmp,mx_mix,my_mix,mz_mix
      REAL    rho11,rho22
      COMPLEX rho21
C     ..
C     .. Local Arrays ..
      REAL b_xc_h(jmtd),b_xav(ntypd)

      fpi = 4*pimach()

c---> calculated the comp. of the local moment vector
      mx = 2*real(qa21(itype))
      my = 2*aimag(qa21(itype))
      mz = chmom(itype,1) - chmom(itype,2)
      WRITE  (6,8025) mx,my
      WRITE (16,8025) mx,my
c---> determine the polar angles of the moment vector in the local frame
      CALL pol_angle(mx,my,mz,betah,alphh)
      WRITE  (6,8026) betah,alphh
      WRITE (16,8026) betah,alphh
 8025 FORMAT(2x,'--> local frame: ','mx=',f9.5,' my=',f9.5)
 8026 FORMAT(2x,'-->',10x,' delta beta=',f9.5,
     +                   '  delta alpha=',f9.5)

      IF (l_relax(itype)) THEN
c--->    rotate the (total (integrated) density matrix to obtain
c--->    it in the global spin coordinate frame
         rho11 = chmom(itype,1)
         rho22 = chmom(itype,2)
         rho21 = qa21(itype)
         CALL rot_den_mat(alph(itype),beta(itype),
     X                       rho11,rho22,rho21)
c--->    determine the polar angles of the mom. vec. in the global frame
         mx = 2*real(rho21)
         my = 2*aimag(rho21)
         mz = rho11 - rho22
         CALL pol_angle(mx,my,mz,betah,alphh)
         WRITE  (6,8027) beta(itype),alph(itype)-alphdiff(itype)
         WRITE (16,8027) beta(itype),alph(itype)-alphdiff(itype)
         WRITE  (6,8028) betah,alphh-alphdiff(itype)
         WRITE (16,8028) betah,alphh-alphdiff(itype)
 8027    FORMAT(2x,'-->',10x,' input beta=',f9.5,
     +                      '  input alpha=',f9.5)
 8028    FORMAT(2x,'-->',10x,'output beta=',f9.5,
     +                      ' output alpha=',f9.5)

c  ff    do the same for mixed density: rho21 = mix_b * rho21
         rho11 = chmom(itype,1)
         rho22 = chmom(itype,2)
         rho21 = qa21(itype)
         rho21 = mix_b * rho21
         CALL rot_den_mat(alph(itype),beta(itype),
     X                       rho11,rho22,rho21)
c--->    determine the polar angles of the mom. vec. in the global frame
         mx_mix = 2*real(rho21)
         my_mix = 2*aimag(rho21)
         mz_mix = rho11 - rho22
         WRITE  (6,8031) mx_mix,my_mix
         WRITE (16,8031) mx_mix,my_mix 
 8031 FORMAT(2x,'--> global frame: ','mixed mx=',f9.5,' mixed my=',f9.5)
c if magnetic moment (in local frame!) is negative, direction of quantization
c has to be antiparallel! 
         mz_tmp = chmom(itype,1) - chmom(itype,2) 
         IF ( mz_tmp .LT. 0.0 ) THEN
            mx_mix = (-1.0) * mx_mix
            my_mix = (-1.0) * my_mix
            mz_mix = (-1.0) * mz_mix
         ENDIF
c calculate angles alpha and beta in global frame
         CALL pol_angle(mx_mix,my_mix,mz_mix,betah,alphh)
         WRITE  (6,8029) betah,alphh-alphdiff(itype)
         WRITE (16,8029) betah,alphh-alphdiff(itype)
 8029    FORMAT(2x,'-->',10x,' new beta  =',f9.5,
     +                      '  new alpha  =',f9.5)
         alph(itype) = alphh
         beta(itype) = betah
      ENDIF

      IF (l_constr) THEN
c--->    calculate the average value of B_xc (<B_xc>)
         DO iri = 1,jri(itype)
            b_xc_h(iri) = (  vr0(iri,itype,1)
     +                     - vr0(iri,itype,2) )*rmsh(iri,itype)
         ENDDO
         CALL intgr3(b_xc_h,rmsh(1,itype),dx(itype),jri(itype),b_xavh)
         b_xav(itype) = fpi*b_xavh/volmts(itype)
c--->    calculate the output constraint B-field (B_con)
!        take negative of absolute value! gb`05
         scale = -abs(b_xav(itype)/(chmom(itype,1)-chmom(itype,2)))
         b_con_outx = scale*mx
         b_con_outy = scale*my
c--->    mix input and output constraint fields
         WRITE  (6,8100) b_con(1,itype),b_con(2,itype)
         WRITE (16,8100) b_con(1,itype),b_con(2,itype)
         WRITE  (6,8200) b_con_outx,b_con_outy
         WRITE (16,8200) b_con_outx,b_con_outy
         b_con(1,itype) = b_con(1,itype) + mix_b*b_con_outx
         b_con(2,itype) = b_con(2,itype) + mix_b*b_con_outy
      ENDIF

 8100 FORMAT (2x,'-->',10x,' input B_con_x=',f12.6,
     +                    '  input B_con_y=',f12.6,
     +                    ' B_xc average=',f12.6)
 8200 FORMAT (2x,'-->',10x,' delta B_con_x=',f12.6,
     +                    ' delta B_con_y=',f12.6)

      END
