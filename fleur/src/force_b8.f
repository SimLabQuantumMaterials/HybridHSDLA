      MODULE m_forceb8
c-------------------
c Implements the surface contribution to the force. Madsen Eq.(B8)
c
c FZJ 15/3-01 GMadsen
c---------------------------------
      CONTAINS
      SUBROUTINE force_b8(
     >                    ntype,ntypd,natd,jspd,neq,ecwk,nq3_fft,
     >                    tau,rmt,kv3,n3d,bmat,lmax,lmaxd,taual,
     >                    nop,mrot,bbmat,jspin,
     X                    force,f_b8)

      USE m_constants, ONLY : pimach
      USE m_sphbes
      USE m_stern

      IMPLICIT NONE
C     ..
C     .. Arguments
      INTEGER, INTENT (IN) :: n3d,nq3_fft,ntypd,ntype,natd,nop,jspin
      INTEGER, INTENT (IN) :: lmaxd,jspd
      INTEGER, INTENT (IN) :: mrot(3,3,nop)
      INTEGER, INTENT (IN) :: kv3(3,n3d),lmax(ntypd),neq(ntypd)
      REAL,    INTENT (IN) :: tau(3,nop),bmat(3,3),bbmat(3,3)
      REAL,    INTENT (IN) :: taual(3,natd)
      COMPLEX, INTENT (IN) :: ecwk(n3d)
      COMPLEX, INTENT (INOUT) :: f_b8(3,ntypd)
      REAL,    INTENT (INOUT) :: force(3,ntypd,jspd)
C     ..
C     .. Local Variables
      INTEGER g(3),nst,stg(3,nop),ia,istr,i,j,jj,jneq
      REAL    fj(0:lmaxd),rotkzz(3),rstg(3,nop)
      REAL    frmt,rmt(ntypd),gl,tpi,pha,s
      COMPLEX taup(nop),factor,fact,fstar(3),fsur2(3)
c
      tpi = 2 * pimach()

      ia=1
      DO jneq=1,ntype
         frmt = 2.0*tpi*rmt(jneq)**2
         fsur2(1:3) = cmplx(0.0,0.0)         

c skip G=(0,0,0) no contribution to ekin
         DO istr=2,nq3_fft

            DO i=1,3
               g(i)     = kv3(i,istr)
               fstar(i) = cmplx(0.0,0.0)
            ENDDO
            CALL stern(
     >                 nop,tau,g,mrot,bbmat,bmat,
     <                 nst,stg,taup,gl,rstg)

            CALL sphbes(lmax(jneq),rmt(jneq)*gl,fj)
            fact = ecwk(istr) * fj(1) / gl

            DO jj=1,nst
               pha=(taual(1,ia)*stg(1,jj)+taual(2,ia)*stg(2,jj)
     $             +taual(3,ia)*stg(3,jj))*tpi
c
ccc   swaped sin and cos because there's an i in the equation
c
               factor = fact * cmplx(-sin(pha),cos(pha)) * taup(jj)
               DO i=1,3
                  fstar(i) = fstar(i) + factor*rstg(i,jj)
               ENDDO
            ENDDO
            DO i=1,3
               fsur2(i)=fsur2(i)+fstar(i)*frmt
            ENDDO
         ENDDO
         DO i=1,3
            f_b8(i,jneq) = f_b8(i,jneq) + fsur2(i)
            force(i,jneq,jspin) = force(i,jneq,jspin) + real(fsur2(i))
         ENDDO
         ia=ia+neq(jneq)
      ENDDO

      END SUBROUTINE force_b8
      END MODULE m_forceb8
