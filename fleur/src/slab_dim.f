      MODULE m_slabdim
      CONTAINS
      SUBROUTINE slab_dim(   
     >                      natd,ntypd,
     >                      ntype,taual,neq,
     <                      nsld)
c***********************************************************************
c     This subroutine calculates  the number of layers in the slab 
c     
c                                   Yury Koroteev 2003-09-30
c***********************************************************************
C                     ABBREVIATIONS
C
C natd                  : in, the number of atoms in the film
C taual(3,natd)         : in, the coordinates of atoms in the film
C ntypd,ntype           : in, the number of mt-sphere types
C neq(ntypd)            : in, the number of mt-spheres of the same type
c-----------------------------------------------------------------------
C nsld                  : out, the number of layers in the film
c-----------------------------------------------------------------------
C znz(nsl)              : work, the z-ordinate of mt-spheres in 
C                               the nsl-layer 
C-----------------------------------------------------------------------
C
	IMPLICIT NONE
C	..
C	..Scalar Argument
	INTEGER, INTENT  (IN) :: natd,ntypd,ntype
	INTEGER, INTENT (OUT) :: nsld
C	..
C	..Array Arguments
	INTEGER, INTENT  (IN) :: neq(ntypd)
	REAL,    INTENT  (IN) :: taual(3,natd)
C	..
C	..Local Scalars 
	INTEGER nz,iz,i,j,na
	REAL    epsz,zs
C	..
C	..Local Arrays 
	REAL    znz(natd)
C	..
C	..Intrinsic Functions 
	INTRINSIC abs
C    ----------------------------------------------
	DATA epsz/1.e-3/ 
C    ----------------------------------------------
C
C --->  Calculate the number of the film layers (nsld)
C
	znz(1) = taual(3,1)
	nz = 1
	na = 0
	DO 1 i=1,ntype
	   DO 2 j=1,neq(i)
	        na = na + 1
	        zs = taual(3,na)
	        DO iz=1,nz
	           IF(abs(zs-znz(iz)).LT.epsz) GO TO 2
	        ENDDO
 	        nz = nz+1
	        znz(nz) = zs
 2         CONTINUE
 1	CONTINUE
	nsld = nz
	IF(nsld.GT.natd)  STOP ' in  slab_dim.F   nsld.GT.natd '  
C
      END SUBROUTINE slab_dim
      END MODULE m_slabdim

