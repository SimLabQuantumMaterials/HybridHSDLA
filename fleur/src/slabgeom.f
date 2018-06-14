      MODULE m_slabgeom
      CONTAINS
      SUBROUTINE slabgeom(   
     >                      natd,ntypd,nsld,
     >                      ntype,pos,z1,neq,area,volmts,
     <                      nsl,zsl,nmtsl,nslat,volsl,volintsl)
c***********************************************************************
c     This subroutine calculates  z-coordinates of film layers, 
c     a number of mt-pheres in each layer, and they typs.
c                                   Yury Koroteev 2003-09-30
c***********************************************************************
C                     ABBREVIATIONS
C
C natd                  : in, the number of atoms in the film
C pos(3,natd)           : in, the coordinates of atoms in the film
C ntypd,ntype           : in, the number of mt-sphere types
C z1                    : in, half the film thickness (0.5*D_tilde)
C neq(ntypd)            : in, the number of mt-spheres of the same type
C area                  : in, the area of the surface unit cell
C volmts(ntypd)         : in, the volume of mt-spheres
C nsld                  : in, the number of layers in the film
c-----------------------------------------------------------------------
C nsl                   : in, the number of layers in the film
C zsl(2,nsld)           : out, z-coordinates of the layers
C nmtsl(ntypd,nsld)     : out, the number of mt-spheres of the ntypd-
C                                type in the nsl-layer of the film
C nslat(natd,nsld)      : out, 
C volsl(nsld)           : out, the volume of film layers  
C volintsl(nsld)        : out, the volume of mt-spheres
C
c-----------------------------------------------------------------------
C znz(nsl)              : work, the z-ordinate of mt-spheres in 
C                               the nsl-layer 
C-----------------------------------------------------------------------
C
	IMPLICIT NONE
C	..
C	..Scalar Argument
	INTEGER, INTENT  (IN) :: natd,ntypd,ntype,nsld
	REAL,    INTENT  (IN) :: z1,area
	INTEGER, INTENT (OUT) :: nsl
C	..
C	..Array Arguments
	INTEGER, INTENT  (IN) :: neq(ntypd)
	REAL,    INTENT  (IN) :: pos(3,natd),volmts(ntypd)
	INTEGER, INTENT (OUT) :: nmtsl(ntypd,nsld),nslat(natd,nsld)
	REAL,    INTENT (OUT) :: zsl(2,nsld),volsl(nsld)  
	REAL,    INTENT (OUT) :: volintsl(nsld)
C	..
C	..Local Scalars 
	INTEGER nz,iz,i,j,na,isum,mt,n
	REAL    epsz,half,zs,w,del,vmt
C	..
C	..Local Arrays 
	REAL    znz(nsld)
C	..
C	..Intrinsic Functions 
	INTRINSIC abs
C    ------------------------------------------------------------------
	DATA epsz/1.e-3/ half/0.5/
C    ----------------------------------------------
C
C --->  Calculate the number of the film layers (nsl)
C
        epsz = epsz * 2 * z1 ! here pos() is used instead of taual() !
	znz(1) = pos(3,1)
	nz = 1
	na = 0
	DO 1 i=1,ntype
	   DO 2 j=1,neq(i)
	        na = na + 1
	        zs = pos(3,na)
	        DO iz=1,nz
	           IF(abs(zs-znz(iz)).LT.epsz) GO TO 2
	        ENDDO
 	        nz = nz+1
	        znz(nz) = zs
 2         CONTINUE
 1	CONTINUE
	nsl = nz
        IF (nsl.GT.nsld) THEN
            WRITE(*,*) 'nsl =',nsl,' > nsld =',nsld
            STOP ' in  slabgeom.f  nsl.GT.nsld '
        ENDIF
C
C ---> Order the film layers
C
	DO  i=1,nsl
	  DO  j=i,nsl
	     IF(znz(j).LT.znz(i)) THEN
	       w      = znz(i)
	       znz(i) = znz(j)
	       znz(j) = w
	     ENDIF
 	  ENDDO
	ENDDO
C
C ---> Construct the z-coordinates of the film layers ( zsl(2,nsl) )
C
	zsl(1,1) = -z1
	DO i=1,nsl-1
	   zsl(2,i) = (znz(i) + znz(i+1)) * half
	   zsl(1,i+1) = zsl(2,i)
 	ENDDO
 	zsl(2,nsl) = z1
C 
C ---> Calculate a number of mt-spheres of the same type
C ---> (nmtsl) in aech layer of the film
C
	DO i=1,nsl
	   del = abs( zsl(2,i) - zsl(1,i) )
           volsl(i) = del*area
	   n = 0
	   vmt = 0.0
	   DO j=1,ntype
	      isum = 0
              DO mt=1,neq(j)
	         n = n + 1
	         zs = pos(3,n)
		 IF((zsl(1,i).LT.zs).AND.(zs.LT.zsl(2,i)))  THEN
                     isum=isum+1
                     nslat(n,i)=1
                   ELSE
                     nslat(n,i)=0
                 ENDIF
              ENDDO
 	      nmtsl(j,i) = isum
	      vmt = vmt + isum*volmts(j)
 	   ENDDO
	   volintsl(i) = volsl(i) - vmt
 	ENDDO
C
      END SUBROUTINE slabgeom
      END MODULE m_slabgeom

