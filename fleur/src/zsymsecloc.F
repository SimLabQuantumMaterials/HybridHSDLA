      MODULE m_zsymsecloc
c*******************************************************
c  Solve the generalized secular equation. 
c  For film-systems exhibiting
c  z-reflexion symmetry, the basis is transformed to
c  even and odd functions and the even-even and odd-odd 
c  blocks are diagonalized separately.
c  If local orbitals are present in a film with z-reflection,
c  locrectify is used to construct linear combinations of
c  the local orbitals that are eigenfunctions of the z-
c  reflexion operation.
c  Frank Freimuth, January 2006
c*******************************************************
      CONTAINS
      SUBROUTINE zsymsecloc(
     >     film,nvd,k1,k2,k3,bkpt,nlotot,
     >     kveclo,ntypd,ntype,nlod,nlo,neq,llo,natd,taual,
     >     ngopr,nop,invtab,mrot,bmat,nv,
     >                  nbasfcn,matsize,neigd,
     >                  nsize,lb,ub,
     >                  l_zref,l_J,matind,nred,
     >                  a,b,
     >                  invsat,invsatnr,
     <                  z,eig,ne)

#include"cpp_double.h"

      USE m_locrectify
      USE m_geneigprobl

      IMPLICIT NONE
      logical,intent(in)::film
      integer,intent(in)::natd
      integer,intent(in)::invsat(natd)
      integer,intent(in)::invsatnr(natd)
      integer,intent(in)::nvd,nv
      integer,intent(in)::k1(nvd)
      integer,intent(in)::k2(nvd)
      integer,intent(in)::k3(nvd)
      real,intent(in)   ::bkpt(3)
      integer,intent(in)::nlotot
      integer,intent(in)::kveclo(nlotot)
      integer,intent(in)::ntypd,ntype,nlod
      integer,intent(in)::nlo(ntypd),neq(ntypd)
      integer,intent(in)::llo(nlod,ntypd)
      real,intent(in)   ::taual(3,natd)
      integer,intent(in)::ngopr(natd)
      integer,intent(in)::nop
      integer,intent(in)::invtab(nop)
      integer,intent(in)::mrot(3,3,nop)
      real,intent(in)   ::bmat(3,3)
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN)    :: nbasfcn,matsize,neigd
      INTEGER, INTENT (INOUT) :: nsize
      INTEGER, INTENT (OUT)   :: ne
      INTEGER, INTENT (IN)    :: nred
      REAL,    INTENT (IN)    :: lb,ub
      LOGICAL, INTENT (IN)    :: l_zref,l_J
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: matind(nbasfcn,2)
      REAL,    INTENT (OUT) :: eig(neigd)
#ifdef CPP_F90

#ifdef CPP_INVERSION
      REAL,  INTENT (INOUT) :: a(:),b(:)
      REAL,  INTENT (INOUT) :: z(:,:)
#else
      COMPLEX, INTENT (INOUT)::a(:),b(:)
      COMPLEX, INTENT (INOUT) :: z(:,:)
#endif

#else

#ifdef CPP_INVERSION
      REAL, ALLOCATABLE, INTENT (INOUT) :: a(:),b(:)
      REAL, ALLOCATABLE, INTENT (INOUT) :: z(:,:)
#else
      COMPLEX, ALLOCATABLE, INTENT (INOUT)::a(:),b(:)
      COMPLEX, ALLOCATABLE, INTENT (INOUT) :: z(:,:)
#endif

#endif

#ifdef CPP_INVERSION
      real locrec(nlotot,nlotot)
      real,allocatable::atemp(:,:)
      real,allocatable::btemp(:,:)
      real,allocatable::aa(:)
      real,allocatable::bb(:)
      real,allocatable::z1(:,:)
      real,allocatable::z2(:,:)
      real,allocatable::ztemp(:,:)
      real recsqrtwo,dobzero
#else
      complex locrec(nlotot,nlotot)
      complex,allocatable:: atemp(:,:)
      complex,allocatable:: btemp(:,:)
      complex,allocatable::aa(:)
      complex,allocatable::bb(:)
      complex,allocatable::z1(:,:)
      complex,allocatable::z2(:,:)
      complex,allocatable::ztemp(:,:)
      complex recsqrtwo,dobzero
#endif
      logical l_verbose
      logical kindlocrec(nlotot)
      real,allocatable::etemp1(:),etemp2(:)
      integer iind,jind,ii,info,iu,pos1,pos2,pos3,i2,j2
      integer ne1,ne2,i1,j1,i,j
      logical,allocatable:: evensort(:)
      integer evenlocs,oddlocs
      integer evenlocindex(nlotot)
      integer oddlocindex(nlotot)


c      print*,"in zsymsecloc"
#ifndef CPP_F90
      deallocate(z)
#endif

!******************************************
! l_zref=.false. => simply call eigensolver
!******************************************
      if(.not.l_zref)then
         call geneigprobl(
     >                    nbasfcn,
     >                    nsize,neigd,lb,ub,l_J,a,b,
     <                    z,eig,ne)

#ifndef CPP_F90
         allocate(a(nbasfcn*(nbasfcn+1)/2))
         allocate(b(nbasfcn*(nbasfcn+1)/2))
#endif
         return
!******************************************
! l_zref=.true. => blockdiagonalize
! hamiltonian and overlap-matrix by
! transforming to even/odd basisfunctions
!******************************************
      else
       inquire(file='verbose',exist=l_verbose)
       if(.not.l_verbose)inquire(file='veryverbose',exist=l_verbose)
       iu=neigd/2
#ifdef CPP_INVERSION
       recsqrtwo=1.0/sqrt(2.0)
       dobzero=0.0
#else
       recsqrtwo=cmplx(1.0/sqrt(2.0),0.0)
       dobzero=cmplx(0.0,0.0)
#endif


      if(nlotot.gt.0)then
!*********************************************
! Find even/odd linear combinations of locs
!*********************************************
      if(l_verbose)then
      print*,"find lincos of locs that are eigenfunctions
     &                                    of zreflection"
         print*,"apws=",nv
         print*,"nlotot=",nlotot
         print*,"basis-size=",nsize
      endif
      if(nsize.ne.(nv+nlotot))stop 'nsize /= nv + nlotot'
      call locrectify(
     >     film,nvd,k1,k2,k3,bkpt,nlotot,
     >     kveclo,ntypd,ntype,nlod,nlo,neq,llo,natd,taual,
     >     ngopr,nop,invtab,mrot,bmat,
     >     invsat,
     >     invsatnr,
     <     locrec,kindlocrec,evenlocs,oddlocs,
     <     evenlocindex,oddlocindex)


!*********************************************
! Perform basis-transformation of Hamiltonian
! and Overlap matrix. The new basis is the basis
! of even / odd (with respect to z-reflection)
! local orbitals.
!*********************************************
      allocate(atemp(nv,nlotot))
      allocate(btemp(nlotot,nv))
      do iind=1,nlotot
         pos1=((nv+iind)*(nv+iind-1))/2
         atemp(1:nv,iind)=a(pos1+1:pos1+nv)
      enddo
#ifdef CPP_INVERSION
      call CPP_BLAS_sgemm('T','T',nlotot,nv,nlotot,real(1.0),
     &    locrec,nlotot,atemp,nv,real(0.0),btemp,nlotot)
#else
      call CPP_BLAS_cgemm('C','T',nlotot,nv,nlotot,cmplx(1.0,0.0),
     &    locrec,nlotot,atemp,nv,cmplx(0.0,0.0),btemp,nlotot)
#endif
      atemp=transpose(btemp)
      do iind=1,nlotot
         pos1=((nv+iind)*(nv+iind-1))/2
         a(pos1+1:pos1+nv)=atemp(1:nv,iind)
      enddo

      do iind=1,nlotot
         pos1=((nv+iind)*(nv+iind-1))/2
         atemp(1:nv,iind)=b(pos1+1:pos1+nv)
      enddo
#ifdef CPP_INVERSION
      call CPP_BLAS_sgemm('T','T',nlotot,nv,nlotot,real(1.0),
     &    locrec,nlotot,atemp,nv,real(0.0),btemp,nlotot)
#else
      call CPP_BLAS_cgemm('C','T',nlotot,nv,nlotot,cmplx(1.0,0.0),
     &    locrec,nlotot,atemp,nv,cmplx(0.0,0.0),btemp,nlotot)
#endif
      atemp=transpose(btemp)
      do iind=1,nlotot
         pos1=((nv+iind)*(nv+iind-1))/2
         b(pos1+1:pos1+nv)=atemp(1:nv,iind)
      enddo
      deallocate(atemp)
      deallocate(btemp)

      allocate(atemp(nlotot,nlotot))
      allocate(btemp(nlotot,nlotot))
      do iind=1,nlotot
       pos1=((nv+iind)*(nv+iind-1))/2
       do jind=1,iind-1
        atemp(iind,jind)=a(pos1+nv+jind)
#ifdef CPP_INVERSION
        atemp(jind,iind)=a(pos1+nv+jind)
#else
        atemp(jind,iind)=conjg(a(pos1+nv+jind))
#endif
       enddo
       atemp(iind,iind)=a(pos1+nv+iind)
      enddo
#ifdef CPP_INVERSION
      call CPP_BLAS_sgemm('T','N',nlotot,nlotot,nlotot,real(1.0),
     &     locrec,nlotot,atemp,nlotot,real(0.0),btemp,nlotot)
      call CPP_BLAS_sgemm('N','N',nlotot,nlotot,nlotot,real(1.0),
     &     btemp,nlotot,locrec,nlotot,real(0.0),atemp,nlotot)
#else
      call CPP_BLAS_cgemm('C','N',nlotot,nlotot,nlotot,cmplx(1.0,0.0),
     &     locrec,nlotot,atemp,nlotot,cmplx(0.0,0.0),btemp,nlotot)
      call CPP_BLAS_cgemm('N','N',nlotot,nlotot,nlotot,cmplx(1.0,0.0),
     &     btemp,nlotot,locrec,nlotot,cmplx(0.0,0.0),atemp,nlotot)
#endif
      do iind=1,nlotot
       pos1=((nv+iind)*(nv+iind-1))/2
       do jind=1,iind
        a(pos1+nv+jind)=atemp(iind,jind)
       enddo
      enddo

      do iind=1,nlotot
       pos1=((nv+iind)*(nv+iind-1))/2
       do jind=1,iind-1
        atemp(iind,jind)=b(pos1+nv+jind)
#ifdef CPP_INVERSION
        atemp(jind,iind)=b(pos1+nv+jind)
#else
        atemp(jind,iind)=conjg(b(pos1+nv+jind))
#endif
       enddo
       atemp(iind,iind)=b(pos1+nv+iind)
      enddo
#ifdef CPP_INVERSION
      call CPP_BLAS_sgemm('T','N',nlotot,nlotot,nlotot,real(1.0),
     &     locrec,nlotot,atemp,nlotot,real(0.0),btemp,nlotot)
      call CPP_BLAS_sgemm('N','N',nlotot,nlotot,nlotot,real(1.0),
     &     btemp,nlotot,locrec,nlotot,real(0.0),atemp,nlotot)
#else
      call CPP_BLAS_cgemm('C','N',nlotot,nlotot,nlotot,cmplx(1.0,0.0),
     &     locrec,nlotot,atemp,nlotot,cmplx(0.0,0.0),btemp,nlotot)
      call CPP_BLAS_cgemm('N','N',nlotot,nlotot,nlotot,cmplx(1.0,0.0),
     &     btemp,nlotot,locrec,nlotot,cmplx(0.0,0.0),atemp,nlotot)
#endif
      do iind=1,nlotot
       pos1=((nv+iind)*(nv+iind-1))/2
       do jind=1,iind
        b(pos1+nv+jind)=atemp(iind,jind)
       enddo
      enddo

      deallocate(atemp)
      deallocate(btemp)

      else
         evenlocs=0
         oddlocs=0
      endif ! nlotot.gt.0

!*********************************************
! Test matind.
!*********************************************
      ii=0
      jind=0
      pos1=0
      do iind=1,nred
         if(matind(iind,1).lt.matind(iind,2)) stop 'mat1mat2'
         if(matind(iind,1).lt.pos1)stop 'matpos1'
         pos1=matind(iind,1)
         if(matind(iind,1).ne.matind(iind,2))then
            ii=ii+1
         else
            jind=jind+1
         endif
      enddo
      if(2*ii+jind.ne.nv)stop 'matind'

!*****************************************************************
!Transform into representation with even-even- and odd-odd-blocks.
!First step: Transform the lapw-lapw-part.
!*****************************************************************
      allocate(    aa(  ((nred+evenlocs)*(nred+evenlocs+1))/2   )  )
      allocate(    bb(  ((nred+evenlocs)*(nred+evenlocs+1))/2   )  )
      i2=0
      j2=0
      DO i=1,nred
          DO j=1,i
            i1=(matind(i,1)-1)*matind(i,1)/2+matind(j,1)
            j1=(matind(i,1)-1)*matind(i,1)/2+matind(j,2)
            i2=i2+1
            aa(i2)=a(i1)+a(j1)
            bb(i2)=b(i1)+b(j1)
            IF ((matind(i,1).NE.matind(i,2)).AND.
     +          (matind(j,1).NE.matind(j,2))) THEN
               j2=j2+1
               a(j2)=a(i1)-a(j1)
               b(j2)=b(i1)-b(j1)
            ENDIF
          ENDDO
      ENDDO

      if(nlotot.gt.0)then
!******************************************************
!the lapw-lo- and lo-lo-parts for the even-even-block.
!******************************************************
      allocate(atemp(nv,nlotot))
      allocate(btemp(nred,evenlocs))
!hamiltonian
      do iind=1,nlotot
         pos1=((nv+iind)*(nv+iind-1))/2
         atemp(1:nv,iind)=a(pos1+1:pos1+nv)
      enddo
      do jind=1,evenlocs
       do iind=1,nred
            btemp(iind,jind)=(atemp(matind(iind,1),evenlocindex(jind))
     &        +atemp(matind(iind,2),evenlocindex(jind)))*recsqrtwo
       enddo
      enddo
      do iind=1,evenlocs
         pos1=((nred+iind)*(nred+iind-1))/2
         aa(pos1+1:pos1+nred)=btemp(1:nred,iind)
      enddo
!metric
      do iind=1,nlotot
         pos1=((nv+iind)*(nv+iind-1))/2
         atemp(1:nv,iind)=b(pos1+1:pos1+nv)
      enddo
      do jind=1,evenlocs
       do iind=1,nred
            btemp(iind,jind)=(atemp(matind(iind,1),evenlocindex(jind))
     &        +atemp(matind(iind,2),evenlocindex(jind)))*recsqrtwo
       enddo
      enddo
      do iind=1,evenlocs
         pos1=((nred+iind)*(nred+iind-1))/2
         bb(pos1+1:pos1+nred)=btemp(1:nred,iind)
      enddo
      deallocate(btemp)
      deallocate(atemp)
!hamiltonian
      do iind=1,evenlocs
         do jind=1,iind
            pos1=((nv+evenlocindex(iind))*(nv+evenlocindex(iind)-1))/2
            pos2=((nred+iind)*(nred+iind-1))/2
            aa(pos2+nred+jind)=a(pos1+nv+evenlocindex(jind))
         enddo
      enddo
!metric
      do iind=1,evenlocs
         do jind=1,iind
            pos1=((nv+evenlocindex(iind))*(nv+evenlocindex(iind)-1))/2
            pos2=((nred+iind)*(nred+iind-1))/2
            bb(pos2+nred+jind)=b(pos1+nv+evenlocindex(jind))
         enddo
      enddo


!******************************************************
!the lapw-lo- and lo-lo-parts for the odd-odd-block
!******************************************************
      allocate(atemp(nv,nlotot))
      allocate(btemp(nv-nred,oddlocs))
!hamiltonian
      do iind=1,nlotot
         pos1=((nv+iind)*(nv+iind-1))/2
         atemp(1:nv,iind)=a(pos1+1:pos1+nv)
      enddo
      do jind=1,oddlocs
       ii=0
       do iind=1,nred
         if(matind(iind,1).ne.matind(iind,2))then
            ii=ii+1
            btemp(ii,jind)=(atemp(matind(iind,1),oddlocindex(jind))
     &        -atemp(matind(iind,2),oddlocindex(jind)))*recsqrtwo
         endif
       enddo
      enddo
      do iind=1,oddlocs
         pos1=((nv-nred+iind)*(nv-nred+iind-1))/2
         a(pos1+1:pos1+nv-nred)=btemp(1:nv-nred,iind)
      enddo
!metric
      do iind=1,nlotot
         pos1=((nv+iind)*(nv+iind-1))/2
         atemp(1:nv,iind)=b(pos1+1:pos1+nv)
      enddo
      do jind=1,oddlocs
       ii=0
       do iind=1,nred
         if(matind(iind,1).ne.matind(iind,2))then
            ii=ii+1
            btemp(ii,jind)=(atemp(matind(iind,1),oddlocindex(jind))
     &        -atemp(matind(iind,2),oddlocindex(jind)))*recsqrtwo
         endif
       enddo
      enddo
      do iind=1,oddlocs
         pos1=((nv-nred+iind)*(nv-nred+iind-1))/2
         b(pos1+1:pos1+nv-nred)=btemp(1:nv-nred,iind)
      enddo
      deallocate(btemp)
      deallocate(atemp)
!hamiltonian
      do iind=1,oddlocs
         do jind=1,iind
            pos1=((nv+oddlocindex(iind))*(nv+oddlocindex(iind)-1))/2
            pos2=((nv-nred+iind)*(nv-nred+iind-1))/2
            a(pos2+nv-nred+jind)=a(pos1+nv+oddlocindex(jind))
         enddo
      enddo
!metric
      do iind=1,oddlocs
         do jind=1,iind
            pos1=((nv+oddlocindex(iind))*(nv+oddlocindex(iind)-1))/2
            pos2=((nv-nred+iind)*(nv-nred+iind-1))/2
            b(pos2+nv-nred+jind)=b(pos1+nv+oddlocindex(jind))
         enddo
      enddo



      endif !nlotot.gt.0

!******************************************************
!Solve the eigenvalue problem for the odd-odd block
!******************************************************

      allocate(etemp2(iu+1))
      call geneigprobl(
     >                    nv-nred+oddlocs,
     >                    nv-nred+oddlocs,iu,lb,ub,l_J,a,b,
     <                    z2,etemp2,ne2)

      if(l_verbose)then   
        print*,"odd block diagonalized"
      endif

!***********************************************************
!Solve the eigenvalue problem for the even-even block
!***********************************************************
      allocate(etemp1(iu+1))
      call geneigprobl(
     >                    nred+evenlocs,
     >                    nred+evenlocs,iu,lb,ub,l_J,aa,bb,
     <                    z1,etemp1,ne1)
      if(l_verbose)then
        print*,"even block diagonalized"
      endif


      ne=ne1+ne2
!********************************************************************
!  Recover eigenvectors of original eigenvalue problem.
!  Sort eigenvalues and eigenvectors according to increasing eigenvalue.
!  etemp1 holds eigenvalues of even block.
!  etemp2 holds eigenvalues of odd block.
!  z1 holds eigenvectors of even block.
!  z2 holds eigenvectors of odd block.
!********************************************************************
#ifndef CPP_F90
      allocate(z(nbasfcn,neigd))
#endif
      allocate(evensort(ne))
      etemp1(ne1+1)=99.9e9
      etemp2(ne2+1)=99.9e9
      jind=1
      iind=1
!evensort(ii)=.true.  => eigenvalue ii belongs to even spectrum
!evensort(ii)=.false. => eigenvalue ii belongs to odd spectrum
      do ii=1,ne
         if(etemp1(iind).lt.etemp2(jind)) then
            evensort(ii)=.true.
            iind=iind+1
         else
            evensort(ii)=.false.
            jind=jind+1
         endif
      enddo
      iind=1 !position in the even-arrays
      jind=1 !position in the odd-arrays
      do ii=1,ne
       if(evensort(ii))then
        eig(ii)=etemp1(iind)
        z(1:nv+nlotot,ii)=dobzero
!Recover the eigenvectors of the original problem for the even block
        do i=1,nred
           z(matind(i,1),ii)=z1(i,iind)*recsqrtwo
           z(matind(i,2),ii)=z(matind(i,2),ii)+z1(i,iind)*recsqrtwo
        enddo !i
        if(nlotot.gt.0) z(nv+1:nv+nlotot,ii)=dobzero
        do pos1=1,evenlocs
           do i=1,nlotot
#ifdef CPP_INVERSION
              z(nv+i,ii)=z(nv+i,ii)+
     &          z1(nred+pos1,iind)*locrec(i,evenlocindex(pos1))
#else
              z(nv+i,ii)=z(nv+i,ii)+
     &          z1(nred+pos1,iind)*conjg(locrec(i,evenlocindex(pos1)))
#endif
           enddo
        enddo
        iind=iind+1
       else
!Recover the eigenvectors of the original problem for the odd block
        eig(ii)=etemp2(jind)
        j1=0
        do i=1,nred
          if(matind(i,1).ne.matind(i,2))then
           j1=j1+1
           z(matind(i,1),ii)=z2(j1,jind)*recsqrtwo
           z(matind(i,2),ii)=-z2(j1,jind)*recsqrtwo
          else
           z(matind(i,1),ii)=dobzero
          endif
        enddo !i
        if(nlotot.gt.0) z(nv+1:nv+nlotot,ii)=dobzero
        do pos1=1,oddlocs
           do i=1,nlotot
#ifdef CPP_INVERSION
              z(nv+i,ii)=z(nv+i,ii)+
     &          z2(nv-nred+pos1,jind)*locrec(i,oddlocindex(pos1))
#else
              z(nv+i,ii)=z(nv+i,ii)+
     &          z2(nv-nred+pos1,jind)*conjg(locrec(i,oddlocindex(pos1)))
#endif
           enddo
        enddo
        jind=jind+1
       endif !evensort
      enddo !ii

#ifndef CPP_F90
      allocate(a(nbasfcn*(nbasfcn+1)/2))
      allocate(b(nbasfcn*(nbasfcn+1)/2))
#endif
      endif !l_zref
      END SUBROUTINE
      END MODULE
