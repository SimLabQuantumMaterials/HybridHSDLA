c*************************************************************
c     Find linear combinations of local orbitals that are
c     eigenfunctions of the z-reflection operation.
c     Used to exploit z-reflection symmetry when solving
c     the eigenvalue problem (film calculations).
c     Frank Freimuth, January 2007
c*************************************************************
#include"cpp_double.h"
      MODULE m_locrectify
      CONTAINS
      SUBROUTINE locrectify(
     >     film,nvd,k1,k2,k3,bkpt,nlotot,
     >     kveclo,ntypd,ntype,nlod,nlo,neq,llo,natd,taual,
     >     ngopr,nop,invtab,mrot,bmat,
     >     invsat,
     >     invsatnr,
     <     locrec,kindlocrec,evenlocs,oddlocs,evenlocindex,
     <     oddlocindex)

      USE m_ylm
      USE m_constants,only:pimach
      USE m_cotra,only:cotra3
      USE m_loccoeff

      IMPLICIT NONE
      logical,intent(in)::film
      integer,intent(in)::natd
      integer,intent(in)::invsat(natd)
      integer,intent(in)::invsatnr(natd)
      integer,intent(in)::nvd
      integer,intent(in)::k1(nvd)
      integer,intent(in)::k2(nvd)
      integer,intent(in)::k3(nvd)
      real,intent(in)::bkpt(3)
      integer,intent(in)::nlotot
      integer,intent(in)::kveclo(nlotot)
      integer,intent(in)::ntypd,ntype,nlod
      integer,intent(in)::nlo(ntypd),neq(ntypd)
      integer,intent(in)::llo(nlod,ntypd)
      real,intent(in)::taual(3,natd)
      integer,intent(in)::ngopr(natd)
      integer,intent(in)::nop
      integer,intent(in)::invtab(nop)
      integer,intent(in)::mrot(3,3,nop)
      real,intent(in)::bmat(3,3)

#ifdef CPP_INVERSION
      real,intent(out)::locrec(nlotot,nlotot)
#else
      complex,intent(out)::locrec(nlotot,nlotot)
#endif
      logical,intent(out)::kindlocrec(nlotot)
      integer,intent(out)::evenlocs
      integer,intent(out)::oddlocs
      integer,intent(out)::evenlocindex(nlotot)
      integer,intent(out)::oddlocindex(nlotot)

      integer numorec
      integer invatom,zrefatom
      logical l_planetwo,l_spacetwo,l_spacefour,l_spacetwoinv
      logical l_crosstwo
      integer ind,loc,nn,angmom,bas,locm,loopi,j,i
      integer lm,locmm,lmm
      integer pos,vec,natom,nap,inap,locnum,ll1
      integer ipiv(28),info
      complex loccoffs(28,28),reccoffs(28)
      complex ylm(16),phase
      real fk(3),tmk,tpi,fkp(3)
!      complex coeffi(1:28,1:28,0:13)     transferred to mod_loccoeff.f
!      logical coffkind(1:28,0:13)        import with USE m_loccoeff
      integer jump(natd)
      integer mqnum,nalo,naup,basindex(natd,nlod)
      integer zrefnap,zrefinap,coffindex
      logical l_verbose,l_veryverbose,l_invsatswap

!      ci=cmplx(0.0,1.0) transferred to mod_loccoeff.f


!diagnostic tools
      inquire(file='verbose',exist=l_verbose)
      inquire(file='veryverbose',exist=l_veryverbose)
      if(l_veryverbose)l_verbose=.true.

      jump(1:natd)=0
      locrec=0.0
      tpi=2*pimach()
      numorec=0

c***********************************************
c   basindex holds first index of subset of locs
c***********************************************
      bas=1
      nalo=1
      do ind=1,ntype
         naup=nalo+neq(ind)-1
         do natom=nalo,naup
            if(invsat(natom).eq.2)cycle
            do loc=1,nlo(ind)
               basindex(natom,loc)=bas
               locnum=(invsat(natom)+1)*(2*llo(loc,ind)+1)
               bas=bas+locnum
            enddo !loc
         enddo !natom
         nalo=naup+1
      enddo !ind

      bas=0
      nalo=1
      do ind=1,ntype
       naup=nalo+neq(ind)-1
       do natom=nalo,naup
        if(l_verbose)then
          print*,"invsat=",invsat(natom)
          print*,"taual=",taual(:,natom)
        endif
        if(invsat(natom).eq.2)cycle
        if(jump(natom).ne.0)then
           bas=bas+jump(natom)
           cycle
        endif

!*******************************************************
!       find out by what kind of symmetry the atoms
!       are interrelated
!*******************************************************
!two atoms related by inversion symmetry are in xy-plane
        l_planetwo=.false.    
!two atoms are related by z-reflection, but not by inversion
        l_spacetwo=.false.
!two atoms are related by z-reflection and by inversion
        l_spacetwoinv=.false.
!four atoms are entangled, by both inversion and z-reflection
        l_spacefour=.false.
!two atoms are related by inversion, but not by z-reflection
        l_crosstwo=.false.
!order of locs corresponds to setup of coeffi-array => l_invsatswap=.false.
!order of locs is reversed with respect to coeffi-array => l_invsatswap=.true.
        l_invsatswap=.false.
        if((abs(abs(taual(3,natom))-0.5).lt.1e-6).and..not.film)then !atom on face
           if(invsat(natom).eq.1)then ! two atoms
              invatom=invsatnr(natom)
              l_planetwo=.true.
           endif
        elseif(abs(taual(3,natom)).gt.1e-6) then !atom not in xy-plane
           do zrefatom=nalo,naup !find z-reflection image
            if(all(abs(taual(1:2,natom)-taual(1:2,zrefatom))
     &               .lt.1e-6))then
             if(abs(taual(3,natom)+taual(3,zrefatom)).lt.1e-6) goto 543
            endif
           enddo !zrefatom
#ifdef CPP_INVERSION
!if there is no z-reflected counterpart there must be an inversion 
!counterpart instead
           if(invsat(natom).eq.1)then
              zrefatom=0 !this switches on l_crosstwo later on
           else
              stop 'nozref-invsat'
           endif
#else
!in the complex version there must be a z-reflected counterpart
           stop 'zrefatom not found' 
#endif
 543       continue

           if(invsat(natom).eq.0) then !no inversion-image
              l_spacetwo=.true.
#ifdef CPP_INVERSION
              stop 'spacetwo & INVERSION'
#endif
           else                        !inversion-image
              invatom=invsatnr(natom)
              if(zrefatom.eq.0)then
                 l_crosstwo=.true.
              elseif(invatom.eq.zrefatom)then    !inversion = z-reflection
                 l_spacetwoinv=.true.
              else                           !inversion /= z-reflection
#ifdef CPP_INVERSION
                 l_spacefour=.true. !four entangled locs
#else
                 stop 'spacefour & no INVERSION'
#endif
                 if(invsat(zrefatom).eq.1)then
                    if(l_verbose)then
                       print*,"spacefour-atoms:"
                       print*,"natom=",natom
                       print*,"invatom=",invatom
                       print*,"zrefatom=",zrefatom
                       print*,"zrefinvatom=",invsatnr(zrefatom)
                    endif
                 elseif(invsat(zrefatom).eq.2)then
                    l_invsatswap=.true.
                    do info=1,natd
                       if(invsatnr(info).eq.zrefatom)goto 723
                    enddo
                    stop 'no invsat for zref'
 723                continue
                    zrefatom=info
                    if(invsat(zrefatom).ne.1)stop 'locrectify: invsat'
                    if(.not.all(abs(taual(1:2,invatom)-
     -                 taual(1:2,zrefatom)).lt.1e-6))stop 'invsat-zref2'
                    if(.not.abs(taual(3,invatom)+taual(3,zrefatom))
     &                  .lt.1e-6) stop 'invsat-zref3'
                    if(l_verbose)then
                       print*,"spacefour-atoms:"
                       print*,"natom=",natom
                       print*,"invatom=",invsat(natom)
                       print*,"SWAPPED:"
                       print*,"zrefatom mapped on invatom by zref"
                       print*,"zrefatom=",zrefatom
                       print*,"zrefinvatom=",invsat(zrefatom)
                    endif
                 elseif(invsat(zrefatom).eq.0)then
                    stop 'spacefour: invsat-zref'
                 endif
              endif
           endif
        else !atom lies in the xy-plane
           if(invsat(natom).eq.1)then ! two atoms
              invatom=invsatnr(natom)
              l_planetwo=.true.
           endif
        endif !symmetry partners

        do loc=1,nlo(ind)
         angmom=llo(loc,ind)
         if(l_verbose)then
           print*,"angmom=",angmom
           print*,"bas=",bas
           print*,"numorec=",numorec
           if(l_planetwo)then
              print*,"planetwo"
           elseif(l_crosstwo)then
              print*,"crosstwo"
           elseif(l_spacetwoinv)then
              print*,"spacetwoinv"
           elseif(l_spacetwo)then
              print*,"spacetwo"
           elseif(l_spacefour)then
              print*,"spacefour"
           else
              print*,"single atom"
           endif
         endif !verbose
         if(angmom.eq.0)then
          if(l_planetwo)then
            kindlocrec(numorec+1)=.true.
            kindlocrec(numorec+2)=.true.
            locrec(bas+1,numorec+1)=1.0
            locrec(bas+2,numorec+2)=1.0
            bas=bas+2
            numorec=numorec+2
            cycle
          elseif(l_spacetwoinv.or.l_crosstwo)then
            coffindex=0
            locnum=2
          elseif(l_spacetwo) then
            locnum=1
            coffindex=0
          elseif(l_spacefour)then
            locnum=2
            coffindex=6
          else
            kindlocrec(numorec+1)=.true.
            locrec(bas+1,numorec+1)=1.0
            bas=bas+1
            numorec=numorec+1
            cycle
          endif !l_planetwo,l_crosstwo,....
         else
          if(l_planetwo)then
             locnum=2*(2*angmom+1)
             if(angmom.eq.1)then
                coffindex=3
             elseif(angmom.eq.2)then
                coffindex=9
             elseif(angmom.eq.3)then
                coffindex=12
             else
                stop 'angmom1'
             endif
          elseif(l_spacetwoinv.or.l_crosstwo)then
             locnum=2*(2*angmom+1)
             if(angmom.eq.1)then
               coffindex=4
             elseif(angmom.eq.2)then
               coffindex=8
             elseif(angmom.eq.3)then
               coffindex=11
             else
                stop 'angmom5'
             endif
          elseif(l_spacetwo) then
             locnum=2*angmom+1
             if(angmom.eq.1)then
                coffindex=4
             elseif(angmom.eq.2)then
                coffindex=8
             elseif(angmom.eq.3)then
                coffindex=11
             else
                stop 'angmom2'
             endif
          elseif(l_spacefour) then
             locnum=2*(2*angmom+1)
             if(angmom.eq.1)then
                coffindex=7
             elseif(angmom.eq.2)then
                coffindex=10
             elseif(angmom.eq.3)then
                coffindex=13
             else
                stop 'angmom3'
             endif
          else
             locnum=2*angmom+1
             if(angmom.eq.1)then
                coffindex=1
             elseif(angmom.eq.2)then
                coffindex=2
             elseif(angmom.eq.3)then
                coffindex=5
             else
                stop 'angmom'
             endif
          endif! l_planetwo
         endif !angmom

         mqnum=2*angmom+1
         ll1=angmom*angmom-1
         loccoffs(:,:)=cmplx(0.0,0.0)
c**************************************************
c        Write the expansion of the locs
c        in terms of spherical harmonics into the
c        loccoffs-array.
c        First index of loccoffs: m-quantum-number.
c        Second index of loccoffs: index of loc.
c**************************************************
         if(l_verbose)print*,"locnum=",locnum
         do locm=1,locnum
          pos=bas+locm
          vec=kveclo(pos)
!The vector of the planewave onto which the loc is matched.
          fk(1)=bkpt(1)+k1(vec)
          fk(2)=bkpt(2)+k2(vec)
          fk(3)=bkpt(3)+k3(vec)

          tmk=tpi*(k1(vec)*taual(1,natom)+
     +             k2(vec)*taual(2,natom)+
     +             k3(vec)*taual(3,natom))

          phase = cmplx(cos(tmk),sin(tmk))
          CALL cotra3(fk,fkp,bmat)
          CALL ylm4(3,fkp,ylm)
          do locmm=1,mqnum
             lmm=ll1+locmm
             loccoffs(locmm,locm)=ci**angmom*phase*conjg(ylm(lmm+1))
          enddo
         enddo !locm
c*************************************************************
c        If locs are entangled by symmetry, the number of locs
c        that have to be treated at the same time is larger.
c        => l_spacetwo: two locs, which are NOT in the
c        xy-plane and which are not related by inversion, are
c        entangled by z-reflexion.
c        => l_spacefour: two pairs the partners of which are
c        connected by inversion are interrelated by z-reflection.
c*************************************************************
         if(l_spacetwo.or.l_spacefour)then
          do locm=1,locnum
           pos=basindex(zrefatom,loc)-1+locm
           vec=kveclo(pos)
           fk(1)=bkpt(1)+k1(vec)
           fk(2)=bkpt(2)+k2(vec)
           fk(3)=bkpt(3)+k3(vec)
           tmk = tpi*(k1(vec)*taual(1,zrefatom)+
     +                k2(vec)*taual(2,zrefatom)+
     +                k3(vec)*taual(3,zrefatom))
           phase = cmplx(cos(tmk),sin(tmk))
           CALL cotra3(fk,fkp,bmat)
           CALL ylm4(3,fkp,ylm)
           do locmm=1,mqnum
             lmm=ll1+locmm
             loccoffs(locnum+locmm,locnum+locm)=
     &                 ci**angmom*phase*conjg(ylm(lmm+1))
           enddo
          enddo !locm
          jump(zrefatom)=jump(zrefatom)+locnum
          if(l_spacetwo)locnum=locnum*2
         endif !l_spacetwo
c*************************************************************
c        If locs are entangled by symmetry, loccoffs has a larger
c        number of components for a given loc.
c        l_planetwo => two locs in xy-plane are entangled
c        by inversion symmetry.
c        l_spacetwoinv => two z-reflexion-symmetric atoms
c        are entangled also by inversion.
c        l_crosstwo => two inversion-symmetric atoms
c        are not directly connected by z-reflexion.
c        l_spacefour =>  two pairs the partners of which are
c        connected by inversion are interrelated by z-reflection.
c*************************************************************
         if(l_planetwo.or.l_crosstwo.or.l_spacetwoinv
     &                            .or.l_spacefour) then
          do locm=1,locnum
            do locmm=1,mqnum
               loccoffs(mqnum+locmm,locm)=(-1)**(locmm-1)*
     &             conjg(loccoffs(mqnum+1-locmm,locm))
            enddo! locmm
          enddo ! locm
         endif ! l_planetwo

c******************************************************************
c        If four atoms are entangled, something remains to be added
c        here:
c******************************************************************
         if(l_spacefour)then
          do locm=1,locnum
            do locmm=1,mqnum
              loccoffs(locnum+mqnum+locmm,locnum+locm)=(-1)**(locmm-1)*
     &             conjg(loccoffs(locnum+mqnum+1-locmm,locnum+locm))
            enddo! locmm
          enddo ! locm
          locnum=locnum*2
         endif
c******************************************************************
c        Find linear combinations of locs that are eigenfunctions
c        of the z-reflection operation by solving linear system
c        of equations. Write the transformation into locrec-matrix.
c******************************************************************
         call CPP_LAPACK_cgetrf(locnum,locnum,loccoffs,28,ipiv,info)
         if(info/=0) stop 'locrectify: cgetrf'
         do loopi=1,locnum
            numorec=numorec+1
            kindlocrec(numorec)=coffkind(loopi,coffindex)
            if(l_verbose)print*,"coffkind=",coffkind(loopi,coffindex)
            if(l_invsatswap)then
              if(l_verbose)print*,"invsatswap" 
              reccoffs(1:locnum/2)=coeffi(1:locnum/2,loopi,coffindex)
              reccoffs(locnum/2+1:locnum/4*3)=
     =           coeffi(locnum/4*3+1:locnum,loopi,coffindex)
              reccoffs(locnum/4*3+1:locnum)=
     =           coeffi(locnum/2+1:locnum/4*3,loopi,coffindex)
            else
              reccoffs(1:locnum)=coeffi(1:locnum,loopi,coffindex)
            endif
            call CPP_LAPACK_cgetrs('N',locnum,1,loccoffs,28,ipiv,
     &              reccoffs,28,info)
            if(info/=0) stop 'locrectify: cgetrs'
#ifdef CPP_INVERSION
            if(any(abs(imag(reccoffs(1:locnum))).gt.1e-6)) then
!In the real version of the code, the hamiltonian and overlap
!matrices are real. Hence the transformation matrices that make
!the locs become eigenfunctions of the z-reflection operation
!have to be real also. Check this here. 
               print*,"Sorry!"
               print*,"The transformation is not purely real."
               print*,"=> Something is wrong here."
               stop 'locrectify: unfortunately complex'
            endif
            if(l_spacefour)then
               if(l_verbose)print*,"spacefour => shifts in bas"
               do locmm=1,locnum/2
                locrec(bas+locmm,numorec)=real(reccoffs(locmm))
               enddo
               if(l_verbose)print*,"basindex=",basindex(zrefatom,loc)
               do locmm=1,locnum/2
                 if(l_verbose)print*,basindex(zrefatom,loc)-1+locmm
                 locrec(basindex(zrefatom,loc)-1+locmm,numorec)=
     &                  real(reccoffs(locmm+locnum/2))
               enddo
            else
              do locmm=1,locnum
                locrec(bas+locmm,numorec)=real(reccoffs(locmm))
              enddo
            endif
#else
            if(l_spacetwo)then
               if(l_verbose)then
                  print*,"spacetwo => shifts in bas"
                  print*,"bas=",bas
                  print*,"locnum=",locnum
                  print*,"basindex=",basindex(zrefatom,loc)
               endif
               locrec(bas+1:bas+locnum/2,numorec)=reccoffs(1:locnum/2)
               locrec(basindex(zrefatom,loc):
     &                 basindex(zrefatom,loc)+locnum/2-1,numorec)=
     &                   reccoffs(1+locnum/2:locnum)
            else
              locrec(bas+1:bas+locnum,numorec)=reccoffs(1:locnum)
            endif
#endif
         enddo !loopi
         if(l_spacetwo) locnum=locnum/2
         if(l_spacefour) locnum=locnum/2
         bas=bas+locnum
        enddo !loc
       enddo !natom
       nalo=naup+1
      enddo !ind


c*********************************************************
c     Determine number of even and odd locs.
c     Prepare sort-arrays.
c*********************************************************
      evenlocs=0
      oddlocs=0
      do loopi=1,nlotot
         if(kindlocrec(loopi)) then
            evenlocs=evenlocs+1
            evenlocindex(evenlocs)=loopi
         else
            oddlocs=oddlocs+1
            oddlocindex(oddlocs)=loopi
         endif
      enddo
      if(l_verbose)then
        print*,"evenlocs=",evenlocs
        print*,"oddlocs=",oddlocs
        print*,"evenlocindex"
        print*,evenlocindex(1:evenlocs)
        print*,"oddlocindex"
        print*,oddlocindex(1:oddlocs)
      endif
      if(l_veryverbose)then
        do info=1,nlotot 
          if(kindlocrec(info))then
            print*,"even orbital"
          else
            print*,"odd orbital"
          endif
          print*,locrec(1:nlotot,info)
        enddo
      endif
      end subroutine
      end module

