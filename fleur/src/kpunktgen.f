      program kpunktgen
      implicit none
      integer num,numkpt,dim
      integer c1,c2,c3,limit
      real ilen,shift
      real i1,i2,i3
      logical l_shift


      print*,"specify dimension [2/3]"
      read(*,*)dim
      print*,"symmetric to origin [t/f]?"
      read(*,*)l_shift
      if(dim==3)then
      print*,"create three dimensional k-point set"   
      print*,"Zahl der Unterliederungen eingeben [integer]"
      read(*,*)num
      print*,"Zahl der Untergliederungen: ",num
      numkpt=num**3
      print*,"Zahl der k-Punkte: ",numkpt
      ilen=1.0/num
      print*,"Intervallaenge: ",ilen
      open(100,file='kpts',form='formatted',status='unknown')
      open(200,file='wannierkpts',form='formatted',status='unknown')
      write(100,'(2x,i3,8x,f7.5)')numkpt,1.0
      limit=num-1
      if (l_shift) then
         shift=limit*ilen/2.0        
      endif

      do c1=0,limit
         do c2=0,limit
            do c3=0,limit
               i1=ilen*c1-shift
               i2=ilen*c2-shift
               i3=ilen*c3-shift
               write(100,'(3x,f7.5,3x,f7.5,3x,f7.5,3x,f7.5)')i1,i2,i3,1.00000
               write(200,'(3x,f7.5,3x,f7.5,3x,f7.5)')i1,i2,i3
            enddo
         enddo
      enddo
      close(100)
      close(200)

      elseif(dim==2)then
      print*,"create two dimensional k-point set"
         print*,"Zahl der Unterliederungen eingeben [integer]"
      read(*,*)num
      print*,"Zahl der Untergliederungen: ",num
      numkpt=num**2
      print*,"Zahl der k-Punkte: ",numkpt
      ilen=1.0/num
      print*,"Intervallaege: ",ilen
      open(100,file='kpts',form='formatted',status='unknown')
      open(200,file='wannierkpts',form='formatted',status='unknown')
      write(100,'(2x,i3,8x,f7.5,8x,1a)')numkpt,1.0,"F"
      limit=num-1
      do c1=0,limit
         do c2=0,limit
               i1=ilen*c1
               i2=ilen*c2
               write(100,'(3x,f7.5,3x,f7.5,3x,f7.5)')i1,i2,1.00000
               write(200,'(3x,f7.5,3x,f7.5,3x,f6.5)')i1,i2,0.00000
         enddo
      enddo
      close(100)
      close(200)
      endif
      end program kpunktgen
