      program test

      real t,f,tt
      integer i
      
      t = 2.0/3.0
      tt = -2.0/3.0

      do i = 1,9

       f = abs(i*t)-int(abs(i*t)) + abs(i*tt)-int(abs(i*tt)) 
       write(*,*) f
      enddo

      end
