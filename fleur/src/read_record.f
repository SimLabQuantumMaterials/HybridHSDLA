      MODULE m_readrecord
!***********************************************************************
!     reads in the next input section which is not an empty line and not
!     a comment. 
!     Section can be
!     
!     either
!     ---> a single line of input
!     or 
!     ---> a fortran name list ( & ... / ) 
!          & Must be the first non-space charcter in a new line. 
!          & ... / can extend over many lines, can contain comments
!          and empty lines. Maximum length is xl_buffer.
!
!     End of line comments in a section are removed.    
!     A comment is everything to the right of the first ! in a line.
!     
!***********************************************************************
      CONTAINS
      SUBROUTINE read_record(
     >                       infh,xl_buffer,
     X                       nline, 
     <                       nbuffer,buffer,ios )

      IMPLICIT NONE

      INTEGER, INTENT (IN)    :: infh            ! input filehandle (5)
      INTEGER, INTENT (IN)    :: xl_buffer       ! maximum length of read record
      INTEGER, INTENT (INOUT) :: nline           ! in: last line read ; on output new read lines added
      INTEGER, INTENT (OUT)   :: nbuffer, ios    ! read buffer & I/O status
      CHARACTER(len=xl_buffer), INTENT (OUT)   :: buffer

      INTEGER           :: i,j,l,n,nw
      LOGICAL           :: building, complete
      CHARACTER(len=xl_buffer)  :: line

      INTEGER, PARAMETER  :: dbgfh=6, errfh=6, bfh=93
!---> initialize some variables

      building = .false.
      complete = .false.

!===> read input

      loop: DO 

        nline = nline + 1
        READ( infh,'(a)',ERR=911,END=999,IOSTAT=ios) line
        LINE = adjustl(line)
        WRITE(dbgfh,'("line:",i2,">",a71)') nline,line(1:71)

        n = SCAN(line,'!')                 ! remove end of line comments

        IF ( n>0 ) THEN
            line = line(1:n-1)
        ENDIF

        n = LEN_TRIM( line )               ! length of line without trailing blanks
        IF ( n == 0 ) CYCLE loop

        IF ( line(1:1)=='&' ) THEN         ! check if beginning of namelist
          IF (building) THEN
            WRITE (errfh,*) 
     &      'missing end of namelist marker / in  or before line', nline 
            STOP 'read_record: ERROR'
          ENDIF
          building = .true.
          buffer = line
          nbuffer = n
          if( line(n:n)=='/' ) complete = .true.

        ELSEIF ( line(n:n)=='/' ) THEN     ! check if end of namelist
          IF (building) THEN
            complete = .true.
            buffer = buffer(1:nbuffer)//' '//line
            nbuffer = nbuffer + 1 + n
          ELSE
            WRITE (errfh,*) 
     &           'out of place end of namelist marker / in line', nline 
            STOP 'read_record: ERROR'
          ENDIF

        ELSEIF ( building ) THEN           ! add line to buffer
          buffer = buffer(1:nbuffer)//' '//line
          nbuffer = nbuffer + 1 + n

        ELSEIF ( n > 0 ) THEN              ! check for non empty lines outside of namelists
          buffer = line
          nbuffer = n
          complete = .true.
        ENDIF

        IF ( complete ) THEN
!dbg      WRITE (dbgfh,'("buffer=>",a71)') buffer(1:71)
          EXIT
        ENDIF
!===> 
      END DO loop

! internal file / namelist fix
      REWIND ( bfh )
      WRITE (bfh,*) buffer
      REWIND ( bfh )
! internal file / namelist fix

      ios = 0
      RETURN

 911  CONTINUE
      WRITE (errfh,*) 'lapw_input: ERROR reading input. ios  =',ios,
     &               ', line =',nline
      STOP 'read_record: ERROR reading input'

 999  CONTINUE

      ios = 1 
      IF ( building ) THEN
        ios = 2
      ELSE
        buffer = '&end /'
      ENDIF
      RETURN

      END SUBROUTINE read_record
      END MODULE m_readrecord
