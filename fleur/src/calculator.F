      MODULE m_calculator
      !
      !  This module implements a parser able to evaluate expressions in
      !  input files
      !

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: evaluate,ASSIGN_var,delete_vars,evaluatefirst
     $     ,makenumberstring


      CHARACTER(len = 10),SAVE,ALLOCATABLE :: var_names(:)
      REAL,ALLOCATABLE,SAVE                :: var_values(:)
      INTEGER,SAVE                         :: n_vars

      CONTAINS
      !<-- S: priv_error(m)

      SUBROUTINE priv_error(m)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
#ifdef CPP_MPI
      include 'mpif.h'
#endif
      !<--Arguments
      CHARACTER(len =*), INTENT(in) :: m
      !>
      !<-- Locals
      INTEGER             :: ierr,n
      !>
      WRITE(*,*) m
      WRITE(6,*) m

      WRITE(*,*) "Defined variables:"

      DO n = 1,n_vars
         WRITE(*,"(a20,'=',f0.10)") var_names(n),var_values(n)
      ENDDO


#ifdef CPP_MPI
      CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#endif

      STOP "Error in expression"
      END SUBROUTINE

      !> 
      !<-- S: priv_increase_storage()

      SUBROUTINE priv_increase_storage()
!-----------------------------------------------
!    increase the storage for the variables if needed
!    This is a very slow procedure, so be careful to
!    adjust the buffer size if it is called very often
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<-- Locals

      CHARACTER(len = 10),ALLOCATABLE :: tmp_names(:)
      REAL   ,ALLOCATABLE            :: tmp_values(:)

      INTEGER,PARAMETER :: min_buffer = 5 

      !>
      IF (ALLOCATED(var_names)) THEN
         IF (n_vars+1 < SIZE(var_names)) RETURN  !nothing to be done
      ENDIF
      !<-- copy old data
      IF (ALLOCATED(var_names)) THEN
         ALLOCATE(tmp_names(SIZE(var_names)))
         tmp_names = var_names
         ALLOCATE(tmp_values(SIZE(var_values)))
         tmp_values = var_values
         DEALLOCATE(var_values,var_names)
      ENDIF
      !>
      
      ALLOCATE(var_values(n_vars+min_buffer))
      ALLOCATE(var_names(n_vars+min_buffer))
      !<-- Copy data back
      IF (ALLOCATED(tmp_names)) THEN
         var_names(:SIZE(tmp_names)) = tmp_names
         var_values(:SIZE(tmp_values)) = tmp_values
         DEALLOCATE(tmp_names,tmp_values)
      ENDIF
      !>


      END SUBROUTINE

      !> 
      !<-- S: delete_vars()
      SUBROUTINE delete_vars()
!-----------------------------------------------
!    
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      IF (ALLOCATED(var_names)) THEN
         DEALLOCATE(var_names)
         DEALLOCATE(var_values)
         n_vars = 0
      ENDIF
      END SUBROUTINE
      !> 
      !<-- S: assign_var(var,value)

      SUBROUTINE ASSIGN_var(var,value)
!-----------------------------------------------
!  assign a value to a variable
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len = *),INTENT(IN) :: var
      REAL   ,INTENT(IN)            :: value
      !>
      !<-- Locals
      INTEGER             :: n
      CHARACTER(len = 10) :: s
      !>
      s = TRIM(ADJUSTL(var))
      DO n = 1,n_vars
         IF (s == var_names(n)) THEN
            !variable exists, new value assigned
            var_values(n) = value
            RETURN
         ENDIF
      ENDDO
      !this is a new variable
      CALL priv_increase_storage()
      n_vars = n_vars+1
      var_names(n_vars) = s
      var_values(n_vars) = value

      END SUBROUTINE

      !> 
      !<-- F: priv_number(string) result(number)
      FUNCTION priv_number(string) result(number)
!-----------------------------------------------
!    read the first part of string as a number
!             (last modified: 06-04-11) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<-- Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      REAL                             :: number
      !>
      !<-- Locals
      INTEGER             :: pos
      LOGICAL             :: dot
      
      !>
      
    
      number = 0.0
      pos = 0
      dot = .false.
      loop:DO
         pos = pos+1
         IF (pos>LEN(string)) THEN
            EXIT
         ENDIF
         SELECT CASE(string(pos:pos))
             CASE ('0':'9') 
                CYCLE loop
             CASE ('+':'-')
                IF(pos /= 1) EXIT loop
             CASE ('.') 
                IF(dot) THEN
                   EXIT loop
                ELSE
                   dot = .TRUE.
                ENDIF
             CASE default
                EXIT loop
         END SELECT
      ENDDO loop
      
      IF (pos == 0) RETURN 
      READ(string(:pos-1),*)number
      IF (pos>LEN(string)) THEN
         string=' '
      ELSE
         string = string(pos:)
      ENDIF
      END function
      !>
      !<-- S: priv_text(string, command, number)

      SUBROUTINE  priv_text(string, command, number)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      CHARACTER(len = 10), INTENT(out) :: command
      REAL, INTENT(out)                :: number
      !>
      !<-- Locals
      INTEGER  :: n,l
      CHARACTER(len = 5),DIMENSION(13), PARAMETER :: commands = (/
     $     'cos( ','sin( ','tan( ','exp( ','log( ','abs( ','sqrt('
     $     ,'acos(','asin(','atan(','cosh(','sinh(','tanh('/)
      !>
      command = ' '
      number = 0.0

      !<--check if this is a function call
      DO n = 1,SIZE(commands)
         l = len_TRIM(commands(n))
         IF (len_TRIM(string) <= l) CYCLE
         IF (string(:l) == commands(n)(:l)) THEN
            command = commands(n)(:l)
            string=string(l:)
            RETURN
         ENDIF
      ENDDO
      !>
      !<-- it must be a variable
      DO n = 1,n_vars
         l = len_TRIM(var_names(n))
         IF (len_TRIM(string)<l) CYCLE
         IF (string(:l) == var_names(n)(:l)) THEN
            command ="variable"
            number  = var_values(n)
            IF (len_TRIM(string)>l) THEN
               string = string(l+1:)
            ELSE
               string=' '
            ENDIF
            RETURN
         ENDIF
      ENDDO
      !>
      CALL priv_error("Unknown character string found: "//TRIM(string
     $     ))
      END SUBROUTINE

      !> 
      !<-- F: priv_operator(string)

      FUNCTION priv_operator(string)result(command)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      CHARACTER(len = 10)              :: command
      !>
      !<-- Locals
      CHARACTER(len = 2),DIMENSION(5), PARAMETER ::  ops = (/ '**', '* '
     $     ,'/ ', '+ ','- ' /)
      INTEGER             :: n,l
      !>

      command = ' '
      IF (len_TRIM(string)<2) CALL priv_error("Parsing error (operator
     $     ):"//trim(string))
      DO n = 1,SIZE(ops)
         l = len_TRIM(ops(n))
         IF (string(:l) == ops(n)(:l)) THEN
            command = ops(n)(:l)
         ENDIF
      ENDDO
      
      IF (command == ' ') CALL priv_error('Unknown operator: '
     $     //TRIM(string))
      
      string = string(len_TRIM(command)+1:)
      END FUNCTION

      !> 
      !<-- F: priv_bracket(string)

      FUNCTION priv_bracket(string)RESULT(substring)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      CHARACTER(len = LEN(string)) :: substring
      !>
      !<-- Locals
      INTEGER             :: n,pos,count
      !>
      pos = 0
      count = 0
      loop:DO
         pos = pos+1
         IF (pos>len_TRIM(string)) CALL priv_error("Unbalanced brackets"
     $        )
         SELECT CASE(string(pos:pos))
            CASE('(')
              count=count+1
            CASE(')')
              count = count-1
         END SELECT
         IF (count == 0) EXIT loop
      ENDDO loop
      substring = TRIM(string(2:pos-1))
      IF (len_TRIM(substring) == 0) CALL priv_error
     $     ("Empty brackets found")
      IF (len_TRIM(string) < pos+1) THEN
         string = ' '
      ELSE
         string = string(pos+1:)
      END IF
      END FUNCTION

      !> 
      !<-- F: priv_calc(string,command,value)

      RECURSIVE FUNCTION priv_calc(string,value,command)
     $     RESULT(number)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      CHARACTER(len = 10), INTENT(in) :: command
      REAL, INTENT(in)                :: value
      REAL                            :: number
      !>
      !<-- Locals
      CHARACTER(len = 10) :: nextcommand
      REAL                :: nextnumber
      !>
      SELECT CASE(command)
         CASE('number','variable','end')
            number = value
         CASE('bracket')
            number = evaluate(priv_bracket(string))
         CASE('+')
            number = value + priv_evaluateblock(string, command)
         CASE('-')
            number = value - priv_evaluateblock(string, command)
         CASE('*')
            number = value * priv_evaluateblock(string, command)
         CASE('/')
            number = priv_evaluateblock(string, command)
            IF (number == 0.0) CALL priv_error("Divide by zero")
            number = value/number
         CASE('**')
            number = priv_evaluateblock(string, command)
            IF (number > 0.0) THEN
               number = value ** number
            ELSEIF (value < 0.0) THEN
               IF (INT(number) == number) THEN
                  number = value ** INT(number)
               ELSE
                  CALL priv_error('x**y, x<0 and y not integer')
               END IF
            ELSEIF (number /= 0.0) THEN
               number = 0.0
            ELSE             
               CALL priv_error('Undefined result 0**0')
            END IF     
         CASE('cos(','sin(','exp(','log(','abs(','sqrt(','acos(','
     $           asin(','atan(','cosh(','tanh(','tan(')
            call priv_getnextatom(string, nextnumber, nextcommand)
            number = priv_calc(string,  nextnumber, nextcommand)
            SELECT CASE (command)
               CASE('sin(')
                  number = SIN(number)
               CASE('cos(')
                  number = COS(number)
               CASE('tan(')
                  number = TAN(number)
               CASE('exp(')
                  number = EXP(number)
               CASE('log(')
                  IF (number <= 0) CALL priv_error("log(x),x <= 0 ")
                  number = LOG(number)
               CASE('abs(')
                  number = ABS(number)
               CASE('sqrt(')
                  IF (number < 0) CALL priv_error("sqrt(x),x <0 ")
                  number = SQRT(number)
               CASE('acos(')
                  IF (ABS(number)>1) CALL priv_error("acos(x), |x|>1")
                  number = ACOS(number)
               CASE('asin(')
                  IF (ABS(number)>1) CALL priv_error("asin(x), |x|>1")
                  number = ASIN(number)
               CASE('atan(')
                  IF (ABS(number)>1) CALL priv_error("atan(x), |x|>1")
                  number = ATAN(number)
               CASE('cosh(')
                  number = COSH(number)
               CASE('sinh(')
                  number = SINH(number)
               CASE('tanh(' )
                  number = TANH(number)
            END SELECT
         CASE default
            CALL priv_error("Parsing error:"//command)
      END SELECT

      END FUNCTION

      !> 
      !<-- S: priv_getnextatom(func, number, command)

      SUBROUTINE priv_getnextatom(string, number, command)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      REAL, INTENT(inout)              :: number
      CHARACTER(len = 10), INTENT(out) :: command
      !>
      
      SELECT CASE(string(1:1))
      CASE('0':'9', '.')
         number = priv_number(string)
         command = 'number'
      CASE('+', '-', '/', '*')
         command = priv_OPERATOR(string)
      CASE('a':'z','A':'Z')
         CALL priv_text(string, command, number)
      CASE('(')
         command = 'bracket'
      CASE default
         IF (len_TRIM(string) > 1) CALL
     $        priv_error('Unknowninput:'//TRIM(string))
         command = 'end'
      END SELECT
      
      END SUBROUTINE

      !> 
      !<-- F: priv_order(command) 
      FUNCTION priv_order(command) result(order)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len =*), INTENT(in) :: command
      INTEGER                       :: order
      !>
      order = 0
      
      SELECT CASE(TRIM(command))
         CASE('+','-')
            order = 10
         CASE('*','/')
            order = 100
         CASE('**')
            order = 1000
      END SELECT

      END FUNCTION
      !> 
      !<-- S: priv_peeknextatom(string, number, command)
      SUBROUTINE priv_peeknextatom(string, number, command)
!-----------------------------------------------
!
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len =*), INTENT(in) :: string
      REAL, INTENT(inout)           :: number
      CHARACTER(len = 10), INTENT(inout) :: command
      !>
      !<-- Locals
      CHARACTER(len   = LEN(string)) :: s
      !>
      s=string
      CALL priv_getnextatom(s,number,command)
      END SUBROUTINE
      !> 
      !<-- F: priv_evaluateblock(func, blockcommand) result(number)
      RECURSIVE FUNCTION priv_evaluateblock(string, blockcommand)
     $     RESULT(number)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: string
      CHARACTER(len=*), INTENT(in)     :: blockcommand
      REAL                             :: number
      !>
      !<-- Locals
      CHARACTER(len = 10) :: command, nextcommand
      REAL                :: nextnumber
      !>

      CALL priv_peeknextatom(string, nextnumber, nextcommand)
      
      IF (TRIM(nextcommand) == 'end') CALL priv_error
     $     ('Premature end of function')
      
      DO WHILE((priv_order(nextcommand) == 0 .OR.
     $     (priv_order(nextcommand) >priv_order(blockcommand))) .AND.
     $     TRIM(nextcommand) /= 'end')
         CALL priv_getnextatom(string, number, command)
         number = priv_calc(string, number, command)
         CALL priv_peeknextatom(string, nextnumber, nextcommand)
      END DO
      END FUNCTION
      !> 
      !<-- F: evaluate(s) 
      RECURSIVE FUNCTION evaluate(s) RESULT(number)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len =*), INTENT(in  ) :: s
      REAL                            :: number
      !>
      !<-- Locals
      CHARACTER(len = 10) :: command
      CHARACTER(len = len_trim(s)) :: tmp_s
      
      !>
      tmp_s  = TRIM(ADJUSTL(s))
      number = 0
      command = ' '
      DO WHILE(command /= 'end')
         CALL priv_getnextatom(tmp_s, number, command)
         number = priv_calc(tmp_s, number, command)
      END DO
      END FUNCTION
      !> 
      !<-- F: makenumberstring(x)
      FUNCTION makenumberstring(x)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      REAL   ,INTENT(IN)     :: x
      CHARACTER(len=20)      :: makenumberstring
      !>
      !<-- Locals
      INTEGER             :: n,xx

      !>
      
      !<-- check in no is simple fraction
      IF (FRACTION(x)>0) THEN
         IF (FRACTION(x) == 1) THEN
            WRITE(makenumberstring,"(i0)") INT(x)
         ELSE
            WRITE(makenumberstring,"(i0,'/',i0)") INT(x*FRACTION(x))
     $           ,FRACTION(x)
         ENDIF
         RETURN
      ENDIF
      !>
      !<-- check if it is the sqrt-root of something
      xx = x**2
       IF (FRACTION(x)>0) THEN
          IF (FRACTION(x) == 1) THEN
             WRITE(makenumberstring,"('sqrt(',i0,')')") x
         ELSE
            WRITE(makenumberstring,"('sqrt(',i0,'/',i0,')')") INT(x
     $           *FRACTION(x)),FRACTION(x)
         ENDIF
         RETURN
      ENDIF
      !>
      !<-- check if it is the tan of something
      xx = atan(x)
       IF (FRACTION(x)>0) THEN
          IF (FRACTION(x) == 1) THEN
             WRITE(makenumberstring,"('tan(',i0,')')") x
         ELSE
            WRITE(makenumberstring,"('tan(',i0,'/',i0,')')") INT(x
     $           *FRACTION(x)),FRACTION(x)
         ENDIF
         RETURN
      ENDIF
      !>

      !ok nothing found
      WRITE(makenumberstring,"(f15.9)") x

      CONTAINS 
      FUNCTION FRACTION(x)
      IMPLICIT NONE
      REAL,INTENT(IN) ::x
      INTEGER         :: fraction,n
      fraction=0
      !<-- Check if number can be represented by a simple fraction
      DO n = 1,25
         IF (ABS((x*n)-INT(x*n))<1E-6) THEN
            fraction = n
            RETURN
         ENDIF
      ENDDO
      !>
      END FUNCTION
      END FUNCTION
      !> 
      !<-- F: evaluateFirst(string)
      FUNCTION evaluateFirst(s,n)result(number)
!-----------------------------------------------
!
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      CHARACTER(len =*), INTENT(inout) :: s
      INTEGER,OPTIONAL                 :: n
      REAL                             :: number
      !>
      !<-- Locals
      INTEGER             :: pos
      !>
      s = TRIM(ADJUSTL(s))
      IF (len_TRIM(s) == 0) THEN
         number = 0
         RETURN
      ENDIF
      pos = INDEX(s," ")
      IF (pos == 0) pos = LEN(s)
      IF (PRESENT(n)) pos = MAX(pos,n)
      number  = evaluate(s(:pos))
      IF (pos<LEN_TRIM(s)) THEN
         s = s(pos:)
      ELSE
         s =" "
      ENDIF


      END FUNCTION
      !> 


      END MODULE m_calculator
