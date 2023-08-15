!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/parse.f90,v 1.4 2022/04/05 17:40:19 malcubi Exp $

! ****************************
! ***   PARAMETER PARSER   ***
! ****************************

  subroutine parse(rank,parfile)

  implicit none

  logical equal
  !logical space,comma

  integer rank
  integer i,j
  integer lword,lvar,lvalue
  integer nvalues
  integer(kind=2) nline

! Arrays for location of commas and names of multiple values.

  integer, allocatable, dimension (:) :: commas
  character(50), allocatable, dimension (:) :: values

! The strings (a,first) have always length 1.

  character(1) a,first

! The string "var" is the name of a parameter, I assume
! a length of 50 to be on the safe side.

  character(50) var

! The strings (line,word,value) are allowed to be very long, since
! we can have a long list of arrays to be output.

  character(1000) line,word,value

! The string "parfile" is passed as an argument to the subroutine,
! so it has assumed length.

  character(len=*) parfile


! *********************
! ***   OPEN FILE   ***
! *********************

! Open parfile.  On error jump ahead.

  open(10,FILE=parfile,STATUS='old',ERR=300)


! **********************
! ***   INITIALIZE   ***
! **********************

! Message.

  if (rank==0) then
     print *, 'Reading parameters ...'
  end if

! Start line counter.

  nline = 0


! *********************
! ***   READ LINE   ***
! *********************

! Read a line.  On end of file jump ahead.

  100 read(10,"(A)",END=200) line

! Eliminate possible spaces at start of line and find
! first character to see if the line is a comment.

  line = adjustl(line)
  first = line(1:1)

! Increment line counter.

  nline = nline + 1


! ****************************
! ***   NOW PROCESS LINE   ***
! ****************************

! Ignore comments or blank lines.

  if ((line.ne." ").and.(first.ne."#")) then

!    Eliminate possible comments at end of line
!    and save whatever is left in the string "word".

     j = index(line,'#')

     if (j.gt.0) then
        word = line(1:j-1)
     else
        word = line
     end if

!    Find length of "word" (without trailing spaces).

     lword = len_trim(word)


!    **********************************
!    ***   FIND NAME OF PARAMETER   ***
!    **********************************

!    Find name of parameter and its length by finding an
!    equal sign and taking whatever is in the left of it
!    (minus trailing spaces).

     equal = .false.

     do i=1,lword
        a = word(i:i)
        if (a=="=") then
           equal = .true.
           var = trim(adjustl(word(1:i-1)))
           exit
        end if
     end do

     lvar = len_trim(var)

!    If there was no equal sign anywhere stop.

     if (.not.equal) then
        if (rank==0) then
           print *
           print *, 'Parfile error, missing equal sign in line:',nline
           print *
           print *, 'Aborting! (subroutine parse.f90)'
           print *
        end if
        call die
     end if

!    We don't allow blank spaces in the middle of parameter names.

     do i=1,lvar
        a = var(i:i)
        if (a==" ") then
           if (rank==0) then
              print *
              print *, 'Parfile error, nonsensical parameter name in line:',nline
              print *
              print *, 'Aborting! (subroutine parse.f90)'
              print *
           end if
           call die
        end if
     end do


!    ***********************************
!    ***   FIND VALUE OF PARAMETER   ***
!    ***********************************

!    Find value of parameter and its length by taking
!    whatever was to the right of the equal sign.

     j = index(word,'=')
     value = trim(adjustl(word(j+1:lword)))

     lvalue = len_trim(value)

!    Eliminate possible trailing tabs.  This is necessary
!    since the tabs in the parfile can get stuck to
!    string type variables and make a mess.  Here I use
!    the ASCII value 9 to check for tabs, but I don't know
!    how portable this number is.

     j=lvalue

     do i=lvalue,1,-1
        a = value(i:i)
        if (a==achar(9)) then
           j = i-1
        end if
     end do

     value = value(1:j)
     lvalue = len_trim(value)

!    If value is empty stop.

     if (value==" ") then
        if (rank==0) then
           print *
           print *, 'Parfile error, missing value of parameter ''',trim(var),''' in line:',nline
           print *
           print *, 'Aborting! (subroutine parse.f90)'
           print *
        end if
        call die
     end if


!    ***************************************
!    ***   DEAL WITH COMMAS AND SPACES   ***
!    ***************************************

!    The first and last characters are not allowed to be commas.

     if (value(1:1)==",") then
        if (rank==0) then
           print *
           print *, 'Parfile error, comma at start of parameter ''',trim(var),''' in line:',nline
           print *
           print *, 'Aborting! (subroutine parse.f90)'
           print *
        end if
        call die
     end if

     if (value(lvalue:lvalue)==",") then
        if (rank==0) then
           print *
           print *, 'Parfile error, comma at end of parameter ''',trim(var),''' in line:',nline
           print *
           print *, 'Aborting! (subroutine parse.f90)'
           print *
        end if
        call die
     end if

!    Find out how many values we have.

     nvalues = 1

     do i=1,lvalue
        a = value(i:i)
        if (a==",") then
           nvalues = nvalues + 1
        end if
     end do

!    Identify comma positions.

     allocate(commas(0:nvalues))

     commas(0) = 0
     commas(nvalues) = lvalue+1

     j = 0

     do i=1,lvalue
        a = value(i:i)
        if (a==",") then
           j = j+1
           commas(j)=i
        end if
     end do

!    Now read values and eliminate spaces.

     allocate(values(1:nvalues))

     do i=1,nvalues
        values(i) = value(commas(i-1)+1:commas(i)-1)
        values(i) = trim(adjustl(values(i)))
     end do

!    Reconstruct the value with no spaces.

     value = values(1)

     do i=2,nvalues
        value = trim(value)//","//trim(values(i))
     end do

     lvalue = len_trim(value)

!    Deal with spaces in the middle of values.

     do i=1,nvalues
        do j=1,len_trim(values(i))
           a = values(i)(j:j)
           if (a==" ") then
              if (rank==0) then
                 print *
                 print *, 'Parfile error, non-sensical value for parameter ''',trim(var),''': ',values(i)
                 print *
                 print *, 'Aborting! (subroutine parse.f90)'
                 print *
              end if
              call die
           end if
        end do
     end do

!    Deallocate array "commas".

     deallocate(commas)


!    ************************
!    ***   ASSIGN VALUE   ***
!    ************************

!    If we got here we gave a parameter name and its value
!    in the strings (var,value).  Now we call assign routine.

     !print *, var,value
     call assign(nline,nvalues,var,value,values)

!    Now deallocate array "values".

     deallocate(values)

  end if


! **************************
! ***   READ NEXT LINE   ***
! **************************

! Jump back to read next line.

  goto 100


! ***********************
! ***   END OF FILE   ***
! ***********************

! End of file has been found, close file and return.

  200 continue

  close(10)

  if (rank==0) then
     print *, 'Finished reading parameters'
     print *
  end if

  return


! *****************
! ***   ERROR   ***
! *****************

! If we get here there was an error opening the file.

  300 continue

  if (rank==0) then
     print *
     print *, 'There was an error opening the parfile'
     print *, 'Are you sure it exists?'
     print *
     print *, 'Aborting! (subroutine parse.f90) '
     print *
  end if

  call die


! ***************
! ***   END   ***
! ***************

  end subroutine parse




