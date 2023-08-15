!$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/src/base/functions.f90,v 1.8 2022/05/10 16:41:39 malcubi Exp $

! ************************************
! ***   VARIOUS USEFUL FUNCTIONS   ***
! ************************************

! Here I define various useful functions:
!
! contains        Checks if a string contains a specific pattern (in a list).
! message         Created a properly formatted message to screen.


! *****************************
! ***   FUNCTION CONTAINS   ***
! *****************************

  logical function contains(string,pattern)

! This function checks if a given string contains a specific pattern
! in a comma or space separated list, that is, the pattern must be
! separated by commas and/or spaces from other elements in the list.

  character(len=*) string,pattern

! Pattern is not in the string.

  if (index(string,pattern)==0) then

     contains = .false.
     return

! The pattern is equal to the string.

  else if (string==pattern) then

     contains = .true.
     !print *,pattern,' EQUAL ',trim(string)
     return

! The pattern is inside the string but not equal to it.
! Then we need to check more carefully.

  else

!    Pattern in between commas.

     if (index(string,","//pattern//",")/=0) then

        contains = .true.
        !print *,pattern,' BETWEEN COMMAS ',trim(string)
        return

!    Pattern at beginning of the string with comma after it.
!    Notice that sometimes the string can start with "(",
!    as in the range for some parameters.

     else if ((index(string,pattern//",")==1).or.(index(string,"("//pattern//",")==1)) then

        contains = .true.
        !print *,pattern,' AT BEGINNING ',trim(string)
        return

!    Pattern at end of the string with comma before it.
!    Notice that sometimes the string can end with ")",
!    as in the range for some parameters.

     else if ((index(string,","//pattern//" ")/=0).or.(index(string,","//pattern//")")/=0)) then

        contains = .true.
        !print *,pattern,' AT END ',trim(string)
        return

!    If we get here the pattern is inside the string, but it is part of a
!    larger pattern that it does not match fully, so we ignore it.
!
!    For example, if we have "mattertype=complexproca" in the parameter
!    file we don't want "contains(mattertype,proca)" to return a true value.

     else

        contains = .false.
        !print *,pattern,' SUBPATTERN ',trim(string)
        return

     end if

  end if

! End.

  end function contains







! ****************************
! ***   FUNCTION MESSAGE   ***
! ****************************

  character(len=*) function message(string1,string2)

! This function creates a properly formatted message
! to be sent out to the screen taking the two input
! strings (string1,string2).

  implicit none

  logical space
  integer i,j,l

  character(len=*) string1,string2
  character(50) work,work2

! Flag for extra space.

  space = .false.

! Find out length of full concatenated string.

  if (string1==" ") then
     l = len(trim(string2))
  else if (string2==" ") then
     l = len(trim(string1))
  else
     l = len(trim(string1)//" "//trim(string2))
  end if

! Length is even, string is OK.

  if (real(l/2)==real(l)/2.d0) then

     j = l/2

     if (string1==" ") then
        work = trim(string2)
     else if (string2==" ") then
        work = trim(string1)
     else
        work = trim(string1)//" "//trim(string2)
     end if

! Length is odd, string needs extra space.

  else

     j = l/2 + 1

     if (string1==" ") then
        space = .true.
        work = trim(string2)//" "
     else if (string2==" ") then
        space = .true.
        work = trim(string1)//" "
     else
        work = trim(string1)//"  "//trim(string2)
     end if

  end if

! Construct message. Notice that the message is constructed
! from right to left

  work2 = "***"

  do i=1,18-j
     work2 = " "//work2
  end do

  if (space) then
     work2 = trim(work)//" "//work2
  else
     work2 = trim(work)//work2
  end if

  do i=1,18-j
     work2 = " "//work2
  end do

  message = "***"//work2

! End.

  end function message


