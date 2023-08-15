#!/usr/bin/env perl

#$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/prl/assign.pl,v 1.2 2018/08/31 16:16:19 malcubi Exp $

print "CREATING FILE FOR PARAMETER ASSIGNMENT ...\n\n";

# Open outfile assign.f90

open(OUTFILE,">src/auto/assign.f90") or die "Can't open assign.f90: $!";

# Write beginning of file assign.f90

print OUTFILE "! Automatically generated file.  Do not edit!\n\n";
print OUTFILE "  subroutine assign(nline,nvalues,var,value,values)\n\n";
print OUTFILE "  use param\n\n";
print OUTFILE "  use procinfo\n\n";
print OUTFILE "  implicit none\n\n";
print OUTFILE "  logical contains\n";
print OUTFILE "  integer i,nvalues\n";
print OUTFILE "  integer(kind=2) nline\n";
print OUTFILE "  character(50) type\n";
print OUTFILE "  character(1000) range\n";
print OUTFILE "  character(len=*) var,value\n";
print OUTFILE "  character(len=*) values(1:nvalues)\n\n";

# Open infile param.f90

open(INFILE,"src/base/param.f90") or die "Can't open param.f90: $!";

# Parse file param.f90 to identify declared parameters

my $line = " ";
my $nline = 0;

while ($line=<INFILE>) {

   $nline = $nline+1;

#  Ignore comments

   if ($line !~ /^!/) {
  
#     Logical parameters

      if ($line =~ /logical/i) {

         if (($line =~ /::/i) && ($line =~ /=/i)) {

            if ($line =~ /::\s*(\S*)\s*=/) {
	        $var = $1;
            }

            if ($line =~ /multiple/i) {
	       print "\n\n";
	       die "Error in param.f90: logical parameter defined as multiple in line ",$nline,"\n\n";
            }

            if ($line =~ /range/i) {
	       print "\n\n";
	       die "Error in param.f90: logical parameter defined with range in line ",$nline,"\n\n";
            }

            print OUTFILE  "  if (var=='",$var,"') then\n";
            print OUTFILE  "     if (nvalues>1) goto 200\n";
            print OUTFILE  "     type = 'LOGICAL'\n";
            print OUTFILE  "     read(value,*,ERR=100) ",$var,"\n";
            print OUTFILE  "     if (nvalues>1) goto 200\n";
            print OUTFILE  "     return\n";
            print OUTFILE  "  end if\n\n";

         } else {
	    print "\n\n";
	    die "Bad parameter assignment in line ",$nline," of file param.f90\n\n";
         }
      }  

#     Integer parameters

      if ($line =~ /integer/i) {

         if (($line =~ /::/) && ($line =~ /=/)) {

            if ($line =~ /::\s*(\S*)\s*=/) {
	       $var = $1;
            }

            if ($line =~ /multiple/i) {
	       print "\n\n";
	       die "Error in param.f90: integer parameter '",$var,"' defined as multiple in line ",$nline,"\n\n";
            }

            if ($line =~ /range/i) {
	       print "\n\n";
	       die "Error in param.f90: integer parameter defined with range in line ",$nline,"\n\n";
            }

            print OUTFILE  "  if (var=='",$var,"') then\n";
            print OUTFILE  "     if (nvalues>1) goto 200\n";
            print OUTFILE  "     type = 'INTEGER'\n";
            print OUTFILE  "     read(value,*,ERR=100) ",$var,"\n";
            print OUTFILE  "     return\n";
            print OUTFILE  "  end if\n\n";

         } else {
	    print "\n\n";
            die "Bad parameter assignment in line ",$nline," of file param.f90\n\n";
         }

      }

#     Real parameters

      if ($line =~ /real/i) {

         if (($line =~ /::/i) && ($line =~ /=/i)) {

            if ($line =~ /::\s*(\S*)\s*=/) {
	       $var = $1;
            }

            if ($line =~ /multiple/i) {
	       print "\n\n";
	       die "Error in param.f90: real parameter '",$var,"' defined as multiple in line ",$nline,"\n\n";
            }

            if ($line =~ /range/i) {
	       print "\n\n";
	       die "Error in param.f90: real parameter defined with range in line ",$nline,"\n\n";
            }

            print OUTFILE  "  if (var=='",$var,"') then\n";
            print OUTFILE  "     if (nvalues>1) goto 200\n";
            print OUTFILE  "     type = 'REAL'\n";
            print OUTFILE  "     read(value,*,ERR=100) ",$var,"\n";
            print OUTFILE  "     return\n";
            print OUTFILE  "  end if\n\n";

         } else {
	    print "\n\n";
	    die "Bad parameter assignment in line ",$nline," of file param.f90\n\n";
         }
      }

#     Character parameters

      if ($line =~ /character/i) {

         if (($line =~ /::/i) && ($line =~ /=/i)) {

            if ($line =~ /::\s*(\S*)\s*=/) {
	       $var = $1;
            }

            print OUTFILE  "  if (var=='",$var,"') then\n";
            print OUTFILE  "     type = 'CHARACTER'\n";
            print OUTFILE  "     $var = value\n";

#           Only string parameters defined as "multiple" can have more than one value.

            if ($line !~ /multiple/i) {
               print OUTFILE  "     if (nvalues>1) goto 200\n";
            }

#           Check range if necessary.

            if ($line =~ /range/i) {
               if ($line =~ /range\s*=\s*\((\S*)\)\s*$/i) {
	          $range = "(".$1.")";
                  $lrange = int(length($range)/80);
                  print OUTFILE  "     range = &\n";
                  $nn = 0;
                  while ($nn < $lrange+1) {
                     print OUTFILE  "     '",substr($range,$nn*80,80),"'";
		     $nn = $nn + 1;
		     if ($nn < $lrange+1) {
			 print OUTFILE  " // &\n";
                     } else {
                         print OUTFILE  "\n";
                     }
		  }
                  print OUTFILE  "     do i=1,nvalues\n";
                  print OUTFILE  "        if ((.not.contains(range,'('//trim(values(i))//')')).and. &\n";
                  print OUTFILE  "            (.not.contains(range,'('//trim(values(i)))).and. &\n";
                  print OUTFILE  "            (.not.contains(range,trim(values(i))//')')).and. &\n";
                  print OUTFILE  "            (.not.contains(range,trim(values(i))))) then\n";
                  print OUTFILE  "           if (rank==0) then\n";
                  print OUTFILE  "              print *\n";
                  print OUTFILE  "              print *, 'Parfile error.'\n";
                  print OUTFILE  "              print *, 'Out of range value ''',trim(values(i)),''' for parameter ''',trim(var),''' in line:',nline\n";
                  print OUTFILE  "              print *\n";
                  print OUTFILE  "              print *, 'Aborting! (subroutine assign.f90)'\n";
                  print OUTFILE  "              print *\n";
                  print OUTFILE  "           end if\n";
                  print OUTFILE  "           call die\n";
                  print OUTFILE  "        end if\n";
                  print OUTFILE  "     end do\n";
	       } else {
	          print "\n\n";
                  print "Bad range assignment in line ",$nline," of file param.f90\n\n";
                  print "Make sure there is an equal sign, parenthesis, and NO spaces\n\n";
	          die
               }
            }

            print OUTFILE  "     return\n";
            print OUTFILE  "  end if\n\n";


         } else {
	    print "\n\n";
	    die "Bad parameter assignment in line ",$nline," of file param.f90\n\n";
         }
      }
   }
}

# Close infile

close(INFILE);

# Write ending of file assign.f90

print OUTFILE "  print *\n";
print OUTFILE "  if (rank==0) then\n";
print OUTFILE "     print *, 'Parfile error, non-existent parameter in line:',nline\n";
print OUTFILE "     print *\n";
print OUTFILE "     print *, 'Aborting! (subroutine assign.f90)'\n";
print OUTFILE "     print *\n";
print OUTFILE "  end if\n";
print OUTFILE "  call die\n";

print OUTFILE "  100 continue\n";
print OUTFILE "  if (rank==0) then\n";
print OUTFILE "     print *\n";
print OUTFILE "     print *, 'There was an error assigning the variable ''',trim(var),''' in line:',nline\n";
print OUTFILE "     print *, 'Are you sure you gave a ',trim(type),' value?'\n";
print OUTFILE "     print *\n";
print OUTFILE "     print *, 'Aborting! (subroutine assign.f90)'\n";
print OUTFILE "     print *\n";
print OUTFILE "  end if\n";
print OUTFILE "  call die\n\n";

print OUTFILE "  200 continue\n";
print OUTFILE "  if (rank==0) then\n";
print OUTFILE "     print *\n";
print OUTFILE "     print *, 'Multiple values not allowed for variable ''',trim(var),''' in line:',nline\n";
print OUTFILE "     print *\n";
print OUTFILE "     print *, 'Aborting! (subroutine assign.f90)'\n";
print OUTFILE "     print *\n";
print OUTFILE "  end if\n";
print OUTFILE "  call die\n\n";

print OUTFILE "  end subroutine assign\n";

# Close outfile.

close(OUTFILE);

